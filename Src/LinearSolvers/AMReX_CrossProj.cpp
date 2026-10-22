#include <AMReX_CrossProj.H>

#if (AMREX_SPACEDIM == 2)

#include <AMReX_GpuContainers.H>
#include <AMReX_ParReduce.H>
#ifdef AMREX_USE_FFT
#include <AMReX_FFT_LocalR2C.H>
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace amrex::GMG {

// Solves G^T G phi = G^T u (= -D u) by PCG, then u -= G phi. G is the
// cell-centered gradient of nodal phi.
//
// Every operator in the hierarchy is K_x (x) M_y + M_x (x) K_y with 1D
// three-point stencils K and M per direction: the fine one has K = second
// difference and M = [1/4,1/2,1/4], and Galerkin coarsening with
// tensor-product interpolation keeps the form. On a periodic grid all of
// them are constant, so each level is a constant 9-point stencil.
//
// Node data are stored on the cell-centered BoxArray: box [lo,hi) holds
// the nodes it owns (a periodic grid has as many nodes as cells), so no
// node is duplicated and ghost nodes come from FillBoundary.
//
// dx == dy: the fine stencil couples only diagonal neighbors, so the nodes
// with i+j even (red) and odd (black) decouple. Eliminating the (odd,odd)
// nodes of red and the (even,odd) nodes of black is exact and leaves the
// same 9-point problem for both colors on a grid coarsened by 2 (black is
// red shifted by one in x). PCG runs on that reduced problem, stored as two
// components of one MultiFab, preconditioned by a V-cycle.
//
// Otherwise PCG runs on the fine grid, preconditioned by a multigrid cycle
// with a coarse space for the smooth modes (linear interpolation) and one
// for the near-null modes, the checkerboard (-1)^(i+j) times smooth modes.
struct CrossProj::Level
{
    Box domain;
    Periodicity period;
    bool modulate = false; // first checkerboard level: P has a factor (-1)^(i+j)
    int ncomp = 1;
    int cheby_degree = 3;  // 2 on the reduced hierarchy
    BoxArray ba;
    DistributionMapping dm;
    BoxArray cba;      // parent BoxArray coarsened, if ba was regridded
    Array<Array<Real,3>,AMREX_SPACEDIM> K, M;  // [dir][offset+1]
    GpuArray<Real,9> sten{};  // [(b+1)*3+(a+1)] for offset (a,b)
    GpuArray<Real,9> rw{};    // restriction weights from the parent, same layout
    GpuArray<Real,4> pw{};    // prolongation weights [b*2+a] for fine offset (a,b)
    MultiFab x, x2, b, r, t;
    // bottom level: gathered data on the device, inverse symbol (including
    // the normalization) and the transform
    bool binit = false;
    Gpu::DeviceVector<Real> bdev;
    Vector<Real> bg;
#ifdef AMREX_USE_FFT
    // real-to-complex, half spectrum along x; one per component so that
    // the plans are never rebuilt for a different pointer
    Vector<std::unique_ptr<FFT::LocalR2C<Real>>> bfft;
    Gpu::DeviceVector<GpuComplex<Real>> bspec;
    Gpu::DeviceVector<Real> binv;
#else
    Array<Vector<Real>,2> bcos, bsin;           // naive DFT tables
    Vector<Real> bre, bim, binv;
#endif
    Real lmax = Real(0.0); // Chebyshev upper bound
};

namespace {

constexpr int regrid_max_size = 32;
constexpr int bottom_max_nodes = 1024; // exact bottom solve at or below this size
constexpr int nhier = 2; // coarse hierarchies: smooth modes and the checkerboard

// Whether a level with this domain gets a coarser level (by 2), or is the
// bottom: small enough for the exact solve, or not coarsenable.
bool coarsenable (Box const& domain)
{
    if (domain.numPts() <= bottom_max_nodes) { return false; }
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        int const len = domain.length(idim);
        if (len % 2 != 0 || len < 8) { return false; }
    }
    return true;
}

// Galerkin P^T A P of a periodic 1D three-point stencil for linear
// interpolation by a factor of 2, optionally with P modulated by (-1)^m.
Array<Real,3> galerkin1d (Array<Real,3> const& a, bool modulate)
{
    auto w = [] (int p) { return (p == 0) ? Real(1.0) : Real(0.5); };
    Array<Real,3> c{};
    for (int oc = -1; oc <= 1; ++oc) {
        for (int p = -1; p <= 1; ++p) {
            for (int q = -1; q <= 1; ++q) {
                int const p2 = p + q - 2*oc;
                if (p2 < -1 || p2 > 1) { continue; }
                Real const sgn = (modulate && (q & 1)) ? Real(-1.0) : Real(1.0);
                c[oc+1] += w(p) * w(p2) * sgn * a[q+1];
            }
        }
    }
    return c;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real stencil_sum (GpuArray<Real,9> const& s, Array4<Real const> const& x,
                  int i, int j, int k, int n) noexcept
{
    return s[0]*x(i-1,j-1,k,n) + s[1]*x(i,j-1,k,n) + s[2]*x(i+1,j-1,k,n)
         + s[3]*x(i-1,j  ,k,n) + s[4]*x(i,j  ,k,n) + s[5]*x(i+1,j  ,k,n)
         + s[6]*x(i-1,j+1,k,n) + s[7]*x(i,j+1,k,n) + s[8]*x(i+1,j+1,k,n);
}

// r = b - A x
void residual (CrossProj::Level& L)
{
    BL_PROFILE("GMG::CrossProj::residual");
    L.x.FillBoundary(L.period);
    auto const sten = L.sten;
    auto const& ra = L.r.arrays();
    auto const& xa = L.x.const_arrays();
    auto const& ba = L.b.const_arrays();
    ParallelFor(L.r, IntVect(0), L.ncomp, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
    {
        ra[bno](i,j,k,n) = ba[bno](i,j,k,n) - stencil_sum(sten, xa[bno], i, j, k, n);
    });
}

// One step of the three-term Chebyshev recurrence
//   x_new = x + alpha (x - x_prev) + beta D^-1 (b - A x),
// with x in L.x and x_prev in L.x2; x_new overwrites L.x2 and the two are
// swapped. ZERO_X: x is zero (not read). FIRST: no alpha term. ZERO_PREV:
// x_prev is zero.
template <bool ZERO_X, bool FIRST, bool ZERO_PREV>
void cheby_step (CrossProj::Level& L, Real alpha, Real beta)
{
    BL_PROFILE("GMG::CrossProj::cheby_step");
    if (!ZERO_X) { L.x.FillBoundary(L.period); }
    auto const sten = L.sten;
    Real const bd = beta / L.sten[4];
    auto const& xa = L.x.const_arrays();
    auto const& ba = L.b.const_arrays();
    auto const& oa = L.x2.arrays();
    ParallelFor(L.x2, IntVect(0), L.ncomp, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
    {
        amrex::ignore_unused(sten, bd, alpha, xa, ba, oa); // nvcc: capture outside if constexpr
        Real v;
        if constexpr (ZERO_X) {
            v = bd * ba[bno](i,j,k,n);
        } else {
            Real const xc = xa[bno](i,j,k,n);
            v = xc + bd * (ba[bno](i,j,k,n) - stencil_sum(sten, xa[bno], i, j, k, n));
            if constexpr (!FIRST) {
                if constexpr (ZERO_PREV) {
                    v += alpha * xc;
                } else {
                    v += alpha * (xc - oa[bno](i,j,k,n));
                }
            }
        }
        oa[bno](i,j,k,n) = v;
    });
    std::swap(L.x, L.x2);
}

// Chebyshev smoothing (4th kind, degree L.cheby_degree) of x with rhs b.
void cheby (CrossProj::Level& L, int nsweeps, bool zero_init)
{
    for (int sweep = 0; sweep < nsweeps; ++sweep) {
        bool const z0 = zero_init && sweep == 0;
        for (int kk = 1; kk <= L.cheby_degree; ++kk) {
            Real const alpha = Real(2*kk-3) / Real(2*kk+1);
            Real const beta = Real(8*kk-4) / (Real(2*kk+1) * L.lmax);
            if (kk == 1) {
                if (z0) { cheby_step<true,true,false>(L, alpha, beta); }
                else    { cheby_step<false,true,false>(L, alpha, beta); }
            } else if (kk == 2 && z0) {
                cheby_step<false,false,true>(L, alpha, beta);
            } else {
                cheby_step<false,false,false>(L, alpha, beta);
            }
        }
    }
}

// Upper end of the spectrum of D^-1 A by power iteration, with a margin.
void estimate_lmax (CrossProj::Level& L)
{
    auto const& va = L.x.arrays();
    ParallelFor(L.x, IntVect(0), L.ncomp, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
    {
        auto h = static_cast<unsigned int>(i)*73856093U ^ static_cast<unsigned int>(j)*19349663U
               ^ static_cast<unsigned int>(k)*83492791U ^ static_cast<unsigned int>(n)*2654435761U;
        h ^= h >> 13; h *= 0x5bd1e995U; h ^= h >> 15;
        va[bno](i,j,k,n) = Real(h % 1000U) / Real(1000.0) - Real(0.5);
    });
    L.b.setVal(Real(0.0));
    Real const dinv = Real(1.0) / L.sten[4];
    Real lam = Real(0.0);
    for (int it = 0; it < 30; ++it) {
        Real const nrm = std::sqrt(MultiFab::Dot(L.x, 0, L.x, 0, L.ncomp, 0));
        L.x.mult(Real(1.0)/nrm, 0, L.ncomp, 0);
        residual(L); // r = -A x
        lam = dinv * std::sqrt(MultiFab::Dot(L.r, 0, L.r, 0, L.ncomp, 0));
        MultiFab::Copy(L.x, L.r, 0, 0, L.ncomp, 0);
    }
    L.lmax = Real(1.1) * lam;
}

// C.b = P^T F.r (weights C.rw)
void restrict_to (CrossProj::Level& F, CrossProj::Level& C)
{
    BL_PROFILE("GMG::CrossProj::restrict");
    F.r.FillBoundary(F.period);
    MultiFab* dst = C.cba.empty() ? &C.b : &C.t;
    auto const w = C.rw;
    auto const& fa = F.r.const_arrays();
    auto const& ca = dst->arrays();
    ParallelFor(*dst, IntVect(0), C.ncomp, [=] AMREX_GPU_DEVICE (int bno, int I, int J, int k, int n) noexcept
    {
        auto const& f = fa[bno];
        Real s = Real(0.0);
        for (int b = -1; b <= 1; ++b) {
        for (int a = -1; a <= 1; ++a) {
            s += w[(b+1)*3+(a+1)] * f(2*I+a, 2*J+b, k, n);
        }}
        ca[bno](I,J,k,n) = s;
    });
    if (dst != &C.b) {
        C.b.ParallelCopy(C.t, 0, 0, C.ncomp);
    }
}

// F.x += P C.x (weights C.pw). Each coarse node writes the fine nodes
// (2I+a, 2J+b), a,b < 2, from the coarse nodes I..I+1, J..J+1.
void prolong_add (CrossProj::Level& F, CrossProj::Level& C)
{
    BL_PROFILE("GMG::CrossProj::prolong");
    MultiFab const* src = &C.x;
    if (C.cba.empty()) {
        C.x.FillBoundary(C.period);
    } else {
        C.t.ParallelCopy(C.x, 0, 0, C.ncomp, IntVect(0), IntVect(1), C.period);
        src = &C.t;
    }
    auto const w = C.pw;
    auto const& fa = F.x.arrays();
    auto const& ca = src->const_arrays();
    ParallelFor(*src, IntVect(0), C.ncomp, [=] AMREX_GPU_DEVICE (int bno, int I, int J, int k, int n) noexcept
    {
        auto const& cx = ca[bno];
        auto const& fx = fa[bno];
        Real const c00 = cx(I  ,J  ,k,n);
        Real const c10 = cx(I+1,J  ,k,n);
        Real const c01 = cx(I  ,J+1,k,n);
        Real const c11 = cx(I+1,J+1,k,n);
        int const i = 2*I;
        int const j = 2*J;
        fx(i  ,j  ,k,n) += w[0]*c00;
        fx(i+1,j  ,k,n) += w[1]*(c00+c10);
        fx(i  ,j+1,k,n) += w[2]*(c00+c01);
        fx(i+1,j+1,k,n) += w[3]*(c00+c10+c01+c11);
    });
}

void set_stencil (CrossProj::Level& L)
{
    for (int b = 0; b < 3; ++b) {
    for (int a = 0; a < 3; ++a) {
        L.sten[b*3+a] = L.K[0][a]*L.M[1][b] + L.M[0][a]*L.K[1][b];
    }}
}

// Transfer weights of a level from its parent: linear interpolation times
// (-1)^(i+j) on the first checkerboard level.
void set_transfer_weights (CrossProj::Level& C)
{
    for (int b = -1; b <= 1; ++b) {
    for (int a = -1; a <= 1; ++a) {
        Real w = ((a == 0) ? Real(1.0) : Real(0.5)) * ((b == 0) ? Real(1.0) : Real(0.5));
        if (C.modulate && ((a + b) & 1)) { w = -w; }
        C.rw[(b+1)*3+(a+1)] = w;
    }}
    for (int b = 0; b < 2; ++b) {
    for (int a = 0; a < 2; ++a) {
        Real w = Real(1.0) / Real((a+1)*(b+1));
        if (C.modulate && ((a + b) & 1)) { w = -w; }
        C.pw[b*2+a] = w;
    }}
}

// Exact solve of the bottom level, A x = b, with the periodic constant
// stencil: transform the right-hand side, divide by the stencil symbol
// (zero for the null modes), transform back. With AMReX FFT the transform
// is FFT::LocalR2C on the device; otherwise a naive DFT on the host. The
// bottom is one box on one rank (no communication), except when the top
// level itself is the bottom: then the right-hand side is gathered by an
// allreduce and every rank solves.
void bottom_solve (CrossProj::Level& L)
{
    BL_PROFILE("GMG::CrossProj::bottom");
    IntVect const dlo = L.domain.smallEnd();
    int const nn[2] = {L.domain.length(0), L.domain.length(1)};
    int const nx = nn[0], ny = nn[1];
    Long const N = Long(nx)*ny;
    int const nc = L.ncomp;
    bool const local = (L.ba.size() == 1);
    if (local && L.b.local_size() == 0) { return; }

    // The stencil symbol at mode (i,j), in double so that the null modes
    // cancel to roundoff in single precision too.
    auto symbol = [&] (int i, int j) {
        double sym = 0.0;
        for (int b = -1; b <= 1; ++b) {
        for (int a = -1; a <= 1; ++a) {
            double const arg = 2.0*amrex::Math::pi<double>()
                * (double(a*i)/nx + double(b*j)/ny);
            sym += double(L.sten[(b+1)*3+(a+1)]) * std::cos(arg);
        }}
        return sym;
    };
    // The null modes cancel only to the roundoff of the Real stencil.
    auto invert = [] (Vector<Real>& v, Real scale) {
        constexpr Real tol = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-10);
        Real smax = Real(0.0);
        for (auto x : v) { smax = std::max(smax, std::abs(x)); }
        for (auto& x : v) {
            x = (std::abs(x) > tol*smax) ? scale/x : Real(0.0);
        }
    };

    if (!L.binit) {
        L.binit = true;
        L.bdev.resize(N*nc);
        L.bg.resize(N*nc);
        Vector<Real> hinv;
#ifdef AMREX_USE_FFT
        for (int n = 0; n < nc; ++n) {
            L.bfft.push_back(std::make_unique<FFT::LocalR2C<Real>>(L.domain.length()));
        }
        int const sx = L.bfft[0]->spectralSize()[0];
        Long const NS = Long(sx)*ny;
        hinv.resize(NS);
        for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < sx; ++i) {
            hinv[Long(j)*sx + i] = Real(symbol(i,j));
        }}
        invert(hinv, L.bfft[0]->scalingFactor());
        L.bspec.resize(NS);
        L.binv.resize(hinv.size());
        Gpu::copy(Gpu::hostToDevice, hinv.begin(), hinv.end(), L.binv.begin());
        Gpu::streamSynchronize();
#else
        for (int d = 0; d < 2; ++d) {
            int const n = nn[d];
            L.bcos[d].resize(Long(n)*n);
            L.bsin[d].resize(Long(n)*n);
            for (int kk = 0; kk < n; ++kk) {
            for (int m = 0; m < n; ++m) {
                double const arg = 2.0*amrex::Math::pi<double>() * double(Long(kk)*m % n) / double(n);
                L.bcos[d][Long(kk)*n+m] = Real(std::cos(arg));
                L.bsin[d][Long(kk)*n+m] = Real(std::sin(arg));
            }}
        }
        hinv.resize(N);
        for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            hinv[Long(j)*nx + i] = Real(symbol(i,j));
        }}
        invert(hinv, Real(1.0)); // the inverse DFT below normalizes
        L.bre.resize(N);
        L.bim.resize(N);
        L.binv = hinv;
#endif
    }

    // gather the right-hand side into the global array on the device
    int const dlo0 = dlo[0], dlo1 = dlo[1];
    Real* const pd = L.bdev.data();
    {
        Gpu::fillAsync(L.bdev.begin(), L.bdev.end(), [=] AMREX_GPU_DEVICE (Real& v, Long) noexcept { v = Real(0.0); });
        auto const& ba = L.b.const_arrays();
        ParallelFor(L.b, IntVect(0), nc, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
        {
            pd[n*N + Long(j-dlo1)*nx + (i-dlo0)] = ba[bno](i,j,k,n);
        });
    }
    if (!local) {
        Gpu::copy(Gpu::deviceToHost, L.bdev.begin(), L.bdev.end(), L.bg.begin());
        ParallelAllReduce::Sum(L.bg.data(), int(N*nc), ParallelContext::CommunicatorSub());
        Gpu::copy(Gpu::hostToDevice, L.bg.begin(), L.bg.end(), L.bdev.begin());
    }

#ifdef AMREX_USE_FFT
    {
        Real const* const pinv = L.binv.data();
        GpuComplex<Real>* const ps = L.bspec.data();
        Long const NS = Long(L.bspec.size());
        for (int n = 0; n < nc; ++n) {
            L.bfft[n]->forward(pd + n*N, ps);
            ParallelFor(NS, [=] AMREX_GPU_DEVICE (Long m) noexcept { ps[m] *= pinv[m]; });
            L.bfft[n]->backward(ps, pd + n*N);
        }
    }
#else
    Gpu::copy(Gpu::deviceToHost, L.bdev.begin(), L.bdev.end(), L.bg.begin());
    // DFT along direction d of (re, im), in place; inverse divides by n
    auto dft = [&] (int d, bool inverse) {
        int const n = nn[d];
        if (n == 1) { return; }
        Long const stride = (d == 0) ? 1 : nx;
        Long const nlines = N / n;
        Real const* ct = L.bcos[d].data();
        Real const* st = L.bsin[d].data();
        Real const sgn = inverse ? Real(-1.0) : Real(1.0);
        Real const scale = inverse ? Real(1.0)/Real(n) : Real(1.0);
        Vector<Real> tr(n), ti(n);
        for (Long l = 0; l < nlines; ++l) {
            Long const base = (d == 0) ? l*nx : l;
            for (int kk = 0; kk < n; ++kk) {
                Real sr = Real(0.0), si = Real(0.0);
                for (int m = 0; m < n; ++m) {
                    Real const c = ct[Long(kk)*n+m];
                    Real const sn = sgn*st[Long(kk)*n+m];
                    Real const xr = L.bre[base + m*stride];
                    Real const xi = L.bim[base + m*stride];
                    sr += xr*c + xi*sn;
                    si += xi*c - xr*sn;
                }
                tr[kk] = sr*scale;
                ti[kk] = si*scale;
            }
            for (int kk = 0; kk < n; ++kk) {
                L.bre[base + kk*stride] = tr[kk];
                L.bim[base + kk*stride] = ti[kk];
            }
        }
    };
    Real const* const hinv = L.binv.data();
    for (int n = 0; n < nc; ++n) {
        std::copy(L.bg.begin() + n*N, L.bg.begin() + (n+1)*N, L.bre.begin());
        std::fill(L.bim.begin(), L.bim.end(), Real(0.0));
        for (int d = 0; d < 2; ++d) { dft(d, false); }
        for (Long m = 0; m < N; ++m) { L.bre[m] *= hinv[m]; L.bim[m] *= hinv[m]; }
        for (int d = 0; d < 2; ++d) { dft(d, true); }
        std::copy(L.bre.begin(), L.bre.end(), L.bg.begin() + n*N);
    }
    Gpu::copy(Gpu::hostToDevice, L.bg.begin(), L.bg.end(), L.bdev.begin());
#endif

    // scatter the solution into x
    {
        Real const* const ps = L.bdev.data();
        auto const& xa = L.x.arrays();
        ParallelFor(L.x, IntVect(0), nc, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
        {
            xa[bno](i,j,k,n) = ps[n*N + Long(j-dlo1)*nx + (i-dlo0)];
        });
    }
}

void alloc_level (CrossProj::Level& L)
{
    L.x.define(L.ba, L.dm, L.ncomp, 1);
    L.x2.define(L.ba, L.dm, L.ncomp, 1);
    L.r.define(L.ba, L.dm, L.ncomp, 1);
    L.b.define(L.ba, L.dm, L.ncomp, 0);
    L.x.setVal(Real(0.0));
    L.x2.setVal(Real(0.0));
    L.r.setVal(Real(0.0));
    set_stencil(L);
    set_transfer_weights(L);
}

// Coarse level of P by 2; regrids when the coarsened boxes are too small
// to be coarsened again.
std::unique_ptr<CrossProj::Level>
make_coarse (CrossProj::Level const& P, bool modulate)
{
    auto C = std::make_unique<CrossProj::Level>();
    C->domain = amrex::coarsen(P.domain, 2);
    C->period = Periodicity(C->domain.length());
    C->modulate = modulate;
    C->ncomp = P.ncomp;
    C->cheby_degree = P.cheby_degree;

    // The bottom level is one box on one rank, so that its exact solve is
    // local and the transfers reach only that rank. Other levels keep the
    // parent's boxes unless they got too small to coarsen again.
    bool const bottom = !coarsenable(C->domain);
    BoxArray const cba = amrex::coarsen(P.ba, 2);
    bool ok = true;
    for (int i = 0, N = int(cba.size()); i < N && ok; ++i) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            int const len = cba[i].length(idim);
            if (len % 2 != 0 || len < 4) { ok = false; }
        }
    }
    if (bottom && cba.size() > 1) { ok = false; }
    if (ok) {
        C->ba = cba;
        C->dm = P.dm;
    } else {
        C->ba = BoxArray(C->domain);
        if (!bottom) { C->ba.maxSize(regrid_max_size); }
        C->dm = DistributionMapping(C->ba);
        C->cba = cba;
        C->t.define(C->cba, P.dm, C->ncomp, 1);
        C->t.setVal(Real(0.0));
    }

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        C->K[idim] = galerkin1d(P.K[idim], modulate);
        C->M[idim] = galerkin1d(P.M[idim], modulate);
    }
    alloc_level(*C);
    if (!bottom) { estimate_lmax(*C); }
    return C;
}

// Reduced right-hand side of the red and black problems from the fine b
// (needs 1 ghost node).
void reduce_rhs (CrossProj::Level const& T, CrossProj::Level& R)
{
    Real const w = -T.sten[0] / T.sten[4]; // minus corner over diagonal (1/4)
    auto const& fa = T.b.const_arrays();
    auto const& ra = R.b.arrays();
    ParallelFor(R.b, [=] AMREX_GPU_DEVICE (int bno, int I, int J, int k) noexcept
    {
        auto const& b = fa[bno];
        int const i = 2*I;
        int const j = 2*J;
        ra[bno](I,J,k,0) = b(i,j,k) + w*(b(i-1,j-1,k) + b(i+1,j-1,k) + b(i-1,j+1,k) + b(i+1,j+1,k));
        ra[bno](I,J,k,1) = b(i+1,j,k) + w*(b(i,j-1,k) + b(i+2,j-1,k) + b(i,j+1,k) + b(i+2,j+1,k));
    });
}

// Fine phi from the reduced solution xr (needs 1 ghost node) and the fine b.
void back_substitute (CrossProj::Level& T, MultiFab const& xr)
{
    Real const w = -T.sten[0] / T.sten[4];
    Real const dinv = Real(1.0) / T.sten[4];
    auto const& fa = T.b.const_arrays();
    auto const& ca = xr.const_arrays();
    auto const& pa = T.x.arrays();
    ParallelFor(T.x, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
    {
        auto const& x = ca[bno];
        int const ip = i & 1;
        int const jp = j & 1;
        int const I = (i - ip) / 2;
        int const J = (j - jp) / 2;
        Real v;
        if (ip == 0 && jp == 0) {
            v = x(I,J,k,0);
        } else if (ip == 1 && jp == 1) {
            v = dinv*fa[bno](i,j,k) + w*(x(I,J,k,0) + x(I+1,J,k,0) + x(I,J+1,k,0) + x(I+1,J+1,k,0));
        } else if (ip == 1) {
            v = x(I,J,k,1);
        } else {
            v = dinv*fa[bno](i,j,k) + w*(x(I-1,J,k,1) + x(I,J,k,1) + x(I-1,J+1,k,1) + x(I,J+1,k,1));
        }
        pa[bno](i,j,k) = v;
    });
}

// max |u| (global)
Real maxabs (MultiFab const& u, int ncomp)
{
    Real s = Real(0.0);
    for (int n = 0; n < ncomp; ++n) {
        s = std::max(s, u.norminf(n, 0, true));
    }
    ParallelAllReduce::Max(s, ParallelContext::CommunicatorSub());
    return s;
}

// q = A p (p needs filled ghosts); returns p.q (global)
Real apply_dot (CrossProj::Level const& P, MultiFab& q, MultiFab const& p)
{
    BL_PROFILE("GMG::CrossProj::pcg_apply");
    auto const sten = P.sten;
    auto const& qa = q.arrays();
    auto const& pa = p.const_arrays();
    Real s = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, q, IntVect(0), P.ncomp,
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept -> GpuTuple<Real>
        {
            Real const v = stencil_sum(sten, pa[bno], i, j, k, n);
            qa[bno](i,j,k,n) = v;
            return {v * pa[bno](i,j,k,n)};
        });
    ParallelAllReduce::Sum(s, ParallelContext::CommunicatorSub());
    return s;
}

// x += alpha p, r -= alpha q; returns max |r| (global)
Real update_norm (CrossProj::Level const& P, MultiFab& x, MultiFab& r, MultiFab const& p,
                  MultiFab const& q, Real alpha)
{
    BL_PROFILE("GMG::CrossProj::pcg_update");
    auto const& xa = x.arrays();
    auto const& ra = r.arrays();
    auto const& pa = p.const_arrays();
    auto const& qa = q.const_arrays();
    Real s = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, r, IntVect(0), P.ncomp,
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept -> GpuTuple<Real>
        {
            xa[bno](i,j,k,n) += alpha * pa[bno](i,j,k,n);
            Real const v = ra[bno](i,j,k,n) - alpha * qa[bno](i,j,k,n);
            ra[bno](i,j,k,n) = v;
            return {std::abs(v)};
        });
    ParallelAllReduce::Max(s, ParallelContext::CommunicatorSub());
    return s;
}

}

CrossProj::CrossProj (Geometry const& geom)
    : m_geom(geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.isAllPeriodic(), "GMG::CrossProj: only supports periodic BC");
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.Domain().length(idim) % 2 == 0 &&
                                         geom.Domain().length(idim) >= 4 &&
                                         geom.Domain().smallEnd(idim) == 0,
            "GMG::CrossProj: the domain needs an even number (>= 4) of cells starting at index 0");
    }
}

CrossProj::~CrossProj () = default;

void
CrossProj::setup (BoxArray const& ba, DistributionMapping const& dm)
{
    BL_PROFILE("GMG::CrossProj::setup");

    for (int i = 0, N = int(ba.size()); i < N; ++i) {
        IntVect const len = ba[i].length();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(len[idim] % 2 == 0, "GMG::CrossProj: box lengths must be even");
        }
    }
    m_grids = ba;
    m_dmap = dm;
    m_vel.define(m_grids, m_dmap, AMREX_SPACEDIM, 1);

    m_top = std::make_unique<Level>();
    Level& T = *m_top;
    T.domain = m_geom.Domain();
    T.period = m_geom.periodicity();
    T.ba = ba;
    T.dm = dm;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Real const hinv2 = Real(1.0) / (m_geom.CellSize(idim)*m_geom.CellSize(idim));
        T.K[idim] = {-hinv2, Real(2.0)*hinv2, -hinv2};
        T.M[idim] = {Real(0.25), Real(0.5), Real(0.25)};
    }

    for (auto& h : m_hier) {
        h.clear();
    }

    {
        Real const rdx = m_geom.CellSize(0) / m_geom.CellSize(1);
        constexpr Real tol = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-12);
        m_reduced = std::abs(Real(1.0) - rdx*rdx) <= tol;
    }

    Level* pcglev = &T;
    if (m_reduced) {
        set_stencil(T);
        T.b.define(ba, dm, 1, 1);
        T.x.define(ba, dm, 1, 1);
        auto R = std::make_unique<Level>();
        R->domain = amrex::coarsen(T.domain, 2);
        R->period = Periodicity(R->domain.length());
        R->ncomp = 2;
        R->cheby_degree = 2;
        R->ba = amrex::coarsen(ba, 2);
        R->dm = dm;
        Real const a0 = T.sten[4];
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            R->K[idim] = {-a0/Real(4.0), a0/Real(2.0), -a0/Real(4.0)};
            R->M[idim] = {Real(0.125), Real(0.75), Real(0.125)};
        }
        alloc_level(*R);
        estimate_lmax(*R);
        m_hier[0].push_back(std::move(R));
        pcglev = m_hier[0][0].get();
    } else {
        alloc_level(T);
        estimate_lmax(T);
    }

    // Hierarchy 0: smooth modes; 1: the checkerboard.
    for (int h = 0; h < (m_reduced ? 1 : nhier); ++h) {
        // Coarsen while possible, until the level is small enough for the
        // exact bottom solve. The fine grid always gets a coarse level.
        Level const* parent = m_hier[h].empty() ? &T : m_hier[h].back().get();
        bool more = m_hier[h].empty() || coarsenable(parent->domain);
        while (more) {
            m_hier[h].push_back(make_coarse(*parent, h == 1 && m_hier[h].empty()));
            parent = m_hier[h].back().get();
            more = coarsenable(parent->domain);
        }
    }

    m_x.define(pcglev->ba, pcglev->dm, pcglev->ncomp, 1);
    m_p.define(pcglev->ba, pcglev->dm, pcglev->ncomp, 1);
    m_q.define(pcglev->ba, pcglev->dm, pcglev->ncomp, 0);
    if (!m_reduced) {
        m_rold.define(pcglev->ba, pcglev->dm, pcglev->ncomp, 0);
    }
}

void
CrossProj::vcycle (int h, int lev)
{
    auto& hier = m_hier[h];
    Level& L = *hier[lev];
    if (lev+1 == int(hier.size())) {
        bottom_solve(L);
        return;
    }
    Level& C = *hier[lev+1];
    cheby(L, 1, true);
    residual(L);
    restrict_to(L, C);
    vcycle(h, lev+1);
    prolong_add(L, C);
    cheby(L, 1, false);
}

// Preconditioner: x = B^-1 b on the PCG level. Reduced: one V-cycle. Fine
// grid: smoothing, a correction from each hierarchy in turn, smoothing.
// Not symmetric, hence the flexible CG in solve.
void
CrossProj::precond ()
{
    BL_PROFILE("GMG::CrossProj::precond");
    if (m_reduced) {
        vcycle(0, 0);
        return;
    }
    Level& T = *m_top;
    cheby(T, 1, true);
    for (int h = 0; h < nhier; ++h) {
        Level& C = *m_hier[h][0];
        residual(T);
        restrict_to(T, C);
        vcycle(h, 0);
        prolong_add(T, C);
    }
    cheby(T, 1, false);
}

void
CrossProj::solve (MultiFab& a_soln, MultiFab const& a_rhs)
{
    BL_PROFILE("GMG::CrossProj::solve");

    AMREX_ALWAYS_ASSERT(a_soln.is_cell_centered() && a_rhs.is_cell_centered() &&
                        a_soln.nComp() >= AMREX_SPACEDIM && a_rhs.nComp() >= AMREX_SPACEDIM);

    if (!m_top || m_grids != a_rhs.boxArray() || m_dmap != a_rhs.DistributionMap()) {
        setup(a_rhs.boxArray(), a_rhs.DistributionMap());
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a_soln.boxArray() == m_grids && a_soln.DistributionMap() == m_dmap,
        "GMG::CrossProj: a_soln and a_rhs must have the same BoxArray and DistributionMapping");
    Level& T = *m_top;
    Level& P = m_reduced ? *m_hier[0][0] : T; // PCG level: P.b is r, P.x is z
    int const nc = P.ncomp;

    MultiFab::Copy(m_vel, a_rhs, 0, 0, AMREX_SPACEDIM, 0);
    m_vel.FillBoundary(m_geom.periodicity());

    auto const dxinv = m_geom.InvCellSizeArray();
    // G averages the differences over the 2 edges of a cell.
    GpuArray<Real,2> const fac{Real(0.5)*dxinv[0], Real(0.5)*dxinv[1]};

    // b = G^T u = -D u at the node (i,j), from the cells (i-1..i, j-1..j)
    {
        BL_PROFILE("GMG::CrossProj::rhs");
        auto const& ba = T.b.arrays();
        auto const& va = m_vel.const_arrays();
        ParallelFor(T.b, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
        {
            auto const& v = va[bno];
            Real s = Real(0.0);
            for (int bb = 0; bb < 2; ++bb) {
            for (int a = 0; a < 2; ++a) {
                int const ii = i-1+a, jj = j-1+bb;
                s += (a  ? -fac[0] : fac[0]) * v(ii,jj,k,0)
                   + (bb ? -fac[1] : fac[1]) * v(ii,jj,k,1);
            }}
            ba[bno](i,j,k) = s;
        });
    }

    Real const bnorm = maxabs(T.b, 1);
    // Below the roundoff in G^T u the residual cannot decrease, e.g., when
    // projecting a field that is already divergence free.
    Real const floor = Real(100.0) * std::numeric_limits<Real>::epsilon()
        * maxabs(m_vel, AMREX_SPACEDIM) * (dxinv[0] + dxinv[1]);
    Real const target = std::max({m_reltol*bnorm, m_abstol, floor});

    if (m_reduced) {
        BL_PROFILE("GMG::CrossProj::reduce");
        T.b.FillBoundary(T.period);
        reduce_rhs(T, P);
    }

    m_x.setVal(Real(0.0));
    m_niters = 0;
    Real rnorm = bnorm;
    if (rnorm > target) {
        // PCG; with the non-symmetric fine-grid preconditioner the flexible
        // variant (beta from z.(r_new - r_old)) is used.
        bool const flexible = !m_reduced;
        precond();
        MultiFab::Copy(m_p, P.x, 0, 0, nc, 0);
        Real rz = MultiFab::Dot(P.b, 0, P.x, 0, nc, 0);
        for (int iter = 1; iter <= m_maxiter; ++iter) {
            m_p.FillBoundary(P.period);
            Real const pq = apply_dot(P, m_q, m_p);
            Real const alpha = rz / pq;
            if (flexible) { MultiFab::Copy(m_rold, P.b, 0, 0, nc, 0); }
            rnorm = update_norm(P, m_x, P.b, m_p, m_q, alpha);
            m_niters = iter;
            if (m_verbose > 1) {
                amrex::Print() << "GMG::CrossProj: iteration " << iter << " residual "
                               << rnorm << " (" << rnorm/bnorm << ")\n";
            }
            if (rnorm <= target) { break; }
            precond();
            BL_PROFILE_VAR("GMG::CrossProj::pcg_p", pcgp);
            Real const rz1 = MultiFab::Dot(P.b, 0, P.x, 0, nc, 0);
            Real const beta = (flexible ? rz1 - MultiFab::Dot(m_rold, 0, P.x, 0, nc, 0) : rz1) / rz;
            MultiFab::Xpay(m_p, beta, P.x, 0, 0, nc, 0);
            BL_PROFILE_VAR_STOP(pcgp);
            rz = rz1;
        }
    }
    if (m_verbose > 0) {
        amrex::Print() << "GMG::CrossProj: " << m_niters << " iterations, residual "
                       << rnorm << " (" << ((bnorm > Real(0.0)) ? rnorm/bnorm : Real(0.0)) << ")\n";
    }
    if (rnorm > target) {
        amrex::Warning("GMG::CrossProj: not converged");
    }

    MultiFab& phi = m_reduced ? T.x : m_x;
    if (m_reduced) {
        BL_PROFILE("GMG::CrossProj::backsub");
        m_x.FillBoundary(P.period);
        back_substitute(T, m_x);
    }
    phi.FillBoundary(T.period);

    // u -= G phi
    {
        BL_PROFILE("GMG::CrossProj::grad");
        auto const& sa = a_soln.arrays();
        auto const& va = a_rhs.const_arrays();
        auto const& pa = phi.const_arrays();
        ParallelFor(a_soln, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
        {
            auto const& ph = pa[bno];
            Real const gx = ph(i+1,j,k) + ph(i+1,j+1,k) - ph(i,j,k) - ph(i,j+1,k);
            Real const gy = ph(i,j+1,k) + ph(i+1,j+1,k) - ph(i,j,k) - ph(i+1,j,k);
            sa[bno](i,j,k,0) = va[bno](i,j,k,0) - fac[0] * gx;
            sa[bno](i,j,k,1) = va[bno](i,j,k,1) - fac[1] * gy;
        });
    }
    Gpu::streamSynchronize();

    IntVect const ng = amrex::elemwiseMin(a_soln.nGrowVect(), IntVect(1));
    if (ng.max() > 0) {
        a_soln.FillBoundary(0, AMREX_SPACEDIM, ng, m_geom.periodicity());
    }
}

}

#endif
