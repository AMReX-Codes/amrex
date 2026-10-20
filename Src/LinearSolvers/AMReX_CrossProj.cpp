#include <AMReX_CrossProj.H>

#if (AMREX_SPACEDIM >= 2)

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
// Every operator in the hierarchy is sum_d K_d (x) M_e (x) M_f with 1D
// three-point stencils K and M per direction: the fine one has K = second
// difference and M = [1/4,1/2,1/4], and Galerkin coarsening with
// tensor-product interpolation keeps the form. On a periodic grid all of
// them are constant, so each level is a constant 9- or 27-point stencil.
//
// Node data are stored on the cell-centered BoxArray: box [lo,hi) holds
// the nodes it owns (a periodic grid has as many nodes as cells), so no
// node is duplicated and ghost nodes come from FillBoundary.
//
// 2D with dx == dy: the fine stencil couples only diagonal neighbors, so the
// nodes with i+j even (red) and odd (black) decouple. Eliminating the
// (odd,odd) nodes of red and the (even,odd) nodes of black is exact and
// leaves the same 9-point problem for both colors on a grid coarsened by 2
// (black is red shifted by one in x). PCG runs on that reduced problem,
// stored as two components of one MultiFab, preconditioned by a V-cycle.
//
// Otherwise PCG runs on the fine grid, preconditioned by a multigrid cycle
// with a coarse space for the smooth modes (linear interpolation) and
// coarse spaces for the near-null modes: in 2D the checkerboard (-1)^(i+j)
// times smooth modes, in 3D the families (-1)^(j+k) f(i), (-1)^(i+k) f(j),
// (-1)^(i+j) f(k), coarsened only transverse to the free direction and
// smoothed with line solves along it.
struct CrossProj::Level
{
    Box domain;
    Periodicity period;
    IntVect ratio{1};      // coarsening ratio from the parent level
    int line_dir = -1;     // family levels: line smoother along this direction
    IntVect modulate{0};   // first near-null level: P has (-1)^(sum of these indices)
    int ncomp = 1;
    int cheby_degree = 3;  // 2 on the reduced 2D hierarchy
    BoxArray cba;
    DistributionMapping dm;
    BoxArray aligned;      // parent BoxArray coarsened, if cba was regridded
    Array<Array<Real,3>,AMREX_SPACEDIM> K, M;  // [dir][offset+1]
    GpuArray<Real,27> sten{};  // [(c+1)*9+(b+1)*3+(a+1)] for offset (a,b,c); c = 0 in 2D
    GpuArray<Real,27> rw{};    // restriction weights from the parent, same layout
    GpuArray<Real,8> pw{};     // prolongation weights [c*4+b*2+a] for fine offset (a,b,c)
    Gpu::DeviceVector<Real> tcp, tdi; // Thomas factors of the line solver by position
    Gpu::DeviceVector<int> lines;     // GPU line solver: (box, i, j, k) of every line's first node
    MultiFab x, x2, b, r, y, talign;
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
    Array<Vector<Real>,3> bcos, bsin;           // naive DFT tables
    Vector<Real> bre, bim, binv;
#endif
    Real lmax = Real(0.0); // Chebyshev upper bound
};

namespace {

constexpr int regrid_max_size = 32;
constexpr int bottom_max_nodes = 1024; // exact bottom solve at or below this size
constexpr int kz = (AMREX_SPACEDIM == 3) ? 1 : 0;
constexpr int nhier = (AMREX_SPACEDIM == 3) ? 4 : 2;


// Whether a level with this domain gets a coarser level (by `rr`), or is
// the bottom: small enough for the exact solve, or not coarsenable.
bool coarsenable (Box const& domain, IntVect const& rr)
{
    if (domain.numPts() <= bottom_max_nodes) { return false; }
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        int const len = domain.length(e);
        if (rr[e] == 2 && (len % 2 != 0 || len < 8)) { return false; }
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
Real stencil_sum (GpuArray<Real,27> const& s, Array4<Real const> const& x,
                  int i, int j, int k, int n) noexcept
{
#if (AMREX_SPACEDIM == 2)
    return s[ 9]*x(i-1,j-1,k,n) + s[10]*x(i,j-1,k,n) + s[11]*x(i+1,j-1,k,n)
         + s[12]*x(i-1,j  ,k,n) + s[13]*x(i,j  ,k,n) + s[14]*x(i+1,j  ,k,n)
         + s[15]*x(i-1,j+1,k,n) + s[16]*x(i,j+1,k,n) + s[17]*x(i+1,j+1,k,n);
#else
    return s[ 0]*x(i-1,j-1,k-1,n)
         + s[ 1]*x(i,j-1,k-1,n)
         + s[ 2]*x(i+1,j-1,k-1,n)
         + s[ 3]*x(i-1,j,k-1,n)
         + s[ 4]*x(i,j,k-1,n)
         + s[ 5]*x(i+1,j,k-1,n)
         + s[ 6]*x(i-1,j+1,k-1,n)
         + s[ 7]*x(i,j+1,k-1,n)
         + s[ 8]*x(i+1,j+1,k-1,n)
         + s[ 9]*x(i-1,j-1,k,n)
         + s[10]*x(i,j-1,k,n)
         + s[11]*x(i+1,j-1,k,n)
         + s[12]*x(i-1,j,k,n)
         + s[13]*x(i,j,k,n)
         + s[14]*x(i+1,j,k,n)
         + s[15]*x(i-1,j+1,k,n)
         + s[16]*x(i,j+1,k,n)
         + s[17]*x(i+1,j+1,k,n)
         + s[18]*x(i-1,j-1,k+1,n)
         + s[19]*x(i,j-1,k+1,n)
         + s[20]*x(i+1,j-1,k+1,n)
         + s[21]*x(i-1,j,k+1,n)
         + s[22]*x(i,j,k+1,n)
         + s[23]*x(i+1,j,k+1,n)
         + s[24]*x(i-1,j+1,k+1,n)
         + s[25]*x(i,j+1,k+1,n)
         + s[26]*x(i+1,j+1,k+1,n);
#endif
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

// y = B^-1 r, with B the diagonal of A or, on family levels, the part of A
// coupling the nodes of a box along the line direction (Thomas algorithm
// with the factors precomputed by position, since the tridiagonal is
// constant).
void apply_prec (CrossProj::Level& L, MultiFab& y, MultiFab const& r)
{
    BL_PROFILE("GMG::CrossProj::prec");
    if (L.line_dir < 0) {
        Real const dinv = Real(1.0) / L.sten[13];
        auto const& ya = y.arrays();
        auto const& ra = r.const_arrays();
        ParallelFor(y, IntVect(0), L.ncomp, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
        {
            ya[bno](i,j,k,n) = ra[bno](i,j,k,n) * dinv;
        });
        return;
    }

    AMREX_ASSERT(L.ncomp == 1);
    int const f = L.line_dir;
    int const stride = (f == 0) ? 1 : ((f == 1) ? 3 : 9);
    Real const tlo = L.sten[13-stride]; // line part of the stencil: (tlo, tdiag, tlo)
    Real const* AMREX_RESTRICT cp = L.tcp.data();
    Real const* AMREX_RESTRICT di = L.tdi.data();
#ifdef AMREX_USE_GPU
    // one launch for all lines of all boxes
    int const* AMREX_RESTRICT ln = L.lines.data();
    int const nlines = int(L.lines.size() / 4);
    auto const& ya4 = y.arrays();
    auto const& ra4 = r.const_arrays();
    ParallelFor(nlines, [=] AMREX_GPU_DEVICE (int l) noexcept
    {
        int const bno = ln[4*l];
        IntVect iv(AMREX_D_DECL(ln[4*l+1], ln[4*l+2], ln[4*l+3]));
        auto const& ya = ya4[bno];
        auto const& ra = ra4[bno];
        int const mlo = iv[f];
        int const mhi = mlo + ya.end[f] - ya.begin[f] - 1; // y has no ghost nodes
        Real dprev = Real(0.0);
        for (int m = mlo; m <= mhi; ++m) {
            iv[f] = m;
            dprev = (ra(iv) - tlo*dprev) * di[m-mlo];
            ya(iv) = dprev;
        }
        for (int m = mhi-1; m >= mlo; --m) {
            IntVect jv = iv;
            iv[f] = m;
            jv[f] = m+1;
            ya(iv) -= cp[m-mlo] * ya(jv);
        }
    });
#else
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(y); mfi.isValid(); ++mfi) {
        Box const& vbx = mfi.validbox();
        int const mlo = vbx.smallEnd(f);
        int const mhi = vbx.bigEnd(f);
        auto const& ya = y.array(mfi);
        auto const& ra = r.const_array(mfi);
        // Line index outer, transverse plane inner (vectorizes for y and z lines).
        auto const lo = amrex::lbound(vbx);
        auto const hi = amrex::ubound(vbx);
        Dim3 const off{.x = (f == 0), .y = (f == 1), .z = (f == 2)};
        auto plane = [&] (int m) {
            Dim3 plo = lo, phi = hi;
            if (f == 0) { plo.x = phi.x = m; } else if (f == 1) { plo.y = phi.y = m; } else { plo.z = phi.z = m; }
            return std::make_pair(plo, phi);
        };
        for (int m = mlo; m <= mhi; ++m) {
            auto const [plo, phi] = plane(m);
            Real const d = di[m-mlo];
            if (m == mlo) {
                LoopConcurrentOnCpu(plo, phi, [&] (int i, int j, int k) noexcept
                {
                    ya(i,j,k) = ra(i,j,k) * d;
                });
            } else {
                LoopConcurrentOnCpu(plo, phi, [&] (int i, int j, int k) noexcept
                {
                    ya(i,j,k) = (ra(i,j,k) - tlo*ya(i-off.x,j-off.y,k-off.z)) * d;
                });
            }
        }
        for (int m = mhi-1; m >= mlo; --m) {
            auto const [plo, phi] = plane(m);
            Real const c = cp[m-mlo];
            LoopConcurrentOnCpu(plo, phi, [&] (int i, int j, int k) noexcept
            {
                ya(i,j,k) -= c * ya(i+off.x,j+off.y,k+off.z);
            });
        }
    }
#endif
}

// One step of the three-term Chebyshev recurrence
//   x_new = x + alpha (x - x_prev) + beta B^-1 (b - A x),
// with x in L.x and x_prev in L.x2; x_new overwrites L.x2 and the two are
// swapped. ZERO_X: x is zero (not read). FIRST: no alpha term. ZERO_PREV:
// x_prev is zero. Point levels do it in one kernel.
template <bool ZERO_X, bool FIRST, bool ZERO_PREV>
void cheby_step (CrossProj::Level& L, Real alpha, Real beta)
{
    BL_PROFILE("GMG::CrossProj::cheby_step");
    if (L.line_dir < 0) {
        if (!ZERO_X) { L.x.FillBoundary(L.period); }
        auto const sten = L.sten;
        Real const bd = beta / L.sten[13];
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
    } else {
        if constexpr (ZERO_X) {
            MultiFab::Copy(L.r, L.b, 0, 0, L.ncomp, 0);
        } else {
            residual(L);
        }
        apply_prec(L, L.y, L.r);
        auto const& ya = L.y.const_arrays();
        auto const& xa = L.x.const_arrays();
        auto const& oa = L.x2.arrays();
        ParallelFor(L.x2, IntVect(0), L.ncomp, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
        {
            amrex::ignore_unused(alpha, beta, xa, ya, oa);
            Real v = beta * ya[bno](i,j,k,n);
            if constexpr (!ZERO_X) {
                Real const xc = xa[bno](i,j,k,n);
                v += xc;
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
    }
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

// Upper end of the spectrum of B^-1 A by power iteration, with a margin.
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
    Real lam = Real(0.0);
    for (int it = 0; it < 30; ++it) {
        Real const nrm = std::sqrt(MultiFab::Dot(L.x, 0, L.x, 0, L.ncomp, 0));
        L.x.mult(Real(1.0)/nrm, 0, L.ncomp, 0);
        residual(L); // r = -A x
        apply_prec(L, L.y, L.r);
        lam = std::sqrt(MultiFab::Dot(L.y, 0, L.y, 0, L.ncomp, 0));
        MultiFab::Copy(L.x, L.y, 0, 0, L.ncomp, 0);
    }
    L.lmax = Real(1.1) * lam;
}

// C.b = P^T F.r (weights C.rw); F.r needs filled ghosts
template <int RX, int RY, int RZ>
void restrict_impl (CrossProj::Level& F, CrossProj::Level& C)
{
    MultiFab* dst = C.aligned.empty() ? &C.b : &C.talign;
    auto const w = C.rw;
    auto const& fa = F.r.const_arrays();
    auto const& ca = dst->arrays();
    ParallelFor(*dst, IntVect(0), C.ncomp, [=] AMREX_GPU_DEVICE (int bno, int I, int J, int K, int n) noexcept
    {
        auto const& f = fa[bno];
        Real s = Real(0.0);
        for (int c = -(RZ-1); c <= RZ-1; ++c) {
        for (int b = -(RY-1); b <= RY-1; ++b) {
        for (int a = -(RX-1); a <= RX-1; ++a) {
            s += w[(c+1)*9+(b+1)*3+(a+1)] * f(RX*I+a, RY*J+b, RZ*K+c, n);
        }}}
        ca[bno](I,J,K,n) = s;
    });
    if (dst != &C.b) {
        C.b.ParallelCopy(C.talign, 0, 0, C.ncomp);
    }
}

void restrict_to (CrossProj::Level& F, CrossProj::Level& C)
{
    BL_PROFILE("GMG::CrossProj::restrict");
    F.r.FillBoundary(F.period);
#if (AMREX_SPACEDIM == 2)
    restrict_impl<2,2,1>(F, C);
#else
    if (C.ratio[0] == 1)      { restrict_impl<1,2,2>(F, C); }
    else if (C.ratio[1] == 1) { restrict_impl<2,1,2>(F, C); }
    else if (C.ratio[2] == 1) { restrict_impl<2,2,1>(F, C); }
    else                      { restrict_impl<2,2,2>(F, C); }
#endif
}

// F.x += P C.x (weights C.pw). Each coarse node writes the fine nodes
// (RX I + a, ...), a < RX, from the coarse nodes I and I+1.
template <int RX, int RY, int RZ>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void prolong_node (Array4<Real const> const& cx, Array4<Real> const& fx,
                   GpuArray<Real,8> const& w, int I, int J, int K, int n) noexcept
{
    Real const c000 = cx(I,J,K,n);
    Real c100 = Real(0.0), c010 = Real(0.0), c110 = Real(0.0);
    Real c001 = Real(0.0), c101 = Real(0.0), c011 = Real(0.0), c111 = Real(0.0);
    if constexpr (RX == 2) { c100 = cx(I+1,J,K,n); }
    if constexpr (RY == 2) { c010 = cx(I,J+1,K,n); }
    if constexpr (RX == 2 && RY == 2) { c110 = cx(I+1,J+1,K,n); }
    if constexpr (RZ == 2) {
        c001 = cx(I,J,K+1,n);
        if constexpr (RX == 2) { c101 = cx(I+1,J,K+1,n); }
        if constexpr (RY == 2) { c011 = cx(I,J+1,K+1,n); }
        if constexpr (RX == 2 && RY == 2) { c111 = cx(I+1,J+1,K+1,n); }
    }
    int const i = RX*I;
    int const j = RY*J;
    int const k = RZ*K;
    fx(i,j,k,n) += w[0]*c000;
    if constexpr (RX == 2) { fx(i+1,j,k,n) += w[1]*(c000+c100); }
    if constexpr (RY == 2) { fx(i,j+1,k,n) += w[2]*(c000+c010); }
    if constexpr (RX == 2 && RY == 2) { fx(i+1,j+1,k,n) += w[3]*(c000+c100+c010+c110); }
    if constexpr (RZ == 2) {
        fx(i,j,k+1,n) += w[4]*(c000+c001);
        if constexpr (RX == 2) { fx(i+1,j,k+1,n) += w[5]*(c000+c100+c001+c101); }
        if constexpr (RY == 2) { fx(i,j+1,k+1,n) += w[6]*(c000+c010+c001+c011); }
        if constexpr (RX == 2 && RY == 2) {
            fx(i+1,j+1,k+1,n) += w[7]*(c000+c100+c010+c110+c001+c101+c011+c111);
        }
    }
}

template <int RX, int RY, int RZ>
void prolong_impl (CrossProj::Level& F, CrossProj::Level& C)
{
    MultiFab const* src = &C.x;
    if (C.aligned.empty()) {
        C.x.FillBoundary(C.period);
    } else {
        C.talign.ParallelCopy(C.x, 0, 0, C.ncomp, IntVect(0), IntVect(1), C.period);
        src = &C.talign;
    }
    auto const w = C.pw;
    auto const& fa = F.x.arrays();
    auto const& ca = src->const_arrays();
    ParallelFor(*src, IntVect(0), C.ncomp, [=] AMREX_GPU_DEVICE (int bno, int I, int J, int K, int n) noexcept
    {
        prolong_node<RX,RY,RZ>(ca[bno], fa[bno], w, I, J, K, n);
    });
}

void prolong_add (CrossProj::Level& F, CrossProj::Level& C)
{
    BL_PROFILE("GMG::CrossProj::prolong");
#if (AMREX_SPACEDIM == 2)
    prolong_impl<2,2,1>(F, C);
#else
    if (C.ratio[0] == 1)      { prolong_impl<1,2,2>(F, C); }
    else if (C.ratio[1] == 1) { prolong_impl<2,1,2>(F, C); }
    else if (C.ratio[2] == 1) { prolong_impl<2,2,1>(F, C); }
    else                      { prolong_impl<2,2,2>(F, C); }
#endif
}

void set_stencil (CrossProj::Level& L)
{
    for (int b = 0; b < 3; ++b) {
    for (int a = 0; a < 3; ++a) {
#if (AMREX_SPACEDIM == 2)
        L.sten[9+b*3+a] = L.K[0][a]*L.M[1][b] + L.M[0][a]*L.K[1][b];
#else
        for (int c = 0; c < 3; ++c) {
            L.sten[c*9+b*3+a] = L.K[0][a]*L.M[1][b]*L.M[2][c] + L.M[0][a]*L.K[1][b]*L.M[2][c]
                + L.M[0][a]*L.M[1][b]*L.K[2][c];
        }
#endif
    }}
}

// Transfer weights of a level from its parent: linear interpolation in the
// coarsened directions, times (-1)^(modulate . fine index).
void set_transfer_weights (CrossProj::Level& C)
{
    int rr[3] = {1,1,1};
    int md[3] = {0,0,0};
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        rr[e] = C.ratio[e];
        md[e] = C.modulate[e];
    }
    for (int c = -1; c <= 1; ++c) {
    for (int b = -1; b <= 1; ++b) {
    for (int a = -1; a <= 1; ++a) {
        int const off[3] = {a, b, c};
        Real w = Real(1.0);
        for (int e = 0; e < 3; ++e) {
            if (off[e] == 0) { continue; }
            w *= (rr[e] == 2) ? Real(0.5) : Real(0.0);
        }
        if ((md[0]*a + md[1]*b + md[2]*c) & 1) { w = -w; }
        C.rw[(c+1)*9+(b+1)*3+(a+1)] = w;
    }}}
    for (int c = 0; c < 2; ++c) {
    for (int b = 0; b < 2; ++b) {
    for (int a = 0; a < 2; ++a) {
        Real w = Real(1.0) / Real((a+1)*(b+1)*(c+1));
        if ((md[0]*a + md[1]*b + md[2]*c) & 1) { w = -w; }
        C.pw[c*4+b*2+a] = w;
    }}}
}

// Thomas factors of the constant tridiagonal (tlo, tdiag, tlo) by position
// in a segment: tdi[m] = 1/den_m, tcp[m] = tlo/den_m.
void set_thomas (CrossProj::Level& L)
{
    if (L.line_dir < 0) { return; }
    int const f = L.line_dir;
    int const stride = (f == 0) ? 1 : ((f == 1) ? 3 : 9);
    Real const tlo = L.sten[13-stride];
    Real const tdiag = L.sten[13];
    int maxlen = 0;
    for (int i = 0, N = int(L.cba.size()); i < N; ++i) {
        maxlen = std::max(maxlen, L.cba[i].length(f));
    }
    Vector<Real> cp(maxlen), di(maxlen);
    Real den = tdiag;
    for (int m = 0; m < maxlen; ++m) {
        di[m] = Real(1.0) / den;
        cp[m] = tlo * di[m];
        den = tdiag - tlo * cp[m];
    }
    L.tcp.resize(maxlen);
    L.tdi.resize(maxlen);
    Gpu::copy(Gpu::hostToDevice, cp.begin(), cp.end(), L.tcp.begin());
    Gpu::copy(Gpu::hostToDevice, di.begin(), di.end(), L.tdi.begin());
#ifdef AMREX_USE_GPU
    Vector<int> lines;
    for (MFIter mfi(L.y); mfi.isValid(); ++mfi) {
        Box lbx = mfi.validbox();
        lbx.setBig(f, lbx.smallEnd(f));
        int const bno = mfi.LocalIndex();
        amrex::LoopOnCpu(lbx, [&] (int i, int j, int k)
        {
            lines.push_back(bno); lines.push_back(i); lines.push_back(j); lines.push_back(k);
        });
    }
    L.lines.resize(lines.size());
    Gpu::copy(Gpu::hostToDevice, lines.begin(), lines.end(), L.lines.begin());
#endif
    Gpu::streamSynchronize();
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
    int const nn[3] = {L.domain.length(0), L.domain.length(1),
                       (AMREX_SPACEDIM == 3) ? L.domain.length(AMREX_SPACEDIM-1) : 1};
    int const nx = nn[0], ny = nn[1], nz = nn[2];
    Long const N = Long(nx)*ny*nz;
    int const nc = L.ncomp;
    bool const local = (L.cba.size() == 1);
    if (local && L.b.local_size() == 0) { return; }

    // The stencil symbol at mode (i,j,k), in double so that the null modes
    // cancel to roundoff in single precision too.
    auto symbol = [&] (int i, int j, int k) {
        double sym = 0.0;
        for (int c = -kz; c <= kz; ++c) {
        for (int b = -1; b <= 1; ++b) {
        for (int a = -1; a <= 1; ++a) {
            double const arg = 2.0*amrex::Math::pi<double>()
                * (double(a*i)/nx + double(b*j)/ny + double(c*k)/nz);
            sym += double(L.sten[(c+1)*9+(b+1)*3+(a+1)]) * std::cos(arg);
        }}}
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
        IntVect const ss = L.bfft[0]->spectralSize();
        int const sx = ss[0];
        Long const NS = Long(sx)*ny*nz;
        hinv.resize(NS);
        for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < sx; ++i) {
            hinv[(Long(k)*ny + j)*sx + i] = Real(symbol(i,j,k));
        }}}
        invert(hinv, L.bfft[0]->scalingFactor());
        L.bspec.resize(NS);
#else
        for (int d = 0; d < 3; ++d) {
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
        for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            hinv[(Long(k)*ny + j)*nx + i] = Real(symbol(i,j,k));
        }}}
        invert(hinv, Real(1.0)); // the inverse DFT below normalizes
        L.bre.resize(N);
        L.bim.resize(N);
        L.binv = hinv;
#endif
#ifdef AMREX_USE_FFT
        L.binv.resize(hinv.size());
        Gpu::copy(Gpu::hostToDevice, hinv.begin(), hinv.end(), L.binv.begin());
        Gpu::streamSynchronize();
#endif
    }

    // gather the right-hand side into the global array on the device
    int const dlo0 = dlo[0], dlo1 = dlo[1], dlo2 = (AMREX_SPACEDIM == 3) ? dlo[AMREX_SPACEDIM-1] : 0;
    Real* const pd = L.bdev.data();
    {
        Gpu::fillAsync(L.bdev.begin(), L.bdev.end(), [=] AMREX_GPU_DEVICE (Real& v, Long) noexcept { v = Real(0.0); });
        auto const& ba = L.b.const_arrays();
        ParallelFor(L.b, IntVect(0), nc, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
        {
            pd[n*N + (Long(k-dlo2)*ny + (j-dlo1))*nx + (i-dlo0)] = ba[bno](i,j,k,n);
        });
    }
    if (!local) {
        Gpu::copy(Gpu::deviceToHost, L.bdev.begin(), L.bdev.end(), L.bg.begin());
        ParallelDescriptor::ReduceRealSum(L.bg.data(), int(N*nc));
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
        Long const stride = (d == 0) ? 1 : ((d == 1) ? nx : Long(nx)*ny);
        Long const nlines = N / n;
        Real const* ct = L.bcos[d].data();
        Real const* st = L.bsin[d].data();
        Real const sgn = inverse ? Real(-1.0) : Real(1.0);
        Real const scale = inverse ? Real(1.0)/Real(n) : Real(1.0);
        Vector<Real> tr(n), ti(n);
        for (Long l = 0; l < nlines; ++l) {
            Long base;
            if (d == 0)      { base = l*nx; }
            else if (d == 1) { base = (l / nx)*(Long(nx)*ny) + (l % nx); }
            else             { base = l; }
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
        for (int d = 0; d < 3; ++d) { dft(d, false); }
        for (Long m = 0; m < N; ++m) { L.bre[m] *= hinv[m]; L.bim[m] *= hinv[m]; }
        for (int d = 0; d < 3; ++d) { dft(d, true); }
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
            xa[bno](i,j,k,n) = ps[n*N + (Long(k-dlo2)*ny + (j-dlo1))*nx + (i-dlo0)];
        });
    }
}

void alloc_level (CrossProj::Level& L)
{
    L.x.define(L.cba, L.dm, L.ncomp, 1);
    L.x2.define(L.cba, L.dm, L.ncomp, 1);
    L.r.define(L.cba, L.dm, L.ncomp, 1);
    L.b.define(L.cba, L.dm, L.ncomp, 0);
    L.y.define(L.cba, L.dm, L.ncomp, 0);
    L.x.setVal(Real(0.0));
    L.x2.setVal(Real(0.0));
    L.r.setVal(Real(0.0));
    set_stencil(L);
    set_transfer_weights(L);
    set_thomas(L);
}

// Coarse level of P by `ratio`; regrids when the aligned boxes are too small
// to be coarsened again.
std::unique_ptr<CrossProj::Level>
make_coarse (CrossProj::Level const& P, IntVect const& ratio, int line_dir, IntVect const& modulate)
{
    auto C = std::make_unique<CrossProj::Level>();
    C->domain = amrex::coarsen(P.domain, ratio);
    C->period = Periodicity(C->domain.length());
    C->ratio = ratio;
    C->line_dir = line_dir;
    C->modulate = modulate;
    C->ncomp = P.ncomp;
    C->cheby_degree = P.cheby_degree;

    // The bottom level is one box on one rank, so that its exact solve is
    // local and the transfers reach only that rank. Other levels keep the
    // parent's boxes unless they got too small to coarsen again; regridded
    // family levels get boxes spanning the free direction, since the box
    // local line solves are a poor smoother once the transverse coupling
    // is weak (256^3 took 10 iterations instead of 7 with 32-cell segments).
    bool const bottom = !coarsenable(C->domain, ratio);
    BoxArray const aligned = amrex::coarsen(P.cba, ratio);
    bool ok = true;
    for (int i = 0, N = int(aligned.size()); i < N && ok; ++i) {
        for (int e = 0; e < AMREX_SPACEDIM; ++e) {
            int const len = aligned[i].length(e);
            if (ratio[e] == 2 && (len % 2 != 0 || len < 4)) { ok = false; }
        }
    }
    if (bottom && aligned.size() > 1) { ok = false; }
    if (ok) {
        C->cba = aligned;
        C->dm = P.dm;
    } else {
        C->cba = BoxArray(C->domain);
        if (!bottom) {
            IntVect ms(regrid_max_size);
            if (line_dir >= 0) { ms[line_dir] = C->domain.length(line_dir); }
            C->cba.maxSize(ms);
        }
        C->dm = DistributionMapping(C->cba);
        C->aligned = aligned;
        C->talign.define(C->aligned, P.dm, C->ncomp, 1);
        C->talign.setVal(Real(0.0));
    }

    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        if (ratio[e] == 1) {
            C->K[e] = P.K[e];
            C->M[e] = P.M[e];
        } else {
            C->K[e] = galerkin1d(P.K[e], modulate[e]);
            C->M[e] = galerkin1d(P.M[e], modulate[e]);
        }
    }
    alloc_level(*C);
    if (!bottom) { estimate_lmax(*C); }
    return C;
}

#if (AMREX_SPACEDIM == 2)
// Reduced right-hand side of the red and black problems from the fine b
// (needs 1 ghost node).
void reduce_rhs (CrossProj::Level const& T, CrossProj::Level& R)
{
    Real const w = -T.sten[9] / T.sten[13]; // minus corner over diagonal (1/4)
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
    Real const w = -T.sten[9] / T.sten[13];
    Real const dinv = Real(1.0) / T.sten[13];
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
#endif

// CPU reductions over the valid boxes with partial accumulators so that the
// loops vectorize. G takes an MFIter and returns the per-point functor
// f(i,j,k,n) of that box.
#ifndef AMREX_USE_GPU
constexpr int nacc = 8;

template <typename G>
Real cpu_sum (MultiFab const& mf, int ncomp, G const& g)
{
    Real s = Real(0.0);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const lo = amrex::lbound(bx);
        auto const hi = amrex::ubound(bx);
        auto const f = g(mfi);
        for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            Real acc[nacc] = {};
            int i = lo.x;
            for (; i + nacc - 1 <= hi.x; i += nacc) {
                for (int l = 0; l < nacc; ++l) { acc[l] += f(i+l, j, k, n); }
            }
            for (; i <= hi.x; ++i) { acc[0] += f(i, j, k, n); }
            for (Real a : acc) { s += a; }
        }}}
    }
    return s;
}

template <typename G>
Real cpu_max (MultiFab const& mf, int ncomp, G const& g)
{
    Real s = Real(0.0);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const lo = amrex::lbound(bx);
        auto const hi = amrex::ubound(bx);
        auto const f = g(mfi);
        for (int n = 0; n < ncomp; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            Real acc[nacc] = {};
            int i = lo.x;
            for (; i + nacc - 1 <= hi.x; i += nacc) {
                for (int l = 0; l < nacc; ++l) { acc[l] = std::max(acc[l], f(i+l, j, k, n)); }
            }
            for (; i <= hi.x; ++i) { acc[0] = std::max(acc[0], f(i, j, k, n)); }
            for (Real a : acc) { s = std::max(s, a); }
        }}}
    }
    return s;
}
#endif

// u.v (global)
Real dot (MultiFab const& u, MultiFab const& v, int ncomp)
{
#ifdef AMREX_USE_GPU
    auto const& ua = u.const_arrays();
    auto const& va = v.const_arrays();
    Real s = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, u, IntVect(0),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            for (int n = 0; n < ncomp; ++n) {
                t += ua[bno](i,j,k,n) * va[bno](i,j,k,n);
            }
            return {t};
        });
#else
    Real s = cpu_sum(u, ncomp, [&] (MFIter const& mfi)
    {
        Array4<Real const> const ua = u.const_array(mfi);
        Array4<Real const> const va = v.const_array(mfi);
        return [=] (int i, int j, int k, int n) { return ua(i,j,k,n) * va(i,j,k,n); };
    });
#endif
    ParallelDescriptor::ReduceRealSum(s);
    return s;
}

// max |u| (global)
Real maxabs (MultiFab const& u, int ncomp)
{
#ifdef AMREX_USE_GPU
    Real s = Real(0.0);
    for (int n = 0; n < ncomp; ++n) {
        s = std::max(s, u.norminf(n, 0, true));
    }
#else
    Real s = cpu_max(u, ncomp, [&] (MFIter const& mfi)
    {
        Array4<Real const> const ua = u.const_array(mfi);
        return [=] (int i, int j, int k, int n) { return std::abs(ua(i,j,k,n)); };
    });
#endif
    ParallelDescriptor::ReduceRealMax(s);
    return s;
}

// q = A p (p needs filled ghosts); returns p.q (global)
Real apply_dot (CrossProj::Level const& P, MultiFab& q, MultiFab const& p)
{
    BL_PROFILE("GMG::CrossProj::pcg_apply");
    auto const sten = P.sten;
    int const nc = P.ncomp;
#ifdef AMREX_USE_GPU
    auto const& qa = q.arrays();
    auto const& pa = p.const_arrays();
    Real s = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, q, IntVect(0),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            for (int n = 0; n < nc; ++n) {
                Real const v = stencil_sum(sten, pa[bno], i, j, k, n);
                qa[bno](i,j,k,n) = v;
                t += v * pa[bno](i,j,k,n);
            }
            return {t};
        });
#else
    // q = A p on each row, then the row's contribution to p.q
    Real s = Real(0.0);
    for (MFIter mfi(q); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const lo = amrex::lbound(bx);
        auto const hi = amrex::ubound(bx);
        Array4<Real> const qq = q.array(mfi);
        Array4<Real const> const pp = p.const_array(mfi);
        for (int n = 0; n < nc; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                qq(i,j,k,n) = stencil_sum(sten, pp, i, j, k, n);
            }
            Real acc[nacc] = {};
            int i = lo.x;
            for (; i + nacc - 1 <= hi.x; i += nacc) {
                for (int l = 0; l < nacc; ++l) { acc[l] += qq(i+l,j,k,n) * pp(i+l,j,k,n); }
            }
            for (; i <= hi.x; ++i) { acc[0] += qq(i,j,k,n) * pp(i,j,k,n); }
            for (Real a : acc) { s += a; }
        }}}
    }
#endif
    ParallelDescriptor::ReduceRealSum(s);
    return s;
}

// x += alpha p, r -= alpha q; returns max |r| (global)
Real update_norm (CrossProj::Level const& P, MultiFab& x, MultiFab& r, MultiFab const& p,
                  MultiFab const& q, Real alpha)
{
    BL_PROFILE("GMG::CrossProj::pcg_update");
    int const nc = P.ncomp;
    auto const& xa = x.arrays();
    auto const& ra = r.arrays();
    auto const& pa = p.const_arrays();
    auto const& qa = q.const_arrays();
#ifdef AMREX_USE_GPU
    Real s = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, r, IntVect(0),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            for (int n = 0; n < nc; ++n) {
                xa[bno](i,j,k,n) += alpha * pa[bno](i,j,k,n);
                Real const v = ra[bno](i,j,k,n) - alpha * qa[bno](i,j,k,n);
                ra[bno](i,j,k,n) = v;
                t = std::max(t, std::abs(v));
            }
            return {t};
        });
#else
    ParallelFor(r, IntVect(0), nc, [=] (int bno, int i, int j, int k, int n) noexcept
    {
        xa[bno](i,j,k,n) += alpha * pa[bno](i,j,k,n);
        ra[bno](i,j,k,n) -= alpha * qa[bno](i,j,k,n);
    });
    Real s = cpu_max(r, nc, [&] (MFIter const& mfi)
    {
        Array4<Real const> const rr = r.const_array(mfi);
        return [=] (int i, int j, int k, int n) { return std::abs(rr(i,j,k,n)); };
    });
#endif
    ParallelDescriptor::ReduceRealMax(s);
    return s;
}

}

CrossProj::CrossProj (Geometry const& geom)
    : m_geom(geom)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.isAllPeriodic(), "GMG::CrossProj: only supports periodic BC");
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(geom.Domain().length(e) % 2 == 0 &&
                                         geom.Domain().length(e) >= 4 &&
                                         geom.Domain().smallEnd(e) == 0,
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
        for (int e = 0; e < AMREX_SPACEDIM; ++e) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(len[e] % 2 == 0, "GMG::CrossProj: box lengths must be even");
        }
    }
    m_grids = ba;
    m_dmap = dm;
    m_vel.define(m_grids, m_dmap, AMREX_SPACEDIM, 1);

    m_top = std::make_unique<Level>();
    Level& T = *m_top;
    T.domain = m_geom.Domain();
    T.period = m_geom.periodicity();
    T.cba = ba;
    T.dm = dm;
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        Real const hinv2 = Real(1.0) / (m_geom.CellSize(e)*m_geom.CellSize(e));
        T.K[e] = {-hinv2, Real(2.0)*hinv2, -hinv2};
        T.M[e] = {Real(0.25), Real(0.5), Real(0.25)};
    }

    for (auto& h : m_hier) { h.clear(); }

#if (AMREX_SPACEDIM == 2)
    {
        Real const rdx = m_geom.CellSize(0) / m_geom.CellSize(1);
        constexpr Real tol = std::is_same_v<Real,float> ? Real(1.e-5) : Real(1.e-12);
        m_reduced = std::abs(Real(1.0) - rdx*rdx) <= tol;
    }
#endif

    Level* pcglev = &T;
    if (m_reduced) {
#if (AMREX_SPACEDIM == 2)
        set_stencil(T);
        T.b.define(ba, dm, 1, 1);
        T.x.define(ba, dm, 1, 1);
        auto R = std::make_unique<Level>();
        R->domain = amrex::coarsen(T.domain, 2);
        R->period = Periodicity(R->domain.length());
        R->ratio = IntVect(2);
        R->ncomp = 2;
        R->cheby_degree = 2;
        R->cba = amrex::coarsen(ba, 2);
        R->dm = dm;
        Real const a0 = T.sten[13];
        for (int e = 0; e < AMREX_SPACEDIM; ++e) {
            R->K[e] = {-a0/Real(4.0), a0/Real(2.0), -a0/Real(4.0)};
            R->M[e] = {Real(0.125), Real(0.75), Real(0.125)};
        }
        alloc_level(*R);
        estimate_lmax(*R);
        m_hier[0].push_back(std::move(R));
        pcglev = m_hier[0][0].get();
#endif
    } else {
        alloc_level(T);
        estimate_lmax(T);
    }

    // Hierarchy 0: smooth modes. 2D: 1 for the checkerboard. 3D: 1-3 for the
    // families with the free direction x, y, z.
    for (int h = 0; h < (m_reduced ? 1 : nhier); ++h) {
        int const free = (AMREX_SPACEDIM == 3) ? h - 1 : -1; // line direction
        IntVect rr(2);
        IntVect mod(h > 0 ? 1 : 0);
        if (free >= 0) {
            rr[free] = 1;
            mod[free] = 0;
        }
        // Coarsen while possible, until the level is small enough for the
        // exact bottom solve. The fine grid always gets a coarse level.
        Level const* parent = m_hier[h].empty() ? &T : m_hier[h].back().get();
        bool more = m_hier[h].empty() || coarsenable(parent->domain, rr);
        while (more) {
            m_hier[h].push_back(make_coarse(*parent, rr, free,
                                            m_hier[h].empty() ? mod : IntVect(0)));
            parent = m_hier[h].back().get();
            more = coarsenable(parent->domain, rr);
        }
    }

    m_x.define(pcglev->cba, pcglev->dm, pcglev->ncomp, 1);
    m_p.define(pcglev->cba, pcglev->dm, pcglev->ncomp, 1);
    m_q.define(pcglev->cba, pcglev->dm, pcglev->ncomp, 0);
    if (!m_reduced) {
        m_rold.define(pcglev->cba, pcglev->dm, pcglev->ncomp, 0);
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
    // G averages the differences over the 2^(dim-1) edges of a cell.
    Real const avg = (AMREX_SPACEDIM == 3) ? Real(0.25) : Real(0.5);
    GpuArray<Real,3> fac{};
    for (int e = 0; e < AMREX_SPACEDIM; ++e) { fac[e] = avg*dxinv[e]; }

    // b = G^T u = -D u at the node (i,j,k), from the cells (i-1..i, ...)
    {
        BL_PROFILE("GMG::CrossProj::rhs");
        auto const& ba = T.b.arrays();
        auto const& va = m_vel.const_arrays();
        ParallelFor(T.b, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
        {
            auto const& v = va[bno];
            Real s = Real(0.0);
            for (int c = 0; c <= kz; ++c) {
            for (int bb = 0; bb < 2; ++bb) {
            for (int a = 0; a < 2; ++a) {
                int const ii = i-1+a, jj = j-1+bb, kk = k-kz+c;
                s += (a  ? -fac[0] : fac[0]) * v(ii,jj,kk,0)
                   + (bb ? -fac[1] : fac[1]) * v(ii,jj,kk,1);
#if (AMREX_SPACEDIM == 3)
                s += (c  ? -fac[2] : fac[2]) * v(ii,jj,kk,2);
#endif
            }}}
            ba[bno](i,j,k) = s;
        });
    }

    Real const bnorm = maxabs(T.b, 1);
    // Below the roundoff in G^T u the residual cannot decrease, e.g., when
    // projecting a field that is already divergence free.
    Real const floor = Real(100.0) * std::numeric_limits<Real>::epsilon()
        * maxabs(m_vel, AMREX_SPACEDIM) * (AMREX_D_TERM(dxinv[0], + dxinv[1], + dxinv[2]));
    Real const target = std::max({m_reltol*bnorm, m_abstol, floor});

#if (AMREX_SPACEDIM == 2)
    if (m_reduced) {
        BL_PROFILE("GMG::CrossProj::reduce");
        T.b.FillBoundary(T.period);
        reduce_rhs(T, P);
    }
#endif

    m_x.setVal(Real(0.0));
    m_niters = 0;
    Real rnorm = bnorm;
    if (rnorm > target) {
        // PCG; with the non-symmetric fine-grid preconditioner the flexible
        // variant (beta from z.(r_new - r_old)) is used.
        bool const flexible = !m_reduced;
        precond();
        MultiFab::Copy(m_p, P.x, 0, 0, nc, 0);
        Real rz = dot(P.b, P.x, nc);
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
            Real const rz1 = dot(P.b, P.x, nc);
            Real const beta = (flexible ? rz1 - dot(m_rold, P.x, nc) : rz1) / rz;
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
#if (AMREX_SPACEDIM == 2)
    if (m_reduced) {
        BL_PROFILE("GMG::CrossProj::backsub");
        m_x.FillBoundary(P.period);
        back_substitute(T, m_x);
    }
#endif
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
            Real g[3] = {Real(0.0), Real(0.0), Real(0.0)};
            for (int c = 0; c <= kz; ++c) {
            for (int bb = 0; bb < 2; ++bb) {
            for (int a = 0; a < 2; ++a) {
                Real const v = ph(i+a,j+bb,k+c);
                g[0] += a  ? v : -v;
                g[1] += bb ? v : -v;
                g[2] += c  ? v : -v;
            }}}
            for (int n = 0; n < AMREX_SPACEDIM; ++n) {
                sa[bno](i,j,k,n) = va[bno](i,j,k,n) - fac[n] * g[n];
            }
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
