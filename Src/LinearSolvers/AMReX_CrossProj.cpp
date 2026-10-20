#include <AMReX_CrossProj.H>

#if (AMREX_SPACEDIM >= 2)

#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParReduce.H>

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
    Box domain;            // cell-centered
    Periodicity period;
    IntVect ratio{1};      // coarsening ratio from the parent level
    int line_dir = -1;     // family levels: line smoother along this direction
    IntVect modulate{0};   // first near-null level: P has (-1)^(sum of these indices)
    int ncomp = 1;
    BoxArray cba;          // cell-centered
    DistributionMapping dm;
    BoxArray aligned;      // nodal parent BoxArray coarsened, if cba was regridded
    Array<Array<Real,3>,AMREX_SPACEDIM> K, M;  // [dir][offset+1]
    GpuArray<Real,27> sten{};  // [(c+1)*9+(b+1)*3+(a+1)] for offset (a,b,c); c = 0 in 2D
    GpuArray<Real,27> rw{};    // restriction weights from the parent, same layout
    GpuArray<Real,8> pw{};     // prolongation weights [c*4+b*2+a] for fine offset (a,b,c)
    std::unique_ptr<iMultiFab> owner;
    MultiFab x, x2, b, r, y, t, talign;
    Real lmax = Real(0.0); // Chebyshev upper bound
    int cheby_degree = 3;  // 2 on the reduced 2D hierarchy
};

namespace {

constexpr int bottom_sweeps = 4;
constexpr int regrid_max_size = 32;
constexpr int kz = (AMREX_SPACEDIM == 3) ? 1 : 0;
constexpr int nhier = (AMREX_SPACEDIM == 3) ? 4 : 2;

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
    Real r = Real(0.0);
    for (int c = -1; c <= 1; ++c) {
    for (int b = -1; b <= 1; ++b) {
    for (int a = -1; a <= 1; ++a) {
        r += s[(c+1)*9+(b+1)*3+(a+1)] * x(i+a,j+b,k+c,n);
    }}}
    return r;
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
// coupling the owned nodes of a box along the line direction.
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

    int const f = L.line_dir;
    int const stride = (f == 0) ? 1 : ((f == 1) ? 3 : 9);
    Real const tlo = L.sten[13-stride]; // line part of the stencil: (tlo, tdiag, tlo)
    Real const tdiag = L.sten[13];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(y); mfi.isValid(); ++mfi) {
        Box const& vbx = mfi.validbox();
        Box lbx = vbx;
        lbx.setBig(f, vbx.smallEnd(f));
        int const mlo = vbx.smallEnd(f);
        int const mhi = vbx.bigEnd(f);
        auto const& ya = y.array(mfi);
        auto const& ra = r.const_array(mfi);
        auto const& ta = L.t.array(mfi);
        auto const& own = L.owner->const_array(mfi);
        ParallelFor(lbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            IntVect iv(AMREX_D_DECL(i,j,k));
            amrex::ignore_unused(k);
            for (int m = mlo; m <= mhi; ++m) {
                iv[f] = m;
                ya(iv) = Real(0.0);
            }
            // The owned nodes of a line in a box are contiguous.
            int s = mlo;
            int e = mhi;
            iv[f] = s;
            while (s <= mhi && !own(iv)) { iv[f] = ++s; }
            iv[f] = e;
            while (e >= s && !own(iv)) { iv[f] = --e; }
            if (s > e) { return; }
            // Thomas algorithm: t holds the modified upper diagonal, y the
            // modified right-hand side, then the solution.
            Real cprev = Real(0.0);
            Real dprev = Real(0.0);
            for (int m = s; m <= e; ++m) {
                iv[f] = m;
                Real const a = (m > s) ? tlo : Real(0.0);
                Real const c = (m < e) ? tlo : Real(0.0);
                Real const den = tdiag - a*cprev;
                cprev = c / den;
                dprev = (ra(iv) - a*dprev) / den;
                ta(iv) = cprev;
                ya(iv) = dprev;
            }
            for (int m = e-1; m >= s; --m) {
                IntVect jv = iv;
                iv[f] = m;
                jv[f] = m+1;
                ya(iv) -= ta(iv) * ya(jv);
            }
        });
    }
    amrex::OverrideSync(y, *L.owner, L.period);
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
        auto h = static_cast<unsigned int>(i*73856093 ^ j*19349663 ^ k*83492791 ^ n*2654435761U);
        h ^= h >> 13; h *= 0x5bd1e995U; h ^= h >> 15;
        va[bno](i,j,k,n) = Real(h % 1000U) / Real(1000.0) - Real(0.5);
    });
    L.b.setVal(Real(0.0));
    Real lam = Real(0.0);
    for (int it = 0; it < 30; ++it) {
        Real const nrm = std::sqrt(MultiFab::Dot(*L.owner, L.x, 0, L.x, 0, L.ncomp, 0));
        L.x.mult(Real(1.0)/nrm, 0, L.ncomp, 0);
        residual(L); // r = -A x
        apply_prec(L, L.y, L.r);
        lam = std::sqrt(MultiFab::Dot(*L.owner, L.y, 0, L.y, 0, L.ncomp, 0));
        MultiFab::Copy(L.x, L.y, 0, 0, L.ncomp, 0);
    }
    L.lmax = Real(1.1) * lam;
}

// C.b = P^T F.r (weights C.rw)
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
// (RX I + a, ...), a < RX; at the high end of a box those are ghost nodes
// of F.x, which the next exchange overwrites.
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
    if (!C.aligned.empty()) {
        C.talign.ParallelCopy(C.x, 0, 0, C.ncomp);
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

void alloc_level (CrossProj::Level& L)
{
    BoxArray const nba = amrex::convert(L.cba, IntVect(1));
    L.x.define(nba, L.dm, L.ncomp, 1);
    L.x2.define(nba, L.dm, L.ncomp, 1);
    L.r.define(nba, L.dm, L.ncomp, 1);
    L.b.define(nba, L.dm, L.ncomp, 0);
    L.y.define(nba, L.dm, L.ncomp, 0);
    if (L.line_dir >= 0) {
        L.t.define(nba, L.dm, L.ncomp, 0);
    }
    L.x.setVal(Real(0.0));
    L.x2.setVal(Real(0.0));
    L.r.setVal(Real(0.0));
    L.owner = amrex::OwnerMask(L.x, L.period);
    set_stencil(L);
    set_transfer_weights(L);
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

    BoxArray const aligned = amrex::coarsen(P.cba, ratio);
    bool ok = true;
    for (int i = 0, N = int(aligned.size()); i < N && ok; ++i) {
        for (int e = 0; e < AMREX_SPACEDIM; ++e) {
            int const len = aligned[i].length(e);
            if (ratio[e] == 2 && (len % 2 != 0 || len < 4)) { ok = false; }
        }
    }
    if (ok) {
        C->cba = aligned;
        C->dm = P.dm;
    } else {
        C->cba = BoxArray(C->domain);
        C->cba.maxSize(regrid_max_size);
        C->dm = DistributionMapping(C->cba);
        C->aligned = amrex::convert(aligned, IntVect(1));
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
    estimate_lmax(*C);
    return C;
}

#if (AMREX_SPACEDIM == 2)
// Reduced right-hand side of the red and black problems from the fine b
// (needs 2 ghost nodes).
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

// Reductions count each node once: on the GPU through the owner mask, on
// the CPU by looping over the cells enclosed by each nodal box, with
// partial accumulators so that the loops vectorize. G takes an MFIter and
// returns the per-point functor f(i,j,k,n) of that box.
#ifndef AMREX_USE_GPU
constexpr int nacc = 8;

template <typename G>
Real cpu_sum (MultiFab const& mf, int ncomp, G const& g)
{
    Real s = Real(0.0);
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        Box const& bx = amrex::enclosedCells(mfi.validbox());
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

// sum over owned nodes of u*v (global)
Real dot_owned (iMultiFab const& owner, MultiFab const& u, MultiFab const& v, int ncomp)
{
#ifdef AMREX_USE_GPU
    auto const& oa = owner.const_arrays();
    auto const& ua = u.const_arrays();
    auto const& va = v.const_arrays();
    Real s = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, u, IntVect(0),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            if (oa[bno](i,j,k)) {
                for (int n = 0; n < ncomp; ++n) {
                    t += ua[bno](i,j,k,n) * va[bno](i,j,k,n);
                }
            }
            return {t};
        });
#else
    amrex::ignore_unused(owner);
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

// q = A p (p needs filled ghosts); returns p.q over owned nodes (global)
Real apply_dot (CrossProj::Level const& P, MultiFab& q, MultiFab const& p)
{
    BL_PROFILE("GMG::CrossProj::pcg_apply");
    auto const sten = P.sten;
    int const nc = P.ncomp;
#ifdef AMREX_USE_GPU
    auto const& qa = q.arrays();
    auto const& pa = p.const_arrays();
    auto const& oa = P.owner->const_arrays();
    Real s = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Real>{}, q, IntVect(0),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept -> GpuTuple<Real>
        {
            Real t = Real(0.0);
            for (int n = 0; n < nc; ++n) {
                Real const v = stencil_sum(sten, pa[bno], i, j, k, n);
                qa[bno](i,j,k,n) = v;
                if (oa[bno](i,j,k)) { t += v * pa[bno](i,j,k,n); }
            }
            return {t};
        });
#else
    // q = A p on each row, then the row's contribution to p.q
    Real s = Real(0.0);
    for (MFIter mfi(q); mfi.isValid(); ++mfi) {
        Box const& nbx = mfi.validbox();
        Box const& cbx = amrex::enclosedCells(nbx);
        auto const lo = amrex::lbound(nbx);
        auto const hi = amrex::ubound(nbx);
        auto const chi = amrex::ubound(cbx);
        Array4<Real> const qq = q.array(mfi);
        Array4<Real const> const pp = p.const_array(mfi);
        for (int n = 0; n < nc; ++n) {
        for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                qq(i,j,k,n) = stencil_sum(sten, pp, i, j, k, n);
            }
            if (j > chi.y || k > chi.z) { continue; }
            Real acc[nacc] = {};
            int i = lo.x;
            for (; i + nacc - 1 <= chi.x; i += nacc) {
                for (int l = 0; l < nacc; ++l) { acc[l] += qq(i+l,j,k,n) * pp(i+l,j,k,n); }
            }
            for (; i <= chi.x; ++i) { acc[0] += qq(i,j,k,n) * pp(i,j,k,n); }
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
                                         geom.Domain().smallEnd(e) % 2 == 0,
            "GMG::CrossProj: the domain needs an even number (>= 4) of cells starting at an even index");
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
        BoxArray const nba = amrex::convert(ba, IntVect(1));
        T.b.define(nba, dm, 1, 2);
        T.x.define(nba, dm, 1, 0);
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
    // families with the free direction x, y, z. Coarsen while the coarse level
    // keeps at least 4 cells in each coarsened direction.
    for (int h = 0; h < (m_reduced ? 1 : nhier); ++h) {
        int const free = (AMREX_SPACEDIM == 3) ? h - 1 : -1; // line direction
        IntVect rr(2);
        IntVect mod(h > 0 ? 1 : 0);
        if (free >= 0) {
            rr[free] = 1;
            mod[free] = 0;
        }
        Level const* parent = m_hier[h].empty() ? &T : m_hier[h].back().get();
        bool more = true; // the fine grid always gets a coarse level
        for (int e = 0; e < AMREX_SPACEDIM && !m_hier[h].empty(); ++e) {
            int const len = parent->domain.length(e);
            if (rr[e] == 2 && (len % 2 != 0 || len < 8)) { more = false; }
        }
        while (more) {
            m_hier[h].push_back(make_coarse(*parent, rr, free,
                                            m_hier[h].empty() ? mod : IntVect(0)));
            parent = m_hier[h].back().get();
            for (int e = 0; e < AMREX_SPACEDIM; ++e) {
                int const len = parent->domain.length(e);
                if (rr[e] == 2 && (len % 2 != 0 || len < 8)) { more = false; }
            }
        }
    }

    BoxArray const pba = amrex::convert(pcglev->cba, IntVect(1));
    m_x.define(pba, pcglev->dm, pcglev->ncomp, 1);
    m_p.define(pba, pcglev->dm, pcglev->ncomp, 1);
    m_q.define(pba, pcglev->dm, pcglev->ncomp, 0);
}

void
CrossProj::vcycle (int h, int lev)
{
    auto& hier = m_hier[h];
    Level& L = *hier[lev];
    if (lev+1 == int(hier.size())) {
        cheby(L, bottom_sweeps, true);
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
// grid: smoothing, a correction from each hierarchy in order and then in
// reverse order (symmetric), smoothing.
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
    for (int pass = 0; pass < 2; ++pass) {
        for (int ih = 0; ih < nhier; ++ih) {
            int const h = (pass == 0) ? ih : nhier-1-ih;
            Level& C = *m_hier[h][0];
            residual(T);
            restrict_to(T, C);
            vcycle(h, 0);
            prolong_add(T, C);
        }
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

    // b = G^T u = -D u
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
        precond();
        MultiFab::Copy(m_p, P.x, 0, 0, nc, 0);
        Real rz = dot_owned(*P.owner, P.b, P.x, nc);
        for (int iter = 1; iter <= m_maxiter; ++iter) {
            m_p.FillBoundary(P.period);
            Real const pq = apply_dot(P, m_q, m_p);
            Real const alpha = rz / pq;
            rnorm = update_norm(P, m_x, P.b, m_p, m_q, alpha);
            m_niters = iter;
            if (m_verbose > 1) {
                amrex::Print() << "GMG::CrossProj: iteration " << iter << " residual "
                               << rnorm << " (" << rnorm/bnorm << ")\n";
            }
            if (rnorm <= target) { break; }
            precond();
            BL_PROFILE_VAR("GMG::CrossProj::pcg_p", pcgp);
            Real const rz1 = dot_owned(*P.owner, P.b, P.x, nc);
            MultiFab::Xpay(m_p, rz1/rz, P.x, 0, 0, nc, 0);
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

#if (AMREX_SPACEDIM == 2)
    if (m_reduced) {
        BL_PROFILE("GMG::CrossProj::backsub");
        m_x.FillBoundary(P.period);
        back_substitute(T, m_x);
    }
#endif
    MultiFab const& phi = m_reduced ? T.x : m_x;

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
