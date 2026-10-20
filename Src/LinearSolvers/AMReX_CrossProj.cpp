#include <AMReX_CrossProj.H>

#if (AMREX_SPACEDIM >= 2)

#include <AMReX_MultiFabUtil.H>

#include <algorithm>
#include <limits>

namespace amrex::GMG {

// Solves G^T G phi = G^T u (= -D u) by PCG, then u -= G phi. G is the
// cell-centered gradient of nodal phi. The preconditioner is a multigrid
// cycle with a coarse space for the smooth modes (linear interpolation) and
// coarse spaces for the near-null modes: in 2D the checkerboard (-1)^(i+j)
// times smooth modes, in 3D the families (-1)^(j+k) f(i), (-1)^(i+k) f(j),
// (-1)^(i+j) f(k), coarsened only transverse to the free direction and
// smoothed with line solves along it.
//
// Every operator in the hierarchy is sum_d K_d (x) M_e (x) M_f with 1D
// three-point stencils K and M per direction: the fine one has K = second
// difference and M = [1/4,1/2,1/4], and Galerkin coarsening with
// tensor-product interpolation keeps the form. On a periodic grid all of
// them are constant, so each level is a constant 9- or 27-point stencil.
struct CrossProj::Level
{
    Box domain;            // cell-centered
    Periodicity period;
    IntVect ratio{1};      // coarsening ratio from the parent level
    int line_dir = -1;     // family levels: line smoother along this direction
    IntVect modulate{0};   // first near-null level: P has (-1)^(sum of these indices)
    BoxArray cba;          // cell-centered
    DistributionMapping dm;
    BoxArray aligned;      // nodal parent BoxArray coarsened, if cba was regridded
    Array<Array<Real,3>,AMREX_SPACEDIM> K, M;  // [dir][offset+1]
    GpuArray<Real,27> sten{};  // [(c+1)*9+(b+1)*3+(a+1)] for offset (a,b,c); c = 0 in 2D
    std::unique_ptr<iMultiFab> owner;
    MultiFab x, b, r, y, d, t;
    Real lmax = Real(0.0); // Chebyshev upper bound
};

namespace {

constexpr int cheby_degree = 3;
constexpr Real cheby_ratio = Real(16.0);
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

void stencil_apply (CrossProj::Level const& L, MultiFab& out, MultiFab& in,
                    MultiFab const* b = nullptr)
{
    BL_PROFILE("GMG::CrossProj::apply");
    in.FillBoundary(L.period);
    auto const sten = L.sten;
    auto const& o = out.arrays();
    auto const& xin = in.const_arrays();
    bool const resid = (b != nullptr);
    auto const& bb = resid ? b->const_arrays() : xin; // r = b - A x if b is given
    ParallelFor(out, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
    {
        auto const& xa = xin[bno];
        Real s = Real(0.0);
        for (int c = -kz; c <= kz; ++c) {
        for (int bj = -1; bj <= 1; ++bj) {
        for (int a = -1; a <= 1; ++a) {
            s += sten[(c+1)*9+(bj+1)*3+(a+1)] * xa(i+a,j+bj,k+c);
        }}}
        o[bno](i,j,k) = resid ? bb[bno](i,j,k) - s : s;
    });
    Gpu::streamSynchronize();
}

// r = b - A x
void residual (CrossProj::Level const& L, MultiFab& r, MultiFab const& b, MultiFab& x)
{
    stencil_apply(L, r, x, &b);
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
        ParallelFor(y, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
        {
            ya[bno](i,j,k) = ra[bno](i,j,k) * dinv;
        });
        Gpu::streamSynchronize();
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

// Chebyshev iteration in B^-1 A on [lmax/ratio, lmax].
void cheby (CrossProj::Level& L, MultiFab& x, MultiFab const& b, int nsweeps, bool zero_init)
{
    BL_PROFILE("GMG::CrossProj::cheby");
    Real const hi = L.lmax;
    Real const lo = hi / cheby_ratio;
    Real const theta = Real(0.5)*(hi + lo);
    Real const delta = Real(0.5)*(hi - lo);
    Real const sig = theta / delta;
    for (int sweep = 0; sweep < nsweeps; ++sweep) {
        if (zero_init && sweep == 0) {
            x.setVal(Real(0.0));
            MultiFab::Copy(L.r, b, 0, 0, 1, 0);
        } else {
            residual(L, L.r, b, x);
        }
        apply_prec(L, L.y, L.r);
        MultiFab::Copy(L.d, L.y, 0, 0, 1, 0);
        L.d.mult(Real(1.0)/theta);
        MultiFab::Add(x, L.d, 0, 0, 1, 0);
        Real rho = Real(1.0) / sig;
        for (int kk = 1; kk < cheby_degree; ++kk) {
            residual(L, L.r, b, x);
            apply_prec(L, L.y, L.r);
            Real const rho1 = Real(1.0) / (Real(2.0)*sig - rho);
            MultiFab::LinComb(L.d, rho1*rho, L.d, 0, Real(2.0)*rho1/delta, L.y, 0, 0, 1, 0);
            MultiFab::Add(x, L.d, 0, 0, 1, 0);
            rho = rho1;
        }
    }
}

// Upper end of the spectrum of B^-1 A by power iteration, with a margin.
void estimate_lmax (CrossProj::Level& L)
{
    auto const& va = L.x.arrays();
    ParallelFor(L.x, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
    {
        auto h = static_cast<unsigned int>(i*73856093 ^ j*19349663 ^ k*83492791);
        h ^= h >> 13; h *= 0x5bd1e995U; h ^= h >> 15;
        va[bno](i,j,k) = Real(h % 1000U) / Real(1000.0) - Real(0.5);
    });
    Gpu::streamSynchronize();
    Real lam = Real(0.0);
    for (int it = 0; it < 30; ++it) {
        Real const nrm = std::sqrt(MultiFab::Dot(*L.owner, L.x, 0, L.x, 0, 1, 0));
        L.x.mult(Real(1.0)/nrm);
        stencil_apply(L, L.r, L.x);
        apply_prec(L, L.y, L.r);
        lam = std::sqrt(MultiFab::Dot(*L.owner, L.y, 0, L.y, 0, 1, 0));
        MultiFab::Copy(L.x, L.y, 0, 0, 1, 0);
    }
    L.lmax = Real(1.1) * lam;
    L.x.setVal(Real(0.0));
}

// cb = P^T fr
void restrict_to (CrossProj::Level const& F, CrossProj::Level const& C, MultiFab& fr, MultiFab& cb)
{
    BL_PROFILE("GMG::CrossProj::restrict");
    fr.FillBoundary(F.period);
    MultiFab tmp;
    MultiFab* dst = &cb;
    if (!C.aligned.empty()) {
        tmp.define(C.aligned, F.dm, 1, 0);
        dst = &tmp;
    }
    GpuArray<int,3> rr{1,1,1};
    GpuArray<int,3> md{0,0,0};
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        rr[e] = C.ratio[e];
        md[e] = C.modulate[e];
    }
    auto const& fa = fr.const_arrays();
    auto const& ca = dst->arrays();
    ParallelFor(*dst, [=] AMREX_GPU_DEVICE (int bno, int I, int J, int K) noexcept
    {
        int const ax = (rr[0] == 2) ? 1 : 0;
        int const ay = (rr[1] == 2) ? 1 : 0;
        int const az = (rr[2] == 2) ? 1 : 0;
        Real s = Real(0.0);
        for (int c = -az; c <= az; ++c) {
        for (int b = -ay; b <= ay; ++b) {
        for (int a = -ax; a <= ax; ++a) {
            int const i = rr[0]*I + a;
            int const j = rr[1]*J + b;
            int const k = rr[2]*K + c;
            Real w = ((a == 0) ? Real(1.0) : Real(0.5)) * ((b == 0) ? Real(1.0) : Real(0.5))
                   * ((c == 0) ? Real(1.0) : Real(0.5));
            if ((md[0]*i + md[1]*j + md[2]*k) & 1) { w = -w; }
            s += w * fa[bno](i,j,k);
        }}}
        ca[bno](I,J,K) = s;
    });
    Gpu::streamSynchronize();
    if (dst != &cb) {
        cb.ParallelCopy(tmp, 0, 0, 1);
    }
}

// fx += P cx
void prolong_add (CrossProj::Level const& F, CrossProj::Level const& C, MultiFab& fx, MultiFab const& cx)
{
    BL_PROFILE("GMG::CrossProj::prolong");
    MultiFab tmp;
    MultiFab const* src = &cx;
    if (!C.aligned.empty()) {
        tmp.define(C.aligned, F.dm, 1, 0);
        tmp.ParallelCopy(cx, 0, 0, 1);
        src = &tmp;
    }
    GpuArray<int,3> rr{1,1,1};
    GpuArray<int,3> md{0,0,0};
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        rr[e] = C.ratio[e];
        md[e] = C.modulate[e];
    }
    auto const& fa = fx.arrays();
    auto const& ca = src->const_arrays();
    ParallelFor(fx, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
    {
        int const iv[3] = {i, j, k};
        int c0[3], n[3];
        for (int e = 0; e < 3; ++e) {
            if (rr[e] == 1) {
                c0[e] = iv[e]; n[e] = 1;
            } else if ((iv[e] & 1) == 0) {
                c0[e] = iv[e] / 2; n[e] = 1;
            } else {
                c0[e] = (iv[e] - 1) / 2; n[e] = 2;
            }
        }
        Real s = Real(0.0);
        for (int c = 0; c < n[2]; ++c) {
        for (int b = 0; b < n[1]; ++b) {
        for (int a = 0; a < n[0]; ++a) {
            s += ca[bno](c0[0]+a, c0[1]+b, c0[2]+c);
        }}}
        s /= Real(n[0]*n[1]*n[2]);
        if ((md[0]*i + md[1]*j + md[2]*k) & 1) { s = -s; }
        fa[bno](i,j,k) += s;
    });
    Gpu::streamSynchronize();
}

void alloc_level (CrossProj::Level& L)
{
    BoxArray const nba = amrex::convert(L.cba, IntVect(1));
    L.x.define(nba, L.dm, 1, 1);
    L.r.define(nba, L.dm, 1, 1);
    L.b.define(nba, L.dm, 1, 0);
    L.y.define(nba, L.dm, 1, 0);
    L.d.define(nba, L.dm, 1, 0);
    L.t.define(nba, L.dm, 1, 0);
    L.x.setVal(Real(0.0));
    L.r.setVal(Real(0.0));
    L.owner = amrex::OwnerMask(L.x, L.period);
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
    alloc_level(T);
    estimate_lmax(T);

    // Hierarchy 0: smooth modes. 2D: 1 for the checkerboard. 3D: 1-3 for the
    // families with the free direction x, y, z. Coarsen while the coarse level
    // keeps at least 4 cells in each coarsened direction.
    for (int h = 0; h < nhier; ++h) {
        m_hier[h].clear();
        int const free = (AMREX_SPACEDIM == 3) ? h - 1 : -1; // line direction
        IntVect rr(2);
        IntVect mod(h > 0 ? 1 : 0);
        if (free >= 0) {
            rr[free] = 1;
            mod[free] = 0;
        }
        Level const* parent = &T;
        bool more = true;
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
}

void
CrossProj::vcycle (int h, int lev)
{
    auto& hier = m_hier[h];
    Level& L = *hier[lev];
    if (lev+1 == int(hier.size())) {
        cheby(L, L.x, L.b, bottom_sweeps, true);
        return;
    }
    Level& C = *hier[lev+1];
    cheby(L, L.x, L.b, 1, true);
    residual(L, L.r, L.b, L.x);
    restrict_to(L, C, L.r, C.b);
    vcycle(h, lev+1);
    prolong_add(L, C, L.x, C.x);
    cheby(L, L.x, L.b, 1, false);
}

// Symmetric multigrid cycle: smoothing, a correction from each hierarchy in
// order and then in reverse order, smoothing.
void
CrossProj::precond (MultiFab& z, MultiFab const& r)
{
    Level& T = *m_top;
    cheby(T, z, r, 1, true);
    for (int pass = 0; pass < 2; ++pass) {
        for (int ih = 0; ih < nhier; ++ih) {
            int const h = (pass == 0) ? ih : nhier-1-ih;
            Level& C = *m_hier[h][0];
            residual(T, T.r, r, z);
            restrict_to(T, C, T.r, C.b);
            vcycle(h, 0);
            prolong_add(T, C, z, C.x);
        }
    }
    cheby(T, z, r, 1, false);
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

    MultiFab vel(m_grids, m_dmap, AMREX_SPACEDIM, 1);
    MultiFab::Copy(vel, a_rhs, 0, 0, AMREX_SPACEDIM, 0);
    vel.FillBoundary(m_geom.periodicity());

    BoxArray const nba = amrex::convert(m_grids, IntVect(1));
    MultiFab b(nba, m_dmap, 1, 0);
    MultiFab x(nba, m_dmap, 1, 1);
    MultiFab r(nba, m_dmap, 1, 1);
    MultiFab z(nba, m_dmap, 1, 1);
    MultiFab p(nba, m_dmap, 1, 1);
    MultiFab q(nba, m_dmap, 1, 0);
    for (auto* mf : {&x, &r, &z, &p}) { mf->setVal(Real(0.0)); }

    auto const dxinv = m_geom.InvCellSizeArray();
    // G averages the differences over the 2^(dim-1) edges of a cell.
    Real const avg = (AMREX_SPACEDIM == 3) ? Real(0.25) : Real(0.5);
    GpuArray<Real,3> fac{};
    for (int e = 0; e < AMREX_SPACEDIM; ++e) { fac[e] = avg*dxinv[e]; }

    // b = G^T u = -D u
    {
        auto const& ba = b.arrays();
        auto const& va = vel.const_arrays();
        ParallelFor(b, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
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
        Gpu::streamSynchronize();
    }

    // PCG; dot products count each node once.
    auto dot = [&] (MultiFab const& u, MultiFab const& v) {
        return MultiFab::Dot(*T.owner, u, 0, v, 0, 1, 0);
    };
    Real const bnorm = b.norminf();
    // Below the roundoff in G^T u the residual cannot decrease, e.g., when
    // projecting a field that is already divergence free.
    Real const floor = Real(100.0) * std::numeric_limits<Real>::epsilon()
        * vel.norminf(0, AMREX_SPACEDIM, IntVect(0)) * (AMREX_D_TERM(dxinv[0], + dxinv[1], + dxinv[2]));
    Real const target = std::max({m_reltol*bnorm, m_abstol, floor});
    MultiFab::Copy(r, b, 0, 0, 1, 0);
    m_niters = 0;
    Real rnorm = bnorm;
    if (rnorm > target) {
        precond(z, r);
        MultiFab::Copy(p, z, 0, 0, 1, 0);
        Real rz = dot(r, z);
        for (int iter = 1; iter <= m_maxiter; ++iter) {
            stencil_apply(T, q, p);
            Real const alpha = rz / dot(p, q);
            MultiFab::Saxpy(x, alpha, p, 0, 0, 1, 0);
            MultiFab::Saxpy(r, -alpha, q, 0, 0, 1, 0);
            rnorm = r.norminf();
            m_niters = iter;
            if (m_verbose > 1) {
                amrex::Print() << "GMG::CrossProj: iteration " << iter << " residual "
                               << rnorm << " (" << rnorm/bnorm << ")\n";
            }
            if (rnorm <= target) { break; }
            precond(z, r);
            Real const rz1 = dot(r, z);
            MultiFab::Xpay(p, rz1/rz, z, 0, 0, 1, 0);
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

    // u -= G phi
    x.FillBoundary(m_geom.periodicity());
    {
        auto const& va = vel.arrays();
        auto const& pa = x.const_arrays();
        ParallelFor(vel, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
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
                va[bno](i,j,k,n) -= fac[n] * g[n];
            }
        });
        Gpu::streamSynchronize();
    }

    MultiFab::Copy(a_soln, vel, 0, 0, AMREX_SPACEDIM, 0);
    IntVect const ng = amrex::elemwiseMin(a_soln.nGrowVect(), IntVect(1));
    if (ng.max() > 0) {
        a_soln.FillBoundary(0, AMREX_SPACEDIM, ng, m_geom.periodicity());
    }
}

}

#endif
