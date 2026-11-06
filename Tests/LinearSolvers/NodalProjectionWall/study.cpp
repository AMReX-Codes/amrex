// The convergence study (study = 1): projections of the Burggraf cavity flow (burggraf), a
// channel with an inflow face and a pressure outflow, and a uniform inflow against no-slip walls
// (channel), and a closed 3D box (box3d, 3D only), with the viscous operators incflo uses
// applied to the projected velocity.  It prints tables (see README.md) and checks nothing.

#include "NodalProjectionWall.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>

using namespace amrex;

namespace {

// Channel flow for the inflow/outflow check, from the stream function psi = H(y) + F(x) G(y)
// with H = 3y^2 - 2y^3, F = x + x^2/2 and G = y^2 (1-y)^2: u = psi_y = 6y(1-y) + F G',
// v = -psi_x = -F' G.  It is divergence-free and zero on the walls y = 0 and 1, and u > 0 at
// x = 0 and 1.  On both walls v_n'' = v_yy = -2 F'(x) = -2 (1+x), and at the x faces
// u_n'' = u_xx = F'' G' = G'(y).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real chG  (Real y) { return y*y*(Real(1)-y)*(Real(1)-y); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real chG1 (Real y) { return Real(2)*y*(Real(1)-y)*(Real(1)-Real(2)*y); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real chF  (Real x) { return x + Real(0.5)*x*x; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real chF1 (Real x) { return Real(1) + x; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real chu (Real x, Real y) { return Real(6)*y*(Real(1)-y) + chF(x)*chG1(y); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real chv (Real x, Real y) { return -chF1(x)*chG(y); }

// The viscous operators tested: MLABecLaplacian then MLTensorOp, each at every maxorder in
// the `maxorders` input (default 2 3).  Operator op uses maxorders[op % nmo], with the
// tensor operator for op >= nmo.
Vector<int> maxorders{2, 3};
int nops () { return 2 * int(maxorders.size()); }
std::string op_name (int op)
{
    int const nmo = int(maxorders.size());
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%-16s maxorder %d",
                  (op < nmo) ? "MLABecLaplacian," : "MLTensorOp,", maxorders[op % nmo]);
    return std::string(buf);
}

struct Result
{
    Real e_max = 0;                       // max |u_projected - u| over all cells
    std::array<Real,4> e_delta{};         // wall-normal error / delta, cells 0-3 from the wall
    std::array<Real,4> e_corner{};        // |u_projected - u| / h^2 in the corner cells: (0,0), (n,0), (0,n), (n,n)
    Vector<Real> c_exact;                 // c for L_h applied to the exact velocity
    Vector<Real> c_proj;                  // c for L_h applied to the projected velocity
    Vector<Real> t_exact;                 // max |L_h u - lap u| of the tangential component in
    Vector<Real> t_proj;                  // the wall cells (as for c), exact and projected u
};

// Max |a(tangential component)| over the cells next to each of the four walls, leaving out
// n/8 cells at each end of every wall (the corners), as wall_fit does
Real wall_tangential_max (Array4<Real const> const& a, int comp0, int n)
{
    int const e = n/8;
    Real m = 0;
    for (int t = e; t < n-e; ++t) {
        for (int w : {0, n-1}) {
            m = std::max(m, std::abs(a(w,t,0,comp0+1)));   // x walls: v
            m = std::max(m, std::abs(a(t,w,0,comp0)));     // y walls: u
        }
    }
    return m;
}

// Least-squares fit of a(normal component) = coef * scale * u_n'' over the cells m cells from
// each of the four walls, leaving out n/8 cells at each end of every wall (the corners).
Real wall_fit (Array4<Real const> const& a, int comp0, int n, int m, Real scale)
{
    Real const h = Real(1)/Real(n);
    int const e = n/8;
    Real num = 0, den = 0;
    for (int t = e; t < n-e; ++t) {
        Real const s = (Real(t)+Real(0.5))*h;
        for (int side = 0; side < 2; ++side) {
            // x walls: normal component u, u_xx = 8 f''(x_wall) g'(y)
            int const i = (side == 0) ? m : n-1-m;
            Real const rx = scale * Real(8) * bf2(Real(side)) * bg1(s);
            num += a(i,t,0,comp0) * rx;
            den += rx * rx;
            // y walls: normal component v, v_yy = -8 f'(x) g''(y_wall)
            int const j = (side == 0) ? m : n-1-m;
            Real const ry = scale * Real(-8) * bf1(s) * bg2(Real(side));
            num += a(t,j,0,comp0+1) * ry;
            den += ry * ry;
        }
    }
    return num/den;
}

Result run (int n, int max_grid_size, Real rtol, int verbose)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), Real(nz)*h)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, 1)};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    // Exact velocity in the valid cells, zero in the ghost cells (incflo calls setBndry(0) on
    // the velocity before the projection)
    MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
    MultiFab lap_exact(ba, dm, AMREX_SPACEDIM, 0);
    vel.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const& v = vel.array(mfi);
        auto const& l = lap_exact.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = (Real(i)+Real(0.5))*h;
            Real const y = (Real(j)+Real(0.5))*h;
            v(i,j,k,0) =  Real(8) * bf(x) * bg1(y);
            v(i,j,k,1) = -Real(8) * bf1(x) * bg(y);
            l(i,j,k,0) =  Real(8) * (bf2(x) * bg1(y) + bf(x) * bg3(y));
            l(i,j,k,1) = -Real(8) * (bf3(x) * bg(y) + bf1(x) * bg2(y));
#if (AMREX_SPACEDIM == 3)
            l(i,j,k,2) = Real(0);
#endif
        });
    }
    vel.FillBoundary(geom.periodicity());
    MultiFab vel_exact(ba, dm, AMREX_SPACEDIM, 1);
    MultiFab::Copy(vel_exact, vel, 0, 0, AMREX_SPACEDIM, 1);

    // 1. Nodal projection, as in AMReX-Hydro's NodalProjector: vel += -sigma grad(phi)
    {
        Array<LinOpBCType,AMREX_SPACEDIM> const bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                                LinOpBCType::Neumann,
                                                                LinOpBCType::Periodic)};
        MLNodeLaplacian linop({geom}, {ba}, {dm});
        linop.setDomainBC(bc, bc);
        MultiFab sigma(ba, dm, 1, 0);
        sigma.setVal(1.0);
        linop.setSigma(0, sigma);
        linop.setWallDivuCorrection(wall_corr);

        BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
        MultiFab rhs(nba, dm, 1, 0);
        MultiFab phi(nba, dm, 1, 1);
        phi.setVal(0.0);
        linop.compRHS({&rhs}, {&vel}, {}, {});

        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        mlmg.solve({&phi}, {&rhs}, rtol, Real(0));

        // Residual of the solve at the (0,n) lid corner node, its two neighbours, and
        // over all nodes: a converged solve leaves them all at the tolerance.  A corner
        // node that the solve does not act on shows up as a residual equal to its rhs.
        {
            MultiFab res(nba, dm, 1, 0);
            mlmg.compResidual({&res}, {&phi}, {&rhs});
            Box const cbox(IntVect(AMREX_D_DECL(0, n-1, 0)), IntVect(AMREX_D_DECL(1, n, 0)));
            BoxArray const cba(cbox);
            DistributionMapping const cdm(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
            MultiFab cres(cba, cdm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
            MultiFab crhs(cba, cdm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
            cres.ParallelCopy(res, 0, 0, 1);
            crhs.ParallelCopy(rhs, 0, 0, 1);
            Gpu::streamSynchronize();
            Real const rmax = res.norminf();
            for (MFIter mfi(cres); mfi.isValid(); ++mfi) {
                auto const& a = cres.const_array(mfi);
                auto const& b = crhs.const_array(mfi);
                char line[256];
                std::snprintf(line, sizeof(line),
                    "  n = %4d  residual: max %9.2e  corner (0,n) %9.2e (rhs %9.2e)  (1,n) %9.2e (rhs %9.2e)  (0,n-1) %9.2e (rhs %9.2e)\n",
                    n, double(rmax), double(a(0,n,0)), double(b(0,n,0)), double(a(1,n,0)), double(b(1,n,0)),
                    double(a(0,n-1,0)), double(b(0,n-1,0)));
                amrex::Print() << line;
            }
        }

        MultiFab flux(ba, dm, AMREX_SPACEDIM, 0);
        mlmg.getFluxes({&flux});
        MultiFab::Add(vel, flux, 0, 0, AMREX_SPACEDIM, 0);
        vel.FillBoundary(geom.periodicity());
    }

    // 2. The viscous operators, as incflo sets them up: L_h u = div(grad u) per component
    //    (MLABecLaplacian), or div(grad u + grad u^T) - 2/3 grad(div u) (MLTensorOp, no bulk
    //    viscosity).  Both reduce to lap u for a divergence-free u.
    Array<LinOpBCType,AMREX_SPACEDIM> const vbc{AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                             LinOpBCType::Dirichlet,
                                                             LinOpBCType::Periodic)};
    LPInfo info;
    info.setMaxCoarseningLevel(0);

    // The ghost cells hold the wall values, which MLMG takes as the Dirichlet values on the
    // wall faces: zero, except u = 16 f(x) above the lid.  That includes the ghost cells past
    // the box edges along the lid, which MLTensorOp's cross terms use.
    auto set_wall_values = [&] (MultiFab& v)
    {
        v.setBndry(0.0);
        for (MFIter mfi(v); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            if (vbx.bigEnd(1) != n-1) { continue; }
            Box const lid = amrex::grow(amrex::adjCellHi(vbx, 1), IntVect(AMREX_D_DECL(1,0,1)));
            auto const& a = v.array(mfi);
            ParallelFor(lid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i >= 0 && i < n) {
                    a(i,j,k,0) = Real(16) * bf((Real(i)+Real(0.5))*h);
                }
            });
        }
        v.FillBoundary(geom.periodicity());
    };

    auto apply_op = [&] (int op, MultiFab& out, MultiFab const& in)
    {
        MultiFab vin(ba, dm, AMREX_SPACEDIM, 1);
        MultiFab::Copy(vin, in, 0, 0, AMREX_SPACEDIM, 0);
        set_wall_values(vin);

        int const nmo = int(maxorders.size());
        int const maxorder = maxorders[op % nmo];
        if (op < nmo) {
            MLABecLaplacian lap({geom}, {ba}, {dm}, info);
            lap.setMaxOrder(maxorder);
            lap.setDomainBC(vbc, vbc);
            lap.setScalars(0.0, -1.0);
            lap.setACoeffs(0, Real(1));
            lap.setBCoeffs(0, Real(1));
            MultiFab ci(ba, dm, 1, 1);
            MultiFab co(ba, dm, 1, 0);
            for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                MultiFab::Copy(ci, vin, comp, 0, 1, 1);
                lap.setLevelBC(0, &ci);
                MLMG mlmg(lap);
                mlmg.apply({&co}, {&ci});
                MultiFab::Copy(out, co, 0, comp, 1, 0);
            }
        } else {
            MLTensorOp tensor({geom}, {ba}, {dm}, info);
            tensor.setMaxOrder(maxorder);
            tensor.setDomainBC(Vector<Array<LinOpBCType,AMREX_SPACEDIM>>(AMREX_SPACEDIM, vbc),
                               Vector<Array<LinOpBCType,AMREX_SPACEDIM>>(AMREX_SPACEDIM, vbc));
            tensor.setScalars(0.0, -1.0);
            tensor.setACoeffs(0, Real(1));
            tensor.setShearViscosity(0, Real(1));
            tensor.setLevelBC(0, &vin);
            MLMG mlmg(tensor);
            mlmg.apply({&out}, {&vin});
        }
        MultiFab::Subtract(out, lap_exact, 0, 0, AMREX_SPACEDIM, 0);
    };

    // All errors in one MultiFab: the projection error, then L_h u - lap u for each operator
    // on the exact and on the projected velocity
    int const ncomp = AMREX_SPACEDIM * (1 + 2*nops());
    MultiFab err(ba, dm, ncomp, 0);
    MultiFab::LinComb(err, Real(1), vel, 0, Real(-1), vel_exact, 0, 0, AMREX_SPACEDIM, 0);
    for (int op = 0; op < nops(); ++op) {
        MultiFab tmp(ba, dm, AMREX_SPACEDIM, 0);
        apply_op(op, tmp, vel_exact);
        MultiFab::Copy(err, tmp, 0, AMREX_SPACEDIM*(1+2*op), AMREX_SPACEDIM, 0);
        apply_op(op, tmp, vel);
        MultiFab::Copy(err, tmp, 0, AMREX_SPACEDIM*(2+2*op), AMREX_SPACEDIM, 0);
    }

    Result r;
    r.c_exact.resize(nops(), Real(0));
    r.c_proj.resize(nops(), Real(0));
    r.t_exact.resize(nops(), Real(0));
    r.t_proj.resize(nops(), Real(0));
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        r.e_max = std::max(r.e_max, err.norminf(comp, 0));
    }

    // Gather onto one box on the I/O rank for the wall fits
    BoxArray const ba1(domain);
    DistributionMapping const dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab all(ba1, dm1, ncomp, 0, MFInfo().SetArena(The_Pinned_Arena()));
    all.ParallelCopy(err, 0, 0, ncomp);
    Gpu::streamSynchronize();
    for (MFIter mfi(all); mfi.isValid(); ++mfi) {
        auto const& a = all.const_array(mfi);
        for (int m = 0; m < 4; ++m) {
            r.e_delta[m] = wall_fit(a, 0, n, m, h*h/Real(8));
        }
        int q = 0;
        for (int cj : {0, n-1}) {
            for (int ci : {0, n-1}) {
                Real e = 0;
                for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
                    e = std::max(e, std::abs(a(ci,cj,0,comp)) / (h*h));
                }
                r.e_corner[q++] = e;
            }
        }
        for (int op = 0; op < nops(); ++op) {
            r.c_exact[op] = -wall_fit(a, AMREX_SPACEDIM*(1+2*op), n, 0, Real(1));
            r.c_proj[op]  = -wall_fit(a, AMREX_SPACEDIM*(2+2*op), n, 0, Real(1));
            r.t_exact[op] = wall_tangential_max(a, AMREX_SPACEDIM*(1+2*op), n);
            r.t_proj[op]  = wall_tangential_max(a, AMREX_SPACEDIM*(2+2*op), n);
        }
    }
    return r;
}

// Where a wall meets an inflow or an outflow face.  The channel flow above (chu, chv) on
// [0,1]^2, with no-slip walls at y = 0 and 1 (Neumann for phi), an inflow face at x = 0
// (LinOpBCType::inflow, the inflow velocity in the ghost cells) and a pressure outflow at x = 1
// (phi = 0, LinOpBCType::Dirichlet), as incflo sets up mass inflow and pressure outflow.  In 3D
// the velocity does not depend on z, w = 0, and z is periodic.
struct ChannelResult
{
    // max |rhs| / h of compRHS for the exact (divergence-free) velocity: O(1) where the
    // divergence is first order, O(h) where it is second order
    Real rhs_int = 0;        // interior nodes
    Real rhs_wall = 0;       // wall nodes, leaving out n/8 nodes at each end
    Real rhs_wall_in = 0;    // wall nodes next to the inflow corner (i = 1)
    Real rhs_corner_in = 0;  // the wall/inflow corner nodes (i = 0)
    Real rhs_inflow = 0;     // inflow face nodes, leaving out n/8 nodes at each end
    Real rhs_wall_out = 0;   // wall nodes next to the outflow face (i = n-1)
    // The projected normal velocity error / delta, delta = (h^2/8) u_n''
    Real e_max = 0;                    // max |u_projected - u| over all cells
    Real e_wall = 0;                   // wall cells (v), fit leaving out n/8 cells at each end
    std::array<Real,3> e_wall_in{};    // wall cells i = 0, 1, 2 (mean of both walls)
    std::array<Real,2> e_wall_out{};   // wall cells i = n-2, n-1 (mean of both walls)
    Real e_inflow = 0;                 // inflow face cells (u), fit leaving out n/8 at each end
};

ChannelResult channel_run (int n, int max_grid_size, Real rtol, int verbose)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), Real(nz)*h)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, 1)};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    // Exact velocity in the valid cells.  Ghost cells: zero at the walls, the inflow velocity
    // on the face at x = 0 (compRHS ignores its tangential part), and first-order extrapolation
    // at the outflow, as incflo fills them (the outflow nodes are Dirichlet, so compRHS does not
    // use those).
    MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
    vel.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        Box const& vbx = mfi.validbox();
        auto const& v = vel.array(mfi);
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = (Real(i)+Real(0.5))*h;
            Real const y = (Real(j)+Real(0.5))*h;
            v(i,j,k,0) = chu(x, y);
            v(i,j,k,1) = chv(x, y);
        });
        // Including the ghost cells at the box corners, along the whole face (j in [0,n)):
        // the nodes at a box edge on the face use them.  The outflow ghost cells get the value of
        // the cell next to them.
        if (vbx.smallEnd(0) == 0) {
            ParallelFor(amrex::grow(amrex::adjCellLo(vbx, 0), IntVect(AMREX_D_DECL(0,1,1))),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < 0 || j >= n) { return; }
                Real const y = (Real(j)+Real(0.5))*h;
                v(i,j,k,0) = chu(Real(0), y);
                v(i,j,k,1) = chv(Real(0), y);
            });
        }
        if (vbx.bigEnd(0) == n-1) {
            ParallelFor(amrex::grow(amrex::adjCellHi(vbx, 0), IntVect(AMREX_D_DECL(0,1,1))),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < 0 || j >= n) { return; }
                Real const x = (Real(n)-Real(0.5))*h;
                Real const y = (Real(j)+Real(0.5))*h;
                v(i,j,k,0) = chu(x, y);
                v(i,j,k,1) = chv(x, y);
            });
        }
    }
    vel.FillBoundary(geom.periodicity());
    MultiFab vel_exact(ba, dm, AMREX_SPACEDIM, 0);
    MultiFab::Copy(vel_exact, vel, 0, 0, AMREX_SPACEDIM, 0);

    ChannelResult r;
    int const e = n/8;
    BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
    DistributionMapping const dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    {
        Array<LinOpBCType,AMREX_SPACEDIM> const lobc{AMREX_D_DECL(LinOpBCType::inflow,
                                                                  LinOpBCType::Neumann,
                                                                  LinOpBCType::Periodic)};
        Array<LinOpBCType,AMREX_SPACEDIM> const hibc{AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                                  LinOpBCType::Neumann,
                                                                  LinOpBCType::Periodic)};
        MLNodeLaplacian linop({geom}, {ba}, {dm});
        linop.setDomainBC(lobc, hibc);
        MultiFab sigma(ba, dm, 1, 0);
        sigma.setVal(1.0);
        linop.setSigma(0, sigma);
        linop.setWallDivuCorrection(wall_corr);

        MultiFab rhs(nba, dm, 1, 0);
        MultiFab phi(nba, dm, 1, 1);
        phi.setVal(0.0);
        linop.compRHS({&rhs}, {&vel}, {}, {});

        // The divergence of the exact velocity, on one box on the I/O rank
        MultiFab rhs1(BoxArray(amrex::surroundingNodes(domain)), dm1, 1, 0,
                      MFInfo().SetArena(The_Pinned_Arena()));
        rhs1.ParallelCopy(rhs, 0, 0, 1);
        Gpu::streamSynchronize();
        for (MFIter mfi(rhs1); mfi.isValid(); ++mfi) {
            auto const& a = rhs1.const_array(mfi);
            auto q = [&] (int i, int j) { return std::abs(a(i,j,0)) / h; };
            for (int i = 1; i < n; ++i) {
                for (int j = 1; j < n; ++j) {
                    r.rhs_int = std::max(r.rhs_int, q(i,j));
                }
            }
            for (int jw : {0, n}) {
                for (int i = e; i <= n-e; ++i) {
                    r.rhs_wall = std::max(r.rhs_wall, q(i,jw));
                }
                r.rhs_wall_in   = std::max(r.rhs_wall_in, q(1,jw));
                r.rhs_corner_in = std::max(r.rhs_corner_in, q(0,jw));
                r.rhs_wall_out  = std::max(r.rhs_wall_out, q(n-1,jw));
            }
            for (int j = e; j <= n-e; ++j) {
                r.rhs_inflow = std::max(r.rhs_inflow, q(0,j));
            }
        }

        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        mlmg.solve({&phi}, {&rhs}, rtol, Real(0));

        MultiFab flux(ba, dm, AMREX_SPACEDIM, 0);
        mlmg.getFluxes({&flux});
        MultiFab::Add(vel, flux, 0, 0, AMREX_SPACEDIM, 0);
    }

    MultiFab err(ba, dm, AMREX_SPACEDIM, 0);
    MultiFab::LinComb(err, Real(1), vel, 0, Real(-1), vel_exact, 0, 0, AMREX_SPACEDIM, 0);
    for (int comp = 0; comp < AMREX_SPACEDIM; ++comp) {
        r.e_max = std::max(r.e_max, err.norminf(comp, 0));
    }
    MultiFab all(BoxArray(domain), dm1, AMREX_SPACEDIM, 0, MFInfo().SetArena(The_Pinned_Arena()));
    all.ParallelCopy(err, 0, 0, AMREX_SPACEDIM);
    Gpu::streamSynchronize();
    Real const d8 = h*h/Real(8);
    for (MFIter mfi(all); mfi.isValid(); ++mfi) {
        auto const& a = all.const_array(mfi);
        // Walls: the normal component v, delta = (h^2/8) v_yy, v_yy = -2 (1+x)
        auto dw = [&] (int i) { return d8 * Real(-2) * chF1((Real(i)+Real(0.5))*h); };
        Real num = 0, den = 0;
        for (int i = e; i < n-e; ++i) {
            for (int jw : {0, n-1}) {
                num += a(i,jw,0,1) * dw(i);
                den += dw(i) * dw(i);
            }
        }
        r.e_wall = num/den;
        for (int m = 0; m < 3; ++m) {
            r.e_wall_in[m] = Real(0.5) * (a(m,0,0,1) + a(m,n-1,0,1)) / dw(m);
        }
        for (int m = 0; m < 2; ++m) {
            int const i = n-2+m;
            r.e_wall_out[m] = Real(0.5) * (a(i,0,0,1) + a(i,n-1,0,1)) / dw(i);
        }
        // The inflow face: the normal component u, delta = (h^2/8) u_xx, u_xx = G'(y).  (The
        // outflow face is left out: incflo's viscous solve is Neumann there, so an offset of the
        // velocity next to it does not become an O(1) force.)
        num = 0;
        den = 0;
        for (int j = e; j < n-e; ++j) {
            Real const df = d8 * chG1((Real(j)+Real(0.5))*h);
            num += a(0,j,0,0) * df;
            den += df * df;
        }
        r.e_inflow = num/den;
    }
    return r;
}

#if (AMREX_SPACEDIM == 3)
// A closed box in 3D, walls on all six faces, so that it has cells on two walls (edges) and on
// three (corners), which the periodic-z cases above do not.  The velocity is u = curl(Phi c) =
// grad(Phi) x c with Phi = F(x) F(y) H(z), F(s) = 16 s^2 (1-s)^2, H(s) = 8 s^2 (1-s) and
// c = (1, 2, 3).  It is divergence-free and zero on five walls (F and F' vanish at 0 and 1, H
// and H' at 0); at z = 1, H = 0 but H' != 0, so that face is a lid moving tangentially with
// velocity proportional to F(x) F(y), like the Burggraf lid.  The normal velocity grows like
// the square of the distance from each wall.  Without the lid, u_n'' would vanish (linearly)
// along every edge, as no-slip on both walls forces; with it, u_n'' != 0 along the lid's
// edges (at the lid's corners it still vanishes, quadratically: the lid velocity has to go
// to zero at the side walls).  The Laplacian is grad(lap Phi) x c.  It is projected and the
// viscous operators applied as in run, with the lid velocity as their Dirichlet value at z = 1.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real box3d_fac (int d, Real s, int der)   // der-th derivative of F (d = 0, 1) or H (d = 2)
{
    if (d < 2) {
        if (der == 0) { return Real(16)*s*s*(Real(1)-s)*(Real(1)-s); }
        if (der == 1) { return Real(16)*(Real(2)*s - Real(6)*s*s + Real(4)*s*s*s); }
        if (der == 2) { return Real(16)*(Real(2) - Real(12)*s + Real(12)*s*s); }
        if (der == 3) { return Real(16)*(Real(-12) + Real(24)*s); }
        return Real(16)*Real(24);
    }
    if (der == 0) { return Real(8)*s*s*(Real(1)-s); }
    if (der == 1) { return Real(16)*s - Real(24)*s*s; }
    if (der == 2) { return Real(16) - Real(48)*s; }
    if (der == 3) { return Real(-48); }
    return Real(0);
}

// The derivative of Phi of order o[d] in each direction d, at p
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real box3d_phi (GpuArray<Real,3> const& p, GpuArray<int,3> const& o)
{
    return box3d_fac(0, p[0], o[0]) * box3d_fac(1, p[1], o[1]) * box3d_fac(2, p[2], o[2]);
}

// Component comp of g x c, c = (1, 2, 3), where g_a is the derivative of Phi of order
// o_extra + 1 in direction a (g = grad Phi for o_extra = 0)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real box3d_cross (GpuArray<Real,3> const& p, GpuArray<int,3> const& o_extra, int comp)
{
    GpuArray<Real,3> g;
    for (int a = 0; a < 3; ++a) {
        GpuArray<int,3> o = o_extra;
        ++o[a];
        g[a] = box3d_phi(p, o);
    }
    Real const c0 = 1, c1 = 2, c2 = 3;
    if (comp == 0) { return g[1]*c2 - g[2]*c1; }
    if (comp == 1) { return g[2]*c0 - g[0]*c2; }
    return g[0]*c1 - g[1]*c0;
}

// The velocity (lap = 0) or its Laplacian (lap = 1), component comp, at (x, y, z)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real box3d_vel (Real x, Real y, Real z, int comp, int lap)
{
    GpuArray<Real,3> const p{x, y, z};
    if (lap == 0) { return box3d_cross(p, GpuArray<int,3>{0, 0, 0}, comp); }
    Real r = 0;
    for (int e = 0; e < 3; ++e) {
        GpuArray<int,3> o{0, 0, 0};
        o[e] = 2;
        r += box3d_cross(p, o, comp);
    }
    return r;
}

// u_n'' at the wall normal to direction d, at the point p on it
Real box3d_unn (int d, GpuArray<Real,3> const& p)
{
    GpuArray<int,3> o{0, 0, 0};
    o[d] = 2;
    return box3d_cross(p, o, d);
}

// Least-squares fit of a(normal component) = coef * scale * u_n'' over the cells m cells from
// each of the six walls, leaving out n/8 cells next to the other walls (the edges).
Real box3d_fit (Array4<Real const> const& a, int comp0, int n, int m, Real scale)
{
    Real const h = Real(1)/Real(n);
    int const e = n/8;
    Real num = 0, den = 0;
    for (int d = 0; d < 3; ++d) {
        int const t1 = (d+1)%3, t2 = (d+2)%3;
        for (int side = 0; side < 2; ++side) {
            for (int b = e; b < n-e; ++b) {
            for (int c = e; c < n-e; ++c) {
                IntVect iv;
                iv[d] = (side == 0) ? m : n-1-m;
                iv[t1] = b;
                iv[t2] = c;
                GpuArray<Real,3> p;
                p[d] = Real(side);
                p[t1] = (Real(b)+Real(0.5))*h;
                p[t2] = (Real(c)+Real(0.5))*h;
                Real const r = scale * box3d_unn(d, p);
                num += a(iv,comp0+d) * r;
                den += r * r;
            }}
        }
    }
    return num/den;
}

struct Box3dResult
{
    Real res = 0;                   // max residual of the solve
    std::array<Real,4> rhs_h{};     // max |rhs| / h of the exact velocity: interior, face, edge
                                    // and corner nodes
    Real e_max = 0;                 // max |u_P - u|
    std::array<Real,2> e_delta{};   // wall-normal error / delta in cells 0 and 1, faces only
    std::array<Real,4> e_h2{};      // max |u_P - u| / h^2: face cells (normal, then tangential
                                    // components), edge cells, corner cells
    Vector<Real> c_exact, c_proj;   // c, faces only
    Vector<std::array<Real,6>> v;   // max |L_h u - lap u|, exact then projected u: tangential
                                    // components in the face cells, edge cells, corner cells
};

Box3dResult box3d_run (int n, int max_grid_size, Real rtol, int verbose)
{
    Real const h = Real(1)/Real(n);
    Box const domain(IntVect(0), IntVect(n-1));
    RealBox const rb({Real(0), Real(0), Real(0)}, {Real(1), Real(1), Real(1)});
    Array<int,AMREX_SPACEDIM> const is_periodic{0, 0, 0};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    // Exact velocity in the valid cells, zero in the ghost cells (the wall velocity)
    MultiFab vel(ba, dm, 3, 1);
    MultiFab lap_exact(ba, dm, 3, 0);
    vel.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        auto const& v = vel.array(mfi);
        auto const& l = lap_exact.array(mfi);
        ParallelFor(mfi.validbox(), 3, [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
        {
            Real const x = (Real(i)+Real(0.5))*h;
            Real const y = (Real(j)+Real(0.5))*h;
            Real const z = (Real(k)+Real(0.5))*h;
            v(i,j,k,c) = box3d_vel(x, y, z, c, 0);
            l(i,j,k,c) = box3d_vel(x, y, z, c, 1);
        });
    }
    vel.FillBoundary();
    MultiFab vel_exact(ba, dm, 3, 1);
    MultiFab::Copy(vel_exact, vel, 0, 0, 3, 1);

    Box3dResult r;
    {
        Array<LinOpBCType,AMREX_SPACEDIM> const bc{LinOpBCType::Neumann, LinOpBCType::Neumann,
                                                   LinOpBCType::Neumann};
        MLNodeLaplacian linop({geom}, {ba}, {dm});
        linop.setDomainBC(bc, bc);
        MultiFab sigma(ba, dm, 1, 0);
        sigma.setVal(1.0);
        linop.setSigma(0, sigma);
        linop.setWallDivuCorrection(wall_corr);
        BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
        MultiFab rhs(nba, dm, 1, 0);
        MultiFab phi(nba, dm, 1, 1);
        phi.setVal(0.0);
        linop.compRHS({&rhs}, {&vel}, {}, {});
        {
            // The exact velocity is divergence-free: |rhs| / h is O(1) where the divergence is
            // first order and O(h) where it is second order
            DistributionMapping const dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
            MultiFab rhs1(BoxArray(amrex::surroundingNodes(domain)), dm1, 1, 0,
                          MFInfo().SetArena(The_Pinned_Arena()));
            rhs1.ParallelCopy(rhs, 0, 0, 1);
            Gpu::streamSynchronize();
            for (MFIter mfi(rhs1); mfi.isValid(); ++mfi) {
                auto const& a = rhs1.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k)
                {
                    int const nw = int(i == 0 || i == n) + int(j == 0 || j == n)
                        + int(k == 0 || k == n);
                    r.rhs_h[nw] = std::max(r.rhs_h[nw], std::abs(a(i,j,k))/h);
                });
            }
        }
        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        mlmg.solve({&phi}, {&rhs}, rtol, Real(0));
        MultiFab res(nba, dm, 1, 0);
        mlmg.compResidual({&res}, {&phi}, {&rhs});
        r.res = res.norminf();
        MultiFab flux(ba, dm, 3, 0);
        mlmg.getFluxes({&flux});
        MultiFab::Add(vel, flux, 0, 0, 3, 0);
        vel.FillBoundary();
    }

    // The viscous operators, Dirichlet with zero wall values on every face
    Array<LinOpBCType,AMREX_SPACEDIM> const vbc{LinOpBCType::Dirichlet, LinOpBCType::Dirichlet,
                                                LinOpBCType::Dirichlet};
    LPInfo info;
    info.setMaxCoarseningLevel(0);
    auto apply_op = [&] (int op, MultiFab& out, MultiFab const& in)
    {
        MultiFab vin(ba, dm, 3, 1);
        MultiFab::Copy(vin, in, 0, 0, 3, 0);
        // The ghost cells hold the wall values: zero, except the lid velocity above z = 1,
        // including past the box edges along the lid (MLTensorOp's cross terms use them)
        vin.setBndry(0.0);
        for (MFIter mfi(vin); mfi.isValid(); ++mfi) {
            Box const& vbx = mfi.validbox();
            if (vbx.bigEnd(2) != n-1) { continue; }
            Box const lid = amrex::grow(amrex::adjCellHi(vbx, 2), IntVect(1,1,0));
            auto const& v = vin.array(mfi);
            ParallelFor(lid, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (i >= 0 && i < n && j >= 0 && j < n) {
                    Real const x = (Real(i)+Real(0.5))*h;
                    Real const y = (Real(j)+Real(0.5))*h;
                    v(i,j,k,0) = box3d_vel(x, y, Real(1), 0, 0);
                    v(i,j,k,1) = box3d_vel(x, y, Real(1), 1, 0);
                }
            });
        }
        vin.FillBoundary();
        int const nmo = int(maxorders.size());
        int const maxorder = maxorders[op % nmo];
        if (op < nmo) {
            MLABecLaplacian lap({geom}, {ba}, {dm}, info);
            lap.setMaxOrder(maxorder);
            lap.setDomainBC(vbc, vbc);
            lap.setScalars(0.0, -1.0);
            lap.setACoeffs(0, Real(1));
            lap.setBCoeffs(0, Real(1));
            MultiFab ci(ba, dm, 1, 1);
            MultiFab co(ba, dm, 1, 0);
            for (int comp = 0; comp < 3; ++comp) {
                MultiFab::Copy(ci, vin, comp, 0, 1, 1);
                lap.setLevelBC(0, &ci);
                MLMG mlmg(lap);
                mlmg.apply({&co}, {&ci});
                MultiFab::Copy(out, co, 0, comp, 1, 0);
            }
        } else {
            MLTensorOp tensor({geom}, {ba}, {dm}, info);
            tensor.setMaxOrder(maxorder);
            tensor.setDomainBC(Vector<Array<LinOpBCType,AMREX_SPACEDIM>>(3, vbc),
                               Vector<Array<LinOpBCType,AMREX_SPACEDIM>>(3, vbc));
            tensor.setScalars(0.0, -1.0);
            tensor.setACoeffs(0, Real(1));
            tensor.setShearViscosity(0, Real(1));
            tensor.setLevelBC(0, &vin);
            MLMG mlmg(tensor);
            mlmg.apply({&out}, {&vin});
        }
        MultiFab::Subtract(out, lap_exact, 0, 0, 3, 0);
    };

    int const ncomp = 3 * (1 + 2*nops());
    MultiFab err(ba, dm, ncomp, 0);
    MultiFab::LinComb(err, Real(1), vel, 0, Real(-1), vel_exact, 0, 0, 3, 0);
    for (int op = 0; op < nops(); ++op) {
        MultiFab tmp(ba, dm, 3, 0);
        apply_op(op, tmp, vel_exact);
        MultiFab::Copy(err, tmp, 0, 3*(1+2*op), 3, 0);
        apply_op(op, tmp, vel);
        MultiFab::Copy(err, tmp, 0, 3*(2+2*op), 3, 0);
    }
    for (int comp = 0; comp < 3; ++comp) {
        r.e_max = std::max(r.e_max, err.norminf(comp, 0));
    }

    // Gather onto one box on the I/O rank
    DistributionMapping const dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab all(BoxArray(domain), dm1, ncomp, 0, MFInfo().SetArena(The_Pinned_Arena()));
    all.ParallelCopy(err, 0, 0, ncomp);
    Gpu::streamSynchronize();
    r.c_exact.resize(nops(), Real(0));
    r.c_proj.resize(nops(), Real(0));
    r.v.resize(nops(), std::array<Real,6>{});
    for (MFIter mfi(all); mfi.isValid(); ++mfi) {
        auto const& a = all.const_array(mfi);
        for (int m = 0; m < 2; ++m) {
            r.e_delta[m] = box3d_fit(a, 0, n, m, h*h/Real(8));
        }
        for (int op = 0; op < nops(); ++op) {
            r.c_exact[op] = -box3d_fit(a, 3*(1+2*op), n, 0, Real(1));
            r.c_proj[op]  = -box3d_fit(a, 3*(2+2*op), n, 0, Real(1));
        }
        // Face cells are those n/8 or more from the other walls, as in the fits
        int const e = n/8;
        amrex::LoopOnCpu(domain, [&] (int i, int j, int k)
        {
            IntVect const iv(i,j,k);
            int nw = 0, dn = 0;
            for (int d = 0; d < 3; ++d) {
                if (iv[d] == 0 || iv[d] == n-1) { ++nw; dn = d; }
            }
            if (nw == 0) { return; }
            if (nw == 1) {
                for (int d = 0; d < 3; ++d) {
                    if (d != dn && (iv[d] < e || iv[d] >= n-e)) { return; }
                }
            }
            for (int c = 0; c < 3; ++c) {
                int const q = (nw == 1) ? ((c == dn) ? 0 : 1) : nw;
                r.e_h2[q] = std::max(r.e_h2[q], std::abs(a(i,j,k,c))/(h*h));
                if (q == 0) { continue; }
                for (int op = 0; op < nops(); ++op) {
                    r.v[op][2*(q-1)]   = std::max(r.v[op][2*(q-1)],   std::abs(a(i,j,k,3*(1+2*op)+c)));
                    r.v[op][2*(q-1)+1] = std::max(r.v[op][2*(q-1)+1], std::abs(a(i,j,k,3*(2+2*op)+c)));
                }
            }
        });
    }
    return r;
}
#endif

// A uniform inflow against no-slip walls: the channel of channel_run, with u = 1 in the inflow
// ghost cells up to the walls, where u = 0.  The data do not match at the wall/inflow corners,
// and the velocity is singular there.  The velocity projected is
//   u = (2/pi)^2 atan2(y, x) atan2(1-y, x),  v = 0,
// which is 1 on the inflow face, 0 on the walls, and changes by O(1) from cell to cell next to
// the corners at every resolution.  The projected velocity carries the inflow flux, 1, through
// every cross-section of the channel.  Returns, at i = n/4, n/2, 3n/4, minus 1:
//   [0..2] the flux h sum_j u(i,j) of the projected velocity, which the volume conservation of
//          the uncorrected stencil makes exact;
//   [3..5] the flux that the volume conservation of the corrected stencil makes exact.  The correction takes the flux through the
//          half-cell faces next to a wall as (5/4) u_1 - (1/4) u_2 instead of u_1 (the first two
//          cells from the wall), which adds (h/2) (u_1 - u_2)/4 at each wall.  It is applied to
//          the right-hand side, i.e. to the velocity before the projection, so these terms are
//          those of that velocity.  Without the correction this is the same as [0..2].
// A node that does not conserve volume adds its error to both, the same at every x downstream.
std::array<Real,6> uniform_inflow_run (int n, int max_grid_size, Real rtol, int verbose)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), Real(nz)*h)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, 1)};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    constexpr Real two_over_pi = Real(2.0)/Real(3.1415926535897932);
    MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
    vel.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        Box const& vbx = mfi.validbox();
        auto const& v = vel.array(mfi);
        ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = (Real(i)+Real(0.5))*h;
            Real const y = (Real(j)+Real(0.5))*h;
            v(i,j,k,0) = two_over_pi*two_over_pi * std::atan2(y, x) * std::atan2(Real(1)-y, x);
        });
        if (vbx.smallEnd(0) == 0) {
            ParallelFor(amrex::grow(amrex::adjCellLo(vbx, 0), IntVect(AMREX_D_DECL(0,1,1))),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j >= 0 && j < n) { v(i,j,k,0) = Real(1); }
            });
        }
        if (vbx.bigEnd(0) == n-1) {
            ParallelFor(amrex::grow(amrex::adjCellHi(vbx, 0), IntVect(AMREX_D_DECL(0,1,1))),
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j >= 0 && j < n) { v(i,j,k,0) = v(n-1,j,k,0); }
            });
        }
    }
    vel.FillBoundary(geom.periodicity());

    Array<LinOpBCType,AMREX_SPACEDIM> const lobc{AMREX_D_DECL(LinOpBCType::inflow,
                                                              LinOpBCType::Neumann,
                                                              LinOpBCType::Periodic)};
    Array<LinOpBCType,AMREX_SPACEDIM> const hibc{AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                              LinOpBCType::Neumann,
                                                              LinOpBCType::Periodic)};
    MLNodeLaplacian linop({geom}, {ba}, {dm});
    linop.setDomainBC(lobc, hibc);
    MultiFab sigma(ba, dm, 1, 0);
    sigma.setVal(1.0);
    linop.setSigma(0, sigma);
    linop.setWallDivuCorrection(wall_corr);
    BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
    MultiFab rhs(nba, dm, 1, 0);
    MultiFab phi(nba, dm, 1, 1);
    phi.setVal(0.0);
    MultiFab vel_in(ba, dm, 1, 0);
    MultiFab::Copy(vel_in, vel, 0, 0, 1, 0);
    linop.compRHS({&rhs}, {&vel}, {}, {});
    MLMG mlmg(linop);
    mlmg.setVerbose(verbose);
    mlmg.solve({&phi}, {&rhs}, rtol, Real(0));
    MultiFab flux(ba, dm, AMREX_SPACEDIM, 0);
    mlmg.getFluxes({&flux});
    MultiFab::Add(vel, flux, 0, 0, AMREX_SPACEDIM, 0);

    // u after the projection (component 0) and before it (component 1), on one box
    DistributionMapping const dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab all(BoxArray(domain), dm1, 2, 0, MFInfo().SetArena(The_Pinned_Arena()));
    all.ParallelCopy(vel, 0, 0, 1);
    all.ParallelCopy(vel_in, 0, 1, 1);
    Gpu::streamSynchronize();
    std::array<Real,6> q{};
    for (MFIter mfi(all); mfi.isValid(); ++mfi) {
        auto const& a = all.const_array(mfi);
        int m = 0;
        for (int i : {n/4, n/2, 3*n/4}) {
            Real s = 0;
            for (int j = 0; j < n; ++j) { s += a(i,j,0); }
            q[m] = s*h - Real(1);
            Real const walls = wall_corr
                ? Real(0.125)*h*((a(i,0,0,1) - a(i,1,0,1)) + (a(i,n-1,0,1) - a(i,n-2,0,1)))
                : Real(0);
            q[m+3] = q[m] + walls;
            ++m;
        }
    }
    return q;
}

}

void run_study (Vector<int> const& n_cell, int max_grid_size, Real rtol, int verbose)
{
    ParmParse pp;
    pp.queryarr("maxorders", maxorders);
    char line[256];
    int burggraf = 1;
    pp.query("burggraf", burggraf);
    Vector<Result> results;
    if (burggraf != 0) {
        for (int n : n_cell) {
            results.push_back(run(n, max_grid_size, rtol, verbose));
        }

        amrex::Print() << "\nProjected Burggraf velocity: max error, and the wall-normal error in the\n"
                       << "cells 0-3 from the wall divided by delta = (h^2/8) u_n''\n"
                       << "     n   max |u_P - u|   cell 0   cell 1   cell 2   cell 3   corner cells / h^2: (0,0) (n,0) (0,n) (n,n)\n";
        for (int s = 0; s < n_cell.size(); ++s) {
            auto const& r = results[s];
            std::snprintf(line, sizeof(line), "  %4d   %13.3e   %6.3f   %6.3f   %6.3f   %6.3f   %8.4f %8.4f %8.4f %8.4f\n",
                          n_cell[s], double(r.e_max), double(r.e_delta[0]), double(r.e_delta[1]),
                          double(r.e_delta[2]), double(r.e_delta[3]), double(r.e_corner[0]),
                          double(r.e_corner[1]), double(r.e_corner[2]), double(r.e_corner[3]));
            amrex::Print() << line;
        }

        amrex::Print() << "\nViscous operator in the wall cells: (L_h u - lap u)_n = -c u_n'' for the normal\n"
                       << "component, and max |L_h u - lap u| of the tangential component\n"
                       << "                                       normal: c                 tangential: max\n"
                       << "  operator                         n   exact u   projected u     exact u   projected u\n";
        for (int op = 0; op < nops(); ++op) {
            for (int s = 0; s < n_cell.size(); ++s) {
                auto const& r = results[s];
                std::snprintf(line, sizeof(line), "  %-28s  %4d   %7.4f   %11.4f   %9.4f   %11.4f\n",
                              (s == 0) ? op_name(op).c_str() : "", n_cell[s],
                              double(r.c_exact[op]), double(r.c_proj[op]),
                              double(r.t_exact[op]), double(r.t_proj[op]));
                amrex::Print() << line;
            }
        }
        amrex::Print() << '\n';
    }

    int channel = 0;
    pp.query("channel", channel);
    if (channel != 0) {
        Vector<ChannelResult> cres;
        for (int n : n_cell) {
            cres.push_back(channel_run(n, max_grid_size, rtol, verbose));
        }
        amrex::Print() << "Channel: no-slip walls at y = 0 and 1, inflow at x = 0, pressure outflow\n"
                       << "(phi = 0) at x = 1.  Nodal divergence of the exact (divergence-free) velocity,\n"
                       << "max |rhs| / h: O(1) where it is first order, O(h) where it is second order\n"
                       << "     n   interior       wall   wall i=1   corner i=0   inflow face   wall i=n-1\n";
        for (int s = 0; s < n_cell.size(); ++s) {
            auto const& c = cres[s];
            std::snprintf(line, sizeof(line), "  %4d   %8.4f   %8.4f   %8.4f     %8.4f      %8.4f     %8.4f\n",
                          n_cell[s], double(c.rhs_int), double(c.rhs_wall), double(c.rhs_wall_in),
                          double(c.rhs_corner_in), double(c.rhs_inflow), double(c.rhs_wall_out));
            amrex::Print() << line;
        }
        amrex::Print() << "\nChannel, projected velocity: max error, and the normal error divided by\n"
                       << "delta = (h^2/8) u_n'' in the wall cells (v; fit away from the ends, and at\n"
                       << "i = 0, 1, 2 next to the inflow and i = n-2, n-1 next to the outflow) and in\n"
                       << "the cells next to the inflow face (u; fit away from the walls)\n"
                       << "     n   max |u_P - u|    walls    i=0      i=1      i=2    i=n-2    i=n-1   inflow\n";
        for (int s = 0; s < n_cell.size(); ++s) {
            auto const& c = cres[s];
            std::snprintf(line, sizeof(line),
                          "  %4d   %13.3e   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f   %6.3f\n",
                          n_cell[s], double(c.e_max), double(c.e_wall), double(c.e_wall_in[0]),
                          double(c.e_wall_in[1]), double(c.e_wall_in[2]), double(c.e_wall_out[0]),
                          double(c.e_wall_out[1]), double(c.e_inflow));
            amrex::Print() << line;
        }
        amrex::Print() << "\nUniform inflow (u = 1 in the inflow ghost cells up to the no-slip walls):\n"
                       << "flux of the projected velocity through x = 1/4, 1/2, 3/4, minus the inflow\n"
                       << "flux 1: h sum_j u (the flux the uncorrected stencil carries exactly), and the\n"
                       << "flux the corrected stencil carries exactly (the same without the correction)\n"
                       << "                   h sum_j u                         carried by the stencil\n"
                       << "     n     x = 1/4    x = 1/2    x = 3/4         x = 1/4    x = 1/2    x = 3/4\n";
        for (int n : n_cell) {
            auto const q = uniform_inflow_run(n, max_grid_size, rtol, verbose);
            std::snprintf(line, sizeof(line),
                          "  %4d  %10.3e %10.3e %10.3e      %10.3e %10.3e %10.3e\n", n,
                          double(q[0]), double(q[1]), double(q[2]),
                          double(q[3]), double(q[4]), double(q[5]));
            amrex::Print() << line;
        }
        amrex::Print() << '\n';
    }
    // The closed box is 3D only; the inputs are read in 2D too, so that they are not
    // reported as unused
    int box3d = 0;
    pp.query("box3d", box3d);
    Vector<int> b3n;
    if (!pp.queryarr("box3d_n_cell", b3n)) {
        b3n = {16, 32, 64};
    }
#if (AMREX_SPACEDIM == 3)
    if (box3d != 0) {
        Vector<Box3dResult> b3;
        for (int n : b3n) {
            b3.push_back(box3d_run(n, max_grid_size, rtol, verbose));
        }
        amrex::Print() << "Closed box, walls on all six faces: u = grad(Phi) x (1,2,3), Phi = F(x)F(y)H(z),\n"
                       << "F = 16 s^2 (1-s)^2, H = 8 s^2 (1-s): no-slip, with a lid at z = 1 moving\n"
                       << "tangentially.  Projected velocity: residual of the solve, max error, the\n"
                       << "wall-normal error / delta = (h^2/8) u_n'' in cells 0 and 1 of the faces, and\n"
                       << "max |u_P - u| / h^2 in the face cells (normal and tangential components), the\n"
                       << "edge cells (on two walls) and the corner cells (on three).  Face cells are n/8\n"
                       << "or more from the other walls.  Bounded / h^2 is second order\n"
                       << "     n    residual   max |u_P - u|   cell 0   cell 1    face n.   face t.      edge    corner\n";
        for (int s = 0; s < b3n.size(); ++s) {
            auto const& q = b3[s];
            std::snprintf(line, sizeof(line),
                          "  %4d   %9.2e   %13.3e   %6.3f   %6.3f   %8.3f  %8.3f  %8.3f  %8.3f\n",
                          b3n[s], double(q.res), double(q.e_max), double(q.e_delta[0]),
                          double(q.e_delta[1]), double(q.e_h2[0]), double(q.e_h2[1]),
                          double(q.e_h2[2]), double(q.e_h2[3]));
            amrex::Print() << line;
        }
        amrex::Print() << "\nClosed box, nodal divergence of the exact (divergence-free) velocity, max |rhs| / h:\n"
                       << "O(1) where it is first order, O(h) where it is second order\n"
                       << "     n    interior        face        edge      corner\n";
        for (int s = 0; s < b3n.size(); ++s) {
            auto const& q = b3[s];
            std::snprintf(line, sizeof(line), "  %4d   %9.4f   %9.4f   %9.4f   %9.4f\n", b3n[s],
                          double(q.rhs_h[0]), double(q.rhs_h[1]), double(q.rhs_h[2]),
                          double(q.rhs_h[3]));
            amrex::Print() << line;
        }
        amrex::Print() << "\nClosed box, viscous operator: c in (L_h u - lap u)_n = -c u_n'' in the face cells\n"
                       << "(normal component), and max |L_h u - lap u| of the tangential components in the\n"
                       << "face cells, and of all components in the edge and the corner cells, for the\n"
                       << "exact (ex.) and the projected (pr.) velocity\n"
                       << "  operator                         n  c, ex.  c, pr.   face t., ex.  pr.          edge, ex.  pr.           corner, ex. pr.\n";
        for (int op = 0; op < nops(); ++op) {
            for (int s = 0; s < b3n.size(); ++s) {
                auto const& q = b3[s];
                std::snprintf(line, sizeof(line),
                              "  %-28s  %4d  %6.3f  %6.3f   %9.2e %9.2e   %9.2e %9.2e   %9.2e %9.2e\n",
                              (s == 0) ? op_name(op).c_str() : "", b3n[s],
                              double(q.c_exact[op]), double(q.c_proj[op]),
                              double(q.v[op][0]), double(q.v[op][1]), double(q.v[op][2]),
                              double(q.v[op][3]), double(q.v[op][4]), double(q.v[op][5]));
                amrex::Print() << line;
            }
        }
        amrex::Print() << '\n';
    }
#endif
}
