// Reproducer: after a nodal projection, a viscous operator applied to the velocity has an O(1)
// error in the wall-normal component in the cells next to a no-slip wall.  See README.md.
//
// The velocity is the Burggraf flow on [0,1]^2: divergence-free, zero on the walls except for
// the tangential velocity of the lid at y = 1, u = 16 f(x).  In 3D it does not depend on z,
// w = 0, and z is periodic.  The program
//
//   1. projects it once with MLNodeLaplacian, the way incflo calls AMReX-Hydro's
//      NodalProjector (velocity ghost cells set to 0, sigma = 1), and
//   2. applies the viscous operators incflo uses (MLABecLaplacian per component, and
//      MLTensorOp), with Dirichlet walls and maxorder 2 or 3, to the exact and to the
//      projected velocity.
//
// It prints, for each resolution, the projected wall-normal velocity error divided by
// delta = (h^2/8) u_n'', and the c in (L_h u - lap u)_n = -c u_n'' in the wall cells, where
// u_n'' is the second normal derivative of the normal velocity at the wall.

#include <AMReX.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <array>
#include <cstdio>
#include <string>

using namespace amrex;

namespace {

// Burggraf flow: u = 8 f(x) g'(y), v = -8 f'(x) g(y), f = x^4 - 2x^3 + x^2, g = y^4 - y^2
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bf  (Real x) { return x*x*x*x - Real(2)*x*x*x + x*x; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bf1 (Real x) { return Real(4)*x*x*x - Real(6)*x*x + Real(2)*x; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bf2 (Real x) { return Real(12)*x*x - Real(12)*x + Real(2); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bf3 (Real x) { return Real(24)*x - Real(12); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bg  (Real y) { return y*y*y*y - y*y; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bg1 (Real y) { return Real(4)*y*y*y - Real(2)*y; }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bg2 (Real y) { return Real(12)*y*y - Real(2); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real bg3 (Real y) { return Real(24)*y; }

constexpr int nops = 4;
const std::array<std::string,nops> op_names{"MLABecLaplacian, maxorder 2",
                                            "MLABecLaplacian, maxorder 3",
                                            "MLTensorOp,      maxorder 2",
                                            "MLTensorOp,      maxorder 3"};

struct Result
{
    Real e_max = 0;                       // max |u_projected - u| over all cells
    std::array<Real,4> e_delta{};         // wall-normal error / delta, cells 0-3 from the wall
    std::array<Real,nops> c_exact{};      // c for L_h applied to the exact velocity
    std::array<Real,nops> c_proj{};       // c for L_h applied to the projected velocity
};

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

        BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
        MultiFab rhs(nba, dm, 1, 0);
        MultiFab phi(nba, dm, 1, 1);
        phi.setVal(0.0);
        linop.compRHS({&rhs}, {&vel}, {}, {});

        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        mlmg.solve({&phi}, {&rhs}, rtol, Real(0));

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

        int const maxorder = (op % 2 == 0) ? 2 : 3;
        if (op < 2) {
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
    int const ncomp = AMREX_SPACEDIM * (1 + 2*nops);
    MultiFab err(ba, dm, ncomp, 0);
    MultiFab::LinComb(err, Real(1), vel, 0, Real(-1), vel_exact, 0, 0, AMREX_SPACEDIM, 0);
    for (int op = 0; op < nops; ++op) {
        MultiFab tmp(ba, dm, AMREX_SPACEDIM, 0);
        apply_op(op, tmp, vel_exact);
        MultiFab::Copy(err, tmp, 0, AMREX_SPACEDIM*(1+2*op), AMREX_SPACEDIM, 0);
        apply_op(op, tmp, vel);
        MultiFab::Copy(err, tmp, 0, AMREX_SPACEDIM*(2+2*op), AMREX_SPACEDIM, 0);
    }

    Result r;
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
        for (int op = 0; op < nops; ++op) {
            r.c_exact[op] = -wall_fit(a, AMREX_SPACEDIM*(1+2*op), n, 0, Real(1));
            r.c_proj[op]  = -wall_fit(a, AMREX_SPACEDIM*(2+2*op), n, 0, Real(1));
        }
    }
    return r;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        // queryarr does not shrink the vector, so the defaults are set only if n_cell is absent
        Vector<int> n_cell;
        if (!pp.queryarr("n_cell", n_cell)) {
            n_cell = {32, 64, 128, 256};
        }
        int max_grid_size = 32;
        pp.query("max_grid_size", max_grid_size);
        Real rtol = 1.e-12;
        pp.query("rtol", rtol);
        int verbose = 0;
        pp.query("verbose", verbose);

        Vector<Result> results;
        for (int n : n_cell) {
            results.push_back(run(n, max_grid_size, rtol, verbose));
        }

        char line[256];
        amrex::Print() << "\nProjected Burggraf velocity: max error, and the wall-normal error in the\n"
                       << "cells 0-3 from the wall divided by delta = (h^2/8) u_n''\n"
                       << "     n   max |u_P - u|   cell 0   cell 1   cell 2   cell 3\n";
        for (int s = 0; s < n_cell.size(); ++s) {
            auto const& r = results[s];
            std::snprintf(line, sizeof(line), "  %4d   %13.3e   %6.3f   %6.3f   %6.3f   %6.3f\n",
                          n_cell[s], double(r.e_max), double(r.e_delta[0]), double(r.e_delta[1]),
                          double(r.e_delta[2]), double(r.e_delta[3]));
            amrex::Print() << line;
        }

        amrex::Print() << "\nViscous operator in the wall cells: (L_h u - lap u)_n = -c u_n''\n"
                       << "  operator                         n   c, exact u   c, projected u\n";
        for (int op = 0; op < nops; ++op) {
            for (int s = 0; s < n_cell.size(); ++s) {
                auto const& r = results[s];
                std::snprintf(line, sizeof(line), "  %-28s  %4d   %10.4f   %14.4f\n",
                              (s == 0) ? op_names[op].c_str() : "", n_cell[s],
                              double(r.c_exact[op]), double(r.c_proj[op]));
                amrex::Print() << line;
            }
        }
        amrex::Print() << '\n';
    }
    amrex::Finalize();
}
