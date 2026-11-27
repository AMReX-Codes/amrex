//
// Checks on MLEBNodeFDLaplacian:
//
//   * reusing an operator must not lose the EB Dirichlet values supplied by
//     the callable setEBDirichlet;
//   * solving del dot (sigma grad phi) = rhs with the AlgMG and hypre bottom
//     solvers must reproduce the native bottom solver.  That solve does not
//     coarsen, so the bottom solve does all of the work, which is the sharpest
//     test of the matrix assembled by MLEBNodeFDLaplacian.  The hypre part is
//     skipped when hypre is not available;
//   * multigrid must stop coarsening before a level that can no longer see
//     part of the EB, even when another part of the EB is still visible;
//   * the Helmholtz term of setACoeffs, in Cartesian and (2D) RZ coordinates,
//     including a new coefficient on an existing operator.
//

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

using namespace amrex;

// A solve that never calls setEBDirichlet makes prepareForSolve pick
// homogeneous Dirichlet on the EB.  A setEBDirichlet call after that must
// still be honored.  postSolve stores the EB value in the covered nodes, so
// those nodes tell us which value was used.
void test_eb_dirichlet_reuse (Geometry const& geom, BoxArray const& grids,
                              DistributionMapping const& dmap,
                              EBFArrayBoxFactory const& factory,
                              Array<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                              Array<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                              int max_coarsening_level, Real reltol, int verbose)
{
    BoxArray const& nba = amrex::convert(grids, IntVect(1));

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);
    MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info, {&factory});
    linop.setDomainBC(lobc, hibc);
    linop.setSigma({AMREX_D_DECL(Real(1.0), Real(1.0), Real(1.0))});

    MultiFab rhs(nba, dmap, 1, 0);
    rhs.setVal(Real(1.0));
    MultiFab sol(nba, dmap, 1, 1);

    auto do_solve = [&] () {
        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        MultiFab rhs_copy(nba, dmap, 1, 0);
        MultiFab::Copy(rhs_copy, rhs, 0, 0, 1, 0);
        sol.setVal(0.0);
        mlmg.solve({&sol}, {&rhs_copy}, reltol, Real(0.0));
    };

    do_solve(); // no setEBDirichlet yet

    Real const phi_eb = 5.0;
    linop.setEBDirichlet([=] AMREX_GPU_HOST_DEVICE (AMREX_D_DECL(Real,Real,Real)) -> Real
                         { return phi_eb; });
    do_solve();

    auto const& levset = factory.getLevelSet();
    ReduceOps<ReduceOpMax, ReduceOpSum> reduce_op;
    ReduceData<Real, Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MFIter mfi(sol); mfi.isValid(); ++mfi) {
        Array4<Real const> const& s = sol.const_array(mfi);
        Array4<Real const> const& ls = levset.const_array(mfi);
        reduce_op.eval(mfi.validbox(), reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            bool const covered = ls(i,j,k) >= Real(0.0);
            return { covered ? std::abs(s(i,j,k)-phi_eb) : Real(0.0),
                     covered ? Long(1) : Long(0) };
        });
    }
    auto const hv = reduce_data.value(reduce_op);
    Real maxdiff = amrex::get<0>(hv);
    Long ncovered = amrex::get<1>(hv);
    ParallelDescriptor::ReduceRealMax(maxdiff);
    ParallelDescriptor::ReduceLongSum(ncovered);

    amrex::Print() << "covered nodes      = " << ncovered << "\n"
                   << "max EB value error = " << maxdiff << '\n';

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ncovered > 0,
        "No covered nodes: the EB Dirichlet value is not being tested");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(maxdiff <= Real(1.e-12),
        "setEBDirichlet was ignored after a solve without it");
}

// The AlgMG and hypre bottom solvers must reproduce the native one.
void test_bottom_solvers (Geometry const& geom, BoxArray const& grids,
                           DistributionMapping const& dmap,
                           EBFArrayBoxFactory const& factory,
                           Array<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                           Array<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                           Real reltol, int verbose)
{
    int bottom_verbose = 0;
    int use_sigma_mf = 1;
    int plot = 0;
    Real phi_eb = 1.0;
    // Single precision cannot reach the double precision tolerances.
    Real bottom_reltol = std::is_same_v<Real,float> ? Real(1.e-6) : Real(1.e-9);
    Real max_rel_diff = std::is_same_v<Real,float> ? Real(1.e-3) : Real(1.e-12);
    {
        ParmParse pp;
        pp.query("bottom_verbose", bottom_verbose);
        pp.query("use_sigma_mf", use_sigma_mf);
        pp.query("phi_eb", phi_eb);
        pp.query("bottom_reltol", bottom_reltol);
        pp.query("plot", plot);
        pp.query("max_rel_diff", max_rel_diff);
    }
    // Agreement to max_rel_diff needs a residual well below the default
    // reltol, which a solve may just barely meet.
    if constexpr (std::is_same_v<Real,double>) {
        reltol = std::min(reltol, Real(1.e-13));
    }

    BoxArray const& nba = amrex::convert(grids, IntVect(1));

    // A smooth, strictly positive sigma so that the variable coefficient
    // path is exercised.
    MultiFab sigma(grids, dmap, 1, 1, MFInfo(), factory);
    // Use whole wavenumbers so that sigma and the source stay single valued
    // when a direction is periodic.
    Real const twopi = Real(2.0)*Math::pi<Real>();
    auto const dx = geom.CellSizeArray();
    auto const problo = geom.ProbLoArray();
    for (MFIter mfi(sigma); mfi.isValid(); ++mfi) {
        Array4<Real> const& s = sigma.array(mfi);
        amrex::ParallelFor(mfi.growntilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            AMREX_D_TERM(Real const x = problo[0] + (i+Real(0.5))*dx[0];,
                         Real const y = problo[1] + (j+Real(0.5))*dx[1];,
                         Real const z = problo[2] + (k+Real(0.5))*dx[2];)
            s(i,j,k) = Real(1.0) + Real(0.5)*std::sin(twopi*x)
                * AMREX_D_TERM(Real(1.0), *std::cos(Real(2.0)*twopi*y), *std::sin(twopi*z));
        });
    }
    sigma.FillBoundary(geom.periodicity());

    MultiFab rhs(nba, dmap, 1, 0);
    for (MFIter mfi(rhs); mfi.isValid(); ++mfi) {
        Array4<Real> const& r = rhs.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            AMREX_D_TERM(Real const x = problo[0] + i*dx[0];,
                         Real const y = problo[1] + j*dx[1];,
                         Real const z = problo[2] + k*dx[2];)
            r(i,j,k) = std::sin(twopi*x)
                * AMREX_D_TERM(Real(1.0), *std::cos(Real(3.0)*twopi*y), *std::cos(Real(2.0)*twopi*z));
        });
    }

    auto do_solve = [&] (BottomSolver bottom_solver, MultiFab& sol)
    {
        // No coarsening, so that the bottom solver does all of the work.
        LPInfo info;
        info.setMaxCoarseningLevel(0);

        MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info,
                                  {&factory});
        linop.setDomainBC(lobc, hibc);
        linop.setEBDirichlet(phi_eb);
        if (use_sigma_mf) {
            linop.setSigma(0, sigma);
        } else {
            linop.setSigma({AMREX_D_DECL(Real(1.0), Real(1.5), Real(0.7))});
        }

        MLMG mlmg(linop);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);
        mlmg.setBottomSolver(bottom_solver);
        mlmg.setBottomTolerance(bottom_reltol);
        mlmg.setBottomMaxIter(1000);

        MultiFab rhs_copy(nba, dmap, 1, 0);
        MultiFab::Copy(rhs_copy, rhs, 0, 0, 1, 0);

        sol.define(nba, dmap, 1, 1);
        sol.setVal(0.0);
        return mlmg.solve({&sol}, {&rhs_copy}, reltol, Real(0.0));
    };

    MultiFab sol_native;

    amrex::Print() << "\n==== native bottom solver ====\n";
    Real const err_native = do_solve(BottomSolver::bicgstab, sol_native);
    Real const smax = sol_native.norminf(0, 0);

    auto compare = [&] (std::string const& name, BottomSolver bottom_solver)
    {
        amrex::Print() << "\n==== " << name << " bottom solver ====\n";
        MultiFab sol;
        Real const err = do_solve(bottom_solver, sol);

        MultiFab diff(nba, dmap, 1, 0);
        MultiFab::Copy(diff, sol, 0, 0, 1, 0);
        MultiFab::Subtract(diff, sol_native, 0, 0, 1, 0);
        Real const dmax = diff.norminf();

        amrex::Print() << "\nfinal residual: native = " << err_native
                       << ", " << name << " = " << err << "\n"
                       << "max |phi|              = " << smax << "\n"
                       << "max |phi_" << name << " - phi|  = " << dmax << "\n"
                       << "relative difference    = " << dmax/smax << '\n';

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dmax <= max_rel_diff*smax,
            "The " + name + " bottom solver did not reproduce the native solution");

        if (plot) {
            MultiFab plotmf(nba, dmap, 3, 0);
            MultiFab::Copy(plotmf, sol_native, 0, 0, 1, 0);
            MultiFab::Copy(plotmf, sol       , 0, 1, 1, 0);
            MultiFab::Copy(plotmf, diff      , 0, 2, 1, 0);
            WriteSingleLevelPlotfile("plot_"+name, plotmf, {"phi_native","phi_"+name,"diff"},
                                     geom, 0.0, 0);
        }
    };

    compare("algmg", BottomSolver::algmg);
#ifdef AMREX_USE_HYPRE
    compare("hypre", BottomSolver::hypre);
#endif
}

// A sphere covering one node at an odd column and a row of 2 mod 4 turns
// into a blocked edge on MG level 1 and would be lost on level 2, while a
// second sphere keeps level 2 in use.  The coarse correction diverges if
// level 2 is used, so the operator must stop at 2 levels.  EB2 does not
// coarsen here; the operator builds its own coarse EB data.
void test_hidden_feature (Geometry const& geom, BoxArray const& grids,
                          DistributionMapping const& dmap,
                          Array<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                          Array<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                          int n_cell, Real reltol, int verbose)
{
    // Smaller grids stop at 2 levels for lack of open nodes instead.
    if (n_cell % 8 != 0 || n_cell < 32) {
        amrex::Print() << "skipped: n_cell must be a multiple of 8 and at least 32\n";
        return;
    }

    Real const dx = geom.CellSize(0);
    EB2::SphereIF anchor(Real(2.5)*dx,
                         {AMREX_D_DECL(Real(0.25), Real(0.25), Real(0.25))}, false);
    EB2::SphereIF hidden(Real(0.5)*dx,
                         {AMREX_D_DECL(Real(n_cell/2+1)*dx, Real(n_cell/2+2)*dx, // NOLINT(bugprone-integer-division)
                                       Real(n_cell/2)*dx)}, false); // NOLINT(bugprone-integer-division)
    EB2::Build(EB2::makeShop(EB2::makeUnion(anchor, hidden)), geom, 0, 0);
    auto factory = makeEBFabFactory(geom, grids, dmap, {2,2,2}, EBSupport::full);
    auto const& ebfactory = *static_cast<EBFArrayBoxFactory const*>(factory.get());

    LPInfo info;
    info.setMaxCoarseningLevel(30);
    MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info, {&ebfactory});
    linop.setDomainBC(lobc, hibc);
    linop.setSigma({AMREX_D_DECL(Real(1.0), Real(1.0), Real(1.0))});
    linop.setEBDirichlet(Real(1.0));

    BoxArray const& nba = amrex::convert(grids, IntVect(1));
    MultiFab rhs(nba, dmap, 1, 0);
    rhs.setVal(Real(1.0));
    MultiFab sol(nba, dmap, 1, 1);
    sol.setVal(Real(0.0));

    MLMG mlmg(linop);
    mlmg.setVerbose(verbose);
    mlmg.solve({&sol}, {&rhs}, reltol, Real(0.0));

    amrex::Print() << "# of MG levels: " << linop.NMGLevels(0) << "\n";
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(linop.NMGLevels(0) == 2,
        "Coarsening did not stop before the level that hides an EB feature");
}

// Stretched cells on the unit domain with \p ncells.  MG coarsens only the
// short directions until the cells are nearly cubic, and must converge
// quickly.
void test_stretched_cells (Array<int,AMREX_SPACEDIM> const& is_periodic,
                           Array<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                           Array<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                           IntVect const& ncells, Array<Real,AMREX_SPACEDIM> const& sigma,
                           int max_grid_size, Real reltol, int verbose)
{
    Box const domain(IntVect(0), ncells-1);
    RealBox const rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray grids(domain);
    grids.maxSize(max_grid_size);
    DistributionMapping const dmap(grids);

    EB2::SphereIF sphere(Real(0.23), {AMREX_D_DECL(Real(0.5), Real(0.5), Real(0.5))}, false);
    EB2::Build(EB2::makeShop(sphere), geom, 0, 0);
    auto factory = makeEBFabFactory(geom, grids, dmap, {2,2,2}, EBSupport::full);
    auto const& ebfactory = *static_cast<EBFArrayBoxFactory const*>(factory.get());

    LPInfo info;
    info.setMaxCoarseningLevel(30);
    MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info, {&ebfactory});
    linop.setDomainBC(lobc, hibc);
    linop.setSigma(sigma); // after define: the MG levels are built at the solve
    linop.setEBDirichlet(Real(1.0));

    BoxArray const& nba = amrex::convert(grids, IntVect(1));
    MultiFab rhs(nba, dmap, 1, 0);
    rhs.setVal(Real(1.0));
    MultiFab sol(nba, dmap, 1, 1);
    sol.setVal(Real(0.0));

    AMREX_ALWAYS_ASSERT(linop.NMGLevels(0) == 1);

    MLMG mlmg(linop);
    mlmg.setVerbose(verbose);
    mlmg.setMaxIter(20);
    mlmg.solve({&sol}, {&rhs}, reltol, Real(0.0));

    IntVect const len1 = linop.Geom(0,1).Domain().length();
    amrex::Print() << "# of MG levels: " << linop.NMGLevels(0)
                   << ", MG level 1: " << len1 << "\n";
    // Cells with weaker coupling count as longer: dx/sqrt(sigma).
    Array<Real,AMREX_SPACEDIM> h;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        h[idim] = Real(1.0) / (Real(ncells[idim]) * std::sqrt(sigma[idim]));
    }
    Real const hmin = *std::min_element(h.begin(), h.end());
    IntVect expected = ncells;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (h[idim] < Real(1.5)*hmin) { expected[idim] /= 2; }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(len1 == expected,
        "MG level 1 should be coarsened in the short directions only");
}

// The Helmholtz term of setACoeffs.  With phi = c on the EB and on the
// Dirichlet domain boundaries, del dot (grad phi) - a phi = -a c has the exact
// discrete solution phi = c for any a: the Laplacian of a constant vanishes in
// every row, including rows cut by the EB (whose row scale must then multiply
// a as it multiplies the right-hand side), rows on a Neumann boundary and, in
// RZ, the axis.  A wrong sign or a missing scale gives phi != c.  The problem
// is solved with multigrid (the coarse levels use the restricted a), without
// coarsening, and with the default, AlgMG and hypre bottom solvers (the
// latter two use the assembled matrix).  Also checked:
//   * a = inf at the covered nodes, which the operator does not use, must not
//     reach the coarse levels;
//   * a single Box with the default bottom solver and no coarsening, so that
//     the bottom solve does all of the work: on GPUs, the default is a
//     single-Box solver without a Helmholtz term, which must not be used;
//   * a second, non-constant right-hand side, for which the matrix solvers
//     and multigrid must agree;
//   * a new a on an existing operator and MLMG, as in a Newton iteration,
//     must give the solution of a new operator.
void test_helmholtz (Geometry const& geom, BoxArray const& grids,
                     DistributionMapping const& dmap,
                     EBFArrayBoxFactory const& factory,
                     Array<LinOpBCType,AMREX_SPACEDIM> const& lobc,
                     Array<LinOpBCType,AMREX_SPACEDIM> const& hibc,
                     int max_coarsening_level, Real reltol, int verbose, bool rz)
{
    Real const c = 2.0;
    Real const tol = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-9);
    if constexpr (std::is_same_v<Real,double>) {
        reltol = std::min(reltol, Real(1.e-12));
    }

    BoxArray const& nba = amrex::convert(grids, IntVect(1));
    auto const dx = geom.CellSizeArray();
    auto const problo = geom.ProbLoArray();
    Real const twopi = Real(2.0)*Math::pi<Real>();
    // a h^2 from 0.01 to 10: from a weak shift to a diagonally dominant row
    Real const a0 = Real(1.0)/(dx[0]*dx[0]);

    MultiFab acoef(nba, dmap, 1, 0);
    MultiFab acoef_inf(nba, dmap, 1, 0);
    MultiFab rhs_const(nba, dmap, 1, 0);
    MultiFab rhs_wave(nba, dmap, 1, 0);
    auto const& levset = factory.getLevelSet();
    for (MFIter mfi(acoef); mfi.isValid(); ++mfi) {
        Array4<Real> const& a = acoef.array(mfi);
        Array4<Real> const& ainf = acoef_inf.array(mfi);
        Array4<Real const> const& ls = levset.const_array(mfi);
        Array4<Real> const& r1 = rhs_const.array(mfi);
        Array4<Real> const& r2 = rhs_wave.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            AMREX_D_TERM(Real const x = problo[0] + i*dx[0];,
                         Real const y = problo[1] + j*dx[1];,
                         Real const z = problo[2] + k*dx[2];)
            Real const s = Real(0.5)*(Real(1.0) + std::sin(twopi*x)
                * AMREX_D_TERM(Real(1.0), *std::cos(twopi*y), *std::cos(twopi*z)));
            a(i,j,k) = a0 * (Real(0.01) + Real(10.0)*s);
            ainf(i,j,k) = (ls(i,j,k) >= Real(0.0)) ? std::numeric_limits<Real>::infinity()
                                                  : a(i,j,k);
            r1(i,j,k) = -a(i,j,k) * c;
            r2(i,j,k) = std::sin(twopi*x)
                * AMREX_D_TERM(Real(1.0), *std::cos(Real(3.0)*twopi*y), *std::cos(Real(2.0)*twopi*z));
        });
    }

    auto setup_linop = [&] (MLEBNodeFDLaplacian& linop, MultiFab const& a)
    {
        linop.setDomainBC(lobc, hibc);
        linop.setRZ(rz);
        linop.setSigma({AMREX_D_DECL(Real(1.0), Real(1.0), Real(1.0))});
        linop.setEBDirichlet(c);
        linop.setACoeffs(0, a);
    };

    auto setup_mlmg = [&] (MLMG& mlmg, BottomSolver bottom_solver)
    {
        mlmg.setVerbose(verbose);
        mlmg.setBottomSolver(bottom_solver);
        mlmg.setBottomTolerance(std::is_same_v<Real,float> ? Real(1.e-6) : Real(1.e-12));
        mlmg.setBottomMaxIter(2000);
        mlmg.setMaxIter(100);
    };

    // The domain Dirichlet nodes take their value from the initial guess.
    // Start from 0 elsewhere, so that the solver has to find phi = c.
    auto init_sol = [&] (MultiFab& sol, BoxArray const& ba, DistributionMapping const& dm)
    {
        sol.define(ba, dm, 1, 1);
        sol.setVal(c);
        Box const ndomain = amrex::surroundingNodes(geom.Domain());
        for (MFIter mfi(sol); mfi.isValid(); ++mfi) {
            Array4<Real> const& s = sol.array(mfi);
            Box const inner = mfi.validbox() & amrex::grow(ndomain, -1);
            amrex::ParallelFor(inner, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                s(i,j,k) = Real(0.0);
            });
        }
    };

    auto solve = [&] (MLMG& mlmg, MultiFab const& rhs, MultiFab& sol) -> int
    {
        MultiFab rhs_copy(rhs.boxArray(), rhs.DistributionMap(), 1, 0);
        MultiFab::Copy(rhs_copy, rhs, 0, 0, 1, 0);
        init_sol(sol, rhs.boxArray(), rhs.DistributionMap());
        mlmg.solve({&sol}, {&rhs_copy}, reltol, Real(0.0));
        return mlmg.getNumIters();
    };

    auto do_solve = [&] (MultiFab const& rhs, int max_coarsening, BottomSolver bottom_solver,
                         MultiFab& sol, MultiFab const& a) -> int
    {
        LPInfo info;
        info.setMaxCoarseningLevel(max_coarsening);
        MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info, {&factory});
        setup_linop(linop, a);
        MLMG mlmg(linop);
        setup_mlmg(mlmg, bottom_solver);
        return solve(mlmg, rhs, sol);
    };

    auto max_diff = [&] (MultiFab const& x, MultiFab const& y) {
        MultiFab d(x.boxArray(), x.DistributionMap(), 1, 0);
        MultiFab::Copy(d, x, 0, 0, 1, 0);
        MultiFab::Subtract(d, y, 0, 0, 1, 0);
        return d.norminf();
    };

    MultiFab exact(nba, dmap, 1, 0);
    exact.setVal(c);

    std::vector<std::pair<std::string,BottomSolver>> bottoms{
        {"default", BottomSolver::Default}, {"bicgstab", BottomSolver::bicgstab},
        {"algmg", BottomSolver::algmg}};
#ifdef AMREX_USE_HYPRE
    bottoms.emplace_back("hypre", BottomSolver::hypre);
#endif

    // constant solution, with and without coarsening, and with a = inf at the
    // covered nodes
    for (auto const& [name, bottom] : bottoms) {
        for (int mc : {max_coarsening_level, 0}) {
            for (bool inf_covered : {false, true}) {
                MultiFab sol;
                int const niters = do_solve(rhs_const, mc, bottom, sol,
                                            inf_covered ? acoef_inf : acoef);
                Real const err = max_diff(sol, exact);
                amrex::Print() << (rz ? "RZ, " : "") << "bottom " << name
                               << (mc > 0 ? ", multigrid" : ", no coarsening")
                               << (inf_covered ? ", a = inf where covered" : "")
                               << ": max |phi - c| / c = " << err/c
                               << " after " << niters << " iterations\n";
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(niters > 0,
                    "The Helmholtz test started from its solution");
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(err <= tol*c,
                    "The Helmholtz solve did not reproduce the constant solution");
            }
        }
    }

    // a single Box, the default bottom solver and no coarsening
    {
        BoxArray const ba1(geom.Domain());
        DistributionMapping const dm1(ba1);
        auto factory1 = makeEBFabFactory(geom, ba1, dm1, {2,2,2}, EBSupport::full);
        BoxArray const nba1 = amrex::convert(ba1, IntVect(1));
        MultiFab a1(nba1, dm1, 1, 0);
        MultiFab rhs1(nba1, dm1, 1, 0);
        MultiFab exact1(nba1, dm1, 1, 0);
        a1.ParallelCopy(acoef);
        rhs1.ParallelCopy(rhs_const);
        exact1.setVal(c);
        LPInfo info;
        info.setMaxCoarseningLevel(0);
        MLEBNodeFDLaplacian linop({geom}, {ba1}, {dm1}, info, {factory1.get()});
        setup_linop(linop, a1);
        MLMG mlmg(linop);
        setup_mlmg(mlmg, BottomSolver::Default);
        MultiFab sol;
        int const niters = solve(mlmg, rhs1, sol);
        Real const err = max_diff(sol, exact1);
        amrex::Print() << (rz ? "RZ, " : "") << "single Box, default bottom, no coarsening"
                       << ": max |phi - c| / c = " << err/c
                       << " after " << niters << " iterations\n";
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(err <= tol*c,
            "The single-Box Helmholtz solve did not reproduce the constant solution");
    }

    // non-constant solution: multigrid and the matrix-based bottom solvers agree
    MultiFab sol_mg;
    int const niters = do_solve(rhs_wave, max_coarsening_level, BottomSolver::bicgstab,
                                sol_mg, acoef);
    Real const smax = sol_mg.norminf(0, 0);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(niters <= 30,
        "Multigrid with the Helmholtz term converged too slowly");
    for (auto const& [name, bottom] : bottoms) {
        MultiFab sol;
        do_solve(rhs_wave, 0, bottom, sol, acoef);
        Real const d = max_diff(sol, sol_mg);
        amrex::Print() << (rz ? "RZ, " : "") << "wave, bottom " << name
                       << " vs. multigrid (" << niters << " iterations): rel. diff. = "
                       << d/smax << "\n";
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(d <= Real(100.0)*tol*smax,
            "The " + name + " bottom solver and multigrid disagree with a Helmholtz term");
    }

    // a new a on an existing operator and MLMG
    MultiFab acoef2(nba, dmap, 1, 0);
    MultiFab::Copy(acoef2, acoef, 0, 0, 1, 0);
    acoef2.mult(Real(3.0));
    for (auto const& [name, bottom] : bottoms) {
        for (int mc : {max_coarsening_level, 0}) {
            MultiFab sol_new;
            do_solve(rhs_wave, mc, bottom, sol_new, acoef2);

            LPInfo info;
            info.setMaxCoarseningLevel(mc);
            MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info, {&factory});
            setup_linop(linop, acoef);
            MLMG mlmg(linop);
            setup_mlmg(mlmg, bottom);
            MultiFab sol_first;
            solve(mlmg, rhs_wave, sol_first);
            linop.setACoeffs(0, acoef2);
            MultiFab sol_reused;
            solve(mlmg, rhs_wave, sol_reused);

            Real const s2 = sol_new.norminf(0, 0);
            Real const d = max_diff(sol_reused, sol_new);
            Real const d_first = max_diff(sol_first, sol_new);
            amrex::Print() << (rz ? "RZ, " : "") << "new a on an existing operator, bottom "
                           << name << (mc > 0 ? ", multigrid" : ", no coarsening")
                           << ": rel. diff. to a new operator = " << d/s2
                           << " (before the new a: " << d_first/s2 << ")\n";
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(d_first > Real(1.e3)*tol*s2,
                "The new a does not change the solution: the reuse test is void");
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(d <= Real(100.0)*tol*s2,
                "A new a on an existing operator and MLMG gives a wrong solution");
        }
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 64;
        int max_grid_size = 32;
        int max_coarsening_level = 30;
        int verbose = 1;
        // Single precision cannot reach the double precision tolerance.
        Real reltol = std::is_same_v<Real,float> ? Real(1.e-4) : Real(1.e-11);
        Vector<int> periodic{AMREX_D_DECL(0,0,0)};
        {
            ParmParse pp;
            pp.queryarr("periodic", periodic, 0, AMREX_SPACEDIM);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("max_coarsening_level", max_coarsening_level);
            pp.query("verbose", verbose);
            pp.query("reltol", reltol);
        }

        Box const domain(IntVect(0), IntVect(n_cell-1));
        RealBox const rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
        Array<int,AMREX_SPACEDIM> is_periodic{};
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            is_periodic[idim] = periodic[idim];
        }
        Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);

        BoxArray grids(domain);
        grids.maxSize(max_grid_size);
        DistributionMapping const dmap(grids);

        EB2::Build(geom, 0, max_coarsening_level);
        auto factory = makeEBFabFactory(geom, grids, dmap, {2,2,2}, EBSupport::full);
        auto const& ebfactory = *static_cast<EBFArrayBoxFactory const*>(factory.get());

        // Mix the boundary conditions so that both the Dirichlet nodes (which
        // are left out of the linear system) and the Neumann nodes (whose ghost
        // node is folded onto its mirror image) are covered.
        Array<LinOpBCType,AMREX_SPACEDIM> lobc
            {AMREX_D_DECL(LinOpBCType::Neumann, LinOpBCType::Dirichlet, LinOpBCType::Neumann)};
        Array<LinOpBCType,AMREX_SPACEDIM> hibc
            {AMREX_D_DECL(LinOpBCType::Dirichlet, LinOpBCType::Neumann, LinOpBCType::Dirichlet)};
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (geom.isPeriodic(idim)) {
                lobc[idim] = LinOpBCType::Periodic;
                hibc[idim] = LinOpBCType::Periodic;
            }
        }

        amrex::Print() << "\n==== EB Dirichlet value after operator reuse ====\n";
        test_eb_dirichlet_reuse(geom, grids, dmap, ebfactory,
                                lobc, hibc, max_coarsening_level, reltol, verbose);

        test_bottom_solvers(geom, grids, dmap, ebfactory,
                            lobc, hibc, reltol, verbose);

        amrex::Print() << "\n==== Helmholtz term ====\n";
        test_helmholtz(geom, grids, dmap, ebfactory, lobc, hibc,
                       max_coarsening_level, reltol, verbose, false);
#if (AMREX_SPACEDIM == 2)
        if (!geom.isPeriodic(0)) {
            amrex::Print() << "\n==== Helmholtz term, RZ ====\n";
            test_helmholtz(geom, grids, dmap, ebfactory, lobc, hibc,
                           max_coarsening_level, reltol, verbose, true);
        }
#endif

        amrex::Print() << "\n==== hidden EB feature ====\n";
        test_hidden_feature(geom, grids, dmap, lobc, hibc, n_cell, reltol, verbose);

        amrex::Print() << "\n==== stretched cells, one short direction ====\n";
        IntVect ncells(n_cell);
        ncells[AMREX_SPACEDIM-1] *= 8;
        Array<Real,AMREX_SPACEDIM> sigma{AMREX_D_DECL(Real(1.0), Real(1.0), Real(1.0))};
        test_stretched_cells(is_periodic, lobc, hibc, ncells, sigma, max_grid_size, reltol, verbose);

        amrex::Print() << "\n==== stretched cells, sigma makes them isotropic ====\n";
        sigma[AMREX_SPACEDIM-1] = Real(1.0)/Real(64.0);
        test_stretched_cells(is_periodic, lobc, hibc, ncells, sigma, max_grid_size, reltol, verbose);
#if (AMREX_SPACEDIM == 3)
        amrex::Print() << "\n==== stretched cells, two short directions ====\n";
        test_stretched_cells(is_periodic, lobc, hibc, IntVect(2*n_cell, 2*n_cell, n_cell/2),
                             {Real(1.0), Real(1.0), Real(1.0)}, max_grid_size, reltol, verbose);
#endif
    }
    amrex::Finalize();
}
