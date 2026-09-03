//
// Solves del dot (sigma grad phi) = rhs with MLEBNodeFDLaplacian and compares
// the solution obtained with the hypre bottom solver against the one obtained
// with the native bottom solver.  Set max_coarsening_level = 0 to make the
// bottom solve do all of the work, which is the sharpest test of the matrix
// assembled by MLEBNodeFDLaplacian::fillIJMatrix.
//

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 64;
        int max_grid_size = 32;
        int max_coarsening_level = 0;
        int verbose = 1;
        int bottom_verbose = 0;
        int use_sigma_mf = 1;
        int plot = 0;
        Real phi_eb = 1.0;
        Real reltol = 1.e-11;
        Real bottom_reltol = 1.e-9;
        Real max_rel_diff = 1.e-12;
        Vector<int> periodic{AMREX_D_DECL(0,0,0)};
        {
            ParmParse pp;
            pp.queryarr("periodic", periodic, 0, AMREX_SPACEDIM);
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("max_coarsening_level", max_coarsening_level);
            pp.query("verbose", verbose);
            pp.query("bottom_verbose", bottom_verbose);
            pp.query("use_sigma_mf", use_sigma_mf);
            pp.query("phi_eb", phi_eb);
            pp.query("reltol", reltol);
            pp.query("bottom_reltol", bottom_reltol);
            pp.query("plot", plot);
            pp.query("max_rel_diff", max_rel_diff);
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

        BoxArray const& nba = amrex::convert(grids, IntVect(1));

        // A smooth, strictly positive sigma so that the variable coefficient
        // path is exercised.
        MultiFab sigma(grids, dmap, 1, 1, MFInfo(), *factory);
        // Use whole wavenumbers so that sigma and the source stay single valued
        // when a direction is periodic.
        Real const twopi = Real(2.0)*Real(3.14159265358979323846);
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

        auto do_solve = [&] (BottomSolver bottom_solver, MultiFab& sol)
        {
            LPInfo info;
            info.setMaxCoarseningLevel(max_coarsening_level);

            MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info,
                                      {static_cast<EBFArrayBoxFactory const*>(factory.get())});
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
            Real const err = mlmg.solve({&sol}, {&rhs_copy}, reltol, Real(0.0));

            // A failed solve often returns NaNs.  Check for them explicitly,
            // because the max-norm checks below silently drop NaNs.
            if (sol.contains_nan(0, sol.nComp(), 0)) {
                amrex::Abort("do_solve: solution contains NaN");
            }

            return err;
        };

        MultiFab sol_native;
        MultiFab sol_hypre;

        amrex::Print() << "\n==== native bottom solver ====\n";
        Real const err_native = do_solve(BottomSolver::bicgstab, sol_native);

        amrex::Print() << "\n==== hypre bottom solver ====\n";
        Real const err_hypre = do_solve(BottomSolver::hypre, sol_hypre);

        MultiFab diff(nba, dmap, 1, 0);
        MultiFab::Copy(diff, sol_hypre, 0, 0, 1, 0);
        MultiFab::Subtract(diff, sol_native, 0, 0, 1, 0);
        Real const dmax = diff.norminf();
        Real const smax = sol_native.norminf(0, 0);

        amrex::Print() << "\nfinal residual: native = " << err_native
                       << ", hypre = " << err_hypre << "\n"
                       << "max |phi|              = " << smax << "\n"
                       << "max |phi_hypre - phi|  = " << dmax << "\n"
                       << "relative difference    = " << dmax/smax << std::endl;

        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(dmax <= max_rel_diff*smax,
            "The hypre bottom solver did not reproduce the native solution");

        if (plot) {
            MultiFab plotmf(nba, dmap, 3, 0);
            MultiFab::Copy(plotmf, sol_native, 0, 0, 1, 0);
            MultiFab::Copy(plotmf, sol_hypre , 0, 1, 1, 0);
            MultiFab::Copy(plotmf, diff      , 0, 2, 1, 0);
            WriteSingleLevelPlotfile("plot", plotmf, {"phi_native","phi_hypre","diff"}, geom, 0.0, 0);
        }
    }
    amrex::Finalize();
}
