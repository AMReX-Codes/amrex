#include "MyTest.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();
    initData();
}

void
MyTest::solve ()
{
    if (composite_solve)
    {
        MLNodeLaplacian linop(geom, grids, dmap);

        linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)},
                          {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet,
                                        LinOpBCType::Dirichlet)});

        for (int ilev = 0; ilev <= max_level; ++ilev) {
            linop.setSigma(ilev, sigma[ilev]);
        }

        MLMG mlmg(linop);
        mlmg.setMaxIter(max_iter);
        mlmg.setMaxFmgIter(max_fmg_iter);
        mlmg.setVerbose(verbose);
        mlmg.setBottomVerbose(bottom_verbose);

        // solution is passed to MLMG::solve to provide an initial guess.
        // Additionally it also provides boundary conditions for Dirichlet
        // boundaries if there are any.
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            MultiFab::Copy(solution[ilev], exact_solution[ilev], 0, 0, 1, 0);
            const Box& interior = amrex::surroundingNodes(
                amrex::grow(geom[ilev].Domain(), -1));
            // Usually we want the best initial guess.  For testing here,
            // we set the domain boundaries to exact solution and zero out
            // the interior.
            solution[ilev].setVal(0.0, interior, 0, 1, 0);
        }

        mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), reltol, 0.0);
    }
    else // solve level by level
    {
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            MLNodeLaplacian linop({geom[ilev]}, {grids[ilev]}, {dmap[ilev]});

            linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)},
                              {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet,
                                            LinOpBCType::Dirichlet)});

            linop.setSigma(0, sigma[ilev]); // set solver's level 0 sigma.

            MLMG mlmg(linop);
            mlmg.setMaxIter(max_iter);
            mlmg.setMaxFmgIter(max_fmg_iter);
            mlmg.setVerbose(verbose);
            mlmg.setBottomVerbose(bottom_verbose);

            // solution is passed to MLMG::solve to provide an initial guess.
            // Additionally it also provides boundary conditions for Dirichlet
            // boundaries if there are any.
            if (ilev == 0) {
                MultiFab::Copy(solution[ilev], exact_solution[ilev], 0, 0, 1, 0);
                const Box& interior = amrex::surroundingNodes(
                    amrex::grow(geom[ilev].Domain(), -1));
                // Usually we want the best initial guess.  For testing here,
                // we set the domain boundaries to exact solution and zero out
                // the interior.
                solution[ilev].setVal(0.0, interior, 0, 1, 0);
            } else {
                // Coarse/fine boundary is Dirichlet.
                // For fine levels, we interpolate from coarse to fine to set up
                // the coarse/fine Dirichlet boundary.
                PhysBCFunctNoOp bcnoop;
                Vector<BCRec> bcrec(1);
                amrex::InterpFromCoarseLevel(solution[ilev], 0.0,
                                             solution[ilev-1], 0, 0, 1,
                                             geom[ilev-1], geom[ilev],
                                             bcnoop, 0, bcnoop, 0,
                                             IntVect{ref_ratio},
                                             &node_bilinear_interp,
                                             bcrec, 0);
            }

            mlmg.solve({&solution[ilev]}, {&rhs[ilev]}, reltol, 0.0);
        }
    }
}

void
MyTest::compute_norms () const
{
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::Print() << "Level " << ilev << "\n";
        MultiFab error(solution[ilev].boxArray(), solution[ilev].DistributionMap(), 1, 0);
        MultiFab::Copy(error, solution[ilev], 0, 0, 1, 0);
        MultiFab::Subtract(error, exact_solution[ilev], 0, 0, 1, 0);

        auto mask = error.OwnerMask(geom[ilev].periodicity());

        amrex::Print() << "    max-norm: " << error.norm0(*mask, 0, 0) << "\n";
        const Real* dx = geom[ilev].CellSize();
        Real dvol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
        amrex::Print() << "    1-norm  : " << error.norm1(0, geom[ilev].periodicity())*dvol << "\n";
    }
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("ref_ratio", ref_ratio);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("composite_solve", composite_solve);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("reltol", reltol);
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);
    dmap.resize(nlevels);

    solution.resize(nlevels);
    rhs.resize(nlevels);
    exact_solution.resize(nlevels);
    sigma.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        grids[ilev].define(domain);
        grids[ilev].maxSize(max_grid_size);
        domain.grow(-n_cell/4);   // fine level cover the middle of the coarse domain
        domain.refine(ref_ratio);
    }

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const BoxArray& nba = amrex::convert(grids[ilev],IntVect::TheNodeVector());
        // These are nodal
        solution      [ilev].define(nba        , dmap[ilev], 1, 0);
        rhs           [ilev].define(nba        , dmap[ilev], 1, 0);
        exact_solution[ilev].define(nba        , dmap[ilev], 1, 0);
        // sigma is cell-centered.
        sigma         [ilev].define(grids[ilev], dmap[ilev], 1, 0);

        const auto dx = geom[ilev].CellSizeArray();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs[ilev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            Array4<Real> const phi = exact_solution[ilev].array(mfi);
            Array4<Real> const rh  = rhs[ilev].array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                constexpr Real pi = 3.1415926535897932;
                constexpr Real tpi = 2.*pi;
                constexpr Real fpi = 4.*pi;
                constexpr Real fac = tpi*tpi*AMREX_SPACEDIM;

                Real x = i*dx[0];
                Real y = j*dx[1];
                Real z = k*dx[2];

                phi(i,j,k) = (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z))
                    + 0.25 * (std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));

                rh(i,j,k) = -fac * (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z))
                    -        fac * (std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));
            });
        }

        sigma[ilev].setVal(1.0);
    }
}

