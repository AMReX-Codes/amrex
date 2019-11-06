#include "MyTest.H"

#include <AMReX_MLNodeTensorLaplacian.H>
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
    MLNodeTensorLaplacian linop(geom, grids, dmap,
                                LPInfo().setMaxCoarseningLevel(max_coarsening_level));

    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});

    linop.setBeta(beta);

    MLMG mlmg(linop);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
#ifdef AMREX_USE_HYPRE
    if (use_hypre) {
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    } else
#endif
    {
        mlmg.setBottomSolver(MLMG::BottomSolver::cg);
    }

    mlmg.setBottomMaxIter(1000);

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

#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
    if (use_hypre) max_coarsening_level = 0;
#endif

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_coarsening_level", max_coarsening_level);
    pp.query("reltol", reltol);

    Vector<Real> vbeta;
    pp.queryarr("beta", vbeta);
    if (!vbeta.empty()) {
        AMREX_D_TERM(beta[0] = vbeta[0];, beta[1] = vbeta[1];, beta[2] = vbeta[2];);
    }
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

    GpuArray<Real,AMREX_SPACEDIM> lbeta{AMREX_D_DECL(beta[0],beta[1],beta[2])};

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const BoxArray& nba = amrex::convert(grids[ilev],IntVect::TheNodeVector());
        // These are nodal
        solution      [ilev].define(nba        , dmap[ilev], 1, 0);
        rhs           [ilev].define(nba        , dmap[ilev], 1, 0);
        exact_solution[ilev].define(nba        , dmap[ilev], 1, 0);

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
                constexpr Real fac = 4.*pi*pi;

                Real x = i*dx[0];
                Real y = j*dx[1];
                Real z = k*dx[2];

                phi(i,j,k) = (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z))
                    + 0.25 * (std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));

                Real d2phidx2 = -fac * (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                      + std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));
                Real d2phidy2 = d2phidx2;
                Real d2phidxdy = fac * (std::sin(tpi*x) * std::sin(tpi*y) * std::cos(tpi*z)
                                      + std::sin(fpi*x) * std::sin(fpi*y) * std::cos(fpi*z));

#if (AMREX_SPACEDIM == 2)
                rh(i,j,k) = (1.0-lbeta[0]*lbeta[0]) * d2phidx2
                    +       (1.0-lbeta[1]*lbeta[1]) * d2phidy2
                    -         2.*lbeta[0]*lbeta[1]  * d2phidxdy;
#else
                Real d2phidz2 = d2phidx2;
                Real d2phidxdz = fac * (std::sin(tpi*x) * std::sin(tpi*z) * std::cos(tpi*y)
                                      + std::sin(fpi*x) * std::sin(fpi*z) * std::cos(fpi*y));
                Real d2phidydz = fac * (std::sin(tpi*y) * std::sin(tpi*z) * std::cos(tpi*x)
                                      + std::sin(fpi*y) * std::sin(fpi*z) * std::cos(fpi*x));
                rh(i,j,k) = (1.0-lbeta[0]*lbeta[0]) * d2phidx2
                    +       (1.0-lbeta[1]*lbeta[1]) * d2phidy2
                    +       (1.0-lbeta[2]*lbeta[2]) * d2phidz2
                    -         2.*lbeta[0]*lbeta[1]  * d2phidxdy
                    -         2.*lbeta[0]*lbeta[2]  * d2phidxdz
                    -         2.*lbeta[1]*lbeta[2]  * d2phidydz;
#endif
            });
        }
    }
}

