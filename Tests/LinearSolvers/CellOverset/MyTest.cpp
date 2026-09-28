#include "MyTest.H"

#include <AMReX_MLABecLaplacian.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

#include <numbers>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initData();
}

//
// Solve L(phi) = rhs
//
void
MyTest::solve ()
{
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        mlmg_lobc[idim] = LinOpBCType::Dirichlet;
        mlmg_hibc[idim] = LinOpBCType::Dirichlet;
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);

    std::unique_ptr<MLABecLaplacian> mlabec;
    if (do_overset) {
        mlabec = std::make_unique<MLABecLaplacian>(geom, grids, dmap,
                                                   GetVecOfConstPtrs(oversetmask),
                                                   info);
    } else {
        mlabec = std::make_unique<MLABecLaplacian>(geom, grids, dmap, info);
    }

    mlabec->setDomainBC(mlmg_lobc, mlmg_hibc);

    mlabec->setScalars(ascalar, bscalar);

    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        mlabec->setLevelBC(ilev, &exact_phi[ilev]);
        mlabec->setACoeffs(ilev, acoef[ilev]);

        Array<MultiFab,AMREX_SPACEDIM> face_bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const BoxArray& ba = amrex::convert(bcoef[ilev].boxArray(),
                                                IntVect::TheDimensionVector(idim));
            face_bcoef[idim].define(ba, bcoef[ilev].DistributionMap(), 1, 0);
        }
        amrex::average_cellcenter_to_face(GetArrOfPtrs(face_bcoef),
                                          bcoef[ilev], geom[ilev]);
        mlabec->setBCoeffs(ilev, amrex::GetArrOfConstPtrs(face_bcoef));
    }

    MLMG mlmg(*mlabec);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

#ifdef AMREX_USE_HYPRE
    if (use_hypre) {
        mlmg.setBottomSolver(amrex::BottomSolver::hypre);
    }
#endif

    Real tol_rel;
    if constexpr (std::is_same_v<double,Real>) {
        tol_rel = Real(1.0e-11);
    } else {
        tol_rel = Real(1.0e-4);
    }

    // In region with overset mask = 0, phi has valid solution and rhs is zero.
    mlmg.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), tol_rel, Real(0.0));

    // A failed solve often returns NaNs.  Check for them explicitly, because
    // the max-norm checks used by these tests silently drop NaNs.
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        if (phi[ilev].contains_nan(0, phi[ilev].nComp(), 0)) {
            amrex::Abort("MyTest::solve: solution contains NaN on level "
                         + std::to_string(ilev));
        }
    }

    if (do_overset) { checkOversetCells(); }
}

//
// The solver must not change phi where the overset mask is 0.  Cells covered
// by a finer level are skipped, because they hold averaged fine data.
//
void
MyTest::checkOversetCells () const
{
    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        iMultiFab known = (ilev < max_level)
            ? amrex::makeFineMask(grids[ilev], dmap[ilev], grids[ilev+1], IntVect(2), 1, 0)
            : iMultiFab(grids[ilev], dmap[ilev], 1, 0);
        if (ilev == max_level) { known.setVal(1); }

        MultiFab diff(grids[ilev], dmap[ilev], 1, 0);
        MultiFab::Copy(diff, phi[ilev], 0, 0, 1, 0);
        MultiFab::Subtract(diff, exact_phi[ilev], 0, 0, 1, 0);

        auto const& ka = known.arrays();
        auto const& ma = oversetmask[ilev].const_arrays();
        amrex::ParallelFor(known, [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
        {
            if (ma[b](i,j,k) != 0) { ka[b](i,j,k) = 0; }
        });
        Gpu::streamSynchronize();

        Real const drift = diff.norm0(known);
        amrex::Print() << " level " << ilev << " max change in overset cells: " << drift << '\n';
        if (drift > Real(0.0)) {
            amrex::Abort("MyTest: solver changed phi in overset cells on level "
                         + std::to_string(ilev));
        }
    }
}

void
MyTest::writePlotfile ()
{
    Vector<std::string> varname = {"solution", "rhs", "exact_solution", "error", "acoef", "bcoef"};
    const int nlevels = max_level + 1;
    Vector<MultiFab> plotmf(nlevels);
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        plotmf[ilev].define(grids[ilev], dmap[ilev], static_cast<int>(varname.size()), 0);
        MultiFab::Copy(plotmf[ilev], phi[ilev]      , 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], rhs[ilev]      , 0, 1, 1, 0);
        MultiFab::Copy(plotmf[ilev], exact_phi[ilev], 0, 2, 1, 0);
        MultiFab::Copy(plotmf[ilev], phi[ilev]      , 0, 3, 1, 0);
        MultiFab::Subtract(plotmf[ilev], plotmf[ilev], 2, 3, 1, 0); // error = soln - exact
        MultiFab::Copy(plotmf[ilev], acoef[ilev]    , 0, 4, 1, 0);
        MultiFab::Copy(plotmf[ilev], bcoef[ilev]    , 0, 5, 1, 0);
        auto dx = geom[ilev].CellSize();
        Real dvol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
        amrex::Print() << " level " << ilev
                       << " max-norm error: " << plotmf[ilev].norminf(3)
                       << " 1-norm error: " << plotmf[ilev].norm1(3)*dvol << '\n';
    }
    WriteMultiLevelPlotfile(plot_file_name, nlevels, GetVecOfConstPtrs(plotmf), varname,
                            geom, 0.0, Vector<int>(nlevels, 0),
                            Vector<IntVect>(nlevels, IntVect(2)));
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("plot_file", plot_file_name);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_coarsening_level", max_coarsening_level);

    pp.query("do_overset", do_overset);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_level >= 0 && max_level <= 2,
                                     "max_level must be 0, 1 or 2");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(n_cell % 8 == 0, "n_cell must be a multiple of 8");

#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
    if (use_hypre) max_coarsening_level = 0;
#endif
}

void
MyTest::initGrids ()
{
    const int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> isperiodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, isperiodic.data());
    Box domain(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain, rb, CoordSys::cartesian, isperiodic);
        if (ilev == 0) {
            grids[ilev].define(domain);
        } else {
            // Fine levels cover [0.125,0.625] and [0.1875,0.5625], so they
            // contain both overset and regular cells.
            const int n = domain.length(0);
            const int lo = n/8 + (ilev-1)*n/16;
            const int hi = 5*n/8 - (ilev-1)*n/16 - 1;
            grids[ilev].define(Box(IntVect(AMREX_D_DECL(lo,lo,lo)),
                                   IntVect(AMREX_D_DECL(hi,hi,hi))));
        }
        grids[ilev].maxSize(max_grid_size);
        domain.refine(2);
    }
}

void
MyTest::initData ()
{
    const int nlevels = max_level + 1;
    dmap.resize(nlevels);
    phi.resize(nlevels);
    rhs.resize(nlevels);
    exact_phi.resize(nlevels);
    acoef.resize(nlevels);
    bcoef.resize(nlevels);
    oversetmask.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);

        phi[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        exact_phi[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        bcoef[ilev].define(grids[ilev], dmap[ilev], 1, 1);
        oversetmask[ilev].define(grids[ilev], dmap[ilev], 1, 0);

        // Middle of the domain, [0.25,0.75] on every level
        const Box& domain = geom[ilev].Domain();
        Box overset_box = amrex::grow(domain, -domain.length(0)/4);

        const auto prob_lo = geom[ilev].ProbLoArray();
        const auto prob_hi = geom[ilev].ProbHiArray();
        const auto dx      = geom[ilev].CellSizeArray();
        auto a = ascalar;
        auto b = bscalar;
        auto loverset = do_overset;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs[ilev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& vbx = mfi.tilebox();
            const Box& gbx = mfi.growntilebox(1);

            auto phifab = phi[ilev].array(mfi);
            auto rhsfab = rhs[ilev].array(mfi);
            auto exact = exact_phi[ilev].array(mfi);
            auto alpha = acoef[ilev].array(mfi);
            auto beta = bcoef[ilev].array(mfi);
            auto mask = oversetmask[ilev].array(mfi);

            amrex::ParallelFor(gbx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                constexpr auto w = amrex::Real(0.05);
                constexpr auto sigma = amrex::Real(10.);
                const amrex::Real theta = amrex::Real(0.5)*std::log(amrex::Real(3.)) / w;

                constexpr amrex::Real pi = std::numbers::pi_v<amrex::Real>;
                constexpr amrex::Real tpi =  amrex::Real(2.)*pi;
                constexpr amrex::Real fpi =  amrex::Real(4.)*pi;
                constexpr amrex::Real fac = static_cast<amrex::Real>(AMREX_SPACEDIM)*amrex::Real(4.)*pi*pi;

                amrex::Real xc = (prob_hi[0] + prob_lo[0])*amrex::Real(0.5);
                amrex::Real yc = (prob_hi[1] + prob_lo[1])*amrex::Real(0.5);
#if (AMREX_SPACEDIM == 2)
                auto zc = amrex::Real(0.0);
#else
                amrex::Real zc = (prob_hi[2] + prob_lo[2])*amrex::Real(0.5);
#endif

                amrex::Real x = prob_lo[0] + dx[0] * (static_cast<amrex::Real>(i) + amrex::Real(0.5));
                amrex::Real y = prob_lo[1] + dx[1] * (static_cast<amrex::Real>(j) + amrex::Real(0.5));
#if (AMREX_SPACEDIM == 2)
                auto z = amrex::Real(0.0);
#else
                amrex::Real z = prob_lo[2] + dx[2] * (static_cast<amrex::Real>(k) + amrex::Real(0.5));
#endif

                amrex::Real r = std::sqrt((x-xc)*(x-xc) + (y-yc)*(y-yc) + (z-zc)*(z-zc));
                amrex::Real tmp = std::cosh(theta*(r-amrex::Real(0.25)));
                amrex::Real dbdrfac = (sigma-amrex::Real(1.))/amrex::Real(2.)/(tmp*tmp) * theta/r;
                dbdrfac *= b;

                // for domain boundary
                x = amrex::min(prob_hi[0], amrex::max(prob_lo[0], x));
                y = amrex::min(prob_hi[1], amrex::max(prob_lo[1], y));
#if (AMREX_SPACEDIM == 3)
                z = amrex::min(prob_hi[2], amrex::max(prob_lo[2], z));
#endif

                beta(i,j,k) = (sigma-amrex::Real(1.))/amrex::Real(2.)*std::tanh(theta*(r-amrex::Real(0.25)))
                    + (sigma+amrex::Real(1.))/amrex::Real(2.);
                exact(i,j,k) = std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                       + amrex::Real(0.25) * std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z);
                phifab(i,j,k) = amrex::Real(0.0);

                if (vbx.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                    alpha(i,j,k) = amrex::Real(1.);
                    rhsfab(i,j,k) = beta(i,j,k)*b*fac*(std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                                     + std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z))
                                + dbdrfac*((x-xc)*(tpi*std::sin(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                                  + pi*std::sin(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z))
                                         + (y-yc)*(tpi*std::cos(tpi*x) * std::sin(tpi*y) * std::cos(tpi*z)
                                                  + pi*std::cos(fpi*x) * std::sin(fpi*y) * std::cos(fpi*z))
                                         + (z-zc)*(tpi*std::cos(tpi*x) * std::cos(tpi*y) * std::sin(tpi*z)
                                                  + pi*std::cos(fpi*x) * std::cos(fpi*y) * std::sin(fpi*z)))
                                                + a * (std::cos(tpi*x) * std::cos(tpi*y) * std::cos(tpi*z)
                                              + amrex::Real(0.25) * std::cos(fpi*x) * std::cos(fpi*y) * std::cos(fpi*z));
                    if (loverset && overset_box.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        mask(i,j,k) = 0;
                        phifab(i,j,k) = exact(i,j,k);
                    } else {
                        mask(i,j,k) = 1;
                    }
                }
            });
        }
    }
}
