#include "MyTest.H"

#include <AMReX_MLEBABecLap.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_GMRES_MLMG.H>

#include <cmath>

using namespace amrex;

MyTest::MyTest ()
{
    readParameters();

    initGrids();

    initializeEB();

    initData();
}

void
MyTest::solve ()
{
    solve(false);

    baseline_phi.resize(max_level+1);
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        baseline_phi[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        MultiFab::Copy(baseline_phi[ilev], phi[ilev], 0, 0, 1, 1);
    }
    writePlotfile(baseline_phi, plot_file_name + "_nomask", false);

    if (do_overset) {
        prepareOversetSolve();
        solve(true);
        reportOversetError();
        writePlotfile(phi, plot_file_name + "_overset", true);
    }
}

void
MyTest::solve (bool use_overset_mask)
{
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        const MultiFab& vfrc = factory[ilev]->getVolFrac();
        MultiFab v(vfrc.boxArray(), vfrc.DistributionMap(), 1, 0,
                   MFInfo(), *factory[ilev]);
        MultiFab::Copy(v, vfrc, 0, 0, 1, 0);
        amrex::EB_set_covered(v, 1.0);
        amrex::Print() << "vfrc min = " << v.min(0) << '\n';
    }

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (geom[0].isPeriodic(idim)) {
            mlmg_lobc[idim] = LinOpBCType::Periodic;
            mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
            mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            mlmg_hibc[idim] = LinOpBCType::Dirichlet;
        }
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);

    std::unique_ptr<MLEBABecLap> mleb;
    if (use_overset_mask) {
        Vector<iMultiFab const*> mask_ptrs(max_level+1);
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            mask_ptrs[ilev] = &overset_mask[ilev];
        }
        mleb = std::make_unique<MLEBABecLap>(geom, grids, dmap, mask_ptrs,
                                             info, amrex::GetVecOfConstPtrs(factory));
    } else {
        mleb = std::make_unique<MLEBABecLap>(geom, grids, dmap, info,
                                             amrex::GetVecOfConstPtrs(factory));
    }
    if (use_hypre || use_petsc) {
        if (factory[0]->isAllRegular()) {
            linop_maxorder = std::min(3,linop_maxorder);
        } else {
            linop_maxorder = 2;
        }
    }
    mleb->setMaxOrder(linop_maxorder);

    mleb->setDomainBC(mlmg_lobc, mlmg_hibc);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mleb->setLevelBC(ilev, &phi[ilev]);
    }

    mleb->setScalars(scalars[0], scalars[1]);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mleb->setACoeffs(ilev, acoef[ilev]);
        mleb->setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcoef[ilev]));
    }

    if (eb_is_dirichlet) {
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            mleb->setEBDirichlet(ilev, phi[ilev], bcoef_eb[ilev]);
        }
    }

    MLMG mlmg(*mleb);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setBottomMaxIter(max_bottom_iter);
    mlmg.setBottomTolerance(bottom_reltol);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    if (use_hypre) {
        mlmg.setBottomSolver(MLMG::BottomSolver::hypre);
    } else if (use_petsc) {
        mlmg.setBottomSolver(MLMG::BottomSolver::petsc);
    }
    const Real tol_rel = reltol;
    const Real tol_abs = 0.0;

    if (verbose) {
        Vector<MultiFab> res(max_level+1);
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            res[ilev].define(rhs[ilev].boxArray(), rhs[ilev].DistributionMap(), 1, 0);
        }
        mlmg.apply(GetVecOfPtrs(res), GetVecOfPtrs(phi)); // res = L(sol)
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            MultiFab::Subtract(res[ilev], rhs[ilev], 0, 0, 1, 0);
            amrex::Print() << "Initial max, 1 and 2-norm residuals at level " << ilev << " = "
                           << res[ilev].norminf(0) << " " << res[ilev].norm1(0) << " "
                           << res[ilev].norm2(0) << '\n';
        }
    }

    if (use_gmres) {
        GMRESMLMGT<MultiFab> gmsolver(mlmg);
        gmsolver.usePrecond(true);
        gmsolver.setVerbose(verbose);
        gmsolver.solve(GetVecOfPtrs(phi), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    } else {
        mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }

    // A failed solve often returns NaNs.  Check for them explicitly, because
    // the max-norm checks used by these tests silently drop NaNs.
    for (int ilev = 0; ilev < int(phi.size()); ++ilev) {
        if (phi[ilev].contains_nan(0, phi[ilev].nComp(), 0)) {
            amrex::Abort("MyTest::solve: solution contains NaN on level "
                         + std::to_string(ilev));
        }
    }

    if (verbose) {
        Vector<MultiFab> res(max_level+1);
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            res[ilev].define(rhs[ilev].boxArray(), rhs[ilev].DistributionMap(), 1, 0);
        }
        mlmg.apply(GetVecOfPtrs(res), GetVecOfPtrs(phi)); // res = L(sol)
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            MultiFab::Subtract(res[ilev], rhs[ilev], 0, 0, 1, 0);
            amrex::Print() << "Final max, 1 and 2-norm residuals at level " << ilev << " = "
                           << res[ilev].norminf(0) << " " << res[ilev].norm1(0) << " "
                           << res[ilev].norm2(0) << '\n';
        }
    }
}

void
MyTest::writePlotfile (Vector<MultiFab> const& plot_phi, std::string const& plotfile,
                       bool include_overset_mask) const
{
    Vector<MultiFab> plotmf(max_level+1);
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        const MultiFab& vfrc = factory[ilev]->getVolFrac();
        int const ncomp = include_overset_mask ? 3 : 2;
        plotmf[ilev].define(grids[ilev],dmap[ilev],ncomp,0);
        MultiFab::Copy(plotmf[ilev], plot_phi[ilev], 0, 0, 1, 0);
        MultiFab::Copy(plotmf[ilev], vfrc, 0, 1, 1, 0);
        if (include_overset_mask) {
            plotmf[ilev].setVal(-1.0, 2, 1, 0);
            auto const& pmf = plotmf[ilev];
            for (MFIter mfi(pmf); mfi.isValid(); ++mfi) {
                amrex::ParallelFor(mfi.validbox(),
                [plot = plotmf[ilev].array(mfi), mask = overset_mask[ilev].const_array(mfi)]
                AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    plot(i,j,k,2) = static_cast<Real>(mask(i,j,k));
                });
            }
        }
    }
    WriteMultiLevelPlotfile(plotfile, max_level+1,
                            amrex::GetVecOfConstPtrs(plotmf),
                            include_overset_mask
                                ? Vector<std::string>{"phi","vfrac","overset_mask"}
                                : Vector<std::string>{"phi","vfrac"},
                            geom, 0.0, Vector<int>(max_level+1,0),
                            Vector<IntVect>(max_level,IntVect{2}));

}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("max_level", max_level);
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);
    pp.query("is_periodic", is_periodic);
    pp.query("eb_is_dirichlet", eb_is_dirichlet);

    pp.query("plot_file", plot_file_name);
    pp.query("do_overset", do_overset);
    if (do_overset) {
        Vector<Real> mask_lo(AMREX_SPACEDIM);
        Vector<Real> mask_hi(AMREX_SPACEDIM);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            mask_lo[idim] = overset_mask_lo[idim];
            mask_hi[idim] = overset_mask_hi[idim];
        }
        pp.getarr("overset_mask_lo", mask_lo);
        pp.getarr("overset_mask_hi", mask_hi);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            overset_mask_lo[idim] = mask_lo[idim];
            overset_mask_hi[idim] = mask_hi[idim];
        }
    }

    scalars.resize(2);
    if (is_periodic) {
        scalars[0] = 0.0;
        scalars[1] = 1.0;
    } else {
        scalars[0] = 1.0;
        scalars[1] = 1.0;
    }
    pp.queryarr("scalars", scalars);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_bottom_iter", max_bottom_iter);
    pp.query("bottom_reltol", bottom_reltol);
    pp.query("reltol", reltol);
    pp.query("linop_maxorder", linop_maxorder);
    pp.query("max_coarsening_level", max_coarsening_level);
#ifdef AMREX_USE_HYPRE
    pp.query("use_hypre", use_hypre);
#endif
#ifdef AMREX_USE_PETSC
    pp.query("use_petsc",use_petsc);
#endif

    pp.query("use_gmres", use_gmres);
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> isperiodic{AMREX_D_DECL(is_periodic,is_periodic,is_periodic)};
    Geometry::Setup(&rb, 0, isperiodic.data());
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
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    dmap.resize(nlevels);
    factory.resize(nlevels);
    phi.resize(nlevels);
    rhs.resize(nlevels);
    rhs_original.resize(nlevels);
    acoef.resize(nlevels);
    bcoef.resize(nlevels);
    bcoef_eb.resize(nlevels);
    overset_mask.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev] = std::make_unique<EBFArrayBoxFactory>
            (eb_level, geom[ilev], grids[ilev], dmap[ilev], Vector<int>{2,2,2}, EBSupport::full);

        phi[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        rhs_original[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        overset_mask[ilev].define(grids[ilev], dmap[ilev], 1, 0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].define(amrex::convert(grids[ilev],IntVect::TheDimensionVector(idim)),
                                     dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        }
        if (eb_is_dirichlet) {
            bcoef_eb[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
            bcoef_eb[ilev].setVal(1.0);
        }

        phi[ilev].setVal(0.0);
        rhs[ilev].setVal(0.0);
        acoef[ilev].setVal(1.0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].setVal(1.0);
        }

        const auto dx = geom[ilev].CellSizeArray();

        if (is_periodic)
        {
            const Real pi = 4.0*std::atan(1.0);

            for (MFIter mfi(rhs[ilev]); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                Array4<Real> const& fab = rhs[ilev].array(mfi);
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rx = (i+0.5)*dx[0];
                    Real ry = (j+0.5)*dx[1];
                    fab(i,j,k) = std::sin(rx*2.*pi + 43.5)*std::sin(ry*2.*pi + 89.);
                });
            }
        }
        else if (eb_is_dirichlet)
        {
            phi[ilev].setVal(10.0);
            phi[ilev].setVal(0.0, 0, 1, 0); // set interior
        }
        else
        {
            // Initialize Dirichlet boundary
            for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                Array4<Real> const& fab = phi[ilev].array(mfi);
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rx = (i+0.5)*dx[0];
                    Real ry = (j+0.5)*dx[1];
                    fab(i,j,k) = std::sqrt(0.5)*(rx + ry);
                });
            }
            phi[ilev].setVal(0.0, 0, 1, 0); // set interior to zero
        }

    }

    for (int ilev = 0; ilev < nlevels; ++ilev) {
        MultiFab::Copy(rhs_original[ilev], rhs[ilev], 0, 0, 1, 0);
    }

    buildOversetMask();
}

void
MyTest::resetPhi ()
{
    int const nlevels = max_level + 1;
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        phi[ilev].setVal(0.0);
        if (eb_is_dirichlet) {
            phi[ilev].setVal(10.0);
            phi[ilev].setVal(0.0, 0, 1, 0);
        } else if (!is_periodic) {
            const auto dx = geom[ilev].CellSizeArray();
            for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.fabbox();
                Array4<Real> const& fab = phi[ilev].array(mfi);
                amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    Real rx = (i+0.5)*dx[0];
                    Real ry = (j+0.5)*dx[1];
                    fab(i,j,k) = std::sqrt(0.5)*(rx + ry);
                });
            }
            phi[ilev].setVal(0.0, 0, 1, 0);
        }
        MultiFab::Copy(rhs[ilev], rhs_original[ilev], 0, 0, 1, 0);
    }
}

void
MyTest::buildOversetMask ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        overset_mask[ilev].setVal(1);

        const auto problo = geom[ilev].ProbLoArray();
        const auto dx = geom[ilev].CellSizeArray();
        auto const mask_lo = overset_mask_lo;
        auto const mask_hi = overset_mask_hi;

        for (MFIter mfi(overset_mask[ilev]); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.validbox();
            auto const& mask = overset_mask[ilev].array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real x = problo[0] + (i + 0.5_rt) * dx[0];
                Real y = problo[1] + (j + 0.5_rt) * dx[1];
#if (AMREX_SPACEDIM == 2)
                Real z = 0.0_rt;
#else
                Real z = problo[2] + (k + 0.5_rt) * dx[2];
#endif
                if (x >= mask_lo[0] && x <= mask_hi[0]
                    && y >= mask_lo[1] && y <= mask_hi[1]
#if (AMREX_SPACEDIM == 2)
                    )
#else
                    && z >= mask_lo[2] && z <= mask_hi[2])
#endif
                {
                    mask(i,j,k) = 0;
                }
            });
        }
    }
}

void
MyTest::prepareOversetSolve ()
{
    resetPhi();

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.validbox();
            auto const& sol = phi[ilev].array(mfi);
            auto const& rhs_arr = rhs[ilev].array(mfi);
            auto const& ref = baseline_phi[ilev].const_array(mfi);
            auto const& mask = overset_mask[ilev].const_array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (mask(i,j,k) == 0) {
                    sol(i,j,k) = ref(i,j,k);
                    rhs_arr(i,j,k) = 0.0_rt;
                }
            });
        }
    }
}

void
MyTest::reportOversetError () const
{
    Real total_l1 = 0.0;
    Real total_vol = 0.0;
    Real total_max = 0.0;

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        auto const dx = geom[ilev].CellSizeArray();
        Real const dvol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
        const MultiFab& vfrc = factory[ilev]->getVolFrac();
        MultiFab err(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        MultiFab::Copy(err, phi[ilev], 0, 0, 1, 0);
        MultiFab::Subtract(err, baseline_phi[ilev], 0, 0, 1, 0);

        MultiFab weights(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        MultiFab::Copy(weights, vfrc, 0, 0, 1, 0);
        for (MFIter mfi(weights); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.validbox();
            auto const& w = weights.array(mfi);
            auto const& mask = overset_mask[ilev].const_array(mfi);
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (mask(i,j,k) == 0) {
                    w(i,j,k) = 0.0_rt;
                }
            });
        }

        Real const level_max = err.norm0(overset_mask[ilev], 0, 0, false);
        MultiFab::Multiply(err, weights, 0, 0, 1, 0);
        Real const level_l1_sum = err.norm1(0, 0, false) * dvol;
        Real const level_vol = weights.sum(0, false) * dvol;
        Real const level_l1 = (level_vol > 0.0) ? (level_l1_sum / level_vol) : 0.0;

        total_l1 += level_l1_sum;
        total_vol += level_vol;
        total_max = std::max(total_max, level_max);

        amrex::Print() << "Overset comparison level " << ilev
                       << ": L1 = " << level_l1
                       << " max_abs = " << level_max << '\n';
    }

    Real const global_l1 = (total_vol > 0.0) ? (total_l1 / total_vol) : 0.0;
    amrex::Print() << "Overset comparison total: L1 = " << global_l1
                   << " max_abs = " << total_max << '\n';
}
