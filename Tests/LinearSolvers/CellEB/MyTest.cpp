#include "MyTest.H"

#include <AMReX_MLEBABecLap.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_EB2.H>

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
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(factory[ilev]->getVolFrac(), "vfrc-"+std::to_string(ilev));

        const MultiFab& vfrc = factory[ilev]->getVolFrac();
        MultiFab v(vfrc.boxArray(), vfrc.DistributionMap(), 1, 0,
                   MFInfo(), *factory[ilev]);
        MultiFab::Copy(v, vfrc, 0, 0, 1, 0);
        amrex::EB_set_covered(v, 1.0);
        amrex::Print() << "vfrc min = " << v.min(0) << std::endl;
    }

    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_lobc;
    std::array<LinOpBCType,AMREX_SPACEDIM> mlmg_hibc;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            mlmg_lobc[idim] = LinOpBCType::Periodic;
            mlmg_hibc[idim] = LinOpBCType::Periodic;
        } else {
            mlmg_lobc[idim] = LinOpBCType::Dirichlet;
            mlmg_hibc[idim] = LinOpBCType::Dirichlet;
        }
    }

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);

    MLEBABecLap mleb (geom, grids, dmap, info, amrex::GetVecOfConstPtrs(factory));
    mleb.setMaxOrder(linop_maxorder);

    mleb.setDomainBC(mlmg_lobc, mlmg_hibc);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mleb.setLevelBC(ilev, &phi[ilev]);
    }

    mleb.setScalars(scalars[0], scalars[1]);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mleb.setACoeffs(ilev, acoef[ilev]);
        mleb.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcoef[ilev]));
    }

    if (eb_is_dirichlet) {
        for (int ilev = 0; ilev <= max_level; ++ilev) {
            mleb.setEBDirichlet(ilev, phi[ilev], bcoef_eb[ilev]);
        }
    }

    MLMG mlmg(mleb);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setBottomMaxIter(max_bottom_iter);
    mlmg.setBottomTolerance(bottom_reltol);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    if (use_hypre) mlmg.setBottomSolver(MLMG::BottomSolver::hypre);

    const Real tol_rel = reltol;
    const Real tol_abs = 0.0;
    mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        amrex::VisMF::Write(phi[0], "phi-"+std::to_string(ilev));
    }
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
    acoef.resize(nlevels);
    bcoef.resize(nlevels);
    bcoef_eb.resize(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        dmap[ilev].define(grids[ilev]);
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        factory[ilev].reset(new EBFArrayBoxFactory(eb_level, geom[ilev], grids[ilev], dmap[ilev],
                                                   {2,2,2}, EBSupport::full));

        phi[ilev].define(grids[ilev], dmap[ilev], 1, 1, MFInfo(), *factory[ilev]);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
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

        const Real* dx = geom[ilev].CellSize();

        if (is_periodic)
        {
            constexpr Real pi = 4.0*std::atan(1.0);

            for (MFIter mfi(rhs[ilev]); mfi.isValid(); ++mfi)
            {
                FArrayBox& fab = rhs[ilev][mfi];
                const Box& bx = fab.box();
                fab.ForEachIV(bx, 0, 1, [=] (Real& p, const IntVect& iv) {
                        Real rx = (iv[0]+0.5)*dx[0];
                        Real ry = (iv[1]+0.5)*dx[1];
                        p = std::sin(rx*2.*pi + 43.5)*std::sin(ry*2.*pi + 89.);
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
                FArrayBox& fab = phi[ilev][mfi];
                const Box& bx = fab.box();
                fab.ForEachIV(bx, 0, 1, [=] (Real& p, const IntVect& iv) {
                        Real rx = (iv[0]+0.5)*dx[0];
                        Real ry = (iv[1]+0.5)*dx[1];
                        p = std::sqrt(0.5)*(rx + ry);
                    });
            }
            phi[ilev].setVal(0.0, 0, 1, 0); // set interior to zero
        }

    }
}

