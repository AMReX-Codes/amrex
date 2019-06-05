#include "MyTest.H"
#include "MyTest_F.H"

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

    addFineGrids();

    initData();
}

void
MyTest::solve ()
{
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        const MultiFab& vfrc = factory[ilev]->getVolFrac();
        MultiFab v(vfrc.boxArray(), vfrc.DistributionMap(), 1, 0,
                   MFInfo(), *factory[ilev]);
        MultiFab::Copy(v, vfrc, 0, 0, 1, 0);
        amrex::EB_set_covered(v, 1.0);
        amrex::Print() << "Level " << ilev << ": vfrc min = " << v.min(0) << std::endl;
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

    if (composite_solve)
    {
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
        
        if (true) { // In this test we assume EB is Dirichlet.
            for (int ilev = 0; ilev <= max_level; ++ilev) {
                mleb.setEBDirichlet(ilev, phieb[ilev], bcoef_eb[ilev]);
            }
        }

        MLMG mlmg(mleb);
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
        mlmg.solve(amrex::GetVecOfPtrs(phi), amrex::GetVecOfConstPtrs(rhs), tol_rel, tol_abs);
    }
    else
    {
        for (int ilev = 0; ilev <= max_level; ++ilev)
        {
            MLEBABecLap mleb({geom[ilev]}, {grids[ilev]}, {dmap[ilev]}, info, {factory[ilev].get()});
            mleb.setMaxOrder(linop_maxorder);

            mleb.setDomainBC(mlmg_lobc, mlmg_hibc);

            if (ilev > 0) {
                mleb.setCoarseFineBC(&phi[ilev-1], 2);
            }
            mleb.setLevelBC(0, &phi[ilev]);

            mleb.setScalars(scalars[0], scalars[1]);

            mleb.setACoeffs(0, acoef[ilev]);
            mleb.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef[ilev]));

            if (true) { // In this test we assume EB is Dirichlet.
                mleb.setEBDirichlet(0, phieb[ilev], bcoef_eb[ilev]);
            }

            MLMG mlmg(mleb);
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
            mlmg.solve({&phi[ilev]}, {&rhs[ilev]}, tol_rel, tol_abs);
        }
    }

    for (int ilev = 0; ilev <= max_level; ++ilev)
    {
        MultiFab mf(phi[ilev].boxArray(),phi[ilev].DistributionMap(), 1, 0);
        MultiFab::Copy(mf,phi[ilev],0,0,1,0);
        MultiFab::Subtract(mf, phiexact[ilev], 0, 0, 1, 0);

        const MultiFab& vfrc = factory[ilev]->getVolFrac();

        MultiFab::Multiply(mf, vfrc, 0, 0, 1, 0);

        Real norminf = mf.norm0();
        Real norm1 = mf.norm1()*AMREX_D_TERM((1.0/n_cell), *(1.0/n_cell), *(1.0/n_cell));
        amrex::Print() << "Level " << ilev << ": weighted max and 1 norms " << norminf << ", " << norm1 << std::endl;        
    }    
}

void
MyTest::writePlotfile ()
{
    Vector<MultiFab> plotmf(max_level+1);
    for (int ilev = 0; ilev <= max_level; ++ilev) {
        const MultiFab& vfrc = factory[ilev]->getVolFrac();
        plotmf[ilev].define(grids[ilev],dmap[ilev],5,0);

        MultiFab::Copy(plotmf[ilev], phi[ilev], 0, 0, 1, 0);

        MultiFab::Copy(plotmf[ilev], phiexact[ilev], 0, 1, 1, 0);

        MultiFab::Copy(plotmf[ilev], phi[ilev], 0, 2, 1, 0);
        MultiFab::Subtract(plotmf[ilev], phiexact[ilev], 0, 2, 1, 0);

        MultiFab::Copy(plotmf[ilev], phi[ilev], 0, 3, 1, 0);
        MultiFab::Subtract(plotmf[ilev], phiexact[ilev], 0, 3, 1, 0);
        MultiFab::Multiply(plotmf[ilev], vfrc, 0, 3, 1, 0);

        MultiFab::Copy(plotmf[ilev], vfrc, 0, 4, 1, 0);    
    }
    WriteMultiLevelPlotfile(plot_file_name, max_level+1,
                            amrex::GetVecOfConstPtrs(plotmf),
                            {"phi","exact","error","error*vfrac","vfrac"},
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
    pp.query("prob_type", prob_type);

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_level == 0 || max_level == 1,
                                     "max_level must be either 0 or 1");

    pp.query("plot_file", plot_file_name);

    scalars.resize(2);
    scalars[0] = 1.0;
    scalars[1] = 1.0;

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
    pp.query("use_petsc", use_petsc);
#endif

    pp.query("composite_solve", composite_solve);
}

void
MyTest::initGrids ()
{
    int nlevels = max_level + 1;
    geom.resize(nlevels);
    grids.resize(nlevels);

    RealBox rb({AMREX_D_DECL(0.,0.,0.)}, {AMREX_D_DECL(1.,1.,1.)});
    std::array<int,AMREX_SPACEDIM> isperiodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, isperiodic.data());
    Box domain0(IntVect{AMREX_D_DECL(0,0,0)}, IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    Box domain = domain0;
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        geom[ilev].define(domain);
        domain.refine(ref_ratio);
    }

    // Fine levels will be added later
    grids[0].define(domain0);
    grids[0].maxSize(max_grid_size);
}

void
MyTest::addFineGrids ()
{
    for (int ilev = 1; ilev <= max_level; ++ilev)
    {
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geom[ilev]);
        BoxList bl = eb_level.boxArray().boxList();
        const Box& domain = geom[ilev].Domain();
        for (Box& b : bl) {
            b &= domain;
        }
        grids[ilev].define(bl);
    }
}

void
MyTest::initData ()
{
    int nlevels = max_level + 1;
    dmap.resize(nlevels);
    factory.resize(nlevels);
    phi.resize(nlevels);
    phiexact.resize(nlevels);
    phieb.resize(nlevels);
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
        phiexact[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        phieb[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].define(amrex::convert(grids[ilev],IntVect::TheDimensionVector(idim)),
                                     dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        }
        bcoef_eb[ilev].define(grids[ilev], dmap[ilev], 1, 0, MFInfo(), *factory[ilev]);
        bcoef_eb[ilev].setVal(1.0);

        phi[ilev].setVal(0.0);
        rhs[ilev].setVal(0.0);
        acoef[ilev].setVal(0.0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].setVal(1.0);
        }

        const Real* dx = geom[ilev].CellSize();
        const Box& domainbox = geom[ilev].Domain();

        const FabArray<EBCellFlagFab>& flags = factory[ilev]->getMultiEBCellFlagFab();
        const MultiCutFab& bcent = factory[ilev]->getBndryCent();
        const MultiCutFab& cent = factory[ilev]->getCentroid();

        for (MFIter mfi(phiexact[ilev],true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const Box& xbx = mfi.nodaltilebox(0);
            const Box& ybx = mfi.nodaltilebox(1);
#if (AMREX_SPACEDIM == 3)
            const Box& zbx = mfi.nodaltilebox(2);
#endif
            auto fabtyp = flags[mfi].getType(bx);
            if (FabType::covered == fabtyp) {
                phiexact[ilev][mfi].setVal(0.0, bx, 0, 1);
                phieb[ilev][mfi].setVal(0.0, bx, 0, 1);
            } else if (FabType::regular == fabtyp) {
                phieb[ilev][mfi].setVal(0.0, bx, 0, 1);
                mytest_set_phi_reg(BL_TO_FORTRAN_BOX(bx),
                                   AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
                                                BL_TO_FORTRAN_BOX(ybx),
                                                BL_TO_FORTRAN_BOX(zbx)),
                                   BL_TO_FORTRAN_ANYD(phiexact[ilev][mfi]),
                                   BL_TO_FORTRAN_ANYD(rhs[ilev][mfi]),
                                   AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bcoef[ilev][0][mfi]),
                                                BL_TO_FORTRAN_ANYD(bcoef[ilev][1][mfi]),
                                                BL_TO_FORTRAN_ANYD(bcoef[ilev][2][mfi])),
                                   dx, &prob_type);
            } else {
                mytest_set_phi_eb(BL_TO_FORTRAN_BOX(bx),
                                  AMREX_D_DECL(BL_TO_FORTRAN_BOX(xbx),
                                               BL_TO_FORTRAN_BOX(ybx),
                                               BL_TO_FORTRAN_BOX(zbx)),
                                  BL_TO_FORTRAN_ANYD(phiexact[ilev][mfi]),
                                  BL_TO_FORTRAN_ANYD(phieb[ilev][mfi]),
                                  BL_TO_FORTRAN_ANYD(rhs[ilev][mfi]),
                                  AMREX_D_DECL(BL_TO_FORTRAN_ANYD(bcoef[ilev][0][mfi]),
                                               BL_TO_FORTRAN_ANYD(bcoef[ilev][1][mfi]),
                                               BL_TO_FORTRAN_ANYD(bcoef[ilev][2][mfi])),
                                  BL_TO_FORTRAN_ANYD(bcoef_eb[ilev][mfi]),
                                  BL_TO_FORTRAN_ANYD(flags[mfi]),
                                  BL_TO_FORTRAN_ANYD(cent[mfi]),
                                  BL_TO_FORTRAN_ANYD(bcent[mfi]),
                                  dx, &prob_type);
            }

            const Box& gbx = mfi.growntilebox(1);
            if (!domainbox.contains(gbx)) {
                mytest_set_phi_boundary(BL_TO_FORTRAN_BOX(gbx),
                                        BL_TO_FORTRAN_BOX(domainbox),
                                        BL_TO_FORTRAN_ANYD(phi[ilev][mfi]),
                                        dx);
            }
        }
    }
}

