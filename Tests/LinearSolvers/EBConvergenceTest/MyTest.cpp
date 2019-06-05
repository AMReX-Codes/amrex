#include "MyTest.H"
#include "MyTest_F.H"
    
#include <AMReX_MLEBABecLap.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
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

    MLEBABecLap mleb (geom, grids, dmap, info, amrex::GetVecOfConstPtrs(factory));

    mleb.setDomainBC(mlmg_lobc, mlmg_hibc);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mleb.setLevelBC(ilev, &phi[ilev]);
    }

    mleb.setScalars(1.0, 1.0);

    for (int ilev = 0; ilev <= max_level; ++ilev) {
        mleb.setACoeffs(ilev, acoef[ilev]);
        mleb.setBCoeffs(ilev, amrex::GetArrOfConstPtrs(bcoef[ilev]));
    }

    MLMG mlmg(mleb);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setBottomMaxIter(max_bottom_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);
    if (use_hypre) mlmg.setBottomSolver(MLMG::BottomSolver::hypre);

    const Real tol_rel = 1.e-12;
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

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("max_bottom_iter", max_bottom_iter);
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
    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
//    std::array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1,0,0)};
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

        phi[ilev].setVal(0.0);
        rhs[ilev].setVal(0.0);
        acoef[ilev].setVal(1.0);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[ilev][idim].setVal(1.0);
        }

        // Construct Manufactured Solution RHS
        const Real* dx = geom[ilev].CellSize();
#if (AMREX_SPACEDIM == 2)
        for (MFIter mfi(rhs[ilev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& rhsfab = rhs[ilev][mfi]; 
            FArrayBox& afab   = acoef[ilev][mfi]; 
            FArrayBox& bfabx  = bcoef[ilev][0][mfi]; 
            FArrayBox& bfaby  = bcoef[ilev][1][mfi];
            const Box& bx = rhsfab.box();
            build_a_2d(bx.loVect(), bx.hiVect(), 
                       geom[ilev].ProbLo(), 
                       geom[ilev].ProbHi(), 
                       BL_TO_FORTRAN_ANYD(afab), 
                       dx); 
                       
             build_b_2d(bx.loVect(), bx.hiVect(), 
                       geom[ilev].ProbLo(), 
                       geom[ilev].ProbHi(), 
                       BL_TO_FORTRAN_ANYD( bfabx),
                       BL_TO_FORTRAN_ANYD( bfaby),
                       dx); // */
            
            build_rhs_2d(bx.loVect(), bx.hiVect(), 
                      BL_TO_FORTRAN_ANYD(rhsfab), 
                      BL_TO_FORTRAN_ANYD(  afab), 
                      BL_TO_FORTRAN_ANYD( bfabx),
                      BL_TO_FORTRAN_ANYD( bfaby),
                      dx, geom[ilev].ProbLo(), 
                      geom[ilev].ProbHi()); 
        }
#elif (AMREX_SPACEDIM ==3)
        for (MFIter mfi(rhs[ilev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& rhsfab = rhs[ilev][mfi]; 
            FArrayBox& afab   = acoef[ilev][mfi]; 
            FArrayBox& bfabx  = bcoef[ilev][0][mfi]; 
            FArrayBox& bfaby  = bcoef[ilev][1][mfi];
            FArrayBox& bfabz  = bcoef[ilev][2][mfi];
            const Box& bx = rhsfab.box();

            build_a_3d(bx.loVect(), bx.hiVect(), 
                       geom[ilev].ProbLo(), 
                       geom[ilev].ProbHi(), 
                       BL_TO_FORTRAN_ANYD(afab), 
                       dx); 
                       
             build_b_3d(bx.loVect(), bx.hiVect(), 
                       geom[ilev].ProbLo(), 
                       geom[ilev].ProbHi(), 
                       BL_TO_FORTRAN_ANYD( bfabx),
                       BL_TO_FORTRAN_ANYD( bfaby),
                       BL_TO_FORTRAN_ANYD( bfabz),
                       dx); // */

 
            build_rhs_3d(bx.loVect(), bx.hiVect(), 
                      BL_TO_FORTRAN_ANYD(rhsfab), 
                      BL_TO_FORTRAN_ANYD(  afab), 
                      BL_TO_FORTRAN_ANYD( bfabx),
                      BL_TO_FORTRAN_ANYD( bfaby),
                      BL_TO_FORTRAN_ANYD( bfabz), 
                      dx, geom[ilev].ProbLo(), 
                      geom[ilev].ProbHi()); 
        }

#endif
        // Initialize Dirichlet boundary
        for (MFIter mfi(phi[ilev]); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = phi[ilev][mfi];
            const Box& bx = fab.box();
            apply_bc(bx.loVect(), bx.hiVect(), 
                      BL_TO_FORTRAN_ANYD(fab),
                      dx, geom[ilev].ProbLo(), 
                      geom[ilev].ProbHi());// */ 
        }

        phi[ilev].setVal(0.0, 0, 1, 0);
    }
}

