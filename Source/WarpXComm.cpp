
#include <WarpX.H>
#include <WarpX_f.H>

#include <AMReX_FillPatchUtil_F.H>

#include <algorithm>
#include <cstdlib>

using namespace amrex;

void
WarpX::ExchangeWithPmlB (int lev)
{
    if (do_pml && pml[lev]->ok()) {
        pml[lev]->ExchangeB({ Bfield_fp[lev][0].get(),
                              Bfield_fp[lev][1].get(),
                              Bfield_fp[lev][2].get() },
                            { Bfield_cp[lev][0].get(),
                              Bfield_cp[lev][1].get(),
                              Bfield_cp[lev][2].get() });
    }
}

void
WarpX::ExchangeWithPmlE (int lev)
{
    if (do_pml && pml[lev]->ok()) {
        pml[lev]->ExchangeE({ Efield_fp[lev][0].get(),
                              Efield_fp[lev][1].get(),
                              Efield_fp[lev][2].get() },
                            { Efield_cp[lev][0].get(),
                              Efield_cp[lev][1].get(),
                              Efield_cp[lev][2].get() });
    }
}

void
WarpX::UpdateAuxilaryData ()
{
    BL_PROFILE("UpdateAuxilaryData()");

    const int use_limiter = 0;

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& crse_period = Geom(lev-1).periodicity();
        const IntVect& ng = Bfield_cp[lev][0]->nGrowVect();
        const DistributionMapping& dm = Bfield_cp[lev][0]->DistributionMap();

        // B field
        {
            MultiFab dBx(Bfield_cp[lev][0]->boxArray(), dm, 1, ng);
            MultiFab dBy(Bfield_cp[lev][1]->boxArray(), dm, 1, ng);
            MultiFab dBz(Bfield_cp[lev][2]->boxArray(), dm, 1, ng);
            dBx.setVal(0.0);
            dBy.setVal(0.0);
            dBz.setVal(0.0);
            dBx.copy(*Bfield_aux[lev-1][0], 0, 0, 1, ng, ng, crse_period);
            dBy.copy(*Bfield_aux[lev-1][1], 0, 0, 1, ng, ng, crse_period);
            dBz.copy(*Bfield_aux[lev-1][2], 0, 0, 1, ng, ng, crse_period);
            MultiFab::Subtract(dBx, *Bfield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dBy, *Bfield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dBz, *Bfield_cp[lev][2], 0, 0, 1, ng);
            
            const Real* dx = Geom(lev-1).CellSize();
            const int ref_ratio = refRatio(lev-1)[0];
#ifdef _OPENMP
#pragma omp parallel
#endif
            {
                std::array<FArrayBox,3> bfab;
                for (MFIter mfi(*Bfield_aux[lev][0]); mfi.isValid(); ++mfi)
                {
                    Box ccbx = mfi.fabbox();
                    ccbx.enclosedCells();
                    ccbx.coarsen(ref_ratio).refine(ref_ratio); // so that ccbx is coarsenable
                    
                    const FArrayBox& cxfab = dBx[mfi];
                    const FArrayBox& cyfab = dBy[mfi];
                    const FArrayBox& czfab = dBz[mfi];
                    bfab[0].resize(amrex::convert(ccbx,Bx_nodal_flag));
                    bfab[1].resize(amrex::convert(ccbx,By_nodal_flag));
                    bfab[2].resize(amrex::convert(ccbx,Bz_nodal_flag));

#if (AMREX_SPACEDIM == 3)
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(bfab[0]),
                                                 BL_TO_FORTRAN_ANYD(bfab[1]),
                                                 BL_TO_FORTRAN_ANYD(bfab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(cyfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &ref_ratio,&use_limiter);
#else
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(bfab[0]),
                                                 BL_TO_FORTRAN_ANYD(bfab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &ref_ratio,&use_limiter);
                    amrex_interp_cc_bfield(ccbx.loVect(), ccbx.hiVect(),
                                           BL_TO_FORTRAN_ANYD(bfab[1]),
                                           BL_TO_FORTRAN_ANYD(cyfab),
                                           &ref_ratio,&use_limiter);
#endif
                    
                    for (int idim = 0; idim < 3; ++idim)
                    {
                        FArrayBox& aux = (*Bfield_aux[lev][idim])[mfi];
                        FArrayBox& fp  =  (*Bfield_fp[lev][idim])[mfi];
                        const Box& bx = aux.box();
                        aux.copy(fp, bx, 0, bx, 0, 1);
                        aux.plus(bfab[idim], bx, bx, 0, 0, 1);
                    }
                }
            }
        }

        // E field
        {
            MultiFab dEx(Efield_cp[lev][0]->boxArray(), dm, 1, ng);
            MultiFab dEy(Efield_cp[lev][1]->boxArray(), dm, 1, ng);
            MultiFab dEz(Efield_cp[lev][2]->boxArray(), dm, 1, ng);
            dEx.setVal(0.0);
            dEy.setVal(0.0);
            dEz.setVal(0.0);
            dEx.copy(*Efield_aux[lev-1][0], 0, 0, 1, ng, ng, crse_period);
            dEy.copy(*Efield_aux[lev-1][1], 0, 0, 1, ng, ng, crse_period);
            dEz.copy(*Efield_aux[lev-1][2], 0, 0, 1, ng, ng, crse_period);
            MultiFab::Subtract(dEx, *Efield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dEy, *Efield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dEz, *Efield_cp[lev][2], 0, 0, 1, ng);
            
            const int ref_ratio = refRatio(lev-1)[0];
#ifdef _OPEMP
#pragma omp parallel
#endif
            {
                std::array<FArrayBox,3> efab;
                for (MFIter mfi(*Efield_aux[lev][0]); mfi.isValid(); ++mfi)
                {
                    Box ccbx = mfi.fabbox();
                    ccbx.enclosedCells();
                    ccbx.coarsen(ref_ratio).refine(ref_ratio); // so that ccbx is coarsenable
                    
                    const FArrayBox& cxfab = dEx[mfi];
                    const FArrayBox& cyfab = dEy[mfi];
                    const FArrayBox& czfab = dEz[mfi];
                    efab[0].resize(amrex::convert(ccbx,Ex_nodal_flag));
                    efab[1].resize(amrex::convert(ccbx,Ey_nodal_flag));
                    efab[2].resize(amrex::convert(ccbx,Ez_nodal_flag));

#if (AMREX_SPACEDIM == 3)
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(efab[0]),
                                        BL_TO_FORTRAN_ANYD(efab[1]),
                                        BL_TO_FORTRAN_ANYD(efab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(cyfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &ref_ratio,&use_limiter);
#else
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(efab[0]),
                                        BL_TO_FORTRAN_ANYD(efab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &ref_ratio,&use_limiter);
                    amrex_interp_nd_efield(ccbx.loVect(), ccbx.hiVect(),
                                           BL_TO_FORTRAN_ANYD(efab[1]),
                                           BL_TO_FORTRAN_ANYD(cyfab),
                                           &ref_ratio);
#endif
                    
                    for (int idim = 0; idim < 3; ++idim)
                    {
                        FArrayBox& aux = (*Efield_aux[lev][idim])[mfi];
                        FArrayBox& fp  =  (*Efield_fp[lev][idim])[mfi];
                        const Box& bx = aux.box();
                        aux.copy(fp, bx, 0, bx, 0, 1);
                        aux.plus(efab[idim], bx, bx, 0, 0, 1);
                    }
                }
            }
        }
    }
}

void
WarpX::FillBoundaryB ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryB(lev);
    }
}

void
WarpX::FillBoundaryE ()
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        FillBoundaryE(lev);
    }
}

void
WarpX::FillBoundaryE(int lev)
{
    const auto& period = Geom(lev).periodicity();

    if (do_pml && pml[lev]->ok()) {
        ExchangeWithPmlE(lev);
        pml[lev]->FillBoundaryE();
    }

    (*Efield_fp[lev][0]).FillBoundary( period );
    (*Efield_fp[lev][1]).FillBoundary( period );
    (*Efield_fp[lev][2]).FillBoundary( period );

    if (lev > 0)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        (*Efield_cp[lev][0]).FillBoundary(cperiod);
        (*Efield_cp[lev][1]).FillBoundary(cperiod);
        (*Efield_cp[lev][2]).FillBoundary(cperiod);
    }
}

void
WarpX::FillBoundaryB(int lev)
{
    const auto& period = Geom(lev).periodicity();

    if (do_pml && pml[lev]->ok())
    {
        ExchangeWithPmlB(lev);
        pml[lev]->FillBoundaryB();
    }

    (*Bfield_fp[lev][0]).FillBoundary(period);
    (*Bfield_fp[lev][1]).FillBoundary(period);
    (*Bfield_fp[lev][2]).FillBoundary(period);

    if (lev > 0)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        (*Bfield_cp[lev][0]).FillBoundary(cperiod);
        (*Bfield_cp[lev][1]).FillBoundary(cperiod);
        (*Bfield_cp[lev][2]).FillBoundary(cperiod);
    }
}

void
WarpX::SyncCurrent ()
{
    BL_PROFILE("SyncCurrent()");

    // Restrict fine patch current onto the coarse patch, before fine patch SumBoundary
    for (int lev = 1; lev <= finest_level; ++lev) 
    {
        current_cp[lev][0]->setVal(0.0);
        current_cp[lev][1]->setVal(0.0);
        current_cp[lev][2]->setVal(0.0);
      
        const IntVect& ref_ratio = refRatio(lev-1);

        std::array<const MultiFab*,3> fine { current_fp[lev][0].get(),
                                             current_fp[lev][1].get(),
                                             current_fp[lev][2].get() };
        std::array<      MultiFab*,3> crse { current_cp[lev][0].get(),
                                             current_cp[lev][1].get(),
                                             current_cp[lev][2].get() };
        SyncCurrent(fine, crse, ref_ratio[0]);
    }

    Vector<Array<std::unique_ptr<MultiFab>,3> > j_fp(finest_level+1);
    Vector<Array<std::unique_ptr<MultiFab>,3> > j_cp(finest_level+1);

    if (WarpX::use_filter) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            IntVect ng = current_fp[lev][0]->nGrowVect();
            ng += 1;
            j_fp[lev][0].reset(new MultiFab(current_fp[lev][0]->boxArray(),
                                            current_fp[lev][0]->DistributionMap(),
                                            1, ng));
            j_fp[lev][1].reset(new MultiFab(current_fp[lev][1]->boxArray(),
                                            current_fp[lev][1]->DistributionMap(),
                                            1, ng));
            j_fp[lev][2].reset(new MultiFab(current_fp[lev][2]->boxArray(),
                                            current_fp[lev][2]->DistributionMap(),
                                            1, ng));
            applyFilter(*j_fp[lev][0], *current_fp[lev][0]);
            applyFilter(*j_fp[lev][1], *current_fp[lev][1]);
            applyFilter(*j_fp[lev][2], *current_fp[lev][2]);
            std::swap(j_fp[lev][0], current_fp[lev][0]);
            std::swap(j_fp[lev][1], current_fp[lev][1]);
            std::swap(j_fp[lev][2], current_fp[lev][2]);
        }
        for (int lev = 1; lev <= finest_level; ++lev) {
            IntVect ng = current_cp[lev][0]->nGrowVect();
            ng += 1;
            j_cp[lev][0].reset(new MultiFab(current_cp[lev][0]->boxArray(),
                                            current_cp[lev][0]->DistributionMap(),
                                            1, ng));
            j_cp[lev][1].reset(new MultiFab(current_cp[lev][1]->boxArray(),
                                            current_cp[lev][1]->DistributionMap(),
                                            1, ng));
            j_cp[lev][2].reset(new MultiFab(current_cp[lev][2]->boxArray(),
                                            current_cp[lev][2]->DistributionMap(),
                                            1, ng));
            applyFilter(*j_cp[lev][0], *current_cp[lev][0]);
            applyFilter(*j_cp[lev][1], *current_cp[lev][1]);
            applyFilter(*j_cp[lev][2], *current_cp[lev][2]);
            std::swap(j_cp[lev][0], current_cp[lev][0]);
            std::swap(j_cp[lev][1], current_cp[lev][1]);
            std::swap(j_cp[lev][2], current_cp[lev][2]);
        }
    }

    // Sum up fine patch
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        current_fp[lev][0]->SumBoundary(period);
        current_fp[lev][1]->SumBoundary(period);
        current_fp[lev][2]->SumBoundary(period);
    }

    // Add fine level's coarse patch to coarse level's fine patch
    for (int lev = 0; lev < finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        const IntVect& ngsrc = current_cp[lev+1][0]->nGrowVect();
        const IntVect ngdst = IntVect::TheZeroVector();
        current_fp[lev][0]->copy(*current_cp[lev+1][0],0,0,1,ngsrc,ngdst,period,FabArrayBase::ADD);
        current_fp[lev][1]->copy(*current_cp[lev+1][1],0,0,1,ngsrc,ngdst,period,FabArrayBase::ADD);
        current_fp[lev][2]->copy(*current_cp[lev+1][2],0,0,1,ngsrc,ngdst,period,FabArrayBase::ADD);
    }

    // Sum up coarse patch
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        current_cp[lev][0]->SumBoundary(cperiod);
        current_cp[lev][1]->SumBoundary(cperiod);
        current_cp[lev][2]->SumBoundary(cperiod);
    }

    if (WarpX::use_filter) {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            std::swap(j_fp[lev][0], current_fp[lev][0]);
            std::swap(j_fp[lev][1], current_fp[lev][1]);
            std::swap(j_fp[lev][2], current_fp[lev][2]);
            MultiFab::Copy(*current_fp[lev][0], *j_fp[lev][0], 0, 0, 1, 0);
            MultiFab::Copy(*current_fp[lev][1], *j_fp[lev][1], 0, 0, 1, 0);
            MultiFab::Copy(*current_fp[lev][2], *j_fp[lev][2], 0, 0, 1, 0);
        }
        for (int lev = 1; lev <= finest_level; ++lev)
        {
            std::swap(j_cp[lev][0], current_cp[lev][0]);
            std::swap(j_cp[lev][1], current_cp[lev][1]);
            std::swap(j_cp[lev][2], current_cp[lev][2]);
            MultiFab::Copy(*current_cp[lev][0], *j_cp[lev][0], 0, 0, 1, 0);
            MultiFab::Copy(*current_cp[lev][1], *j_cp[lev][1], 0, 0, 1, 0);
            MultiFab::Copy(*current_cp[lev][2], *j_cp[lev][2], 0, 0, 1, 0);
        }
    }

    // sync shared nodal edges
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        current_fp[lev][0]->OverrideSync(period);
        current_fp[lev][1]->OverrideSync(period);
        current_fp[lev][2]->OverrideSync(period);
    }
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        current_cp[lev][0]->OverrideSync(cperiod);
        current_cp[lev][1]->OverrideSync(cperiod);
        current_cp[lev][2]->OverrideSync(cperiod);
    }
}

void
WarpX::SyncCurrent (const std::array<const amrex::MultiFab*,3>& fine,
                    const std::array<      amrex::MultiFab*,3>& crse,
                    int ref_ratio)
{
    BL_ASSERT(ref_ratio == 2);
    const IntVect& ng = fine[0]->nGrowVect()/ref_ratio;

#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        FArrayBox ffab;
        for (int idim = 0; idim < 3; ++idim)
        {
            for (MFIter mfi(*crse[idim],true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(ng);
                Box fbx = amrex::grow(amrex::refine(bx,ref_ratio),1);
                ffab.resize(fbx);
                fbx &= (*fine[idim])[mfi].box();
                ffab.setVal(0.0);
                ffab.copy((*fine[idim])[mfi], fbx, 0, fbx, 0, 1);
                WRPX_SYNC_CURRENT(bx.loVect(), bx.hiVect(),
                                   BL_TO_FORTRAN_ANYD((*crse[idim])[mfi]),
                                   BL_TO_FORTRAN_ANYD(ffab),
                                   &idim);
            }
        }
    }
}

void
WarpX::SyncRho (amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rhof,
                amrex::Vector<std::unique_ptr<amrex::MultiFab> >& rhoc)
{
    if (!rhof[0]) return;

    // Restrict fine patch onto the coarse patch, before fine patch SumBoundary
    for (int lev = 1; lev <= finest_level; ++lev) 
    {
        rhoc[lev]->setVal(0.0);      
        const IntVect& ref_ratio = refRatio(lev-1);
        SyncRho(*rhof[lev], *rhoc[lev], ref_ratio[0]);
    }

    Vector<std::unique_ptr<MultiFab> > rho_f_g(finest_level+1);
    Vector<std::unique_ptr<MultiFab> > rho_c_g(finest_level+1);

    if (WarpX::use_filter) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            const int ncomp = rhof[lev]->nComp();
            IntVect ng = rhof[lev]->nGrowVect();
            ng += 1;
            rho_f_g[lev].reset(new MultiFab(rhof[lev]->boxArray(),
                                            rhof[lev]->DistributionMap(),
                                            ncomp, ng));
            applyFilter(*rho_f_g[lev], *rhof[lev]);
            std::swap(rho_f_g[lev], rhof[lev]);
        }
        for (int lev = 1; lev <= finest_level; ++lev) {
            const int ncomp = rhoc[lev]->nComp();
            IntVect ng = rhoc[lev]->nGrowVect();
            ng += 1;
            rho_c_g[lev].reset(new MultiFab(rhoc[lev]->boxArray(),
                                            rhoc[lev]->DistributionMap(),
                                            ncomp, ng));
            applyFilter(*rho_c_g[lev], *rhoc[lev]);
            std::swap(rho_c_g[lev], rhoc[lev]);
        }
    }

    // Sum up fine patch
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        rhof[lev]->SumBoundary(period);
    }

    // Add fine level's coarse patch to coarse level's fine patch
    for (int lev = 0; lev < finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        const int ncomp = rhoc[lev+1]->nComp();
        const IntVect& ngsrc = rhoc[lev+1]->nGrowVect();
        const IntVect ngdst = IntVect::TheZeroVector();
        rhof[lev]->copy(*rhoc[lev+1],0,0,ncomp,ngsrc,ngdst,period,FabArrayBase::ADD);
    }

    // Sum up coarse patch
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        rhoc[lev]->SumBoundary(cperiod);
    }

    if (WarpX::use_filter) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            std::swap(rho_f_g[lev], rhof[lev]);
            MultiFab::Copy(*rhof[lev], *rho_f_g[lev], 0, 0, rhof[lev]->nComp(), 0);
        }
        for (int lev = 1; lev <= finest_level; ++lev) {
            std::swap(rho_c_g[lev], rhoc[lev]);
            MultiFab::Copy(*rhoc[lev], *rho_c_g[lev], 0, 0, rhoc[lev]->nComp(), 0);
        }
    }

    // sync shared nodal points
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        rhof[lev]->OverrideSync(period);
    }
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& cperiod = Geom(lev-1).periodicity();
        rhoc[lev]->OverrideSync(cperiod);
    }
}

void
WarpX::SyncRho (const MultiFab& fine, MultiFab& crse, int ref_ratio)
{
    BL_ASSERT(ref_ratio == 2);
    const IntVect& ng = fine.nGrowVect()/ref_ratio;
    const int nc = fine.nComp();

#ifdef _OPEMP
#pragma omp parallel
#endif
    {
        FArrayBox ffab;
        for (MFIter mfi(crse,true); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(ng);
            Box fbx = amrex::grow(amrex::refine(bx,ref_ratio),1);
            ffab.resize(fbx, nc);
            fbx &= fine[mfi].box();
            ffab.setVal(0.0);
            ffab.copy(fine[mfi], fbx, 0, fbx, 0, nc);
            WRPX_SYNC_RHO(bx.loVect(), bx.hiVect(),
                          BL_TO_FORTRAN_ANYD(crse[mfi]),
                          BL_TO_FORTRAN_ANYD(ffab),
                          &nc);
        }
    }
}

