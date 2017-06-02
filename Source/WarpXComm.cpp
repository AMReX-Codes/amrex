
#include <WarpX.H>
#include <WarpX_f.H>

#include <AMReX_FillPatchUtil_F.H>

#include <algorithm>
#include <cstdlib>

using namespace amrex;

void
WarpX::ExchangeWithPML ()
{
    if (do_pml) {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            ExchangeWithPML(lev);
        }
    }
}

void
WarpX::ExchangeWithPML (int lev)
{
    if (do_pml) {
        pml[lev]->Exchange({ Efield_fp[lev][0].get(),
                             Efield_fp[lev][1].get(),
                             Efield_fp[lev][2].get() },
                           { Bfield_fp[lev][0].get(),
                             Bfield_fp[lev][1].get(),
                             Bfield_fp[lev][2].get() },
                           { Efield_cp[lev][0].get(),
                             Efield_cp[lev][1].get(),
                             Efield_cp[lev][2].get() },
                           { Bfield_cp[lev][0].get(),
                             Bfield_cp[lev][1].get(),
                             Bfield_cp[lev][2].get() });
    }
}

void
WarpX::UpdateAuxiliaryData ()
{
    if (do_pml) {
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            ExchangeWithPML(lev);
            pml[lev]->FillBoundary();
        }
    }

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        Bfield_fp[lev][0]->FillBoundary(period);
        Bfield_fp[lev][1]->FillBoundary(period);
        Bfield_fp[lev][2]->FillBoundary(period);
        Efield_fp[lev][0]->FillBoundary(period);
        Efield_fp[lev][1]->FillBoundary(period);
        Efield_fp[lev][2]->FillBoundary(period);
    }

    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& crse_period = Geom(lev-1).periodicity();
        const int ng = Bfield_cp[lev][0]->nGrow();
        const DistributionMapping& dm = Bfield_cp[lev][0]->DistributionMap();

        // B field
        {
            MultiFab dBx(Bfield_cp[lev][0]->boxArray(), dm, 1, ng);
            MultiFab dBy(Bfield_cp[lev][1]->boxArray(), dm, 1, ng);
            MultiFab dBz(Bfield_cp[lev][2]->boxArray(), dm, 1, ng);
            dBx.copy(*Bfield_aux[lev-1][0], 0, 0, 1, 0, ng, crse_period);
            dBy.copy(*Bfield_aux[lev-1][1], 0, 0, 1, 0, ng, crse_period);
            dBz.copy(*Bfield_aux[lev-1][2], 0, 0, 1, 0, ng, crse_period);
            Bfield_cp[lev][0]->FillBoundary(crse_period);
            Bfield_cp[lev][1]->FillBoundary(crse_period);
            Bfield_cp[lev][2]->FillBoundary(crse_period);
            MultiFab::Subtract(dBx, *Bfield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dBy, *Bfield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dBz, *Bfield_cp[lev][2], 0, 0, 1, ng);
            
            const Real* dx = Geom(lev-1).CellSize();
            const int ref_ratio = refRatio(lev-1)[0];
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

#if (BL_SPACEDIM == 3)
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(bfab[0]),
                                                 BL_TO_FORTRAN_ANYD(bfab[1]),
                                                 BL_TO_FORTRAN_ANYD(bfab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(cyfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &ref_ratio);
#else
                    amrex_interp_div_free_bfield(ccbx.loVect(), ccbx.hiVect(),
                                                 BL_TO_FORTRAN_ANYD(bfab[0]),
                                                 BL_TO_FORTRAN_ANYD(bfab[2]),
                                                 BL_TO_FORTRAN_ANYD(cxfab),
                                                 BL_TO_FORTRAN_ANYD(czfab),
                                                 dx, &ref_ratio);
                    amrex_interp_cc_bfield(ccbx.loVect(), ccbx.hiVect(),
                                           BL_TO_FORTRAN_ANYD(bfab[1]),
                                           BL_TO_FORTRAN_ANYD(cyfab),
                                           &ref_ratio);
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
            dEx.copy(*Efield_aux[lev-1][0], 0, 0, 1, 0, ng, crse_period);
            dEy.copy(*Efield_aux[lev-1][1], 0, 0, 1, 0, ng, crse_period);
            dEz.copy(*Efield_aux[lev-1][2], 0, 0, 1, 0, ng, crse_period);
            Efield_cp[lev][0]->FillBoundary(crse_period);
            Efield_cp[lev][1]->FillBoundary(crse_period);
            Efield_cp[lev][2]->FillBoundary(crse_period);
            MultiFab::Subtract(dEx, *Efield_cp[lev][0], 0, 0, 1, ng);
            MultiFab::Subtract(dEy, *Efield_cp[lev][1], 0, 0, 1, ng);
            MultiFab::Subtract(dEz, *Efield_cp[lev][2], 0, 0, 1, ng);
            
            const int ref_ratio = refRatio(lev-1)[0];
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

#if (BL_SPACEDIM == 3)
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(efab[0]),
                                        BL_TO_FORTRAN_ANYD(efab[1]),
                                        BL_TO_FORTRAN_ANYD(efab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(cyfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &ref_ratio);
#else
                    amrex_interp_efield(ccbx.loVect(), ccbx.hiVect(),
                                        BL_TO_FORTRAN_ANYD(efab[0]),
                                        BL_TO_FORTRAN_ANYD(efab[2]),
                                        BL_TO_FORTRAN_ANYD(cxfab),
                                        BL_TO_FORTRAN_ANYD(czfab),
                                        &ref_ratio);
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

    if (do_pml) {
        ExchangeWithPML(lev);
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

    if (do_pml)
    {
        ExchangeWithPML(lev);
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
WarpX::AverageDownB ()
{
// xxxxx
#if 0
    for (int lev = finest_level; lev > 0; --lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);
        const Geometry& crse_geom = Geom(lev-1);
#if (BL_SPACEDIM == 3)
        Array<const MultiFab*> fine {Bfield[lev][0].get(), Bfield[lev][1].get(), Bfield[lev][2].get()};
#else
        Array<const MultiFab*> fine {Bfield[lev][0].get(), Bfield[lev][2].get()};
#endif
        Array<std::unique_ptr<MultiFab> > crse(BL_SPACEDIM);
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            BoxArray cba = fine[idim]->boxArray();
            const DistributionMapping& dm = fine[idim]->DistributionMap();
            cba.coarsen(ref_ratio);
            crse[idim].reset(new MultiFab(cba, dm, 1, 0));
        }

        amrex::average_down_faces(fine, amrex::GetArrOfPtrs(crse), ref_ratio);

#if (BL_SPACEDIM == 3)
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            Bfield[lev-1][idim]->copy(*crse[idim], crse_geom.periodicity());
        }
#else
        Bfield[lev-1][0]->copy(*crse[0], crse_geom.periodicity());
        Bfield[lev-1][2]->copy(*crse[1], crse_geom.periodicity());
        amrex::average_down(*Bfield[lev][1], *Bfield[lev-1][1], 0, 1, ref_ratio);
#endif        
    }
#endif
}

void
WarpX::AverageDownE ()
{
//xxxxx    AverageDownEorCurrent(Efield);
}

void
WarpX::AverageDownCurrent ()
{
//xxxxx    AverageDownEorCurrent(current);
}

void
WarpX::AverageDownEorCurrent (Array<std::array<std::unique_ptr<amrex::MultiFab>,3> >& data)
{
    for (int lev = finest_level; lev > 0; --lev)
    {
        const IntVect& ref_ratio = refRatio(lev-1);
        const auto& crse_period = Geom(lev-1).periodicity();
#if (BL_SPACEDIM == 3)
        Array<const MultiFab*> fine {data[lev][0].get(), data[lev][1].get(), data[lev][2].get()};
#else
        Array<const MultiFab*> fine {data[lev][0].get(), data[lev][2].get()};
#endif
        
        Array<std::unique_ptr<MultiFab> > crse(BL_SPACEDIM);
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            BoxArray cba = fine[idim]->boxArray();
            const DistributionMapping& dm = fine[idim]->DistributionMap();
            cba.coarsen(ref_ratio);
            crse[idim].reset(new MultiFab(cba, dm, 1, 0));
        }

        amrex::average_down_edges(fine, amrex::GetArrOfPtrs(crse), ref_ratio);

#if (BL_SPACEDIM == 3)
        for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
            data[lev-1][idim]->copy(*crse[idim], crse_period);
        }
#else
        data[lev-1][0]->copy(*crse[0], crse_period);
        data[lev-1][2]->copy(*crse[1], crse_period);
//xxxxx amrex::average_down_nodes(*data[lev][1], *data[lev-1][1], 0, 1, ref_ratio);
#endif        
    }
}

void
WarpX::SyncCurrent ()
{
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
#if (BL_SPACEDIM == 3)
        SyncCurrent3D(fine, crse, ref_ratio[0]);
#else
        SyncCurrent2D(fine, crse, ref_ratio[0]);
#endif
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
        const int ngsrc = current_cp[lev+1][0]->nGrow();
        const int ngdst = 0;
        current_fp[lev][0]->copy(*current_cp[lev+1][0],0,0,1,ngsrc,ngdst,period,FabArrayBase::ADD);
        current_fp[lev][1]->copy(*current_cp[lev+1][1],0,0,1,ngsrc,ngdst,period,FabArrayBase::ADD);
        current_fp[lev][2]->copy(*current_cp[lev+1][2],0,0,1,ngsrc,ngdst,period,FabArrayBase::ADD);
    }

    // Sum up coarse patch
    for (int lev = 1; lev <= finest_level; ++lev)
    {
        const auto& period = Geom(lev).periodicity();
        current_cp[lev][0]->SumBoundary(period);
        current_cp[lev][1]->SumBoundary(period);
        current_cp[lev][2]->SumBoundary(period);
    }
}

void
WarpX::SyncCurrent3D (const std::array<const amrex::MultiFab*,3>& fine,
                      const std::array<      amrex::MultiFab*,3>& crse,
                      int ref_ratio)
{
    BL_ASSERT(ref_ratio == 2);
    int ng = fine[0]->nGrow()/ref_ratio;

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
                WARPX_SYNC_CURRENT(bx.loVect(), bx.hiVect(),
                                   BL_TO_FORTRAN_ANYD((*crse[idim])[mfi]),
                                   BL_TO_FORTRAN_ANYD(ffab),
                                   &idim);
            }
        }
    }
}

void
WarpX::SyncCurrent2D (const std::array<const amrex::MultiFab*,3>& fine,
                      const std::array<      amrex::MultiFab*,3>& crse,
                      int ref_ratio)
{
    BL_ASSERT(ref_ratio == 2);
    int ng = fine[0]->nGrow()/ref_ratio;

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
                WARPX_SYNC_CURRENT(bx.loVect(), bx.hiVect(),
                                   BL_TO_FORTRAN_ANYD((*crse[idim])[mfi]),
                                   BL_TO_FORTRAN_ANYD(ffab),
                                   &idim);
            }
        }
    }
}

