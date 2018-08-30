
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil_F.H>
#include <AMReX_MultiFabUtil_F.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex
{

void
EB_set_covered (MultiFab& mf, Real val)
{
    Vector<Real> vals(mf.nComp(), val);
    EB_set_covered(mf, 0, mf.nComp(), vals);
}

void
EB_set_covered (MultiFab& mf, int icomp, int ncomp, const Vector<Real>& vals)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];
        const auto& flagfab = amrex::getEBCellFlagFab(fab);
        amrex_eb_set_covered(BL_TO_FORTRAN_BOX(bx),
                             BL_TO_FORTRAN_N_ANYD(fab,icomp),
                             BL_TO_FORTRAN_ANYD(flagfab),
                             vals.data(),&ncomp);
    }
}

void
EB_set_covered_faces (const Array<MultiFab*,AMREX_SPACEDIM>& umac, Real val)
{
    const auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(umac[0]->Factory()));
    if (factory == nullptr) return;

    const auto& area = factory->getAreaFrac();
    const auto& flags = factory->getMultiEBCellFlagFab();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(*umac[0],true); mfi.isValid(); ++mfi)
    {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
        {
            const Box& bx = mfi.tilebox(IntVect::TheDimensionVector(idim));
            auto fabtyp = flags[mfi].getType(amrex::enclosedCells(bx));
            if (fabtyp == FabType::covered) {
                (*umac[idim])[mfi].setVal(0.0, bx, 0, 1);
            } else if (fabtyp != FabType::regular) {
                amrex_eb_set_covered_faces(BL_TO_FORTRAN_BOX(bx),
                                           BL_TO_FORTRAN_ANYD((*umac[idim])[mfi]),
                                           BL_TO_FORTRAN_ANYD((*area[idim])[mfi]));
            }
        }
    }
}

void
EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, const MultiFab& vol_fine,
                 const MultiFab& vfrac_fine, int scomp, int ncomp, const IntVect& ratio)
{
    BL_PROFILE("EB_average_down");

    BL_ASSERT(S_fine.ixType().cellCentered());
    BL_ASSERT(S_crse.ixType().cellCentered());

    const DistributionMapping& fine_dm = S_fine.DistributionMap();
    BoxArray crse_S_fine_BA = S_fine.boxArray();
    crse_S_fine_BA.coarsen(ratio);

    
    MultiFab crse_S_fine(crse_S_fine_BA,fine_dm,ncomp,0,MFInfo(),FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        auto& crse_fab = crse_S_fine[mfi];
        const auto& fine_fab = S_fine[mfi];

        const auto& flag_fab = amrex::getEBCellFlagFab(fine_fab);
        FabType typ = flag_fab.getType(amrex::refine(tbx,ratio));
        
        if (typ == FabType::regular || typ == FabType::covered)
        {
#if (AMREX_SPACEDIM == 3)
            BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(fine_fab,scomp),
                 BL_TO_FORTRAN_N(crse_fab,0),
                 ratio.getVect(),&ncomp);
#else
            BL_FORT_PROC_CALL(BL_AVGDOWN_WITH_VOL,bl_avgdown_with_vol)
                (tbx.loVect(), tbx.hiVect(),
                 BL_TO_FORTRAN_N(fine_fab,scomp),
                 BL_TO_FORTRAN_N(crse_fab,0),
                 BL_TO_FORTRAN(vol_fine[mfi]),
                 ratio.getVect(),&ncomp);
#endif
        }
        else if (typ == FabType::singlevalued)
        {
#if (AMREX_SPACEDIM == 1)
            amrex::Abort("1D EB not supported");
#else
            amrex_eb_avgdown_sv(BL_TO_FORTRAN_BOX(tbx),
                                BL_TO_FORTRAN_N_ANYD(fine_fab,scomp),
                                BL_TO_FORTRAN_N_ANYD(crse_fab,0),
                                BL_TO_FORTRAN_ANYD(vol_fine[mfi]),
                                BL_TO_FORTRAN_ANYD(vfrac_fine[mfi]),
                                ratio.getVect(),&ncomp);
#endif
        }
        else
        {
            amrex::Abort("multi-valued avgdown to be implemented");
        }
    }

    S_crse.copy(crse_S_fine,0,scomp,ncomp);
}


void
EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp, int ratio)
{
    EB_average_down(S_fine, S_crse, scomp, ncomp, IntVect(ratio));
}

void
EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp, const IntVect& ratio)
{
    if (!S_fine.hasEBFabFactory())
    {
        amrex::average_down(S_fine, S_crse, scomp, ncomp, ratio);
    }
    else
    {
        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(S_fine.Factory());
        const auto& vfrac_fine = factory.getVolFrac();

        BL_ASSERT(S_crse.nComp() == S_fine.nComp());
        BL_ASSERT(S_crse.is_cell_centered() && S_fine.is_cell_centered());

        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        if (crse_S_fine_BA == S_crse.boxArray() 
            and S_fine.DistributionMap() == S_crse.DistributionMap())
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(S_crse,true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();
                auto& crse_fab = S_crse[mfi];
                const auto& fine_fab = S_fine[mfi];

                const auto& flag_fab = amrex::getEBCellFlagFab(fine_fab);
                FabType typ = flag_fab.getType(amrex::refine(tbx,ratio));

                if (typ == FabType::regular || typ == FabType::covered)
                {
                    BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                        (tbx.loVect(), tbx.hiVect(),
                         BL_TO_FORTRAN_N(fine_fab,scomp),
                         BL_TO_FORTRAN_N(crse_fab,scomp),
                         ratio.getVect(),&ncomp);
                }
                else
                {
                    amrex_eb_avgdown(BL_TO_FORTRAN_BOX(tbx),
                                     BL_TO_FORTRAN_N_ANYD(fine_fab,scomp),
                                     BL_TO_FORTRAN_N_ANYD(crse_fab,scomp),
                                     BL_TO_FORTRAN_ANYD(vfrac_fine[mfi]),
                                     ratio.getVect(),&ncomp);
                }
            }
        }
        else
        {
            MultiFab crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(),
                                 ncomp, 0, MFInfo(),FArrayBoxFactory());

#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.tilebox();
                auto& crse_fab = crse_S_fine[mfi];
                const auto& fine_fab = S_fine[mfi];
                
                const auto& flag_fab = amrex::getEBCellFlagFab(fine_fab);
                FabType typ = flag_fab.getType(amrex::refine(tbx,ratio));
                
                if (typ == FabType::regular || typ == FabType::covered)
                {
                    BL_FORT_PROC_CALL(BL_AVGDOWN,bl_avgdown)
                        (tbx.loVect(), tbx.hiVect(),
                         BL_TO_FORTRAN_N(fine_fab,scomp),
                         BL_TO_FORTRAN_N(crse_fab,0),
                         ratio.getVect(),&ncomp);
                }
                else if (typ == FabType::singlevalued)
                {
                    amrex_eb_avgdown(BL_TO_FORTRAN_BOX(tbx),
                                     BL_TO_FORTRAN_N_ANYD(fine_fab,scomp),
                                     BL_TO_FORTRAN_N_ANYD(crse_fab,0),
                                     BL_TO_FORTRAN_ANYD(vfrac_fine[mfi]),
                                     ratio.getVect(),&ncomp);
                }
                else
                {
                    amrex::Abort("multi-valued avgdown to be implemented");
                }
            }

            S_crse.copy(crse_S_fine,0,scomp,ncomp);
        }
    }
}


void EB_average_down_faces (const Vector<const MultiFab*>& fine, const Vector< MultiFab*>& crse,
                            const IntVect& ratio, int ngcrse)
{
    BL_ASSERT(crse.size()  == AMREX_SPACEDIM);
    BL_ASSERT(fine.size()  == AMREX_SPACEDIM);
    BL_ASSERT(crse[0]->nComp() == fine[0]->nComp());

    int ncomp = crse[0]->nComp();
    if (!(*fine[0]).hasEBFabFactory())
    {
        amrex::average_down_faces(fine, crse, ratio, ngcrse);
    }
    else 
    {
        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>((*fine[0]).Factory());
        const auto&  aspect = factory.getAreaFrac();

        if (isMFIterSafe(*fine[0], *crse[0]))
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (int n=0; n<AMREX_SPACEDIM; ++n) {
                for (MFIter mfi(*crse[n],true); mfi.isValid(); ++mfi)
                {
                    const auto& flag_fab = amrex::getEBCellFlagFab((*fine[n])[mfi]);
                    const Box& tbx = mfi.growntilebox(ngcrse);
                    FabType typ = flag_fab.getType(amrex::refine(tbx,ratio));
               
                    if(typ == FabType::regular || typ == FabType::covered) 
                    {    

                        BL_FORT_PROC_CALL(BL_AVGDOWN_FACES,bl_avgdown_faces)
                            (tbx.loVect(),tbx.hiVect(),
                             BL_TO_FORTRAN((*fine[n])[mfi]),
                             BL_TO_FORTRAN((*crse[n])[mfi]),
                             ratio.getVect(),n, ncomp);
                    }
                    else
                    {
                       amrex_eb_avgdown_faces(tbx.loVect(), tbx.hiVect(), 
                                              BL_TO_FORTRAN_ANYD((*fine[n])[mfi]), 
                                              BL_TO_FORTRAN_ANYD((*crse[n])[mfi]),
                                              BL_TO_FORTRAN_ANYD((*aspect[n])[mfi]),
                                              ratio.getVect(), &n, &ncomp);
                    }
                }
            }
        }
        else
        {
            std::array<MultiFab,AMREX_SPACEDIM> ctmp;
            Vector<MultiFab*> vctmp(AMREX_SPACEDIM);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                BoxArray cba = fine[idim]->boxArray();
                cba.coarsen(ratio);
                ctmp[idim].define(cba, fine[idim]->DistributionMap(), ncomp, ngcrse, MFInfo(), FArrayBoxFactory());
                vctmp[idim] = &ctmp[idim];
            }
            EB_average_down_faces(fine, vctmp, ratio, ngcrse);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                crse[idim]->ParallelCopy(ctmp[idim],0,0,ncomp,ngcrse,ngcrse);
            }
        }
    }
}


void EB_average_down_boundaries (const MultiFab& fine, MultiFab& crse,
                                 const IntVect& ratio, int ngcrse)
{
    int ncomp = crse.nComp();

    if (!fine.hasEBFabFactory())
    {
        crse.setVal(0.0, 0, ncomp, ngcrse);
    }
    else
    {
        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(fine.Factory());
        const auto& flags = factory.getMultiEBCellFlagFab();
        const auto& barea = factory.getBndryArea();

        if (isMFIterSafe(fine, crse))
        {
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(crse, MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
            {
                const Box& tbx = mfi.growntilebox(ngcrse);
                FabType typ = flags[mfi].getType(amrex::refine(tbx,ratio));

                if (FabType::covered == typ || FabType::regular == typ) {
                    crse[mfi].setVal(0.0, tbx, 0, 1);
                } else {
                    amrex_eb_avgdown_boundaries(tbx.loVect(), tbx.hiVect(),
                                                BL_TO_FORTRAN_ANYD(fine[mfi]),
                                                BL_TO_FORTRAN_ANYD(crse[mfi]),
                                                BL_TO_FORTRAN_ANYD(barea[mfi]),
                                                ratio.getVect(), &ncomp);
                }
            }
        }
        else
        {
            BoxArray cba = fine.boxArray();
            cba.coarsen(ratio);
            MultiFab ctmp(cba, fine.DistributionMap(), ncomp, ngcrse, MFInfo(), FArrayBoxFactory());
            EB_average_down_boundaries(fine, ctmp, ratio, ngcrse);
            crse.ParallelCopy(ctmp, 0, 0, ncomp, ngcrse, ngcrse);
        }
    }    
}


void EB_computeDivergence (MultiFab& divu, const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                           const Geometry& geom)
{
    if (!divu.hasEBFabFactory())
    {
        amrex::computeDivergence(divu, umac, geom);
    }
    else
    {
        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(divu.Factory());
        const auto& flags = factory.getMultiEBCellFlagFab();
        const auto& vfrac = factory.getVolFrac();
        const auto& area = factory.getAreaFrac();
        const auto& fcent = factory.getFaceCent();

        iMultiFab cc_mask(divu.boxArray(), divu.DistributionMap(), 1, 1);
        cc_mask.setVal(0);

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            std::vector< std::pair<int,Box> > isects;
            const std::vector<IntVect>& pshifts = geom.periodicity().shiftIntVect();
            const BoxArray& ba = cc_mask.boxArray();
            for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
            {
                IArrayBox& fab = cc_mask[mfi];
                const Box& bx = fab.box();
                for (const auto& iv : pshifts)
                {
                    ba.intersections(bx+iv, isects);
                    for (const auto& is : isects)
                    {
                        fab.setVal(1, is.second-iv);
                    }
                }
            }
        }

        const Real* dxinv = geom.InvCellSize();
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(divu,MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& flagfab = flags[mfi];
            auto& divufab = divu[mfi];
            AMREX_D_TERM(const FArrayBox& ufab = (*umac[0])[mfi];,
                         const FArrayBox& vfab = (*umac[1])[mfi];,
                         const FArrayBox& wfab = (*umac[2])[mfi];);

            auto fabtyp = flagfab.getType(bx);
            if (fabtyp == FabType::covered) {
                divufab.setVal(0.0, bx, 0, 1);
            } else if (fabtyp == FabType::regular) {
                amrex_compute_divergence (BL_TO_FORTRAN_BOX(bx),
                                          BL_TO_FORTRAN_ANYD(divufab),
                                          AMREX_D_DECL(BL_TO_FORTRAN_ANYD(ufab),
                                                       BL_TO_FORTRAN_ANYD(vfab),
                                                       BL_TO_FORTRAN_ANYD(wfab)),
                                          dxinv);
            } else {
                amrex_compute_eb_divergence(BL_TO_FORTRAN_BOX(bx),
                                            BL_TO_FORTRAN_ANYD(divufab),
                                            AMREX_D_DECL(BL_TO_FORTRAN_ANYD(ufab),
                                                         BL_TO_FORTRAN_ANYD(vfab),
                                                         BL_TO_FORTRAN_ANYD(wfab)),
                                            BL_TO_FORTRAN_ANYD(cc_mask[mfi]),
                                            BL_TO_FORTRAN_ANYD(flagfab),
                                            BL_TO_FORTRAN_ANYD(vfrac[mfi]),
                                            AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*area[0])[mfi]),
                                                         BL_TO_FORTRAN_ANYD((*area[1])[mfi]),
                                                         BL_TO_FORTRAN_ANYD((*area[2])[mfi])),
                                            AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*fcent[0])[mfi]),
                                                         BL_TO_FORTRAN_ANYD((*fcent[1])[mfi]),
                                                         BL_TO_FORTRAN_ANYD((*fcent[2])[mfi])),
                                            dxinv);
            }
        }
    }
}

}
