
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBMultiFabUtil_F.H>
#include <AMReX_MultiFabUtil_F.H>
#include <AMReX_EBCellFlag.H>

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

        const auto& flag_fab = amrex::getEBCellFlagFab(S_fine[mfi]);
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

}
