
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
EB_set_covered (MultiFab& mf)
{
    EB_set_covered(mf, 0, mf.nComp());
}

void
EB_set_covered (MultiFab& mf, int icomp, int ncomp)
{
    BL_PROFILE("EB_set_covered");

    Array<Real> minvals(ncomp);
    for (int i = icomp; i < icomp+ncomp; ++i) {
        minvals[i] = mf.min(i,0,true);
    }
    ParallelDescriptor::ReduceRealMin(minvals.data(), ncomp);

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
                             minvals.data(),&ncomp);
    }
}

void
EB_set_covered (MultiFab& mf, int icomp, int ncomp, Real val)
{
    BL_PROFILE("EB_set_covered_val");

    Array<Real> vals(ncomp, val);

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
EB_set_single_valued_cells (MultiFab& mf, int icomp, int ncomp, Real val)
{
    BL_PROFILE("EB_set_single_valued_cells");

    Array<Real> vals(ncomp, val);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox& fab = mf[mfi];
        const auto& flagfab = amrex::getEBCellFlagFab(fab);
        amrex_eb_set_single_valued_cells(BL_TO_FORTRAN_BOX(bx),
                                         BL_TO_FORTRAN_N_ANYD(fab,icomp),
                                         BL_TO_FORTRAN_ANYD(flagfab),
                                         vals.data(),&ncomp);
    }
}

void
EB_set_volume_fraction (MultiFab& mf)
{
    BL_PROFILE("EB_set_volume_fraction");

    BL_ASSERT(mf.nComp() == 1);

    const Box& domain = mf.getDomain();

// xxxxx
#if 0

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        EBFArrayBox& fab = dynamic_cast<EBFArrayBox&>(mf[mfi]);
        fab.setVal(1.0, bx);
        FabType typ = fab.getType();
        if (typ != FabType::regular)
        {
            if (typ == FabType::covered) {
                fab.setVal(0.0, bx);
            }
            else
            {
                const auto& ebisbox = fab.getEBISBox();

                const Box& bx_sect = bx & domain;
                for (BoxIterator bi(bx_sect); bi.ok(); ++bi)
                {
                    const IntVect& iv = bi();
                    const auto& vofs = ebisbox.getVoFs(iv);
                    Real vtot = 0.0;
                    for (const auto& vi : vofs)
                    {
                        vtot += ebisbox.volFrac(vi);
                    }
                    fab(iv) = vtot;
                }
            }
        }
    }
#endif
}

void
EB_set_bndry_centroid (MultiFab& mf)
{
    BL_PROFILE("EB_set_bndry_centroid");

    BL_ASSERT(mf.nComp() == AMREX_SPACEDIM);

    const Box& domain = mf.getDomain();

// xxxxx
#if 0

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        EBFArrayBox& fab = dynamic_cast<EBFArrayBox&>(mf[mfi]);
        fab.setVal(-1.0, bx);
        FabType typ = fab.getType();
        if (typ != FabType::regular)
        {
            if (typ == FabType::covered) {
                fab.setVal(0.0, bx);
            }
            else
            {
                const auto& ebisbox = fab.getEBISBox();

                const Box& bx_sect = bx & domain;
                for (BoxIterator bi(bx_sect); bi.ok(); ++bi)
                {
                    const IntVect& iv = bi();
                    const auto& vofs = ebisbox.getVoFs(iv);
                    for (const auto& vi : vofs)
                    {
                        const auto& bcent = ebisbox.bndryCentroid(vi);
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                            fab(iv,idim) = bcent[idim];
                        }
                    }
                }
            }
        }
    }
#endif
}

void
EB_set_area_fraction_face_centroid (std::array<MultiFab,AMREX_SPACEDIM>& areafrac,
                                    std::array<MultiFab,AMREX_SPACEDIM>& facecent)
{
    BL_PROFILE("EB_set_area_fraction_face_centroid");

// xxxxx
#if 0

    const Box& domain = areafrac[0].getDomain();
    const auto& cellflagmf = amrex::getMultiEBCellFlagFab(areafrac[0]);
    const EBLevel& eblevel = amrex::getEBLevel(areafrac[0]);
    const EBISLayout& ebisl = eblevel.getEBISL();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(areafrac[0]); mfi.isValid(); ++mfi)
    {
        const Box& bx = amrex::enclosedCells(mfi.fabbox());
        const auto& flagfab = cellflagmf[mfi];
        FabType typ = flagfab.getType(bx);

        if (typ == FabType::covered)
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                areafrac[idim][mfi].setVal(0.0);
            }
        }
        else
        {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                areafrac[idim][mfi].setVal(1.0);
            }
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            facecent[idim][mfi].setVal(0.0);
        }

        if (typ == FabType::singlevalued)
        {
            const auto& ebisbox = ebisl[mfi];
            const Box& bx_sect = bx & domain;
            for (BoxIterator bi(bx_sect); bi.ok(); ++bi)
            {
                const IntVect& iv = bi();
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                {
                    {
                        const auto& lo_faces = ebisbox.getAllFaces(iv, idim, Side::Lo);
                        if (lo_faces.size() == 0) {
                            areafrac[idim][mfi](iv) = 0.0;
                        } else if (lo_faces.size() == 1) {
                            areafrac[idim][mfi](iv) = ebisbox.areaFrac(lo_faces[0]);
                            const RealVect& rv = ebisbox.centroid(lo_faces[0]);
                            int icomp = 0;
                            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                                if (dir != idim) {
                                    facecent[idim][mfi](iv,icomp) = rv[dir];
                                    ++icomp;
                                }
                            }
                        } else {
                            amrex::Abort("EB_set_area_fraction_face_centroid: multi-value not supported");
                        }
                    }

                    if (iv[idim] == bx_sect.bigEnd(idim))
                    {
                        const IntVect& ivhi = iv + IntVect::TheDimensionVector(idim);
                        const auto& hi_faces = ebisbox.getAllFaces(iv, idim, Side::Hi);
                        if (hi_faces.size() == 0) {
                            areafrac[idim][mfi](ivhi) = 0.0;
                        } else if (hi_faces.size() == 1) {
                            areafrac[idim][mfi](ivhi) = ebisbox.areaFrac(hi_faces[0]);
                            const RealVect& rv = ebisbox.centroid(hi_faces[0]);
                            int icomp = 0;
                            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                                if (dir != idim) {
                                    facecent[idim][mfi](ivhi,icomp) = rv[dir];
                                    ++icomp;
                                }
                            }
                        } else {
                            amrex::Abort("EB_set_area_fraction_face_centroid: multi-value not supported");
                        }                        
                    }
                }
            }
        }
    }
#endif
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

    MultiFab crse_S_fine = amrex::makeMultiEBFab(crse_S_fine_BA,fine_dm,ncomp,0,MFInfo(),S_crse);

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(crse_S_fine,true); mfi.isValid(); ++mfi)
    {
        const Box& tbx = mfi.tilebox();
        auto& crse_fab = crse_S_fine[mfi];
        const auto& fine_fab = S_fine[mfi];
        const auto& flag_fab = amrex::getEBCellFlagFab(crse_fab);

        FabType typ = flag_fab.getType(tbx);
        
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
                                BL_TO_FORTRAN_ANYD(flag_fab),
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

MultiFab 
makeMultiEBFab (const BoxArray& ba, const DistributionMapping& dm,
                int ncomp, int ngrow, const MFInfo& info,
                const MultiFab& mold)
{
#if 0
    const auto& fact_mold = dynamic_cast<EBFArrayBoxFactory const&>(mold.Factory());
    const auto& eblevel_mold = fact_mold.getEBLevel();

    if (ba == mold.boxArray() && dm == mold.DistributionMap() && ngrow <= eblevel_mold.getGhost()) 
    {
        return MultiFab(ba, dm, ncomp, ngrow, info, fact_mold);
    }
    else
    {
        EBLevel eblevel(ba, dm, eblevel_mold.getDomain(), ngrow);
        EBFArrayBoxFactory fact(Geometry(), eblevel);
        return MultiFab(ba, dm, ncomp, ngrow, info, fact);
    }
#else
    return MultiFab();
#endif
}

}
