
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MultiFabUtil_C.H>
#include <AMReX_EBMultiFabUtil_C.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>

#include <AMReX_EBMultiFabUtil_F.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex
{

void
EB_set_covered (MultiFab& mf, Real val)
{
    Vector<Real> vals(mf.nComp(), val);
    EB_set_covered(mf, 0, mf.nComp(), 0, vals);
}

void
EB_set_covered (MultiFab& mf, int icomp, int ncomp, const Vector<Real>& vals)
{
    EB_set_covered(mf, icomp, ncomp, 0, vals);
}

void
EB_set_covered (MultiFab& mf, int icomp, int ncomp, int ngrow, Real val)
{
    // xxxxx todo gpu
    Vector<Real> vals(ncomp, val);
    EB_set_covered(mf, icomp, ncomp, ngrow, vals);
}

void
EB_set_covered (MultiFab& mf, int icomp, int ncomp, int ngrow, const Vector<Real>& a_vals)
{
    AMREX_ALWAYS_ASSERT(mf.ixType().cellCentered() || mf.ixType().nodeCentered());
    bool is_cell_centered = mf.ixType().cellCentered();
    int ng = std::min(mf.nGrow(),ngrow);

    Gpu::AsyncArray<Real> vals_aa(a_vals.data(), ncomp);
    Real const* AMREX_RESTRICT vals = vals_aa.data();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(ng);
        FArrayBox& fab = mf[mfi];
        const auto& flagfab = amrex::getEBCellFlagFab(fab);

        Array4<Real> const& arr = mf.array(mfi);
        Array4<EBCellFlag const> const& flagarr = flagfab.array();

        if (is_cell_centered) {
            AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                if (flagarr(i,j,k).isCovered()) {
                    arr(i,j,k,n+icomp) = vals[n];
                }
            });
        } else {
            amrex_eb_set_covered_nodes(BL_TO_FORTRAN_BOX(bx),
                                       BL_TO_FORTRAN_N_ANYD(fab,icomp),
                                       BL_TO_FORTRAN_ANYD(flagfab),
                                       vals, &ncomp);
        }
    }
}

void
EB_set_covered_faces (const Array<MultiFab*,AMREX_SPACEDIM>& umac, Real val)
{
    const auto factory = dynamic_cast<EBFArrayBoxFactory const*>(&(umac[0]->Factory()));
    if (factory == nullptr) return;

    const auto& area = factory->getAreaFrac();
    const auto& flags = factory->getMultiEBCellFlagFab();
    const int ncomp = umac[0]->nComp();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*umac[0],TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        AMREX_D_TERM(const Box& xbx = mfi.tilebox(IntVect::TheDimensionVector(0));,
                     const Box& ybx = mfi.tilebox(IntVect::TheDimensionVector(1));,
                     const Box& zbx = mfi.tilebox(IntVect::TheDimensionVector(2)););
        AMREX_D_TERM(Array4<Real> const& u = umac[0]->array(mfi);,
                     Array4<Real> const& v = umac[1]->array(mfi);,
                     Array4<Real> const& w = umac[2]->array(mfi););

        auto fabtyp = flags[mfi].getType();
        if (fabtyp == FabType::covered)
        {
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                xbx, txbx,
                {
                    const auto lo = amrex::lbound(txbx);
                    const auto hi = amrex::ubound(txbx);
                    for (int n = 0; n < ncomp; ++n) {
                        for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            u(i,j,k,n) = 0.0;
                        }}}
                    }
                }
                ,ybx, tybx,
                {
                    const auto lo = amrex::lbound(tybx);
                    const auto hi = amrex::ubound(tybx);
                    for (int n = 0; n < ncomp; ++n) {
                        for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            v(i,j,k,n) = 0.0;
                        }}}
                    }
                }
#if (AMREX_SPACEDIM == 3)
                ,zbx, tzbx,
                {
                    const auto lo = amrex::lbound(tzbx);
                    const auto hi = amrex::ubound(tzbx);
                    for (int n = 0; n < ncomp; ++n) {
                        for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            w(i,j,k,n) = 0.0;
                        }}}
                    }
                }
#endif
                );
        }
        else if (fabtyp == FabType::singlevalued)
        {
            AMREX_D_TERM(Array4<Real const> const& ax = area[0]->array(mfi);,
                         Array4<Real const> const& ay = area[1]->array(mfi);,
                         Array4<Real const> const& az = area[2]->array(mfi););
            AMREX_LAUNCH_HOST_DEVICE_LAMBDA (
                xbx, txbx,
                {
                    const auto lo = amrex::lbound(txbx);
                    const auto hi = amrex::ubound(txbx);
                    for (int n = 0; n < ncomp; ++n) {
                        for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            if (ax(i,j,k) == 0.0) u(i,j,k,n) = 0.0;
                        }}}
                    }
                }
                ,ybx, tybx,
                {
                    const auto lo = amrex::lbound(tybx);
                    const auto hi = amrex::ubound(tybx);
                    for (int n = 0; n < ncomp; ++n) {
                        for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            if (ay(i,j,k) == 0.0) v(i,j,k,n) = 0.0;
                        }}}
                    }
                }
#if (AMREX_SPACEDIM == 3)
                ,zbx, tzbx,
                {
                    const auto lo = amrex::lbound(tzbx);
                    const auto hi = amrex::ubound(tzbx);
                    for (int n = 0; n < ncomp; ++n) {
                        for (int k = lo.z; k <= hi.z; ++k) {
                        for (int j = lo.y; j <= hi.y; ++j) {
                        for (int i = lo.x; i <= hi.x; ++i) {
                            if (az(i,j,k) == 0.0) w(i,j,k,n) = 0.0;
                        }}}
                    }
                }
#endif
                );
        }
    }
}

void
EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, const MultiFab& vol_fine,
                 const MultiFab& vfrac_fine, int scomp, int ncomp, const IntVect& ratio)
{
    BL_PROFILE("EB_average_down");

    AMREX_ASSERT(S_fine.ixType().cellCentered());
    AMREX_ASSERT(S_crse.ixType().cellCentered());

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
            amrex_avgdown(tbx, crse_fab.array(), fine_fab.array(), 0, scomp, ncomp, ratio);
#else
            amrex_avgdown_with_vol(tbx, crse_fab.array(), fine_fab.array(), vol_fine[mfi].array(),
                                   0, scomp, ncomp, ratio);
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
                    amrex_avgdown(tbx,crse_fab.array(),fine_fab.array(),scomp,scomp,ncomp,ratio);
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
                    amrex_avgdown(tbx,crse_fab.array(),fine_fab.array(),0,scomp,ncomp,ratio);
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


void EB_average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                            const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                            int ratio, int ngcrse)
{
    EB_average_down_faces(fine, crse, IntVect{ratio}, ngcrse);
}

void EB_average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                            const Array<MultiFab*,AMREX_SPACEDIM>& crse,
                            const IntVect& ratio, int ngcrse)
{
    AMREX_ASSERT(crse[0]->nComp() == fine[0]->nComp());

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
                        amrex_avgdown_faces(tbx, (*crse[n])[mfi].array(), (*fine[n])[mfi].array(), 0, 0, ncomp, ratio, n);
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
            Array<MultiFab,AMREX_SPACEDIM> ctmp;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                BoxArray cba = fine[idim]->boxArray();
                cba.coarsen(ratio);
                ctmp[idim].define(cba, fine[idim]->DistributionMap(), ncomp, ngcrse, MFInfo(), FArrayBoxFactory());
            }
            EB_average_down_faces(fine, amrex::GetArrOfPtrs(ctmp), ratio, ngcrse);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                crse[idim]->ParallelCopy(ctmp[idim],0,0,ncomp,ngcrse,ngcrse);
            }
        }
    }
}

void EB_average_down_boundaries (const MultiFab& fine, MultiFab& crse,
                                 int ratio, int ngcrse)
{
    EB_average_down_boundaries(fine,crse,IntVect{ratio},ngcrse);
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
                    crse[mfi].setVal(0.0, tbx, 0, ncomp);
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

        const GpuArray<Real,AMREX_SPACEDIM> dxinv = geom.InvCellSizeArray();
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

            const auto fabtyp = flagfab.getType(bx);
            if (fabtyp == FabType::covered) {
                divufab.setVal(0.0, bx, 0, 1);
            } else if (fabtyp == FabType::regular) {
                amrex_compute_divergence(bx,divufab.array(),AMREX_D_DECL(ufab.array(),vfab.array(),wfab.array()),dxinv);
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
                                            dxinv.data());
            }
        }
    }
}

void
EB_average_face_to_cellcenter (MultiFab& ccmf, int dcomp,
                               const Array<MultiFab const*,AMREX_SPACEDIM>& fmf)
{
    AMREX_ASSERT(ccmf.nComp() >= dcomp + AMREX_SPACEDIM);

    if (!fmf[0]->hasEBFabFactory())
    {
        average_face_to_cellcenter(ccmf, dcomp, fmf);
    }
    else
    {
        const auto& factory = dynamic_cast<EBFArrayBoxFactory const&>(fmf[0]->Factory());
        const auto& flags = factory.getMultiEBCellFlagFab();
        const auto& area = factory.getAreaFrac();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(ccmf,MFItInfo().EnableTiling().SetDynamic(true)); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            const auto& flagfab = flags[mfi];
            auto& ccfab = ccmf[mfi];
            AMREX_D_TERM(const auto& xfab = (*fmf[0])[mfi];,
                         const auto& yfab = (*fmf[1])[mfi];,
                         const auto& zfab = (*fmf[2])[mfi];);
            const auto fabtyp = flagfab.getType(bx);
            if (fabtyp == FabType::covered) {
                ccfab.setVal(0.0, bx, dcomp, 1);
            } else if (fabtyp == FabType::regular) {
                amrex_avg_fc_to_cc(bx,ccfab.array(),AMREX_D_DECL(xfab.array(),yfab.array(),zfab.array()),dcomp);
            } else {
                amrex_eb_avg_fc_to_cc(BL_TO_FORTRAN_BOX(bx),
                                      BL_TO_FORTRAN_N_ANYD(ccfab,dcomp),
                                      AMREX_D_DECL(BL_TO_FORTRAN_ANYD(xfab),
                                                   BL_TO_FORTRAN_ANYD(yfab),
                                                   BL_TO_FORTRAN_ANYD(zfab)),
                                      AMREX_D_DECL(BL_TO_FORTRAN_ANYD((*area[0])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*area[1])[mfi]),
                                                   BL_TO_FORTRAN_ANYD((*area[2])[mfi])),
                                      BL_TO_FORTRAN_ANYD(flagfab));
            }
        }
    }
}

}
