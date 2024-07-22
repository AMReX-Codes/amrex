#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EB2.H>
#include <AMReX_WriteEBSurface.H>
#include <AMReX_EBToPVD.H>

namespace amrex {

void WriteEBSurface (const BoxArray & ba, const DistributionMapping & dmap, const Geometry & geom,
                     const EBFArrayBoxFactory * ebf) {

    EBToPVD eb_to_pvd;

    const Real* dx     = geom.CellSize();
    const Real* problo = geom.ProbLo();

    MultiFab mf_ba(ba, dmap, 1, 0, MFInfo(), *ebf);

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();
        const auto * my_flag_ptr = &my_flag;

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered ||
            my_flag.getType(bx) == FabType::regular) { continue; }

        std::array<const CutFab *, AMREX_SPACEDIM> areafrac;
        const CutFab * bndrycent;

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            areafrac[d]  =   &(*ebf->getAreaFrac()[d])[mfi];
        }
        bndrycent = &(ebf->getBndryCent()[mfi]);

#ifdef AMREX_USE_GPU
        std::unique_ptr<EBCellFlagFab> host_flag;
        if (my_flag.arena()->isManaged() || my_flag.arena()->isDevice()) {
            host_flag = std::make_unique<EBCellFlagFab>(my_flag.box(), my_flag.nComp(),
                                                  The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(host_flag->dataPtr(), my_flag.dataPtr(),
                                   host_flag->nBytes());
            Gpu::streamSynchronize();
            my_flag_ptr = host_flag.get();
        }

        std::array<std::unique_ptr<CutFab>, AMREX_SPACEDIM> areafrac_h;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (areafrac[d]->arena()->isManaged() || areafrac[d]->arena()->isDevice()) {
                areafrac_h[d] = std::make_unique<CutFab>(areafrac[d]->box(), areafrac[d]->nComp(),
                                                           The_Pinned_Arena());
                Gpu::dtoh_memcpy_async(areafrac_h[d]->dataPtr(), areafrac[d]->dataPtr(),
                                       areafrac[d]->size()*sizeof(Real));
                Gpu::streamSynchronize();
                areafrac[d] = areafrac_h[d].get();
            }
        }

        std::unique_ptr<CutFab> bndrycent_h;
        if (bndrycent->arena()->isManaged() || bndrycent->arena()->isDevice()) {
            bndrycent_h = std::make_unique<CutFab>(bndrycent->box(), bndrycent->nComp(),
                                                        The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(bndrycent_h->dataPtr(), bndrycent->dataPtr(),
                                   bndrycent->size()*sizeof(Real));
            Gpu::streamSynchronize();
            bndrycent = bndrycent_h.get();
        }
#endif

        eb_to_pvd.EBToPolygon(
                problo, dx,
                bx, my_flag_ptr->const_array(),
                bndrycent->const_array(),
                areafrac[0]->const_array(),
                areafrac[1]->const_array(),
                areafrac[2]->const_array());
    }

    int cpu = ParallelDescriptor::MyProc();
    int nProcs = ParallelDescriptor::NProcs();

    eb_to_pvd.WriteEBVTP(cpu);

    if(ParallelDescriptor::IOProcessor()) {
        EBToPVD::WritePVTP(nProcs);
    }

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();
        const auto * my_flag_ptr = &my_flag;

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered ||
            my_flag.getType(bx) == FabType::regular) { continue; }

#ifdef AMREX_USE_GPU
        std::unique_ptr<EBCellFlagFab> host_flag;
        if (my_flag.arena()->isManaged() || my_flag.arena()->isDevice()) {
            host_flag = std::make_unique<EBCellFlagFab>(my_flag.box(), my_flag.nComp(),
                                                  The_Pinned_Arena());
            Gpu::dtoh_memcpy_async(host_flag->dataPtr(), my_flag.dataPtr(),
                                   host_flag->nBytes());
            Gpu::streamSynchronize();
            my_flag_ptr = host_flag.get();
        }
#endif

        eb_to_pvd.EBGridCoverage(cpu, problo, dx, bx, my_flag_ptr->const_array());
    }
}

}
