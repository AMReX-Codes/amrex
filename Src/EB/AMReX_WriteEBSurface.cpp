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

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered ||
            my_flag.getType(bx) == FabType::regular) continue;

        std::array<const MultiCutFab *, AMREX_SPACEDIM> areafrac;
        const MultiCutFab * bndrycent;

        areafrac  =   ebf->getAreaFrac();
        bndrycent = &(ebf->getBndryCent());

        eb_to_pvd.EBToPolygon(
                problo, dx,
                bx, my_flag.const_array(),
                bndrycent->const_array(mfi),
                areafrac[0]->const_array(mfi),
                areafrac[1]->const_array(mfi),
                areafrac[2]->const_array(mfi));
    }

    int cpu = ParallelDescriptor::MyProc();
    int nProcs = ParallelDescriptor::NProcs();

    eb_to_pvd.WriteEBVTP(cpu);

    if(ParallelDescriptor::IOProcessor())
        eb_to_pvd.WritePVTP(nProcs);

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered ||
            my_flag.getType(bx) == FabType::regular) continue;

        eb_to_pvd.EBGridCoverage(cpu, problo, dx, bx, my_flag.const_array());
    }
}

}

