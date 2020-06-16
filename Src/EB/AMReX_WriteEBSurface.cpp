#include <AMReX.H>

#include <AMReX_ParmParse.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_DistributionMapping.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Torus.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Difference.H>
#include <AMReX_EB2_GeometryShop.H>

#include <AMReX_EB_LSCore.H>
#include <AMReX_WriteEBSurface.H>
#include <AMReX_WriteEB_F.H>

namespace amrex {

void WriteEBSurface (const BoxArray & ba, const DistributionMapping & dmap, const Geometry & geom,
                     const EBFArrayBoxFactory * ebf) {

    const Real* dx     = geom.CellSize();
    const Real* problo = geom.ProbLo();

    MultiFab mf_ba(ba, dmap, 1, 0, MFInfo(), *ebf);

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered or
            my_flag.getType(bx) == FabType::regular) continue;

        std::array<const MultiCutFab *, AMREX_SPACEDIM> areafrac;
        const MultiCutFab * bndrycent;

        areafrac  =   ebf->getAreaFrac();
        bndrycent = &(ebf->getBndryCent());

        amrex_eb_to_polygon(problo, dx, BL_TO_FORTRAN_BOX(bx),
                            BL_TO_FORTRAN_3D(my_flag),
                            BL_TO_FORTRAN_3D((* bndrycent)[mfi]),
                            BL_TO_FORTRAN_3D((* areafrac[0])[mfi]),
                            BL_TO_FORTRAN_3D((* areafrac[1])[mfi]),
                            BL_TO_FORTRAN_3D((* areafrac[2])[mfi]) );
    }

    int cpu = ParallelDescriptor::MyProc();
    int nProcs = ParallelDescriptor::NProcs();

    amrex_write_eb_vtp(& cpu);

    if(ParallelDescriptor::IOProcessor())
        amrex_write_pvtp(& nProcs);

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const &>(mf_ba[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered or
            my_flag.getType(bx) == FabType::regular) continue;

        amrex_eb_grid_coverage(& cpu, problo, dx, BL_TO_FORTRAN_BOX(bx), BL_TO_FORTRAN_3D(my_flag));
    }
}

}

