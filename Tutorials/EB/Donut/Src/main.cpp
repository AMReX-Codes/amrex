#include <AMReX.H>

#include <AMReX_ParmParse.H>

#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_DistributionMapping.H>

#include "EB_F.H"


using namespace amrex;


void WriteEBSurface (const BoxArray & ba, const DistributionMapping & dmap, const Geometry & geom,
                     const EBFArrayBoxFactory & ebf) {

    const Real * dx = geom.CellSize();

    MultiFab mf_ba(ba, dmap, 1, 0, MFInfo(), ebf);

    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto & sfab    = static_cast<EBFArrayBox const&>((mf_ba)[mfi]);
        const auto & my_flag = sfab.getEBCellFlagFab();

        const Box & bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered or
            my_flag.getType(bx) == FabType::regular) continue;

        std::array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
        const MultiCutFab * bndrycent;

        areafrac  =   ebf->getAreaFrac();
        bndrycent = &(ebf->getBndryCent());

        mfix_eb_to_polygon(dx, bx.loVect(), bx.hiVect(),
                           my_flag.dataPtr(), my_flag.loVect(), my_flag.hiVect(),
                           (* bndrycent)[mfi].dataPtr(),
                           (* bndrycent)[mfi].loVect(), (* bndrycent)[mfi].hiVect(),
                           (* areafrac[0])[mfi].dataPtr(),
                           (* areafrac[0])[mfi].loVect(), (* areafrac[0])[mfi].hiVect(),
                           (* areafrac[1])[mfi].dataPtr(),
                           (* areafrac[1])[mfi].loVect(), (* areafrac[1])[mfi].hiVect(),
                           (* areafrac[2])[mfi].dataPtr(),
                           (* areafrac[2])[mfi].loVect(), (* areafrac[2])[mfi].hiVect());
    }

    int cpu = ParallelDescriptor::MyProc();
    int nProcs = ParallelDescriptor::NProcs();

    mfix_write_eb_vtp(& cpu);

    if(ParallelDescriptor::IOProcessor())
        mfix_write_pvtp(& nProcs);


    for (MFIter mfi(mf_ba); mfi.isValid(); ++mfi) {

        const auto& sfab = static_cast<EBFArrayBox const&>((mf_ba)[mfi]);
        const auto& my_flag = sfab.getEBCellFlagFab();

        const Box& bx = mfi.validbox();

        if (my_flag.getType(bx) == FabType::covered or
            my_flag.getType(bx) == FabType::regular) continue;

        mfix_eb_grid_coverage(& cpu, dx, bx.loVect(), bx.hiVect(),
                              my_flag.dataPtr(), my_flag.loVect(), my_flag.hiVect());
    }
}


int main (int argc, char * argv[]) {
    amrex::Initialize(argc, argv);
    // Issue an error if AMR input file is not given
    if ( argc < 2 )
        amrex::Abort("AMReX input file missing");

    // {...} in order to ensure everything has gone out of scope (and therefore
    // is deallocated) before amrex::Finalize is called.
    {

    }

    amrex::Finalize();
    return 0;
}
