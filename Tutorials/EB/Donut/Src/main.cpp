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
#include <AMReX_EB2_GeometryShop.H>

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

        areafrac  =   ebf.getAreaFrac();
        bndrycent = &(ebf.getBndryCent());

        mfix_eb_to_polygon(dx, BL_TO_FORTRAN_BOX(bx),
                           BL_TO_FORTRAN_3D(my_flag),
                           BL_TO_FORTRAN_3D((* bndrycent)[mfi]),
                           BL_TO_FORTRAN_3D((* areafrac[0])[mfi]),
                           BL_TO_FORTRAN_3D((* areafrac[1])[mfi]),
                           BL_TO_FORTRAN_3D((* areafrac[2])[mfi]) );
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

        mfix_eb_grid_coverage(& cpu, dx, BL_TO_FORTRAN_BOX(bx), BL_TO_FORTRAN_3D(my_flag));
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

        /************************************************************************
         * Load GRID parameters                                                 *
         ***********************************************************************/

        // AMREX_SPACEDIM: number of dimensions
        int n_cell, max_grid_size, max_lev;
        Vector<int> is_periodic(AMREX_SPACEDIM,0);  // non-periodic in all direction by default
        Vector<Real> prob_lo(AMREX_SPACEDIM,0);
        Vector<Real> prob_hi(AMREX_SPACEDIM,0);
        {
            // ParmParse is way of reading inputs from the inputs file
            ParmParse pp;

            // We need to get n_cell from the inputs file - this is the number of
            // cells on each side of a square (or cubic) domain.
            pp.get("n_cell", n_cell);

            // The domain is broken into boxes of size max_grid_size
            pp.get("max_grid_size", max_grid_size);

            pp.getarr("is_periodic", is_periodic);

            pp.getarr("prob_lo", prob_lo);
            pp.getarr("prob_hi", prob_hi);

            pp.get("max_lev", max_lev);
        }

        /************************************************************************
         * BUILD BoxArray, Geometry, and DistributionMapping                    *
         ***********************************************************************/

        BoxArray grids;
        Geometry geom;
        {
            IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
            IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
            Box domain(dom_lo, dom_hi);

            // Initialize the boxarray "grids" from the single box "domain"
            grids.define(domain);
            // Break up boxarray "grids" into chunks no larger than
            // "max_grid_size" along a direction
            grids.maxSize(max_grid_size);

            // This defines the physical box, [-1,1] in each direction.
            RealBox real_box( prob_lo[0], prob_lo[1],
                              prob_lo[2], prob_hi[0],
                              prob_hi[1], prob_hi[2] );

            // This defines a Geometry object
            geom.define(domain, & real_box, CoordSys::cartesian, is_periodic.data());
        }

        DistributionMapping dmap(grids, ParallelDescriptor::NProcs());


        RealArray center{10., 10., 10.};
        EB2::TorusIF donut(10, 5, center, false);

        auto gshop = EB2::makeShop(donut);
        EB2::Build(gshop, geom, max_lev, max_lev);
    }

    amrex::Finalize();
    return 0;
}
