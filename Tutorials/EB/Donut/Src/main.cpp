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

#include <AMReX_WriteEBSurface.H>
#include <AMReX_WriteEB_F.H>

#include <AMReX_EB_LSCore.H>

using namespace amrex;

inline
void set_tooth_pos(RealArray & tooth_pos, Real angle, Real radius, const RealArray & center) {
    tooth_pos[0] = center[0] + radius * std::cos(angle);
    tooth_pos[1] = center[1] + radius * std::sin(angle);
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
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "grids" from the single box "domain"
        grids.define(domain);
        // Break up boxarray "grids" into chunks no larger than "max_grid_size"
        // along a direction
        grids.maxSize(max_grid_size);

        // This defines the physical box, [-1,1] in each direction.
        RealBox real_box( prob_lo[0], prob_lo[1],
                          prob_lo[2], prob_hi[0],
                          prob_hi[1], prob_hi[2] );

        // This defines a Geometry object
        geom.define(domain, & real_box, CoordSys::cartesian, is_periodic.data());

        DistributionMapping dmap(grids, ParallelDescriptor::NProcs());


        /************************************************************************
         * Basic DONUT                                                          *
         ***********************************************************************/

        RealArray donut_center{30., 30., 30.};
        EB2::TorusIF donut(10, 5, donut_center, false);


        /************************************************************************
         * Bite Mark                                                            *
         ***********************************************************************/

        Real radius = 5;
        Real tooth_radius = 1.3;
        RealArray bite_center{15., 30., 30.};
        EB2::CylinderIF bite_0(radius, 2, bite_center, false);

        Real angle_0 = -1.45;
        RealArray tooth_pos = bite_center;

        set_tooth_pos(tooth_pos, 0 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_1(tooth_radius, 2, tooth_pos, false);

        set_tooth_pos(tooth_pos, 0.4 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_2(tooth_radius, 2, tooth_pos, false);

        set_tooth_pos(tooth_pos, 0.8 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_3(tooth_radius, 2, tooth_pos, false);

        set_tooth_pos(tooth_pos, 1.2 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_4(tooth_radius, 2, tooth_pos, false);

        set_tooth_pos(tooth_pos, 1.6 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_5(tooth_radius, 2, tooth_pos, false);

        set_tooth_pos(tooth_pos, 2.0 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_6(tooth_radius, 2, tooth_pos, false);

        set_tooth_pos(tooth_pos, 2.4 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_7(tooth_radius, 2, tooth_pos, false);

        set_tooth_pos(tooth_pos, 2.8 + angle_0, radius - tooth_radius*0.6, bite_center);
        EB2::CylinderIF tooth_8(tooth_radius, 2, tooth_pos, false);

        auto bite = EB2::makeUnion(bite_0,
                                   tooth_1, tooth_2, tooth_3, tooth_4, tooth_5,
                                   tooth_6, tooth_7, tooth_8);


        /************************************************************************
         * Take a bite OUT of DONUT                                             *
         ***********************************************************************/

        auto bite_donut = EB2::makeDifference(donut, bite);


        /************************************************************************
         * Construct EB surface                                                 *
         ***********************************************************************/

        auto gshop = EB2::makeShop(bite_donut);
        EB2::Build(gshop, geom, 0, 0);

        const EB2::IndexSpace & ebis_donut = EB2::IndexSpace::top();

        int eb_grow = 1;
        EBFArrayBoxFactory ebf_donut ( ebis_donut.getLevel(geom), geom, grids, dmap,
                                       {eb_grow, eb_grow, eb_grow}, EBSupport::full );

        Print() << "Writing EB surface" << std::endl;
        WriteEBSurface (grids, dmap, geom, &ebf_donut);


        /************************************************************************
         * Build AMR mesh around donut                                          *
         ***********************************************************************/

        Print() << "Building adaptive mesh" << std::endl;

        Box domain_crse = domain;
        domain_crse.coarsen(4);

        const IntVect & dom_crse_lo = domain_crse.smallEnd();
        const IntVect & dom_crse_hi = domain_crse.bigEnd();
        // Picket-fence principle
        IntVect n_cells_crse= dom_crse_hi - dom_crse_lo + IntVect{1, 1, 1};
        Vector<int> v_cells = {
            AMREX_D_DECL(n_cells_crse[0], n_cells_crse[1], n_cells_crse[2])
        };

        amrex::Print() << "Declaring AMR levelset:" << std::endl
                       << "coarsest level: " << domain_crse << " n_cells: " << n_cells_crse
                       << std::endl;

        LSCore<decltype(bite_donut)> ls_core(gshop, & real_box, max_lev, v_cells);
        ls_core.InitData();
        ls_core.WritePlotFile();

        amrex::Print() << " ... done" <<std::endl;
    }

    amrex::Finalize();
    return 0;
}
