#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#include "make_shapes.H"

using namespace amrex;

int main (int argc, char* argv[]) {

    amrex::Initialize(argc, argv);
    // Issue an error if AMR input file is not given
    if ( argc < 2 )
       amrex::Abort("AMReX input file missing");


    /****************************************************************************
     * Load LEVEL-SET parameters                                                *
     ***************************************************************************/

    int levelset__refinement    = 1;
    int levelset__eb_refinement = 1;
    int levelset__pad           = 1;
    int levelset__eb_pad        = 1;
    {
        ParmParse pp("eb");

        // Parameters used be the level-set algorithm. Refer to LSFactory (or
        // mfix_level.H) for more details:
        //   -> refinement: how well resolved (fine) the (level-set/EB-facet)
        //                  grid needs to be (note: a fine level-set grid means
        //                  that distances and normals are computed accurately)
        //   -> pad:        how many (refined) grid points _outside_ the
        //                  problem domain the grid extends (avoids edge cases
        //                  in physical domain)
        pp.query("levelset__refinement", levelset__refinement);
        pp.query("levelset__eb_refinement", levelset__eb_refinement);
        pp.query("levelset__pad", levelset__pad);
        pp.query("levelset__eb_pad", levelset__eb_pad);
    }


    /****************************************************************************
     * Load GRID parameters                                                     *
     ***************************************************************************/

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, n_lev;
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

        pp.get("n_lev", n_lev);
    }


    /****************************************************************************
     * BUILD BoxArray, Geometry, and DistributionMapping                        *
     ***************************************************************************/

    BoxArray grids;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "grids" from the single box "domain"
        grids.define(domain);
        // Break up boxarray "grids" into chunks no larger than "max_grid_size" along a direction
        grids.maxSize(max_grid_size);

        // This defines the physical box, [-1,1] in each direction.
        RealBox real_box( prob_lo[0], prob_lo[1], prob_lo[2], prob_hi[0], prob_hi[1], prob_hi[2] );

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    DistributionMapping dmap(grids, ParallelDescriptor::NProcs());

    // Level-Set: initialize container for level set
    // level-set MultiFab is defined here, and set to (fortran) huge(amrex_real)
    //            -> use min to intersect new eb boundaries (in update)
    // If you plant to use untions (max), then use level_set->invert() here
    // first...

    int lev = 0;
    std::unique_ptr<LSFactory> level_set = std::unique_ptr<LSFactory>(
        new LSFactory(lev, levelset__refinement, levelset__eb_refinement,
                      levelset__pad, levelset__eb_pad, grids, geom, dmap) );

    LSCoreBase * ls_core;

    // Constructs EB, followed by level-set or ls_core
    make_my_eb2(n_lev, grids, dmap, geom, level_set.get(), ls_core);

    ls_core->InitData();
    ls_core->WritePlotFile();
    ls_core->WriteCheckpointFile();


    const MultiFab * ls_data = level_set->get_data();

    VisMF::Write(* ls_data, "LevelSet");

    level_set.reset();
    delete ls_core;

    amrex::Finalize();
    return 0;
}
