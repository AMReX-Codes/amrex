#include "AMReX_EB_levelset.H"
#include "AMReX_EB_geom.H"

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include "AMReX_BoxIterator.H"
#include <AMReX_EBCellFlag.H>

#include <AMReX_EB_F.H>
#include <AMReX_EB_utils.H>

#include <AMReX_EB2.H>

namespace amrex {

LSFactory::LSFactory(int lev, int ls_ref, int eb_ref, int ls_pad, int eb_pad,
                     const BoxArray & ba, const Geometry & geom, const DistributionMapping & dm,
                     int a_eb_tile_size)
    : amr_lev(lev),
      ls_grid_ref(ls_ref), eb_grid_ref(eb_ref),
      ls_grid_pad(ls_pad), eb_grid_pad(eb_pad),
      eb_tile_size(a_eb_tile_size),
      dx_vect(AMREX_D_DECL(geom.CellSize()[0]/ls_ref,
                           geom.CellSize()[1]/ls_ref,
                           geom.CellSize()[2]/ls_ref)),
      dx_eb_vect(AMREX_D_DECL(geom.CellSize()[0]/eb_ref,
                              geom.CellSize()[1]/eb_ref,
                              geom.CellSize()[2]/eb_ref))
{

    BL_PROFILE("LSFactory::LSFactory()");

    // Init geometry over which the level set and EB are defined
    init_geom(ba, geom, dm);

    // Initialize MultiFab pointers storing level-set data
    //    -> ls_grid:  nodal MultiFab storing signed distance function to the nearest wall
    //    -> ls_valid: cell-centered iMultiFab storing integer flags assuming the following values:
    //         -1 : not all nodes around cell have been initialized
    //          0 : none of the cell's neighbours contain negative vlaues of ls_phi on its nodes
    //          1 : the cell is in the neighbourhood of phi < 0
    ls_grid  = std::unique_ptr<MultiFab>(new MultiFab);
    ls_valid = std::unique_ptr<iMultiFab>(new iMultiFab);

    // Temporary MultiFab used for generating EB factories.
    eb_grid = std::unique_ptr<MultiFab>(new MultiFab);

    // Define ls_grid and ls_valid, growing them by ls_pad
    // Note: box arrays (such as ls_ba) are initialized in init_geom() -> update_ba()
    ls_grid->define(ls_ba, dm, 1, ls_pad);
    ls_valid->define(ls_ba, dm, 1, ls_pad);
    ls_valid->setVal(-1);

    // Define eb_grid, growing it by eb_pad
    eb_grid->define(eb_ba, dm, 1, eb_pad);


    // Initialize by setting all ls_phi = huge(c_real)
    ls_grid->setVal(std::numeric_limits<Real>::max());
}



LSFactory::LSFactory(const LSFactory & other) :
    LSFactory(other.get_amr_level(),
              other.get_ls_ref(), other.get_eb_ref(),
              other.get_ls_pad(), other.get_eb_pad(),
              other.get_ba(), other.get_geom(), other.get_dm() )
{

    BL_PROFILE("LSFactory::LSFactory(LSFactory)");


    //ls_grid  = other.copy_data();
    //ls_valid = other.copy_valid();
}



LSFactory::~LSFactory() {

    BL_PROFILE("LSFactory::~LSFactory()");

    ls_grid.reset();
    ls_valid.reset();
    eb_grid.reset();
}



void LSFactory::update_ba(const BoxArray & new_ba, const DistributionMapping & dm) {

    BL_PROFILE("LSFactory::update_ba()");

    base_ba = new_ba;

    // Refined versions of both the cell-centered and nodal (phi)
    const BoxArray phi_ba = amrex::convert(new_ba, IntVect::TheNodeVector());

    ls_dm = dm;

    ls_ba = phi_ba;
    ls_ba.refine(ls_grid_ref);

    cc_ba = new_ba;
    cc_ba.refine(ls_grid_ref);

    eb_ba = new_ba;
    eb_ba.refine(eb_grid_ref);
}



void LSFactory::init_geom(const BoxArray & ba, const Geometry & geom,
                          const DistributionMapping & dm) {

    BL_PROFILE("LSFactory::init_geom()");

    base_geom = geom;

    // Initialize Geometry objects for the level set and the EB, note that the
    // Geometry objects reflect the refined (and padded) box arrays, preventing
    // periodic fill operations from "spilling over" from refined/padded
    // indices.
    update_ba(ba, dm);

    geom_ls = LSUtility::make_ls_geometry(* this, geom);
    geom_eb = LSUtility::make_eb_geometry(* this, geom);
}



AMREX_GPU_HOST_DEVICE bool LSFactory::neighbour_is_valid(
        Box const & bx, Array4<Real const> const & phi, int i, int j, int k, int n_pad
    ) {

        GpuArray<int, 3> phlo = bx.loVect3d();
        GpuArray<int, 3> phhi = bx.hiVect3d();

        int ilo = std::max(i-n_pad, phlo[0]);
        int ihi = std::min(i+n_pad, phhi[0]);

        int jlo = std::max(j-n_pad, phlo[1]);
        int jhi = std::min(j+n_pad, phhi[1]);

        int klo = std::max(k-n_pad, phlo[2]);
        int khi = std::min(k+n_pad, phhi[2]);

        for (int kk=klo; kk<=khi; ++kk) {
            for (int jj=jlo; jj<=jhi; ++jj){
                AMREX_PRAGMA_SIMD
                for (int ii=ilo; ii<ihi; ++ii){
                    if(phi(ii, jj, kk) <= 0){
                        return true;
                    }
                }
            }
        }

        return false;
}



void LSFactory::fill_valid_kernel(){

    BL_PROFILE("LSFactory::fill_valid_kernel()");

    int search_radius = 1;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi( * ls_grid, true); mfi.isValid(); ++ mfi) {
        Box tile_box                          = mfi.tilebox();
        Array4<Real const> const & ls_tile    = ls_grid->array(mfi);
        Array4<int>        const & valid_tile = ls_valid->array(mfi);

        ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {

                    bool valid_cell = neighbour_is_valid(
                            tile_box, ls_tile, i, j, k, search_radius
                        );

                    if (valid_cell) valid_tile(i, j, k) = 1;
                    else            valid_tile(i, j, k) = 0;
                }
            );
    }
}



void LSFactory::fill_valid(int n){

    BL_PROFILE("LSFactory::fill_valid(n)");

    if(n <= 0) return;

    fill_valid_kernel();
    return fill_valid(n-1);
}



void LSFactory::fill_valid(){

    BL_PROFILE("LSFactory::fill_valid()");

    fill_valid(ls_grid_pad);

   /****************************************************************************
    * Set boundary values of valid_grid                                        *
    ****************************************************************************/

    BL_PROFILE_REGION_START("LSFactory::fill_valid()::bcs");

    // Simulation domain
    Box domain(geom_ls.Domain());

    // Int array flagging periodic directions => no need to fill the periodic
    // ones as they are filled by FillBoundary
    IntVect periodic(
            AMREX_D_DECL(
                geom_ls.isPeriodic(0),
                geom_ls.isPeriodic(1),
                geom_ls.isPeriodic(2)
            )
        );

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi( * ls_grid); mfi.isValid(); ++mfi){

        Box                   tile_box = mfi.growntilebox();
        Array4<int> const & valid_tile = ls_valid->array(mfi);

        int lo, hi;

        //_______________________________________________________________________
        // Apply x-boundary condition

        lo = domain.smallEnd(0);
        hi = domain.bigEnd(0);

        if (tile_box.smallEnd(0) < lo ) {
            if (periodic[0] == 0){
                ParallelFor(tile_box,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                            if (i < lo) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    );
            }
        }

        if (tile_box.bigEnd(0) > hi ) {
            if (periodic[0] == 0) {
                ParallelFor(tile_box,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                            if (i > hi) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    );
            }
        }

        //-----------------------------------------------------------------------

#if (AMREX_SPACEDIM >= 2)

        //_______________________________________________________________________
        // Apply y-boundary condition

        lo = domain.smallEnd(1);
        hi = domain.bigEnd(1);

        if (tile_box.smallEnd(1) < lo ) {
            if (periodic[1] == 0){
                ParallelFor(tile_box,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                            if (j < lo) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    );
            }
        }

        if (tile_box.bigEnd(1) > hi ) {
            if (periodic[1] == 0) {
                ParallelFor(tile_box,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                            if (j > hi) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    );
            }
        }

        //-----------------------------------------------------------------------

#endif


#if (AMREX_SPACEDIM >= 3)

        //_______________________________________________________________________
        // Apply z-boundary condition

        lo = domain.smallEnd(2);
        hi = domain.bigEnd(2);

        if (tile_box.smallEnd(2) < lo ) {
            if (periodic[2] == 0){
                ParallelFor(tile_box,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                            if (k < lo) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    );
            }
        }

        if (tile_box.bigEnd(2) > hi ) {
            if (periodic[2] == 0) {
                ParallelFor(tile_box,
                        [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                            if (k > hi) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    );
            }
        }

        //-----------------------------------------------------------------------
#endif

    }

    BL_PROFILE_REGION_STOP("LSFactory::fill_valid()::bcs");


    // Avoid FillBoundary in recursive steps
    BL_PROFILE_REGION_START("LSFactory::fill_valid()::FillBoundary");
    ls_valid->FillBoundary(geom_ls.periodicity());
    BL_PROFILE_REGION_STOP("LSFactory::fill_valid()::FillBoundary");
}



std::unique_ptr<Vector<Real>> LSFactory::eb_facets(const FArrayBox & norm_tile,
                                                   const CutFab & bcent_tile,
                                                   const EBCellFlagFab & flag_tile,
                                                   const RealVect & dx_eb,
                                                   const Box & eb_search)
{
    BL_PROFILE("LSFactory::eb_facets(tile)");

    // 1-D list of eb-facet data. Format:
    // { px_1, py_1, pz_1, nx_1, ny_1, nz_1, px_2, py_2, ... , nz_N }
    //   ^                 ^
    //   |                 +---- {nx, ny, nz} is the normal vector pointing _towards_ the facet
    //   +-----------------------{px, py, pz} is the position vector of the facet centre
    std::unique_ptr<Vector<Real>> facet_list;


    const Dim3 lo = amrex::lbound(eb_search);
    const Dim3 hi = amrex::ubound(eb_search);

    // TODO: This can probably be parallelized using ReduceSum -- but I think
    // the advantages of parallelizing this rely on parallizing the next kernel
    // also. This is something to revisit in the future

    int n_facets = 0;
    Array4<EBCellFlag const> const & flag_array = flag_tile.array();
    for (int k=lo.z; k<=hi.z; ++k) {
        for (int j=lo.y; j<=hi.y; ++j) {
            for (int i=lo.x; i<=hi.x; ++i) {
                if (flag_array(i, j, k).isSingleValued()) {
                    n_facets ++;
                }
            }
        }
    }


    // TODO: Parallelizing this kenerl would require GPU atomics -- I don't
    // know the status of how well they are implemented in AMReX. This is
    // something to revisit in the future.

    int facet_list_size = 6 * n_facets;
    facet_list = std::unique_ptr<Vector<Real>>(new Vector<Real>(facet_list_size));

    if (n_facets > 0) {
        int i_facet = 0;
        Array4<Real const> const & norm_array  = norm_tile.array();
        Array4<Real const> const & bcent_array = bcent_tile.array();
        for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    if (flag_array(i, j, k).isSingleValued()) {

                        // Compute facet center
                        RealVect eb_cent, cell_corner{(Real) i, (Real) j, (Real) k};
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            eb_cent[d] = ( bcent_array(i, j, k, d)
                                          + cell_corner[d] + 0.5 )*dx_eb[d];
                        }

                        // Add data to eb facet vector
                        for (int d=0; d<AMREX_SPACEDIM; ++d) {
                            (* facet_list)[i_facet + d] = eb_cent[d];
                            (* facet_list)[i_facet + AMREX_SPACEDIM + d] = norm_array(i, j, k, d);
                        }

                        // Increment facet counter by the facet data length
                        i_facet += 2*AMREX_SPACEDIM;
                    }
                }
            }
        }
    }

    return facet_list;
}



void LSFactory::update_intersection(const MultiFab & ls_in, const iMultiFab & valid_in) {

    BL_PROFILE("LSFactory::update_intersection()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi( * ls_grid, true); mfi.isValid(); ++mfi){
        // Operate on growntilebox instead of tilebox => then you don't have to
        // "intersect" any non-periodic boundary values of manually (as was the
        // case with old Fortran code)
        Box tile_box = mfi.growntilebox();

        Array4<int  const> const & valid_in_tile = valid_in.array(mfi);
        Array4<Real const> const & ls_in_tile    = ls_in.array(mfi);
        Array4<Real      > const & ls_tile       = ls_grid->array(mfi);
        Array4<int       > const & valid_tile    = ls_valid->array(mfi);

        ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (valid_in_tile(i, j, k) == 1) {
                        Real in_node = ls_in_tile(i, j, k);
                        Real ls_node = ls_tile(i, j, k);
                        if (in_node < ls_node) {
                            ls_tile(i, j, k) = in_node;
                            if (ls_node <= 0) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    }
                }
            );

    }

// This is (probably) no longer needed -- I don't know why it was ever needed
// (cf comment on grown tile boxes above)
//    /****************************************************************************
//     * Set boundary values of ls_grid                                           *
//     ****************************************************************************/
// 
//     BL_PROFILE_REGION_START("LSFactory::update_intersection()::bcs");
// 
//     // Simulation domain
//     Box domain(geom_ls.Domain());
// 
//     // Int array flagging periodic directions => no need to fill the periodic
//     // ones as they are filled by FillBoundary
//     IntVect periodic(
//             AMREX_D_DECL(
//                 geom_ls.isPeriodic(0),
//                 geom_ls.isPeriodic(1),
//                 geom_ls.isPeriodic(2)
//             )
//         );
// 
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//     for(MFIter mfi( * ls_grid); mfi.isValid(); ++mfi){
//         const auto & valid_in_tile = valid_in[mfi];
//         const auto & ls_in_tile = ls_in[mfi];
//         auto & v_tile = (* ls_valid)[mfi];
//         auto & ls_tile = (* ls_grid)[mfi];
// 
//         amrex_eb_update_levelset_intersection_bcs( BL_TO_FORTRAN_3D(valid_in_tile),
//                                                    BL_TO_FORTRAN_3D(ls_in_tile),
//                                                    BL_TO_FORTRAN_3D(v_tile),
//                                                    BL_TO_FORTRAN_3D(ls_tile),
//                                                    periodic.getVect(),
//                                                    domain.loVect(), domain.hiVect()  );
//     }
// 
//     BL_PROFILE_REGION_STOP("LSFactory::update_intersection()::bcs");

    BL_PROFILE_REGION_START("LSFactory::update_intersection()::FillBoundary");
    ls_grid->FillBoundary(geom_ls.periodicity());
    BL_PROFILE_REGION_STOP("LSFactory::update_intersection()::FillBoundary");

    fill_valid();
}



void LSFactory::update_union(const MultiFab & ls_in, const iMultiFab & valid_in) {

    BL_PROFILE("LSFactory::update_union()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi( * ls_grid, true); mfi.isValid(); ++mfi){
        // Operate on growntilebox instead of tilebox => then you don't have to
        // "intersect" any non-periodic boundary values of manually (as was the
        // case with old Fortran code)
        Box tile_box = mfi.growntilebox();

        Array4<int  const> const & valid_in_tile = valid_in.array(mfi);
        Array4<Real const> const & ls_in_tile    = ls_in.array(mfi);
        Array4<Real      > const & ls_tile       = ls_grid->array(mfi);
        Array4<int       > const & valid_tile    = ls_valid->array(mfi);

        ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (valid_in_tile(i, j, k) == 1) {
                        Real in_node = ls_in_tile(i, j, k);
                        Real ls_node = ls_tile(i, j, k);
                        if (in_node > ls_node) {
                            ls_tile(i, j, k) = in_node;
                            if (ls_node <= 0) {
                                valid_tile(i, j, k) = 1;
                            }
                        }
                    }
                }
            );

    }

// This is (probably) no longer needed -- I don't know why it was ever needed
// (cf comment on grown tile boxes above)
//    /****************************************************************************
//     * Set boundary values of ls_grid                                           *
//     ****************************************************************************/
// 
//     BL_PROFILE_REGION_START("LSFactory::update_union()::bcs");
// 
//     // Simulation domain
//     Box domain(geom_ls.Domain());
// 
//     // Int array flagging periodic directions => no need to fill the periodic
//     // ones as they are filled by FillBoundary
//     IntVect periodic(
//             AMREX_D_DECL(
//                 geom_ls.isPeriodic(0),
//                 geom_ls.isPeriodic(1),
//                 geom_ls.isPeriodic(2)
//             )
//         );
// 
// #ifdef _OPENMP
// #pragma omp parallel
// #endif
//     for(MFIter mfi( * ls_grid); mfi.isValid(); ++mfi){
//         const auto & valid_in_tile = valid_in[mfi];
//         const auto & ls_in_tile = ls_in[mfi];
//         auto & v_tile = (* ls_valid)[mfi];
//         auto & ls_tile = (* ls_grid)[mfi];
// 
//         amrex_eb_update_levelset_union_bcs( BL_TO_FORTRAN_3D(valid_in_tile),
//                                             BL_TO_FORTRAN_3D(ls_in_tile),
//                                             BL_TO_FORTRAN_3D(v_tile),
//                                             BL_TO_FORTRAN_3D(ls_tile),
//                                             periodic.getVect(),
//                                             domain.loVect(), domain.hiVect()  );
//     }
//     BL_PROFILE_REGION_STOP("LSFactory::update_union()::bcs");

    BL_PROFILE_REGION_START("LSFactory::update_union()::FillBoundary");
    ls_grid->FillBoundary(geom_ls.periodicity());
    BL_PROFILE_REGION_START("LSFactory::update_union()::FillBoundary");

    fill_valid();
}



std::unique_ptr<MultiFab> LSFactory::copy_data(const DistributionMapping& dm) const {

    BL_PROFILE("LSFactory::copy_data()");

    std::unique_ptr<MultiFab> cpy(new MultiFab(ls_ba, dm, 1, ls_grid_pad));
    cpy->copy(* ls_grid, 0, 0, 1, ls_grid_pad, ls_grid_pad);

    BL_PROFILE_REGION_START("LSFactory::copy_data()::FillBoundary");
    cpy->FillBoundary(geom_ls.periodicity());
    BL_PROFILE_REGION_STOP("LSFactory::copy_data()::FillBoundary");

    return cpy;
}



std::unique_ptr<iMultiFab> LSFactory::copy_valid(const DistributionMapping& dm) const {

    BL_PROFILE("LSFactory::copy_valid()");

    std::unique_ptr<iMultiFab> cpy(new iMultiFab(ls_ba, dm, 1, ls_grid_pad));
    cpy->copy(* ls_valid, 0, 0, 1, ls_grid_pad, ls_grid_pad);

    BL_PROFILE_REGION_START("LSFactory::copy_valid()::FillBoundary");
    cpy->FillBoundary(geom_ls.periodicity());
    BL_PROFILE_REGION_STOP("LSFactory::copy_valid()::FillBoundary");
    return cpy;
}



std::unique_ptr<MultiFab> LSFactory::coarsen_data() const {

    BL_PROFILE("LSFactory::coarsen_data()");

    // No refinement => do nothing
    if(ls_grid_ref == 1)
        return copy_data(ls_grid.get()->DistributionMap());

    // Target for coarse nodal version of the level-set MultiFab
    std::unique_ptr<MultiFab> ls_crse = std::unique_ptr<MultiFab>(new MultiFab);
    const MultiFab * ls_fine = ls_grid.get(); // Pointer to fine level-set MultiFab

    // Coarse nodal level-set BoxArray (amrex::average_down requires coarse BA)
    BoxArray crse_ba = ls_fine->boxArray();
    crse_ba.coarsen(ls_grid_ref);
    int ng_crse = ls_fine->nGrow() / ls_grid_ref;

    IntVect ratio = IntVect{AMREX_D_DECL(ls_grid_ref, ls_grid_ref, ls_grid_ref)};
    ls_crse->define(crse_ba, ls_fine->DistributionMap(), ls_fine->nComp(), ng_crse);
    // Don't use average_down, average_down_nodal fills ghost cells
    // amrex::average_down(* ls_fine, * ls_crse, 0, 1, ls_grid_ref);
    amrex::average_down_nodal(* ls_fine, * ls_crse, ratio, ng_crse );

    // Do this in the write plot file function
    // Now map the nodal MultiFab to the cell-centered MultiFab:
    //amrex::average_node_to_cellcenter(* ls_crse, 0, * ls_crse, 0, 1);
    return ls_crse;
}



void LSFactory::regrid(const BoxArray & ba, const DistributionMapping & dm)
{

    BL_PROFILE("LSFactory::regrid()");

    // Regrids the level-set data whenever the
    // DistributionMapping has changed:
    //      -> Rebuilds the nodal levelset (ls_ba), cell-centered valid
    //         (cc_ba), and eb (eb_ba) BoxArrays
    //          -> ls_ba, cc_ba, and eb_ba are all inherited
    update_ba(ba, dm);

    int ng = ls_grid_pad;
    std::unique_ptr<MultiFab> ls_grid_new
        = std::unique_ptr<MultiFab>(new MultiFab(ls_ba, dm, 1, ng));

    ls_grid_new->copy(* ls_grid, 0, 0, 1, ng, ng);
    ls_grid_new->FillBoundary(geom_ls.periodicity());
    ls_grid = std::move(ls_grid_new);

    std::unique_ptr<iMultiFab> ls_valid_new
        = std::unique_ptr<iMultiFab>(new iMultiFab(ls_ba, dm, 1,ng));

    ls_valid_new->copy(* ls_valid, 0, 0, 1, ng, ng);
    ls_valid_new->FillBoundary(geom_ls.periodicity());
    ls_valid = std::move(ls_valid_new);
}



void LSFactory::invert() {

   BL_PROFILE("LSFactory::invert()");

    ls_grid->invert(1.0, 0, 1, ls_grid->nGrow());
}



void LSFactory::set_data(const MultiFab & mf_ls){

    BL_PROFILE("LSFactory::set_data()");

    ls_grid->copy(mf_ls, 0, 0, 1, ls_grid_pad, ls_grid_pad);
    ls_grid->FillBoundary(geom_ls.periodicity());

    fill_valid();
}



void LSFactory::fill_data (MultiFab & data, iMultiFab & valid,
                           const EBFArrayBoxFactory & eb_factory,
                           const MultiFab & eb_impfunc,
                           int ebt_size, int ls_ref, int eb_ref,
                           const Geometry & geom, const Geometry & geom_eb) {

    // No profiling in function specialization

    LSFactory::fill_data(data, valid, eb_factory, eb_impfunc,
                         IntVect{AMREX_D_DECL(ebt_size, ebt_size, ebt_size)},
                         ls_ref, eb_ref, geom, geom_eb);
}



void LSFactory::fill_data (MultiFab & data, iMultiFab & valid,
                           const EBFArrayBoxFactory & eb_factory,
                           const MultiFab & eb_impfunc,
                           const IntVect & ebt_size, int ls_ref, int eb_ref,
                           const Geometry & geom, const Geometry & geom_eb) {

    BL_PROFILE("LSFactory::fill_data()");

    /****************************************************************************
     *                                                                          *
     * Sets: `data` MultiFab containing level-set and a `valid` indicating that *
     * the corresponding cell in `data` as been filled by a valid (i.e. the     *
     * value of the level-set was informed by nearby EB facets) level-set       *
     * function                                                                 *
     *                                                                          *
     ***************************************************************************/

    RealVect dx(AMREX_D_DECL(geom.CellSize(0),
                             geom.CellSize(1),
                             geom.CellSize(2)));

    RealVect dx_eb(AMREX_D_DECL(geom_eb.CellSize(0),
                                geom_eb.CellSize(1),
                                geom_eb.CellSize(2)));

    // Don't use the ls_grid_pad for the eb_padding (goes wrong!)
    const int ls_pad = data.nGrow();

    // Doesn't work with iMultiFabs
    //BL_ASSERT(isMFIterSafe(data, valid));

    const BoxArray & ls_ba            = data.boxArray();
    const BoxArray & eb_ba            = eb_factory.boxArray();
    const DistributionMapping & ls_dm = data.DistributionMap();


    /****************************************************************************
     *                                                                          *
     * Access EB Cut-Cell data:                                                 *
     *                                                                          *
     ***************************************************************************/

    const MultiCutFab & bndrycent = eb_factory.getBndryCent();
    const auto & flags = eb_factory.getMultiEBCellFlagFab();

    // make sure to use the EB-factory's ngrow for the eb-padding;
    const int eb_pad = flags.nGrow();

    /****************************************************************************
     *                                                                          *
     * Compute normals data (which are stored on MultiFab over the EB Grid)     *
     *                                                                          *
     ***************************************************************************/

    MultiFab normal(eb_ba, ls_dm, 3, eb_pad); //deliberately use levelset DM
    amrex::FillEBNormals(normal, eb_factory, geom_eb);


    /****************************************************************************
     *                                                                          *
     * Fill local MultiFab with eb_factory's level-set data. Note the role of   *
     * the temporary eb_valid iMultiFab:                                        *
     *  -> eb_valid = 1 if the corresponding eb_ls location could be projected  *
     *                  onto the eb-facets => level-set sign is OK              *
     *  -> eb_valid = 0 if eb_ls is the fall-back (euclidian) distance to the   *
     *                  nearest eb-facet => the sign needs to be checked        *
     *                                                                          *
     ***************************************************************************/

    iMultiFab eb_valid(ls_ba, ls_dm, 1, ls_pad);
    eb_valid.setVal(0);


    /****************************************************************************
     *                                                                          *
     * Identify the shortest length scale (used to compute the minimum distance *
     * to nearest un-detected EB facet)                                         *
     *                                                                          *
     ***************************************************************************/

    const Real min_dx = LSUtility::min_dx(geom_eb);


    /****************************************************************************
     *                                                                          *
     * Loop over EB tile boxes (ebt) and 1. search for EB facets, the 2. find   *
     * least minimum (thresholded) distance to B facets.                        *
     *                                                                          *
     ***************************************************************************/

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(data, ebt_size * std::max(1, ls_ref/eb_ref)); mfi.isValid(); ++mfi)
    {
        //_______________________________________________________________________
        // Fill grown tile box => fill ghost cells as well
        Box tile_box = mfi.growntilebox();


        //_______________________________________________________________________
        // Don't do anything for the current tile if EB facets are ill-defined
        if (! bndrycent.ok(mfi)){

            auto & ls_tile = data[mfi];
            // TODO: check if this still needs to run on the host only?
            ls_tile.setVal<RunOn::Host>( min_dx *( eb_pad + 1), tile_box );

            Array4<Real const> const & if_tile = eb_impfunc.array(mfi);
            Array4<int const>  const &  v_tile = eb_valid.array(mfi);
            Array4<Real      > const &     phi = data.array(mfi);

            ParallelFor(tile_box,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        if (v_tile(i, j, k) == 0) {
                            Real levelset_node = Math::abs( phi(i, j, k) );
                            if (if_tile(i, j, k) <= 0) {
                                phi(i, j, k) = levelset_node;
                            } else {
                                phi(i, j, k) = -levelset_node;
                            }
                        }
                    }
                );

            continue;
        }


        //_______________________________________________________________________
        // Construct search box over which to look for EB facets

        // mfi inherits from ls_grid which might not have the right refinement
        Box eb_search = mfi.tilebox();
        eb_search.coarsen(ls_ref);
        eb_search.refine(eb_ref);

        // Grow search box out from (correctly refined) tile box. NOTE: grow the
        // tile box AT MOST by eb_grid_pad => ensures that EBFactory is defined
        // for the whole search box
        eb_search.enclosedCells(); // search box must be cell-centered
        eb_search.grow(eb_pad);

        const auto & flag       = flags[mfi];
        const auto & if_tile    = eb_impfunc[mfi];
        const auto & norm_tile  = normal[mfi];
        const auto & bcent_tile = bndrycent[mfi];

        auto & v_tile      = eb_valid[mfi];
        auto & ls_tile     = data[mfi];
        auto & region_tile = valid[mfi];


        //_______________________________________________________________________
        // Construct EB facets
        std::unique_ptr<Vector<Real>> facets = eb_facets(norm_tile, bcent_tile, flag,
                                                         dx_eb, eb_search);
        int len_facets = facets->size();


        //_______________________________________________________________________
        // Fill local level-set
        if (len_facets > 0) {

            // HACK: copyingg EB facet data so that the lambda function can
            // capture it by value. Why by value? => In case we want to run
            // this code on device.
            auto facet_vect = *facets;

            Array4<Real> const & ls_array = ls_tile.array();
            Array4<int > const &  v_array = v_tile.array();

            ParallelFor(tile_box,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                        RealVect pos_node {AMREX_D_DECL( (Real) i, (Real) j, (Real) k)};
                        for(int d=0; d<AMREX_SPACEDIM; ++d)
                            pos_node[d] = pos_node[d] * dx[d];


                        Real min_dist=0;
                        bool proj_valid=true;
                        geom::closest_dist (min_dist, proj_valid,
                                            facet_vect, dx_eb, pos_node);

                        ls_array(i, j, k) = min_dist;
                        if (proj_valid) {
                            v_array(i, j, k) = 1;
                        } else {
                            v_array(i, j, k) = 0;
                        }
                    }
                );


            // TODO: check if this still needs to run on the host only?
            region_tile.setVal<RunOn::Host>(1);
        } else {
            // TODO: check if this still needs to run on the host only?
            ls_tile.setVal<RunOn::Host>( min_dx * ( eb_pad + 1 ) , tile_box );
        }


        //_______________________________________________________________________
        // Threshold local level-set
        Real ls_threshold = min_dx * (eb_pad+1); //eb_pad => we know that any EB
                                                 //is _at least_ eb_pad away from
                                                 //the edge of the eb search box

        Array4<Real> const & phi = ls_tile.array();

        Real phi_th = ls_threshold;
        if (phi_th < 0) phi_th = std::numeric_limits<Real>::max();

        ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (phi(i, j, k) >  phi_th) phi(i, j, k) =  phi_th;
                    if (phi(i, j, k) < -phi_th) phi(i, j, k) = -phi_th;
                }
            );


        //_______________________________________________________________________
        // Validate level-set (here so that tile-wise assignment is still validated)
        Array4<Real const> const & if_array = if_tile.array();
        Array4<int  const> const &  v_array = v_tile.array();

        ParallelFor(tile_box,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) {
                    if (v_array(i, j, k) == 0) {
                        Real levelset_node = Math::abs( phi(i, j, k) );
                        if (if_array(i, j, k) <= 0) {
                            phi(i, j, k) = levelset_node;
                        } else {
                            phi(i, j, k) = -levelset_node;
                        }
                    }
                }
            );

    }
}



void LSFactory::fill_data (MultiFab & data, iMultiFab & valid,
                           const MultiFab & mf_impfunc,
                           int eb_pad, const Geometry & eb_geom) {

    /****************************************************************************
     *                                                                          *
     * Set valid (Here: "valid" region defined as all nodes because IF is       *
     * defined everywhere, unless we are using the threshold eb_pad*eb_dx)      *
     *                                                                          *
     ***************************************************************************/

    if (eb_pad <= 0) {
        valid.setVal(1);
    } else {
        valid.setVal(0);
    }


    /****************************************************************************
     *                                                                          *
     * GeometryService convention:                                              *
     *      -- implicit_function(r) < 0 : r in fluid (outside of EB)            *
     *      -- implicit_function(r) > 0 : r not in fluid (inside EB)            *
     *   => If implicit_function is a signed-distance function, we need to      *
     *      invert sign                                                         *
     *                                                                          *
     ***************************************************************************/

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(mf_impfunc, true); mfi.isValid(); ++ mfi) {
        const FArrayBox & in_fab  = mf_impfunc[mfi];
              FArrayBox & out_fab = data[mfi];

        // NOTE: growntilebox => flip also the ghost cells. They don't get
        // flipped twice because we don't call FillBoundary.
        for (BoxIterator bit(mfi.growntilebox()); bit.ok(); ++ bit)
            out_fab(bit(), 0) = - in_fab(bit(), 0);
    }


    /****************************************************************************
     *                                                                          *
     * Threshold local level-set:  eb_pad => we know that any EB is _at least_  *
     * eb_pad away from the edge of the current tile box (consistent with EB-   *
     * based level-set filling).                                                *
     *                                                                          *
     ***************************************************************************/

    if (eb_pad > 0) {
        const Real min_dx = LSUtility::min_dx(eb_geom);
        const Real ls_threshold = min_dx * (eb_pad+1);

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(data, true); mfi.isValid(); ++ mfi) {

            FArrayBox & ls_tile = data[mfi];
            IArrayBox &  v_tile = valid[mfi];

            //___________________________________________________________________
            // Apply threshold to grown tile box => fill ghost cells as well
            Box tile_box = mfi.growntilebox();

            amrex_eb_threshold_levelset(BL_TO_FORTRAN_BOX(tile_box), & ls_threshold,
                                        BL_TO_FORTRAN_3D(ls_tile));


            //___________________________________________________________________
            // Set valid to 1 whenever the box contains values within
            // +/- ls_threshold
            for (BoxIterator bit(mfi.growntilebox()); bit.ok(); ++ bit) {
                if (std::abs(ls_tile(bit(), 0)) <= ls_threshold) {
                    v_tile.setVal<RunOn::Host>(1, tile_box);
                    break;
                }
            }
        }
    }

}



void LSFactory::intersect_data (MultiFab & data, iMultiFab & valid,
                                const MultiFab & data_in, const iMultiFab & valid_in,
                                const Geometry & /*geom_ls*/) {

    BL_PROFILE("LSFactory::intersect_data()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.growntilebox();

        const auto & valid_in_tile = valid_in[mfi];
        const auto & ls_in_tile = data_in[mfi];
        auto & v_tile = valid[mfi];
        auto & ls_tile = data[mfi];

        amrex_eb_update_levelset_intersection(tile_box.loVect(), tile_box.hiVect(),
                                              BL_TO_FORTRAN_3D(valid_in_tile),
                                              BL_TO_FORTRAN_3D(ls_in_tile),
                                              BL_TO_FORTRAN_3D(v_tile),
                                              BL_TO_FORTRAN_3D(ls_tile)              );
    }
}



void LSFactory::union_data (MultiFab & data, iMultiFab & valid,
                            const MultiFab & data_in, const iMultiFab & valid_in,
                            const Geometry & /*geom_ls*/) {

    BL_PROFILE("LSFactory::union_data()");

#ifdef _OPENMP
#pragma omp parallel
#endif
    for(MFIter mfi(data, true); mfi.isValid(); ++mfi){
        Box tile_box = mfi.growntilebox();

        const auto & valid_in_tile = valid_in[mfi];
        const auto & ls_in_tile = data_in[mfi];
        auto & v_tile = valid[mfi];
        auto & ls_tile = data[mfi];

        amrex_eb_update_levelset_intersection(tile_box.loVect(), tile_box.hiVect(),
                                              BL_TO_FORTRAN_3D(valid_in_tile),
                                              BL_TO_FORTRAN_3D(ls_in_tile),
                                              BL_TO_FORTRAN_3D(v_tile),
                                              BL_TO_FORTRAN_3D(ls_tile)              );
    }
}



std::unique_ptr<iMultiFab> LSFactory::Fill(const EBFArrayBoxFactory & eb_factory,
                                           const MultiFab & mf_impfunc) {
    return Fill(eb_factory, mf_impfunc, eb_tile_size);
}



std::unique_ptr<iMultiFab> LSFactory::Fill(const EBFArrayBoxFactory & eb_factory,
                                           const MultiFab & mf_impfunc,
                                           int ebt_size) {
    return Fill(eb_factory, mf_impfunc,
                IntVect{AMREX_D_DECL(ebt_size, ebt_size, ebt_size)});
}



std::unique_ptr<iMultiFab> LSFactory::Fill(const EBFArrayBoxFactory & eb_factory,
                                           const MultiFab & mf_impfunc,
                                           const IntVect & ebt_size) {

    /****************************************************************************
     *                                                                          *
     * Returns: iMultiFab indicating region that has been filled by a valid     *
     * level-set function (i.e. the value of the level-set was informed by      *
     * nearby EB facets)                                                        *
     *                                                                          *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> region_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    region_valid->define(ls_ba, ls_dm, 1, ls_grid_pad);
    region_valid->setVal(0);


    LSFactory::fill_data(* ls_grid, * region_valid, eb_factory, mf_impfunc,
                         ebt_size, ls_grid_ref, eb_grid_ref, geom_ls, geom_eb);


    fill_valid();

    return region_valid;
}



std::unique_ptr<iMultiFab> LSFactory::Fill(const MultiFab & mf_impfunc,
                                           bool apply_threshold) {

    /****************************************************************************
     *                                                                          *
     * Returns: iMultiFab indicating region that has been filled by a valid     *
     * level-set function (i.e. the value of the level-set was informed by      *
     * nearby EB facets)                                                        *
     *                                                                          *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> region_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    region_valid->define(ls_ba, ls_dm, 1, ls_grid_pad);
    region_valid->setVal(0);


    int ls_threshold_pad = -1;
    if (apply_threshold)
        ls_threshold_pad = eb_grid_pad;

    LSFactory::fill_data(* ls_grid, * region_valid, mf_impfunc,
                         ls_threshold_pad, geom_eb);

    fill_valid();

    return region_valid;

}



std::unique_ptr<iMultiFab> LSFactory::Intersect(const EBFArrayBoxFactory & eb_factory,
                                                const MultiFab & mf_impfunc) {

    return Intersect(eb_factory, mf_impfunc, eb_tile_size);

}



std::unique_ptr<iMultiFab> LSFactory::Intersect(const EBFArrayBoxFactory & eb_factory,
                                                const MultiFab & mf_impfunc,
                                                int ebt_size){
    return Intersect(eb_factory, mf_impfunc,
                     IntVect{AMREX_D_DECL(ebt_size, ebt_size, ebt_size)});
}



std::unique_ptr<iMultiFab> LSFactory::Intersect(const EBFArrayBoxFactory & eb_factory,
                                                const MultiFab & mf_impfunc,
                                                const IntVect & ebt_size) {

    /****************************************************************************
     *                                                                          *
     * Returns: iMultiFab indicating region that has been filled by a valid     *
     * level-set function (i.e. the value of the level-set was informed by      *
     * nearby EB facets)                                                        *
     *                                                                          *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> region_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    region_valid->define(ls_ba, ls_dm, 1, ls_grid_pad);
    region_valid->setVal(0);

    // Local MultiFab storing level-set data for this eb_factory
    MultiFab  eb_ls(ls_ba, ls_dm, 1, ls_grid_pad);

    // Fill local MultiFab with EBF level-set
    LSFactory::fill_data(eb_ls, * region_valid, eb_factory, mf_impfunc,
                         ebt_size, ls_grid_ref, eb_grid_ref, geom_ls, geom_eb);


    // Update LSFactory using local eb level-set
    update_intersection(eb_ls, * region_valid);


    return region_valid;

}



std::unique_ptr<iMultiFab> LSFactory::Intersect(const MultiFab & mf_impfunc,
                                                bool apply_threshold) {

    /****************************************************************************
     *                                                                          *
     * Returns: iMultiFab indicating region that has been filled by a valid     *
     * level-set function (i.e. the value of the level-set was informed by      *
     * nearby EB facets)                                                        *
     *                                                                          *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> region_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    region_valid->define(ls_ba, ls_dm, 1, ls_grid_pad);
    region_valid->setVal(0);

    // Local MultiFab storing level-set data for this eb_factory
    MultiFab  if_ls(ls_ba, ls_dm, 1, ls_grid_pad);


    int ls_threshold_pad = -1;
    if (apply_threshold)
        ls_threshold_pad = eb_grid_pad;

    // Fill local MultiFab with implicit-function data
    LSFactory::fill_data(if_ls, * region_valid, mf_impfunc, ls_threshold_pad, geom_eb);


    // Update LSFactory using local (corrected) implicit function MultiFab
    update_intersection(if_ls, * region_valid);

    return region_valid;
}



std::unique_ptr<iMultiFab> LSFactory::Union(const EBFArrayBoxFactory & eb_factory,
                                            const MultiFab & mf_impfunc) {

    return Union(eb_factory, mf_impfunc, eb_tile_size);

}



std::unique_ptr<iMultiFab> LSFactory::Union(const EBFArrayBoxFactory & eb_factory,
                                            const MultiFab & mf_impfunc,
                                            int ebt_size){
    return Union(eb_factory, mf_impfunc,
                     IntVect{AMREX_D_DECL(ebt_size, ebt_size, ebt_size)});
}



std::unique_ptr<iMultiFab> LSFactory::Union(const EBFArrayBoxFactory & eb_factory,
                                            const MultiFab & mf_impfunc,
                                            const IntVect & ebt_size) {

    /****************************************************************************
     *                                                                          *
     * Returns: iMultiFab indicating region that has been filled by a valid     *
     * level-set function (i.e. the value of the level-set was informed by      *
     * nearby EB facets)                                                        *
     *                                                                          *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> region_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    region_valid->define(ls_ba, ls_dm, 1, ls_grid_pad);
    region_valid->setVal(0);

    // Local MultiFab storing level-set data for this eb_factory
    MultiFab  eb_ls(ls_ba, ls_dm, 1, ls_grid_pad);

    // Fill local MultiFab with EBF level-set
    LSFactory::fill_data(eb_ls, * region_valid, eb_factory, mf_impfunc,
                         ebt_size, ls_grid_ref, eb_grid_ref, geom_ls, geom_eb);


    // Update LSFactory using local eb level-set
    update_union(eb_ls, * region_valid);


    return region_valid;

}


std::unique_ptr<iMultiFab> LSFactory::Union(const MultiFab & mf_impfunc,
                                            bool apply_threshold) {

    /****************************************************************************
     *                                                                          *
     * Returns: iMultiFab indicating region that has been filled by a valid     *
     * level-set function (i.e. the value of the level-set was informed by      *
     * nearby EB facets)                                                        *
     *                                                                          *
     ***************************************************************************/

    std::unique_ptr<iMultiFab> region_valid = std::unique_ptr<iMultiFab>(new iMultiFab);
    region_valid->define(ls_ba, ls_dm, 1, ls_grid_pad);
    region_valid->setVal(0);

    // Local MultiFab storing level-set data for this eb_factory
    MultiFab  if_ls(ls_ba, ls_dm, 1, ls_grid_pad);


    int ls_threshold_pad = -1;
    if (apply_threshold)
        ls_threshold_pad = eb_grid_pad;

    // Fill local MultiFab with implicit-function data
    LSFactory::fill_data(if_ls, * region_valid, mf_impfunc, ls_threshold_pad, geom_eb);


    // Update LSFactory using local (corrected) implicit function MultiFab
    update_union(if_ls, * region_valid);

    return region_valid;
}


}
