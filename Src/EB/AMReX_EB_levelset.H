#ifndef AMREX_EB_LEVELSET_H_
#define AMREX_EB_LEVELSET_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_GeometryShop.H>

#include <AMReX_EBFabFactory.H>
#include <AMReX_EBCellFlag.H>


namespace amrex {

class LSFactory {
    private:

        //! MultiFab pointers storing level-set data
        //!    -> ls_grid: nodal MultiFab storing signed distance function to the
        //!                nearest wall
        //!    -> ls_valid: cell-centered iMultiFab storing integer flags
        //!                 assuming the following values:
        //!         -1 : not all nodes around cell have been initialized
        //!          0 : none of the cell's neighbours contain negative vlaues of
        //!              ls_phi on its nodes
        //!          1 : the cell is in the neighbourhood of phi < 0
        std::unique_ptr<MultiFab>  ls_grid;
        std::unique_ptr<iMultiFab> ls_valid;

        //! AMR information (needed for EBISLevel)
        const int amr_lev;

        DistributionMapping ls_dm;

        //! Refined box array (over which the level-set MultiFab is constructed)
        BoxArray ls_ba, cc_ba, eb_ba;
        const int ls_grid_ref, eb_grid_ref;
        //! Data structures encloding (padded) EB index space
        std::unique_ptr<MultiFab> eb_grid;

        //! Geometry objects reflecting refined (and padded) BoxArrays. This is
        //! required to prevent PeriodicFill operations from "spilling" refined
        //! indices over periodic boundaries
        Geometry geom_ls, geom_eb;

        //! Padding of the leve-set and EB grids:
        const int ls_grid_pad, eb_grid_pad;

        //! Tiling for local level-set filling
        int eb_tile_size;

        //! Baseline BoxArray and Geometry from which refined quantities are
        //! derived. These are mainly kept around for the copy constructor.
        BoxArray base_ba;
        Geometry base_geom;

        //! Cell size of the refined grids
        const RealVect dx_vect, dx_eb_vect;

        //! Updates internal BoxArray and DistributionMapping given a
        //! "reference" BoxArray at a coarser level. Uses `ls_grid_ref` and
        //! `eb_grid_ref` to construct refined versions of `ba`.
        void update_ba(const BoxArray & new_ba, const DistributionMapping & dm);

        //! Initialize internal BoxArray, and Geometry objects given reference
        //! BoxArray and Geometry (at a coarser level of refinement)
        void init_geom(const BoxArray & ba, const Geometry & geom, const DistributionMapping & dm);

        void fill_valid_kernel();
        void fill_valid(int n);
        void fill_valid();

    public:
        LSFactory(int lev, int ls_ref, int eb_ref, int ls_pad, int eb_pad,
                  const BoxArray & ba, const Geometry & geom, const DistributionMapping & dm,
                  int eb_tile_size = 32);
        LSFactory(const LSFactory & other);
        ~LSFactory();

        void regrid(const BoxArray & ba, const DistributionMapping & dm);

        void invert();

        void set_data(const MultiFab & mf_ls);


        /************************************************************************
         *                                                                      *
         * Utility to gather local EB facets                                    *
         *                                                                      *
         ***********************************************************************/

        //! Construct eb_facets list from eb_factory. The facet list is a 1-D
        //! list of facet bcent positions, followed by normals => N facets would
        //! be stored in a Vector<Real> of length 6*N. Only facets that fall
        //! within the Box `eb_search` are included in the list. The RealVect
        //! `eb_dx` provides the discretization of the EB data. NOTE:
        //! `be_search` must be contained in the `norm_tile`, `bcent_tile` and
        //! `flag_tile`.
        static std::unique_ptr<Vector<Real>> eb_facets(const FArrayBox & norm_tile,
                                                       const CutFab & bcent_tile,
                                                       const EBCellFlagFab & flag_tile,
                                                       const RealVect & eb_dx,
                                                       const Box & eb_search);


        /************************************************************************
         *                                                                      *
         * Utilities that fill the level-set MultiFab using local EB Factory    *
         * NOTE: returns an iMultiFab flagging "valid" data                     *
         *                                                                      *
         ***********************************************************************/

        //! Fills level-set MultiFab `data` locally from EBFArrayBoxFactory
        //! `eb_factory`. Same as `fill_data` but with cubic EB search box.
        static void fill_data (MultiFab & data, iMultiFab & valid,
                               const EBFArrayBoxFactory & eb_factory,
                               const MultiFab & eb_impfunc,
                               int ebt_size, int ls_ref, int eb_ref,
                               const Geometry & geom, const Geometry & geom_eb);

        //! Fills level-set MultiFab `data` locally from EBFArrayBoxFactory
        //! `eb_factory`. Also fills iMultiFab tagging cells which are nearby to
        //! EB surface. Only EB facets in a box size of `ebt_size` are
        //! considered. Any EB facets that are outside this box are ignored =>
        //! the min/max value of the level-set are +/- `eb_pad *
        //! min(geom_eb.CellSize(:))`.
        static void fill_data (MultiFab & data, iMultiFab & valid,
                               const EBFArrayBoxFactory & eb_factory,
                               const MultiFab & eb_impfunc,
                               const IntVect & ebt_size, int ls_ref, int eb_ref,
                               const Geometry & geom, const Geometry & geom_eb);

        //! Fills level-set MultiFab `data` from implicit function MultiFab
        //! `mf_impfunc`. Also fills iMultiFab tagging cells whose values are
        //! informed by nearby EB surfaces (in the case of implicit-function
        //! filling without threshold, this covers the whole domain). This works
        //! as long as the implicit function is a signed distance function. The
        //! int `eb_pad` and Geometry object `geom_eb` can be used to specify a
        //! level-set threshold (`eb_pad*min(geom_eb.CellSize(:))`) to make the
        //! result compatible with Intersections and Unions with EB-filled
        //! level-sets. No threshold is applied with `eb_pad < 0`
        static void fill_data (MultiFab & data, iMultiFab & valid,
                               const MultiFab & mf_impfunc,
                               int eb_pad, const Geometry & eb_geom);


        //! Updates the level-set MultiFab `data` with the result of
        //! intersection operation between `data` and another level-set MultiFab
        //! `data_in`, for all cells where `valid_in == 1`. The iMultiFab
        //! `valid` is also updated with the new valid flags corresponding to
        //! the intersection level-set.
        static void intersect_data (MultiFab & data, iMultiFab & valid,
                                    const MultiFab & data_in, const iMultiFab & valid_in,
                                    const Geometry & geom_ls);


        //! Updates the level-set MultiFab `data` with the result of union
        //! operation between `data` and another level-set MultiFab `data_in`,
        //! for all cells where `valid_in == 1`. The iMultiFab `valid` is also
        //! updated with the new valid flags corresponding to the intersection
        //! level-set.
        static void union_data (MultiFab & data, iMultiFab & valid,
                                const MultiFab & data_in, const iMultiFab & valid_in,
                                const Geometry & geom_ls);



        /************************************************************************
         *                                                                      *
         * Public member functions to update level-set data                     *
         *                                                                      *
         ***********************************************************************/

        //! Fills (overwrites) level-set data locally. The local EB search
        //! radius is managed by LSFactory. EB search box is defaulted to
        //! `LSFactory::eb_tile_size`. Level-set values are guaranteed to be
        //! accurate up to `LSFactory::eb_grid_pad` cells. Returns: iMultiFab
        //! indicating region that has been filled by a valid level-set function
        //! (i.e. the value of the level-set was informed by nearby EB facets)
        std::unique_ptr<iMultiFab> Fill(const EBFArrayBoxFactory & eb_factory,
                                        const MultiFab & mf_impfunc);

        //! Fills (overwrites) level-set data locally. The local EB search
        //! radius is specified by the user corresponding to a cube with
        //! edge-length = `ebt_size`. Level-set values are guaranteed to be
        //! accurate up to `LSFactory::eb_grid_pad` cells. Returns: iMultiFab
        //! indicating region that has been filled by a valid level-set function
        //! (i.e. the value of the level-set was informed by nearby EB facets)
        std::unique_ptr<iMultiFab> Fill(const EBFArrayBoxFactory & eb_factory,
                                        const MultiFab & mf_impfunc,
                                        int ebt_size);

        //! Fills (overwrites) level-set data locally. The local EB search
        //! radius is specified by the user corresponding to a box with edge
        //! lengths given by the IntVect `ebt_size`. Level-set values are
        //! guaranteed to be accurate up to `LSFactory::eb_grid_pad` cells.
        //! Returns: iMultiFab indicating region that has been filled by a valid
        //! level-set function (i.e. the value of the level-set was informed by
        //! nearby EB facets)
        std::unique_ptr<iMultiFab> Fill(const EBFArrayBoxFactory & eb_factory,
                                        const MultiFab & mf_impfunc,
                                        const IntVect & ebt_size);

        //! Fills (overwrites) level-set data locally. The level-set is given by
        //! an implicit function which is defined on a MultiFab `mf_impfunc`,
        //! which has the same resolution, and at least as many ghost-cells, as
        //! the LSFactory. This can be guaranteed using the
        //! LSUtility::GShopLSFactory. Requires that the implicit function is a
        //! signed distance function. If `apply_threshold = true` level-set is
        //! truncated to to a length scale corresponding to `eb_grid_pad` cells
        //! (making it compatible with the EBFArrayBoxFactory version). Returns:
        //! iMultiFab indicating region that has been filled by a valid
        //! level-set function (i.e. the value of the level-set was informed by
        //! nearby EB facets)
        std::unique_ptr<iMultiFab> Fill(const MultiFab & mf_impfunc,
                                        bool apply_threshold = false);

        //! Performs intersection operation with the level-set representation of
        //! `eb_factory`. The implicit function (mf_impfunc) is needed to select
        //! the inside/outside edge-cases where the level-set cannot be
        //! determined by projection onto EB facets.
        std::unique_ptr<iMultiFab> Intersect(const EBFArrayBoxFactory & eb_factory,
                                             const MultiFab & mf_impfunc);

        //! Performs intersection operation with the level-set representation of
        //! `eb_factory`. The integer `ebt_size` specifies the (cubic) tile size
        //! used to search of EB facets.
        std::unique_ptr<iMultiFab> Intersect(const EBFArrayBoxFactory & eb_factory,
                                             const MultiFab & mf_impfunc,
                                             int ebt_size);

        //! Performs intersection operation with the level-set representation of
        //! `eb_factory`. The IntVect `ebt_size` specifies the tile size used to
        //! search of EB facets.
        std::unique_ptr<iMultiFab> Intersect(const EBFArrayBoxFactory & eb_factory,
                                             const MultiFab & mf_impfunc,
                                             const IntVect & ebt_size);

        //! Performs intersection operation with the level-set given in the
        //! MultiFab `mf_impfunc`. If the flag `apply_threshold == true`, then
        //! the level-set in `mf_impfunc` will be thresholded.
        std::unique_ptr<iMultiFab> Intersect(const MultiFab & mf_impfunc,
                                             bool apply_threshold = false);

        //! Performs union operation with the level-set representation of
        //! `eb_factory`. The implicit function (mf_impfunc) is needed to select
        //! the inside/outside edge-cases where the level-set cannot be
        //! determined by projection onto EB facets.
        std::unique_ptr<iMultiFab> Union(const EBFArrayBoxFactory & eb_factory,
                                         const MultiFab & mf_impfunc);

        //! Performs union operation with the level-set representation of
        //! `eb_factory`. The integer `ebt_size` specifies the (cubic) tile size
        //! used to search of EB facets.
        std::unique_ptr<iMultiFab> Union(const EBFArrayBoxFactory & eb_factory,
                                         const MultiFab & mf_impfunc,
                                         int ebt_size);

        //! Performs union operation with the level-set representation of
        //! `eb_factory`. The IntVect `ebt_size` specifies the tile size used to
        //! search of EB facets.
        std::unique_ptr<iMultiFab> Union(const EBFArrayBoxFactory & eb_factory,
                                         const MultiFab & mf_impfunc,
                                         const IntVect & ebt_size);

        //! Performs union operation with the level-set given in the MultiFab
        //! `mf_impfunc`. If the flag `apply_threshold == true`, then the
        //! level-set in `mf_impfunc` will be thresholded.
        std::unique_ptr<iMultiFab> Union(const MultiFab & mf_impfunc,
                                         bool apply_threshold = false);

        //! Update internal `ls_grid`` and `ls_valid` taking the min of
        //! `ls_grid` and `ls_in`
        void update_intersection(const MultiFab & ls_in, const iMultiFab & valid_in);

        //! Update internal `ls_grid`` and `ls_valid` taking the max of
        //! `ls_grid` and `ls_in`
        void update_union(const MultiFab & ls_in, const iMultiFab & valid_in);

        //! Return mutable pointer to level-set data
              MultiFab *  get_data()       {return ls_grid.get();};
        //! Return constant pointer to level-set data
        const MultiFab *  get_data() const {return ls_grid.get();};
        //! Return constant pointer to level-set data
        const iMultiFab * get_valid() const {return ls_valid.get();};

        //! Copy level-set data into a new MultiFab with BoxArray matching the
        //! LSFactory `ls_ba`, and DistributionMapping `dm`. Returns a
        //! unique_ptr to the new MultiFab.
        std::unique_ptr<MultiFab>  copy_data(const DistributionMapping & dm) const;

        //! Copy level-set valid tags into a new iMultiFab with BoxArray
        //! matching the LSFactory `ls_ba`, and DistributionMapping `dm`.
        //! Returns a unique_ptr to the new iMultiFab.
        std::unique_ptr<iMultiFab> copy_valid(const DistributionMapping & dm) const;

        //! Coarsen level-set data by a coarsening ration of `ls_grid_ref`. The
        //! coarsened data is saved into a new MultiFab with coarsened version
        //! of the LSFactory BoxArray `ls_ba`, and DistributionMapping `dm_ls`.
        //! Returns a unique_ptr to the new MultiFab.
        std::unique_ptr<MultiFab>  coarsen_data() const;

        //! Return level-set refinement
        int get_ls_ref() const {return ls_grid_ref;};
        //! Return level-set ghost cells (in units of the refined grid)
        int get_ls_pad() const {return ls_grid_pad;};
        //! Return EB grid refinement
        int get_eb_ref() const {return eb_grid_ref;};
        //! Return EB-grid ghost cells (in units of the refined eb-grid)
        int get_eb_pad() const {return eb_grid_pad;};

        //! Return AMR level
        int get_amr_level() const {return amr_lev;};


        //! Return BoxArrays:
        //!    -> ls: levelset (nodal)
        //!    -> cc: cell-centered (dual to ls) used by valid
        //!    -> eb: BoxArray used by EBIS (cell-centered)
        const BoxArray & get_ls_ba() const {return ls_ba;};
        const BoxArray & get_cc_ba() const {return cc_ba;};
        const BoxArray & get_eb_ba() const {return eb_ba;};
        const BoxArray & get_ba()    const {return base_ba;};

        const Geometry & get_geom()          const {return base_geom;};
        const Geometry & get_ls_geom()       const {return geom_ls;};
        const Geometry & get_eb_geom()       const {return geom_eb;};
        const DistributionMapping & get_dm() const {return ls_dm;};
};



//! Template class that handles filling implicit functions from GeometryShop.
//! The implicit function is filled into a MultiFab with the same CellSize,
//! BoxArray, DistributionMapping, and number of ghost cells as the LSFactory.
template <class F>
class GShopLSFactory {
    private:
        EB2::GeometryShop<F> m_gshop;
        Geometry             m_geom;
        BoxArray             m_ls_ba;
        DistributionMapping  m_ls_dm;
        int                  m_ls_pad;

    public:
        GShopLSFactory( const EB2::GeometryShop<F> & gshop, const LSFactory & level_set )
            : m_gshop(gshop),
              m_geom(level_set.get_ls_geom()),
              m_ls_ba(level_set.get_ls_ba()),
              m_ls_dm(level_set.get_dm()),
              m_ls_pad(level_set.get_ls_pad())
        {}

        GShopLSFactory (const EB2::GeometryShop<F> & gshop, const Geometry & geom,
                        const BoxArray & ba, const DistributionMapping & dm, int ls_pad)
            : m_gshop(gshop), m_geom(geom), m_ls_ba(ba), m_ls_dm(dm), m_ls_pad(ls_pad)
        {
            // Ensure that the inputted (refence) BoxArray is nodal
            BoxArray ba_nd = amrex::convert(ba, IntVect::TheNodeVector());
            m_ls_ba = ba_nd;
        }

        GShopLSFactory (const EB2::GeometryShop<F> & gshop, const Geometry & geom_crse,
                        const BoxArray & ba_crse, const DistributionMapping & dm,
                        int ls_pad, int ls_ref)
            : m_gshop(gshop), m_ls_dm(dm), m_ls_pad(ls_pad)
        {
            // Set up refined geometry
            Box dom = geom_crse.Domain();
            dom.refine(ls_ref);
            Geometry geom(dom);
            m_geom = geom;

            // Set up refined BoxArray. NOTE: ensure that the inputted
            // (reference) BoxArray is nodal
            BoxArray ba = amrex::convert(ba_crse, IntVect::TheNodeVector());
            ba.refine(ls_ref);
            m_ls_ba = ba;
        }

        //! Fill MultiFab (with the same CellSize, BoxArray,
        //! DistributionMapping, and number of ghost cells as the LSFactory)
        //! with the implicit function values. The filled MultiFab is returned
        //! as a std::unique_ptr<MultiFab>.
        std::unique_ptr<MultiFab> fill_impfunc() {
            std::unique_ptr<MultiFab> mf_impfunc = std::unique_ptr<MultiFab>(new MultiFab);
            mf_impfunc->define(m_ls_ba, m_ls_dm, 1, m_ls_pad);

#ifdef _OPENMP
#pragma omp parallel
#endif
            for(MFIter mfi(* mf_impfunc); mfi.isValid(); ++ mfi)
                m_gshop.fillFab((* mf_impfunc)[mfi], m_geom);

            mf_impfunc->FillBoundary(m_geom.periodicity());

            return mf_impfunc;
        }
};



//! Collection of static functions that help with the construction of LSFactory
//! objects.
class LSUtility {
    public:


    /****************************************************************************
     *                                                                          *
     * Utility functions for creating geometry objects with the correct Domain  *
     * and CellSize as requires by the LSFactory object.                        *
     *                                                                          *
     ***************************************************************************/

    //! Construct geometry with the Domain (Box) inherited from the Geometry
    //! `geom`, and CellSize matching the LSFactory `level_set` resolution.
    static Geometry make_ls_geometry(const LSFactory & level_set, const Geometry & geom){
        // Set up refined geometry
        Box dom_ls = geom.Domain();
        dom_ls.refine(level_set.get_ls_ref());
        Geometry geom_ls(dom_ls);

        return geom_ls;
    }


    //! Construct geometry with the Domain (Box) inherited from the Geometry
    //! `geom`, and CellSize matching the LSFactory `level_set` EB resolution.
    //! This is particularly useful for building EBFArrayBoxFactory objects
    //! whose Geometry is what's expected by the LSFactory's filling routines.
    static Geometry make_eb_geometry(const LSFactory & level_set, const Geometry & geom){
        // Set up refined geometry
        Box dom_eb = geom.Domain();
        dom_eb.refine(level_set.get_eb_ref());
        Geometry geom_eb(dom_eb);

        return geom_eb;
    }


    //! Returns the min CellSize of the Geometry object `geom`.
    static Real min_dx (const Geometry & geom) {

        const Real * dx = geom.CellSize();

#if (AMREX_SPACEDIM == 1)
        const Real min_dx = dx[0];
#elif (AMREX_SPACEDIM == 2)
        const Real min_dx = std::min(dx[0], dx[1]);
#elif (AMREX_SPACEDIM == 3)
        const Real min_dx = std::min(dx[0], std::min(dx[1], dx[2]));
#endif

        return min_dx;
    }


    //! Returns the minimum eb_pad (EB ghost cells) to ensure that the level-set
    //! is able to resolve a length scale of at least `min_resolution`.
    static int min_eb_pad (Real min_resolution, const Geometry & geom_eb) {

        const Real min_dx = LSUtility::min_dx(geom_eb);

        return min_resolution/min_dx;
    }


    //! Debug utility to help interrogate runtime flag types.
    static void PrintFlagType(const Box & tile_box, const EBCellFlagFab & flag){

        if(flag.getType(tile_box) == FabType::covered)
            amrex::Print() << "flag covered" << std::endl;
        else if(flag.getType(tile_box) == FabType::regular)
            amrex::Print() << "flag regular" << std::endl;
        else if(flag.getType(tile_box) == FabType::singlevalued)
            amrex::Print() << "flag singlevalued" << std::endl;
        else if(flag.getType(tile_box) == FabType::multivalued)
            amrex::Print() << "flag multivalued" << std::endl;
        else if(flag.getType(tile_box) == FabType::undefined)
            amrex::Print() << "flag undefined" << std::endl;
        else
            amrex::Print() << "flag ill-formed" << std::endl;
    }
};



}

#endif
