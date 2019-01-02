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

        // MultiFab pointers storing level-set data
        //    -> ls_phi: nodal MultiFab storing signed distance function to the
        //               nearest wall
        //    -> ls_valid: cell-centered iMultiFab storing integer flags
        //                 assuming the following values:
        //         -1 : not all nodes around cell have been initialized
        //          0 : none of the cell's neighbours contain negative vlaues of
        //              ls_phi on its nodes
        //          1 : the cell is in the neighbourhood of phi < 0

        std::unique_ptr<MultiFab>  ls_grid;
        std::unique_ptr<iMultiFab> ls_valid;

        // AMR information (needed for EBISLevel)
        const int amr_lev;

        DistributionMapping ls_dm;

        // Refined box array (over which the level-set MultiFab is constructed)
        BoxArray ls_ba, cc_ba, eb_ba;
        const int ls_grid_ref, eb_grid_ref;
        // Data structures encloding (padded) EB index space
        std::unique_ptr<MultiFab> eb_grid;

        // Geometry objects reflecting refined (and padded) BoxArrays. This is
        // required to prevent PeriodicFill operations from "spilling" refined
        // indices over periodic boundaries
        Geometry geom_ls, geom_eb;

        // Padding of the leve-set and EB grids:
        const int ls_grid_pad, eb_grid_pad;

        // Tiling for local level-set filling
        int eb_tile_size;

        // Baseline BoxArray and Geometry from which refined quantities are
        // derived. These are mainly kept around for the copy constructor.
        BoxArray base_ba;
        Geometry base_geom;

        // Cell size of the refined grids
        const RealVect dx_vect, dx_eb_vect;

        std::unique_ptr<Vector<Real>> eb_facets(const EBFArrayBoxFactory & eb_factory);
        std::unique_ptr<Vector<Real>> eb_facets(const EBFArrayBoxFactory & eb_factory,
                                                const BoxArray & ba, const DistributionMapping & dm);

        void update_ba(const BoxArray & new_ba, const DistributionMapping & dm);
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

        //! Fills (overwrites) level-set data locally. Same as `fill_data_loc`
        //! but with cubic EB search box.
        static void fill_data (MultiFab & data, iMultiFab & valid,
                               const EBFArrayBoxFactory & eb_factory,
                               const MultiFab & eb_impfunc,
                               int ebt_size, int ls_ref, int eb_ref,
                               const Geometry & geom, const Geometry & geom_eb);

        //! Fills (overwrites) level-set data locally. Also fills iMultiFab
        //! tagging cells which are nearby to EB surface. Only EB facets in a
        //! box size of `ebt_size` are considered. Any EB facets that are
        //! outside this box are ignored => the min/max value of the level-set
        //! are +/- `eb_pad * min(geom_eb.CellSize(:))`.
        static void fill_data (MultiFab & data, iMultiFab & valid,
                               const EBFArrayBoxFactory & eb_factory,
                               const MultiFab & eb_impfunc,
                               const IntVect & ebt_size, int ls_ref, int eb_ref,
                               const Geometry & geom, const Geometry & geom_eb);


        /************************************************************************
         *                                                                      *
         * Public member functions to update level-set data                     *
         *                                                                      *
         ***********************************************************************/

        //! Fills (overwrites) level-set data locally. The local EB search
        //! radius is managed by LSFactory. Returns: iMultiFab indicating region
        //! that has been filled by a valid level-set function (i.e. the value
        //! of the level-set was informed by nearby EB facets)
        std::unique_ptr<iMultiFab> fill(const EBFArrayBoxFactory & eb_factory,
                                        const MultiFab & mf_impfunc);

        //! Fills (overwrites) level-set data locally. The local EB search
        //! radius is specified by the user corresponding to a cube with
        //! edge-length = `ebt_size`. Returns: iMultiFab indicating region
        //! that has been filled by a valid level-set function (i.e. the value
        //! of the level-set was informed by nearby EB facets)
        std::unique_ptr<iMultiFab> fill(const EBFArrayBoxFactory & eb_factory,
                                        const MultiFab & mf_impfunc,
                                        int ebt_size);

        //! Fills (overwrites) level-set data locally. The local EB search
        //! radius is specified by the user corresponding to a box with edge
        //! lengths given by the IntVect `ebt_size`. Returns: iMultiFab indicating region
        //! that has been filled by a valid level-set function (i.e. the value
        //! of the level-set was informed by nearby EB facets)
        std::unique_ptr<iMultiFab> fill(const EBFArrayBoxFactory & eb_factory,
                                        const MultiFab & mf_impfunc,
                                        const IntVect & ebt_size);

        std::unique_ptr<iMultiFab> intersection_ebf(const EBFArrayBoxFactory & eb_factory,
                                                    const MultiFab & mf_impfunc);
        std::unique_ptr<iMultiFab> union_ebf(const EBFArrayBoxFactory & eb_factory,
                                             const MultiFab & mf_impfunc);
        std::unique_ptr<iMultiFab> intersection_impfunc(const MultiFab & mf_impfunc);
        std::unique_ptr<iMultiFab> union_impfunc(const MultiFab & mf_impfunc);

        // Update internal `ls_phi` and `ls_valid` based on the min of `ls_phi`
        // and `ls_in`
        void update_intersection(const MultiFab & ls_in, const iMultiFab & valid_in);
        void update_union(const MultiFab & ls_in, const iMultiFab & valid_in);

        // Return level-set data
        const MultiFab *  get_data() const {return ls_grid.get();};
        const iMultiFab * get_valid() const {return ls_valid.get();};

        // Copy level-set data
        std::unique_ptr<MultiFab>  copy_data(const DistributionMapping & dm) const;
        std::unique_ptr<iMultiFab> copy_valid(const DistributionMapping & dm) const;
        std::unique_ptr<MultiFab>  coarsen_data() const;

        // Return grid parameters
        int get_ls_ref() const {return ls_grid_ref;};
        int get_ls_pad() const {return ls_grid_pad;};
        int get_eb_ref() const {return eb_grid_ref;};
        int get_eb_pad() const {return eb_grid_pad;};

        // Return AMR level
        int get_amr_level() const {return amr_lev;};


        // Return BoxArrays:
        //    -> ls: levelset (nodal)
        //    -> cc: cell-centered (dual to ls) used by valid
        //    -> eb: BoxArray used by EBIS (cell-centered)
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
            for(MFIter mfi(* mf_impfunc, true); mfi.isValid(); ++ mfi)
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



    //! Returns the minimum eb_pad (EB ghost cells) to ensure that the level-set
    //! is able to resolve a length scale of at least `min_resolution`.

    static int min_eb_pad (Real min_resolution, const Geometry & geom_eb) {

        RealVect dx_eb(AMREX_D_DECL(geom_eb.CellSize()[0],
                                    geom_eb.CellSize()[1],
                                    geom_eb.CellSize()[2]));

#if (AMREX_SPACEDIM == 1)
        const Real min_dx = dx_eb[0];
#elif (AMREX_SPACEDIM == 2)
        const Real min_dx = std::min(dx_eb[0], dx_eb[1]);
#elif (AMREX_SPACEDIM == 3)
        const Real min_dx = std::min(dx_eb[0], std::min(dx_eb[1], dx_eb[2]));
#endif

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
