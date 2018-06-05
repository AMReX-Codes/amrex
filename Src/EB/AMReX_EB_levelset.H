#ifndef AMREX_EB_LEVELSET_H_
#define AMREX_EB_LEVELSET_H_

#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_RealVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_PolynomialIF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

class LSFactory {
    private:

        // MultiFab pointers storing level-set data
        //    -> ls_phi:   nodal MultiFab storing signed distance function to the nearest wall
        //    -> ls_valid: cell-centered iMultiFab storing integer flags assuming the following values:
        //         -1 : not all nodes around cell have been initialized
        //          0 : none of the cell's neighbours contain negative vlaues of ls_phi on its nodes
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

        // Baseline BoxArray and Geometry from which refined quantities are
        // derived. These are mainly kept around for the copy constructor.
        BoxArray base_ba;
        Geometry base_geom;

        // Cell size of the refined grids
        const RealVect dx_vect, dx_eb_vect;

        // Construct eb_facets list from eb_factory
        std::unique_ptr<Vector<Real>> eb_facets(const EBFArrayBoxFactory & eb_factory);

        // Construct nodal MultiFab containing EBIS Implicit function velues
        std::unique_ptr<MultiFab> ebis_impfunc(const EBIndexSpace & eb_is);

        void update_ba(const BoxArray& new_ba, const DistributionMapping & dm);
        void init_geom(const BoxArray& ba, const Geometry& geom, const DistributionMapping & dm);

    public:
        LSFactory(int lev, int ls_ref, int eb_ref, int ls_pad, int eb_pad,
                  const BoxArray& ba, const Geometry& geom, const DistributionMapping& dm);
        LSFactory(const LSFactory & other);
        ~LSFactory();

        void regrid(const BoxArray& ba, const DistributionMapping& dm);

        void invert();

        std::unique_ptr<iMultiFab> intersection_ebf(const EBFArrayBoxFactory & eb_factory, const EBIndexSpace & eb_is);
        std::unique_ptr<iMultiFab> union_ebf(const EBFArrayBoxFactory & eb_factory, const EBIndexSpace & eb_is);
        std::unique_ptr<iMultiFab> intersection_ebis(const EBIndexSpace & eb_is);
        std::unique_ptr<iMultiFab> union_ebis(const EBIndexSpace & eb_is);

        // Unpdate internal `ls_phi` and `ls_valid` based on the min of
        // `ls_phi` and `ls_in`
        void update_intersection(const MultiFab & ls_in, const iMultiFab & valid_in);
        void update_union(const MultiFab & ls_in, const iMultiFab & valid_in);

        // Return level-set data
        const MultiFab *  get_data() const {return ls_grid.get();};
        const iMultiFab * get_valid() const {return ls_valid.get();};

        // Copy level-set data
        std::unique_ptr<MultiFab>  copy_data(const DistributionMapping& dm) const;
        std::unique_ptr<iMultiFab> copy_valid(const DistributionMapping& dm) const;
        std::unique_ptr<MultiFab>  coarsen_data() const;

        // Return grid parameters
        const int get_ls_ref() const {return ls_grid_ref;};
        const int get_ls_pad() const {return ls_grid_pad;};
        const int get_eb_ref() const {return eb_grid_ref;};
        const int get_eb_pad() const {return eb_grid_pad;};

        // Return AMR level
        const int get_amr_level() const {return amr_lev;};


        // Return BoxArrays:
        //    -> ls: levelset (nodal)
        //    -> cc: cell-centered (dual to ls) used by valid
        //    -> eb: BoxArray used by EBIS (cell-centered)
        const BoxArray & get_ls_ba() const {return ls_ba;};
        const BoxArray & get_cc_ba() const {return cc_ba;};
        const BoxArray & get_eb_ba() const {return eb_ba;};
        const BoxArray & get_ba()    const {return base_ba;};

        const Geometry & get_geom()          const {return base_geom;};
        const DistributionMapping & get_dm() const {return ls_dm;};
};



class LSUtility {
    public:
        static Geometry make_ls_geometry(const LSFactory & level_set, const Geometry& geom){
            // Set up refined geometry
            Box dom_ls = geom.Domain();
            dom_ls.refine(level_set.get_ls_ref());
            Geometry geom_ls(dom_ls);

            return geom_ls;
        }

        static Geometry make_eb_geometry(const LSFactory & level_set, const Geometry& geom){
            // Set up refined geometry
            Box dom_eb = geom.Domain();
            dom_eb.refine(level_set.get_eb_ref());
            Geometry geom_eb(dom_eb);

            return geom_eb;
        }


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



class PolynomialDF : public PolynomialIF {
    private:
        int order;

    public:
        PolynomialDF(const Vector<PolyTerm> & a_polynomial,
                     const bool             & a_inside);


        virtual Real value(const RealVect         & a_point,
                           const Vector<PolyTerm> & a_polynomial) const;


        virtual Real value(const RealVect & a_point) const;


        virtual BaseIF * newImplicitFunction() const;

};

}

#endif
