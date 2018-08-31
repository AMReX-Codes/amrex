
#include <AMReX_EBFabFactory.H>
#include <AMReX_EB2_IF_Polynomial.H>
#include <AMReX_EB_levelset.H>

#include <AMReX_ParmParse.H>

#include <algorithm>

#include <make_shapes.H>

using namespace amrex;


void
make_my_eb2(int lev, const BoxArray & grids, const DistributionMapping & dmap,
            const Geometry & geom, LSFactory * level_set)
{
    if (geom.isAllPeriodic()) return;

    ParmParse pp("eb");

    bool make_poly     = false;
    bool make_cylinder = false;

    pp.query("make_poly", make_poly);
    pp.query("make_cylinder", make_cylinder);

    if (make_poly && make_cylinder)
        amrex::Abort("Must specify only one of make_poly OR make_cylinder");
    if (!make_poly && !make_cylinder)
        amrex::Abort("Must specify one of make_poly or make_cylinder");


    bool eb_verbosity = true;
    int max_level = 0;
    int grid_size = 16;


   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/
   /* POLY GEOMETRY                                                            */
   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/

    if (make_poly) {

      amrex::Print() << "Using poly2 geometry" << std::endl;

      // Construct Polynomial implicit function, and GeometryShop object
      // GeometryShop's PolynomialIF is not a signed distance function...
      //      => it's easier to use PolynomialIF and to build an
      //         EBFArrayBoxFactory which defines our EB surface now
      //          => define the level set as the (signed) distance to the
      //             closest point on the EB-facets

      std::unique_ptr<CylinderIF> impfunc = make_poly_eb2_geom(lev, SpaceDim, "poly2");

      Geometry geom_eb = LSUtility::make_eb_geometry(* level_set, geom);
      EB2::GeometryShop<CylinderIF> gshop(* impfunc);

      GShopLSFactory<CylinderIF>    cylinder_ls_gshop(gshop, * level_set);

      // Implicit function used by LSFactory
      //  -- returned MF has the same DM as LSFactory
      std::unique_ptr<MultiFab> cylinder_mf_impfunc = cylinder_ls_gshop.fill_impfunc();


      EB2::Build(gshop, geom_eb, max_level, max_level);

      const EB2::IndexSpace & cylinder_ebis = EB2::IndexSpace::top();
      const EB2::Level &      cylinder_lev  = cylinder_ebis.getLevel(geom);


      int eb_pad = level_set->get_eb_pad();
      EBFArrayBoxFactory eb_factory(cylinder_lev, geom_eb, level_set->get_eb_ba(), dmap,
                                   {eb_pad, eb_pad, eb_pad}, EBSupport::full);

      //Fill Level-set using eb_factory
      level_set->intersection_ebf(eb_factory, * cylinder_mf_impfunc);



   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/
   /* CYLINDER GEOMETRY                                                        */
   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/

    } else if (make_cylinder) {

      amrex::Print() << "Using cylinder geometry" << std::endl;

      int cyl_dir;
      Real cyl_radius;
      Real cyl_length;
      Vector<Real> transvec(3);

      {
         ParmParse pp("cyl");
         pp.get("cyl_dir",cyl_dir);
         pp.get("cyl_radius",cyl_radius);
         pp.get("cyl_length",cyl_length);
         pp.getarr("cyl_translate", transvec,  0, 3);
      }

      RealVect translation = RealVect(transvec);

      // The make_cylinder_geom function unions the cylinder level-set function
      // "onto" the level-set factory input argument (`level_set`). As
      // `level_set` is initialized to fortran `huge`, we have to invert it
      // first. Otherwise the uninion operation (which takes a max) ignores the
      // new level-set.
      level_set->invert();

      int lev  = 0;
      make_cylinder_eb2_geom(cyl_dir, cyl_radius, cyl_length, translation,
                             lev, geom, dmap, level_set);
    }
}
