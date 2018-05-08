#include <AMReX_GeometryShop.H>
#include <AMReX_SphereIF.H>
#include <AMReX_PlaneIF.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_FlatPlateGeom.H>
#include <AMReX_EBISLayout.H>
#include <AMReX_EBGraph.H>
#include <AMReX_EBDebugOut.H>
#include <AMReX_EBCellFAB.H>
#include <AMReX_EBCellFactory.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBIndexSpace.H>
#include <AMReX_UnionIF.H>
#include <AMReX_TransformIF.H>
#include <AMReX_ComplementIF.H>
#include <AMReX_IntersectionIF.H>
#include <AMReX_LatheIF.H>
#include <AMReX_PolynomialIF.H>
#include <AMReX_AnisotropicDxPlaneIF.H>
#include <AMReX_AnisotropicIF.H>

#include <AMReX_ParmParse.H>

#include <AMReX_EBTower.H>

#include <algorithm>
#include <AMReX_EB_levelset.H>

#include <make_shapes.H>

void
make_my_eb(int lev, const BoxArray & grids, const DistributionMapping & dmap,
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

      std::unique_ptr<BaseIF> impfunc = make_poly_geom(lev, SpaceDim, "poly2");

      Geometry geom_eb = LSUtility::make_eb_geometry(* level_set, geom);
      GeometryShop gshop(* impfunc, eb_verbosity);
      AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                     RealVect::Zero,  // ......... origin of EBIndexSpace
                                     geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace
                                     gshop,  // ............ GeometryShop object
                                     grid_size, max_level);

      EBTower::Build();
      int eb_pad = level_set->get_eb_pad();
      EBFArrayBoxFactory eb_factory(geom_eb, level_set->get_eb_ba(), dmap,
                                   {eb_pad, eb_pad, eb_pad}, EBSupport::full);

      //Fill Level-set using eb_factory
      level_set->intersection_ebf(eb_factory, * AMReX_EBIS::instance());



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
      make_cylinder_geom(cyl_dir, cyl_radius, cyl_length, translation,
                         lev, geom, dmap, level_set);
    }
}
