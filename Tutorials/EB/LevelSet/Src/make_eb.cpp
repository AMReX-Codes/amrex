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
#include <eb_levelset.H>

#include <make_shapes.H>

void
make_my_eb(int lev, const BoxArray& grids, const DistributionMapping& dmap,
           const Geometry& geom, LSFactory * level_set)
{
    if (geom.isAllPeriodic()) return;

    // Implicit functions for:
    //    * impfunc       -> all EBs in the domain
    //    * impfunc_poly2 -> EBs belonging to the (polynomial) walls
    std::unique_ptr<BaseIF> impfunc;
    std::unique_ptr<BaseIF> impfunc_poly2;
    std::unique_ptr<BaseIF> impfunc_cyl;

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

    GeometryShop gshop;

   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/
   /* POLY GEOMETRY                                                            */
   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/

    if (make_poly) {

      amrex::Print() << "Using poly2 geometry" << std::endl;
      impfunc_poly2 = make_poly_geom(lev, SpaceDim, "poly2");
      impfunc.reset(impfunc_poly2->newImplicitFunction());

      // Define components of the GeometryShop separately:
      gshop.define(* impfunc_poly2, eb_verbosity);

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
      RealVect translation;

      {
         ParmParse pp("eb");
         pp.get("cyl_dir",cyl_dir);
         pp.get("cyl_radius",cyl_radius);
         pp.get("cyl_length",cyl_length);
         pp.getarr("cyl_translate", transvec,  0, 3);
      }

      translation = RealVect(transvec);

      int         lev  = 0;
      bool water_tight = false;

      impfunc_cyl = make_cylinder_geom(cyl_dir, cyl_radius, cyl_length, translation, lev, water_tight);

      impfunc.reset(impfunc_cyl->newImplicitFunction());

      // Define components of the GeometryShop separately:
      gshop.define(* impfunc_cyl, eb_verbosity);
    }

   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/
   /* FOR BOTH                                                                 */
   /****************************************************************************/
   /****************************************************************************/
   /****************************************************************************/

   /****************************************************************************
    *                                                                          *
    * Fill Level-set using:                                                    *
    *      -> Planes (where the GeometryShop's implicit function is a signed   *
    *         distance): implicit function's value                             *
    *      -> Poly2 (where GeometryShop's implicit function is singed but not  *
    *         a distance): min distance to EB facets                           *
    * Note: this requires building and destroying the EBTower (twice), so any  *
    * EBTower data built before this will be lost...                           *
    *                                                                          *
    ****************************************************************************/

    Geometry geom_eb = LSUtility::make_eb_geometry(* level_set, geom);

    // Note GeometryShop's behaviour wrt anisotropic cells: * use x-component of dx as reference length-scale
    //                                                      * rescale y, z- components wrt to dx[0] (dx(1))
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1, above]
                                   gshop,  // ............ GeometryShop object
                                   grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an
    //         EBFArrayBoxFactory which defines our EB surface now
    //          => define the level set as the (signed) distance to the
    //             closest point on the EB-facets
    int eb_pad = level_set->get_eb_pad();
    EBFArrayBoxFactory eb_factory(geom_eb, level_set->get_eb_ba(), dmap,
                                  {eb_pad, eb_pad, eb_pad}, EBSupport::full);
    level_set->intersection_ebf(eb_factory, * AMReX_EBIS::instance());
    EBTower::Destroy();

    /***************************************************************************
     *                                                                         *
     * Build standard EB Factories                                             *
     *                                                                         *
     ***************************************************************************/

    // dx : cell size used by EBIndexSpace
    Box domain(geom.Domain());
    Real dx = geom.CellSize()[0];

    GeometryShop gshop(*impfunc, eb_verbosity);
    AMReX_EBIS::instance()->define(domain, RealVect::Zero, dx, gshop, grid_size, max_level);

    // set up ebfactory
    int m_eb_basic_grow_cells = 2;
    int m_eb_volume_grow_cells = 2;
    int m_eb_full_grow_cells = 2;
    EBSupport m_eb_support_level = EBSupport::full;

    EBTower::Build();

    std::unique_ptr<amrex::EBFArrayBoxFactory> ebfactory =
                         std::unique_ptr<EBFArrayBoxFactory>(
                         new EBFArrayBoxFactory(geom, grids, dmap,
                         {m_eb_basic_grow_cells, m_eb_volume_grow_cells, m_eb_full_grow_cells},
                         m_eb_support_level));
}
