#include <AMReX_PlaneIF.H>
#include <AMReX_AllRegularService.H>
#include <AMReX_TransformIF.H>
#include <AMReX_IntersectionIF.H>

#include <make_shapes.H>

using namespace amrex;

std::unique_ptr<BaseIF> make_cylinder_geom(int dir, Real radius, Real length, const RealVect & translation,
                                           int lev, bool water_tight, Geometry& geom, DistributionMapping& dmap, 
                                           std::unique_ptr<LSFactory> level_set) 
{
    // Construct a cylinder implicit function with finite radius and axis
    // offset (translation). The cylinder can be oriented along any cartesian
    // axis (dir).
    std::unique_ptr<BaseIF> cylinder_IF;

    // Polynomial defining (curved) cylinder walls parallel to a given axis:
    //     IF = a^2 + b^2 - R^2
    // where a, b \in {a, x, z} - {axis} for example, if the cylinder lies on
    // the y-axis => IF = x^2 + z^2 - R^2
    Vector<PolyTerm> poly;
    for(int idir = 0; idir < 3; idir++) {
        // Constucts the coefficient vector describing a cylinder with
        // orientation given by axis `dir`:
        //    *  coefvec[0] = R^2 term
        //    *  coefvec[2] = {x,y,z}^2 term
        Vector<Real> coefvec(3);
        if( idir == dir) coefvec = { - std::pow(radius, 2), 0. ,0.};
        else             coefvec = {0., 0., 1};

        for(int lc = 0; lc < 3; lc++) {
            // x^(lc) term
            IntVect powers = IntVect::Zero;
            powers[idir] = lc;
            PolyTerm mono = {.coef = coefvec[lc], .powers = powers};
            poly.push_back(mono);
        }
    }


    // Internal flow cylinder
    PolynomialIF cylinder0(poly, true);
    TransformIF cylinder1(cylinder0);
    cylinder1.translate(translation);

    // box to clip to correct length
    RealVect normal = RealVect::Zero , center = RealVect::Zero;
    Vector<std::unique_ptr<BaseIF>> planes;

    center[dir] = 0.0;
    normal[dir] = 1.0;
    planes.push_back(std::unique_ptr<BaseIF>(new PlaneIF(normal, center, true)));

    center[dir] = length;
    normal[dir] =-1.0;
    planes.push_back(std::unique_ptr<BaseIF>(new PlaneIF(normal, center, true)));

    // The IntersectionIF constructor requires a vector of pointers to the
    Vector<BaseIF *> plane_ptrs = {planes[0].get(), planes[1].get()};
    IntersectionIF bounding_box(plane_ptrs);
    TransformIF walls(bounding_box);
    walls.translate(translation);

    Vector<BaseIF * > funcs(2);
    funcs[0] = & cylinder0;
    funcs[1] = & bounding_box;

    IntersectionIF cylinder(funcs);

    TransformIF cylinder_trans(cylinder);
    cylinder_trans.translate(translation);

    cylinder_IF.reset(cylinder_trans.newImplicitFunction());

    // IF we are not using the water-tight mode, return now:

    if(! water_tight)
        return cylinder_IF;

    // ELSE construct the level-set by unioning each cylinder (intersected
    // component) of the CLR using the level-set
    //   => corners are much more cleanly resolved, but is much slower.

    int max_level = 0;
    int grid_size = 16;
    bool eb_verbosity = true;

    /**************************************************************************
     *                                                                        *
     * Fill Level-set using:                                                  *
     *      -> Walls (where the GeometryShop's implicit function is a signed  *
     *         distance): implicit function's value                           *
     *      -> Cylinder (where GeometryShop's implicit function is singed but *
     *         not a distance): min distance to EB facets                     *
     *      Note: this requires building and destroying the EBTower (twice),  *
     *      so any EBTower data built before this will be lost...             *
     *                                                                        *
     **************************************************************************/

    // Define both components of the GeometryShop separately:
    GeometryShop gshop_upoly(cylinder1, eb_verbosity);
    GeometryShop gshop_walls(walls, eb_verbosity);

    // Define temporary level-sets used for constructing the cylinder:
    LSFactory ls_cylinder(* level_set);
    LSFactory ls_walls(ls_cylinder);

    Geometry geom_eb = LSUtility::make_eb_geometry(ls_cylinder, geom);

    // Define the EBIS first using only the walls...
    // Note GeometryShop's behaviour wrt anisotropic cells: * use x-component of dx as reference length-scale
    //                                                      * rescale y, z- components wrt to dx[0] (dx(1))
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1]
                                   gshop_walls,  // ............ GeometryShop object
                                   grid_size, max_level);
    // [1]: EBIndexSpace internally assumes an isotropic grid. Any anisotropic
    // implicit function (e.g AnisotrpicPlaneIF) uses dx as a reference, and
    // rescales dy and dz wrt dx. => dx goes here.

    EBTower::Build();
    // GeometryShop's Planes' implicit function is actually a signed distance function
    //      => it's just easier to fill the level-set this way
    ls_cylinder.intersection_ebis(* AMReX_EBIS::instance());
    EBTower::Destroy();

    // Define the EBIS using only the poly (after deleting the walls-only EBTower)...
    // Note GeometryShop's behaviour wrt anisotropic cells: * use x-component of dx as reference length-scale
    //                                                      * rescale y, z- components wrt to dx[0] (dx(1))
    AMReX_EBIS::instance()->define(geom_eb.Domain(),
                                   RealVect::Zero,  // ......... origin of EBIndexSpace
                                   geom_eb.CellSize()[0],  // .. reference cell size of EBIndexSpace [1, above]
                                   gshop_upoly,  // ............ GeometryShop object
                                   grid_size, max_level);

    EBTower::Build();
    // GeometryShop's PolynomialIF is not a signed distance function...
    //      => it's easier to use PolynomialIF to build an EBFArrayBoxFactory
    //         which defines our EB surface
    //          => define the level set as the (signed) distance to the closest
    //             point on the EB-facets
    int eb_grow = level_set->get_eb_pad();
    EBFArrayBoxFactory eb_factory_poly(geom_eb, level_set->get_eb_ba(), dmap,
                                       {eb_grow, eb_grow, eb_grow}, EBSupport::full);

    // Only EB facets that are "in range" (within `n_pad` of the local
    // BoxArray) are considered for filling the EB level-set. flag_valid = {1,0}
    // indicates if a cell's BoxArray contained any EB facets (1), or if the
    // cell's value is invalid (0) because all EB facets where too far away in
    // order to be considered.
    std::unique_ptr<iMultiFab> flag_valid = ls_cylinder.intersection_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    ls_walls.intersection_ebf(eb_factory_poly, * AMReX_EBIS::instance());
    EBTower::Destroy();

    level_set->update_union(* ls_cylinder.get_data(), * flag_valid);

    return cylinder_IF;
}
