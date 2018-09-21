
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>

#include <cmath>
#include <algorithm>

#include "MyTest.H"
#include "MyEB.H"

using namespace amrex;

void
MyTest::initializeEB ()
{
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "rotated_box")
    {
        EB2::BoxIF box({AMREX_D_DECL(0.25,0.25,0.25)},
                       {AMREX_D_DECL(0.75,0.75,0.75)}, true);
        auto gshop = EB2::makeShop(EB2::translate(
                                       EB2::rotate(
                                           EB2::translate(box, {AMREX_D_DECL(-0.5,-0.5,-0.5)}),
                                           std::atan(1.0)*0.3, 2),
                                       {AMREX_D_DECL(0.5,0.5,0.5)}));
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);        
    }
    else if (geom_type == "flower")
    {
        FlowerIF flower(0.3, 0.15, 6, {AMREX_D_DECL(0.5,0.5,0.5)}, true);
#if (AMREX_SPACEDIM == 2)
        auto gshop = EB2::makeShop(flower);
#else
        EB2::PlaneIF planelo({0.,0.,0.1},{0.,0., -1.});
        EB2::PlaneIF planehi({0.,0.,0.9},{0.,0.,  1.});
        auto gshop = EB2::makeShop(EB2::makeUnion(flower,planelo,planehi));
#endif
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
    }
#if (AMREX_SPACEDIM == 3)
    else if (geom_type == "csg")
    {
        EB2::SphereIF sphere(0.5, {0.5,0.5,0.5}, false);
        EB2::BoxIF cube({0.1,0.1,0.1}, {0.9,0.9,0.9}, false);
        auto cubesphere = EB2::makeIntersection(sphere, cube);
        
        EB2::CylinderIF cylinder_x(0.25, 0, {0.5,0.5,0.5}, false);
        EB2::CylinderIF cylinder_y(0.25, 1, {0.5,0.5,0.5}, false);
        EB2::CylinderIF cylinder_z(0.25, 2, {0.5,0.5,0.5}, false);
        auto three_cylindres = EB2::makeUnion(cylinder_x, cylinder_y, cylinder_z);

        auto csg = EB2::makeDifference(cubesphere, three_cylindres);

        auto gshop = EB2::makeShop(csg);
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
    }
#endif
    else
    {
        EB2::Build(geom.back(), max_level, max_level+max_coarsening_level);
    }
}
