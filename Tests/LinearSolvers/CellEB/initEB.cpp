
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

    if (geom_type == "combustor")
    {
        amrex::Abort("initializeEB: todo");
    }
    else if (geom_type == "rotated_box")
    {
        EB2::BoxIF box({AMREX_D_DECL(0.25,0.25,0.25)},
                       {AMREX_D_DECL(0.75,0.75,0.75)}, false);
        auto gshop = EB2::makeShop(EB2::translate(
                                       EB2::rotate(
                                           EB2::translate(box, {AMREX_D_DECL(-0.5,-0.5,-0.5)}),
                                           std::atan(1.0)*0.3, 2),
                                       {AMREX_D_DECL(0.5,0.5,0.5)}));
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);        
    }
    else if (geom_type == "two_spheres")
    {
        EB2::SphereIF sphere1(0.2, {AMREX_D_DECL(0.45, 0.4, 0.58)}, false);
        EB2::SphereIF sphere2(0.2, {AMREX_D_DECL(0.55, 0.42, 0.6)}, false);
        auto twospheres = EB2::makeUnion(sphere1, sphere2);
        auto gshop = EB2::makeShop(twospheres);
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
    }
    else if (geom_type == "two_spheres_one_box")
    {
        EB2::SphereIF sphere1(0.2, {AMREX_D_DECL(0.5, 0.48, 0.5)}, false);
        EB2::SphereIF sphere2(0.2, {AMREX_D_DECL(0.55, 0.58, 0.5)}, false);
        EB2::BoxIF box({AMREX_D_DECL(0.25,0.75,0.5)}, {AMREX_D_DECL(0.75,0.8,0.75)}, false);
        auto twospheres = EB2::makeUnion(sphere1, sphere2, box);
        auto gshop = EB2::makeShop(twospheres);
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
    }
    else if (geom_type == "flower")
    {
        FlowerIF flower(0.2, 0.1, 6, {AMREX_D_DECL(0.5,0.5,0.5)}, false);
        auto gshop = EB2::makeShop(flower);
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
    }
    else
    {
        EB2::Build(geom.back(), max_level, max_level+max_coarsening_level);
    }
}
