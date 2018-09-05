
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
#include <AMReX_EB2_IF_Rotation.H>
#include <AMReX_EB2_IF_Scale.H>
#include <AMReX_EB2_IF_Translation.H>
#include <AMReX_EB2_IF_Lathe.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Plane.H>

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
    else
    {
        EB2::Build(geom.back(), max_level, max_level+max_coarsening_level);
    }
}
