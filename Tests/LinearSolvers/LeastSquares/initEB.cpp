
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
    constexpr Real pi = 3.1415926535897932;

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
    else if (geom_type == "channel") {
        bool fluid_inside = true;
        pp.get("channel_has_fluid_inside", fluid_inside);
#if (AMREX_SPACEDIM == 2)
        Vector<Real> pt_on_top_wall(3);
        Real height;
        Real rotation = 0.0;

        pp.getarr("channel_pt_on_top_wall", pt_on_top_wall, 0, 3);
        pp.get("channel_rotation", rotation);
        pp.get("channel_height", height);
        rotation = (rotation/180.) * pi;

        EB2::PlaneIF left({AMREX_D_DECL(pt_on_top_wall[0],pt_on_top_wall[1],0.0)},
                         {AMREX_D_DECL(-std::sin(rotation),std::cos(rotation),0.0)},
                         fluid_inside);
        EB2::PlaneIF right({AMREX_D_DECL(pt_on_top_wall[0], pt_on_top_wall[1] - (height/std::cos(rotation)),0.0)},
                         {AMREX_D_DECL(std::sin(rotation),-std::cos(rotation),0.0)},
                         fluid_inside);
        auto channel = EB2::makeUnion(left, right);
        auto gshop = EB2::makeShop(channel);
        EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
#else
        bool askew = false;
        pp.get("channel_askew", askew);

        if(askew) {
           Vector<Real> pt_on_top_wall(3);
           Real height;
           Vector<Real> rotation(2);

           pp.getarr("channel_pt_on_top_wall", pt_on_top_wall, 0, 3);
           pp.getarr("channel_rotation", rotation, 0, 2);
           pp.get("channel_height", height);

           Real alpha = (rotation[0]/180.) * pi;
           Real gamma = (rotation[1]/180.) * pi;

           Vector<Real> norm(3);
           Real norm_mag = std::sqrt(std::sin(alpha)*std::sin(alpha) +
                      std::cos(alpha)*std::cos(alpha)*std::cos(gamma)*std::cos(gamma) +
                      std::sin(gamma)*std::sin(gamma));
           norm[0] = -std::sin(gamma) / norm_mag;
           norm[1] = std::cos(alpha)*std::cos(gamma) / norm_mag;
           norm[2] = -std::sin(alpha) / norm_mag;

           EB2::PlaneIF left({AMREX_D_DECL(pt_on_top_wall[0],pt_on_top_wall[1],pt_on_top_wall[2])},
                         {AMREX_D_DECL(norm[0],norm[1],norm[2])},
                         fluid_inside);

           EB2::PlaneIF right({AMREX_D_DECL(pt_on_top_wall[0] - height*norm[0],
                                            pt_on_top_wall[1] - height*norm[1],
                                            pt_on_top_wall[2] - height*norm[2])},
                         {AMREX_D_DECL(-norm[0],-norm[1],-norm[2])},
                         fluid_inside);

           auto channel = EB2::makeUnion(left, right);
           auto gshop = EB2::makeShop(channel);
           EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
        }
        else {
           Vector<Real> lo(3);
           Vector<Real> hi(3);
           pp.getarr("channel_lo", lo, 0, 3);
           pp.getarr("channel_hi", hi, 0, 3);

           EB2::BoxIF channel({AMREX_D_DECL(lo[0],lo[1],lo[2])}, {AMREX_D_DECL(hi[0],hi[1],hi[2])}, fluid_inside);
           auto gshop = EB2::makeShop(channel);
           EB2::Build(gshop, geom.back(), max_level, max_level+max_coarsening_level);
        }
#endif
    }
    else
    {
        EB2::Build(geom.back(), max_level, max_level+max_coarsening_level);
    }
}
