
#include <AMReX_Config.H>
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

    if (geom_type == "channel") { // For more geometries see LeastSquares
        bool fluid_inside = true;
        pp.get("channel_has_fluid_inside", fluid_inside);
#if (AMREX_SPACEDIM == 2)
        Vector<Real> pt_on_top_wall(3);
        Real height;
        Real rotation = 0.0;

        pp.getarr("channel_pt_on_top_wall", pt_on_top_wall, 0, 3);
        pp.get("channel_rotation", rotation);
        pp.get("channel_height", height);
        rotation = (rotation/180.) * M_PI;

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

           Real alpha = (rotation[0]/180.) * M_PI;
           Real gamma = (rotation[1]/180.) * M_PI;

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
