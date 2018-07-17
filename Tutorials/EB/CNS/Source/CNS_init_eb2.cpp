
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Union.H>
#include <AMReX_EB2_IF_Intersection.H>
#include <AMReX_EB2_IF_Complement.H>
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

using namespace amrex;

void
initialize_EB2 (const Geometry& geom, const int required_coarsening_level,
                const int max_coarsening_level)
{
    BL_PROFILE("initializeEB2");

    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);
    
    if (geom_type == "combustor")
    {
        EB2::PlaneIF farwall({AMREX_D_DECL(0.45,0.,0.)},
                             {AMREX_D_DECL(1.  ,0.,0.)});
        auto ramp = EB2::makeIntersection(EB2::PlaneIF({AMREX_D_DECL(0.25, 0.75, 0.)},
                                                       {AMREX_D_DECL(0.  , -1. , 0.)}),
                                          EB2::PlaneIF({AMREX_D_DECL(0.25, 0.75, 0.)},
                                                       {AMREX_D_DECL(.69 , -.19, 0.)}),
                                          EB2::PlaneIF({AMREX_D_DECL(0.06, 0.  , 0.)},
                                                       {AMREX_D_DECL(1.  , 0.  , 0.)}));
        EB2::BoxIF pipe({AMREX_D_DECL(0.06, -1.0, -100.0)},
                        {AMREX_D_DECL(0.08, 0.5, 100.0)}, false);
        EB2::BoxIF flat_corner({AMREX_D_DECL(0.05999, -1.0, -100.0)},
                               {AMREX_D_DECL(1.0, 0.26, 100.0)}, false);
        auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);
        auto pr = EB2::translate(EB2::lathe(polys), {AMREX_D_DECL(0.5,0.5,0.)});

        auto gshop = EB2::makeShop(pr);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level);
    }
    else
    {
        EB2::Build(geom, max_coarsening_level, max_coarsening_level);
    }
}
