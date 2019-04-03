
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>

#include <AMReX_ParmParse.H>

#include <cmath>
#include <algorithm>

using namespace amrex;

void
initialize_EB2 (const Geometry& geom, const int required_coarsening_level,
                const int max_coarsening_level)
{
    BL_PROFILE("initializeEB2");

    ParmParse ppeb2("eb2");
    std::string geom_type;
    ppeb2.get("geom_type", geom_type);

    if (geom_type == "combustor")
    {
        ParmParse pp("combustor");
        
        Real fwl; 
        pp.get("far_wall_loc",fwl);

        EB2::PlaneIF farwall({AMREX_D_DECL(fwl,0.,0.)},
                             {AMREX_D_DECL(1. ,0.,0.)});
        
        Vector<Real> pl1pt, pl2pt, pl2nm, pl3pt; 
        pp.getarr("ramp_plane1_point", pl1pt);
        pp.getarr("ramp_plane2_point", pl2pt);
        pp.getarr("ramp_plane2_normal", pl2nm);
        pp.getarr("ramp_plane3_point", pl3pt);

        auto ramp = EB2::makeIntersection(EB2::PlaneIF({pl1pt[0], pl1pt[1], 0.},
                                                       {      0.,      -1., 0.}),
                                          EB2::PlaneIF({pl2pt[0], pl2pt[1], 0.},
                                                       {pl2nm[0], pl2nm[1], 0.}),
                                          EB2::PlaneIF({pl3pt[0], pl3pt[1], 0.},
                                                       {      1.,       0., 0.}));

        Vector<Real> pipelo, pipehi; 
        pp.getarr("pipe_lo", pipelo);
        pp.getarr("pipe_hi", pipehi);
        
        EB2::BoxIF pipe({pipelo[0], pipelo[1], -1.}, {pipehi[0], pipehi[1], 1.}, false);

        // where does plane 1 and plane 2 intersect?
        Real k2 = std::abs(pl2nm[0]/pl2nm[1]);
        Real secty = pl2pt[1] + k2*(pl3pt[0]-pl2pt[0]);
        // How much do we cut?
        Real dx = geom.CellSize(0);
        Real dycut = 4.*(1.+max_coarsening_level)*std::min(dx, k2*dx);
        EB2::BoxIF flat_corner({pl3pt[0], 0., -1.}, {1.e10, secty+dycut, 1.}, false);
        
        auto polys = EB2::makeUnion(farwall, ramp, pipe, flat_corner);

        Real lenx = Geometry::ProbLength(0);
        Real leny = Geometry::ProbLength(1);
        auto pr = EB2::translate(EB2::lathe(polys), {lenx*0.5, leny*0.5, 0.});
        
        auto gshop = EB2::makeShop(pr);
        EB2::Build(gshop, geom, max_coarsening_level, max_coarsening_level, 4);
    }
    else
    {
        EB2::Build(geom, max_coarsening_level, max_coarsening_level, 4);
    }
}


