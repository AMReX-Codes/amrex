
#include <AMReX_EB2_IF_AllRegular.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Torus.H>
#include <AMReX_EB2_IF_Spline.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>
#include <AMReX.H>

namespace amrex { namespace EB2 {

Vector<std::unique_ptr<IndexSpace> > IndexSpace::m_instance;

int max_grid_size = 64;
bool compare_with_ch_eb = false;

void Initialize ()
{
    ParmParse pp("eb2");
    pp.query("max_grid_size", max_grid_size);
    pp.query("compare_with_ch_eb", compare_with_ch_eb);

    amrex::ExecOnFinalize(Finalize);
}

void Finalize ()
{
    IndexSpace::clear();
}

const IndexSpace* TopIndexSpaceIfPresent() noexcept {
    if (IndexSpace::size() > 0) {
        return &IndexSpace::top();
    }
    return nullptr;
}

void
Build (const Geometry& geom, int required_coarsening_level,
       int max_coarsening_level, int ngrow)
{
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "all_regular")
    {
        EB2::AllRegularIF rif;
        EB2::GeometryShop<EB2::AllRegularIF> gshop(rif);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow);
    }
    else if (geom_type == "box")
    {
        RealArray lo;
        pp.get("box_lo", lo);

        RealArray hi;
        pp.get("box_hi", hi);

        bool has_fluid_inside;
        pp.get("box_has_fluid_inside", has_fluid_inside);

        EB2::BoxIF bf(lo, hi, has_fluid_inside);

        EB2::GeometryShop<EB2::BoxIF> gshop(bf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow);
    }
    else if (geom_type == "cylinder")
    {
        RealArray center;
        pp.get("cylinder_center", center);

        Real radius;
        pp.get("cylinder_radius", radius);

        Real height = -1.0;
        pp.query("cylinder_height", height);

        int direction;
        pp.get("cylinder_direction", direction);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(direction >=0 && direction < AMREX_SPACEDIM,
                                         "eb2.cylinder_direction is invalid");

        bool has_fluid_inside;
        pp.get("cylinder_has_fluid_inside", has_fluid_inside);

        EB2::CylinderIF cf(radius, height, direction, center, has_fluid_inside);

        EB2::GeometryShop<EB2::CylinderIF> gshop(cf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow);
    }
    else if (geom_type == "plane")
    {
        RealArray point;
        pp.get("plane_point", point);

        RealArray normal;
        pp.get("plane_normal", normal);

        EB2::PlaneIF pf(point, normal);

        EB2::GeometryShop<EB2::PlaneIF> gshop(pf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow);
    }
    else if (geom_type == "sphere")
    {
        RealArray center;
        pp.get("sphere_center", center);

        Real radius;
        pp.get("sphere_radius", radius);

        bool has_fluid_inside;
        pp.get("sphere_has_fluid_inside", has_fluid_inside);

        EB2::SphereIF sf(radius, center, has_fluid_inside);

        EB2::GeometryShop<EB2::SphereIF> gshop(sf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow);
    }
    else if (geom_type == "torus")
    {
        RealArray center;
        pp.get("torus_center", center);

        Real small_radius;
        pp.get("torus_small_radius", small_radius);

        Real large_radius;
        pp.get("torus_large_radius", large_radius);

        bool has_fluid_inside = true;

        EB2::TorusIF sf(large_radius, small_radius, center, has_fluid_inside);

        EB2::GeometryShop<EB2::TorusIF> gshop(sf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow);
    }
    else
    {
        amrex::Abort("geom_type "+geom_type+ " not supported");
    }
}

namespace {
static int comp_max_crse_level (Box cdomain, const Box& domain)
{
    int ilev;
    for (ilev = 0; ilev < 30; ++ilev) {
        if (cdomain.contains(domain)) break;
        cdomain.refine(2);
    }
    if (cdomain != domain) ilev = -1;
    return ilev;
}
}

int
maxCoarseningLevel (const Geometry& geom)
{
    const Box& domain = amrex::enclosedCells(geom.Domain());
    const Box& cdomain = IndexSpace::top().coarsestDomain();
    return comp_max_crse_level(cdomain, domain);
}

int
maxCoarseningLevel (IndexSpace const* ebis, const Geometry& geom)
{
    const Box& domain = amrex::enclosedCells(geom.Domain());
    const Box& cdomain = ebis->coarsestDomain();
    return comp_max_crse_level(cdomain,domain);
}

}}
