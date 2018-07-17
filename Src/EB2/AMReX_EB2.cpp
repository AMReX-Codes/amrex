
#include <AMReX_EB2_IF_AllRegular.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Sphere.H>

#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>
#include <AMReX.H>

namespace amrex { namespace EB2 {

Vector<std::unique_ptr<IndexSpace> > IndexSpace::m_instance;

bool use_eb2 = false;
int max_grid_size = 64;
bool compare_with_ch_eb = false;

void Initialize ()
{
    ParmParse pp("eb2");
    pp.query("use_eb2", use_eb2);
    pp.query("max_grid_size", max_grid_size);
    pp.query("compare_with_ch_eb", compare_with_ch_eb);

    amrex::ExecOnFinalize(Finalize);
}

void Finalize ()
{
    IndexSpace::clear();
}

void
Build (const Geometry& geom, int required_coarsening_level, int max_coarsening_level)
{
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "all_regular")
    {
        EB2::AllRegularIF rif;
        EB2::GeometryShop<EB2::AllRegularIF> gshop(rif);
        EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
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
        EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
    }
    else if (geom_type == "cylinder")
    {
        RealArray center;
        pp.get("cylinder_center", center);

        Real radius;
        pp.get("cylinder_radius", radius);

        Real height;
        pp.get("cylinder_height", height);

        int direction;
        pp.get("cylinder_direction", direction);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(direction >=0 && direction <= AMREX_SPACEDIM,
                                         "eb2.cylinder_direction is invalid");

        bool has_fluid_inside;
        pp.get("cylinder_has_fluid_inside", has_fluid_inside);

        EB2::CylinderIF cf(radius, height, direction, center, has_fluid_inside);

        EB2::GeometryShop<EB2::CylinderIF> gshop(cf);
        EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
    }
    else if (geom_type == "plane")
    {
        RealArray point;
        pp.get("plane_point", point);

        RealArray normal;
        pp.get("plane_normal", normal);

        EB2::PlaneIF pf(point, normal);

        EB2::GeometryShop<EB2::PlaneIF> gshop(pf);
        EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
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
        EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
    }
    else
    {
        amrex::Abort("geom_type "+geom_type+ " not supported");
    }
}

const Level&
getLevel (const Geometry& geom)
{
    return IndexSpace::top().getLevel(geom);
}

int
maxCoarseningLevel ()
{
    return IndexSpace::top().maxCoarseningLevel();
}

}}
