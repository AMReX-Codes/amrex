
#include <AMReX_EB2_SphereIF.H>
#include <AMReX_EB2_CylinderIF.H>

#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>

namespace amrex { namespace EB2 {

Vector<std::unique_ptr<IndexSpace> > IndexSpace::m_instance;

void
Initialize (const Geometry& geom, const Info& info)
{
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "sphere")
    {
        std::vector<Real> vc;
        pp.getarr("sphere_center", vc);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vc.size() >= AMREX_SPACEDIM,
                                         "eb2.sphere_center doesn't have enough items");
        RealArray center{AMREX_D_DECL(vc[0],vc[1],vc[2])};

        Real radius;
        pp.get("sphere_radius", radius);

        bool has_fluid_inside;
        pp.get("sphere_has_fluid_inside", has_fluid_inside);

        EB2::SphereIF sf(radius, center, has_fluid_inside);

        EB2::GeometryShop<EB2::SphereIF> gshop(sf);
        EB2::Initialize(gshop, geom, info);
    }
    else if (geom_type == "cylinder")
    {
        std::vector<Real> vc;
        pp.getarr("cylinder_center", vc);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(vc.size() >= AMREX_SPACEDIM,
                                         "eb2.cylinder_center doesn't have enough items");
        RealArray center{AMREX_D_DECL(vc[0],vc[1],vc[2])};

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
        EB2::Initialize(gshop, geom, info);
    }
    else
    {
        amrex::Abort("geom_type "+geom_type+ " not supported");
    }
}

const Level&
getLevel (const Geometry& geom)
{
    IndexSpace::top().getLevel(geom);
}

}}
