
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

void
Build (const Geometry& geom, int max_coarsening_level)
{
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "box")
    {
        std::vector<Real> vlo;
        pp.getarr("box_lo", vlo);
        RealArray lo{AMREX_D_DECL(vlo[0],vlo[1],vlo[2])};

        std::vector<Real> vhi;
        pp.getarr("box_hi", vhi);
        RealArray hi{AMREX_D_DECL(vhi[0],vhi[1],vhi[2])};

        bool has_fluid_inside;
        pp.get("box_has_fluid_inside", has_fluid_inside);
        
        EB2::BoxIF bf(lo, hi, has_fluid_inside);

        EB2::GeometryShop<EB2::BoxIF> gshop(bf);
        EB2::Build(gshop, geom, max_coarsening_level);
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
        EB2::Build(gshop, geom, max_coarsening_level);
    }
    else if (geom_type == "plane")
    {
        std::vector<Real> vpoint;
        pp.getarr("plane_point", vpoint);
        RealArray point{AMREX_D_DECL(vpoint[0],vpoint[1],vpoint[2])};

        std::vector<Real> vnormal;
        pp.getarr("plane_normal", vnormal);
        RealArray normal{AMREX_D_DECL(vnormal[0],vnormal[1],vnormal[2])};

        EB2::PlaneIF pf(point, normal);

        EB2::GeometryShop<EB2::PlaneIF> gshop(pf);
        EB2::Build(gshop, geom, max_coarsening_level);
    }
    else if (geom_type == "sphere")
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
        EB2::Build(gshop, geom, max_coarsening_level);
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

}}
