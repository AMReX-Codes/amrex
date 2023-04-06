
#include <AMReX_EB2_IF_AllRegular.H>
#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Cylinder.H>
#include <AMReX_EB2_IF_Ellipsoid.H>
#include <AMReX_EB2_IF_Plane.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EB2_IF_Torus.H>
#include <AMReX_EB2_IF_Spline.H>
#include <AMReX_EB2_IF_Parser.H>
#include <AMReX_EB2_GeometryShop.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IndexSpace_STL.H>
#include <AMReX_EB2_IndexSpace_chkpt_file.H>
#include <AMReX_ParmParse.H>
#include <AMReX.H>
#include <algorithm>

namespace amrex::EB2 {

AMREX_EXPORT Vector<std::unique_ptr<IndexSpace> > IndexSpace::m_instance;

AMREX_EXPORT int max_grid_size = 64;
AMREX_EXPORT bool extend_domain_face = true;
AMREX_EXPORT int num_coarsen_opt = 0;

void Initialize ()
{
    ParmParse pp("eb2");
    pp.queryAdd("max_grid_size", max_grid_size);
    pp.queryAdd("extend_domain_face", extend_domain_face);
    pp.queryAdd("num_coarsen_opt", num_coarsen_opt);

    amrex::ExecOnFinalize(Finalize);
}

void Finalize ()
{
    IndexSpace::clear();
}

bool ExtendDomainFace ()
{
    return extend_domain_face;
}

int NumCoarsenOpt ()
{
    return num_coarsen_opt;
}

void
IndexSpace::push (IndexSpace* ispace)
{
    auto r = std::find_if(m_instance.begin(), m_instance.end(),
                          [=] (const std::unique_ptr<IndexSpace>& x) -> bool
                          { return x.get() == ispace; });
    if (r == m_instance.end()) {
        m_instance.emplace_back(ispace);
    } else if (r+1 != m_instance.end()) {
        std::rotate(r, r+1, m_instance.end());
    }
}

void
IndexSpace::erase (IndexSpace* ispace)
{
    auto r = std::find_if(m_instance.begin(), m_instance.end(),
                          [=] (const std::unique_ptr<IndexSpace>& x) -> bool
                          { return x.get() == ispace; });
    if (r != m_instance.end()) {
        m_instance.erase(r);
    }
}

const IndexSpace* TopIndexSpaceIfPresent() noexcept {
    if (IndexSpace::size() > 0) {
        return &IndexSpace::top();
    }
    return nullptr;
}

void
Build (const Geometry& geom, int required_coarsening_level,
       int max_coarsening_level, int ngrow, bool build_coarse_level_by_coarsening,
       bool a_extend_domain_face, int a_num_coarsen_opt)
{
    ParmParse pp("eb2");
    std::string geom_type;
    pp.get("geom_type", geom_type);

    if (geom_type == "all_regular")
    {
        EB2::AllRegularIF rif;
        EB2::GeometryShop<EB2::AllRegularIF> gshop(rif);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
                   a_extend_domain_face, a_num_coarsen_opt);
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
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
                   a_extend_domain_face, a_num_coarsen_opt);
    }
    else if (geom_type == "cylinder")
    {
        RealArray center;
        pp.get("cylinder_center", center);

        Real radius;
        pp.get("cylinder_radius", radius);

        Real height = -1.0;
        pp.queryAdd("cylinder_height", height);

        int direction;
        pp.get("cylinder_direction", direction);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(direction >=0 && direction < 3,
                                         "eb2.cylinder_direction is invalid");

        bool has_fluid_inside;
        pp.get("cylinder_has_fluid_inside", has_fluid_inside);

        EB2::CylinderIF cf(radius, height, direction, center, has_fluid_inside);

        EB2::GeometryShop<EB2::CylinderIF> gshop(cf);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
                   a_extend_domain_face, a_num_coarsen_opt);
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
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
                   a_extend_domain_face, a_num_coarsen_opt);
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
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
                   a_extend_domain_face, a_num_coarsen_opt);
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
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
                   a_extend_domain_face, a_num_coarsen_opt);
    }
    else if (geom_type == "parser")
    {
        std::string fn_string;
        pp.get("parser_function", fn_string);
        Parser parser(fn_string);
        parser.registerVariables({"x","y","z"});
        EB2::ParserIF pif(parser.compile<3>());
        EB2::GeometryShop<EB2::ParserIF,Parser> gshop(pif,parser);
        EB2::Build(gshop, geom, required_coarsening_level,
                   max_coarsening_level, ngrow, build_coarse_level_by_coarsening,
                   a_extend_domain_face, a_num_coarsen_opt);
    }
    else if (geom_type == "stl")
    {
        std::string stl_file;
        pp.get("stl_file", stl_file);
        Real stl_scale = 1._rt;
        pp.queryAdd("stl_scale", stl_scale);
        std::vector<Real> stl_center{0.0_rt, 0.0_rt, 0.0_rt};
        pp.queryAdd("stl_center", stl_center);
        int stl_reverse_normal = 0;
        pp.queryAdd("stl_reverse_normal", stl_reverse_normal);
        IndexSpace::push(new IndexSpaceSTL(stl_file, stl_scale,
                                           {stl_center[0], stl_center[1], stl_center[2]},
                                           stl_reverse_normal,
                                           geom, required_coarsening_level,
                                           max_coarsening_level, ngrow,
                                           build_coarse_level_by_coarsening,
                                           a_extend_domain_face,
                                           a_num_coarsen_opt));
    }
    else
    {
        amrex::Abort("geom_type "+geom_type+ " not supported");
    }
}

void addFineLevels (int num_new_fine_levels)
{
    BL_PROFILE("EB2::addFineLevels()");
    auto *p = const_cast<IndexSpace*>(TopIndexSpace());
    if (p) {
        p->addFineLevels(num_new_fine_levels);
    }
}

void
BuildFromChkptFile (std::string const& fname,
                    const Geometry& geom, int required_coarsening_level,
                    int max_coarsening_level, int ngrow, bool build_coarse_level_by_coarsening,
                    bool a_extend_domain_face)
{
    ChkptFile chkpt_file(fname);
    IndexSpace::push(new IndexSpaceChkptFile(chkpt_file,
                     geom, required_coarsening_level,
                     max_coarsening_level, ngrow,
                     build_coarse_level_by_coarsening,
                     a_extend_domain_face));
}

namespace {
int comp_max_crse_level (Box cdomain, const Box& domain)
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

}
