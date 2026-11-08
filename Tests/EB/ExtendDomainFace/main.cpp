// Test of extend_domain_face for EB2::GeometryShop and STL geometries.
//
// A sphere crosses the lower y face of the domain. For each geometry source
// and for extend_domain_face = true and false, the test checks on every
// required level:
//
//  - the number of cells stored outside the domain (nGrowVect);
//  - with extend_domain_face = true, that the cells outside the lower y face
//    are an extrusion of the domain face, i.e., every cell in a column below
//    the face has the same volume fraction and cell flag;
//  - with extend_domain_face = false, that some columns are not uniform,
//    i.e., the actual surface is used outside the domain;
//  - that the volume fractions inside the domain do not depend on
//    extend_domain_face.

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF_Sphere.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParReduce.H>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

using namespace amrex;

namespace {

constexpr Real radius = Real(0.3);
constexpr Real cx = Real(0.5);
constexpr Real cy = Real(0.1);
#if (AMREX_SPACEDIM == 3)
constexpr Real cz = Real(0.5);
#else
constexpr Real cz = Real(0.0); // 2D STL geometries are sliced at z = 0
#endif

// Write a closed, outward-oriented triangulated sphere in ASCII STL format.
void write_sphere_stl (std::string const& fname, int nphi, int ntheta)
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::ofstream ofs(fname);
    ofs << std::setprecision(17) << "solid sphere\n";
    auto point = [&] (int it, int ip) {
        double const pi = 3.14159265358979323846;
        double const th = pi * it / ntheta;
        double const ph = 2.0 * pi * ip / nphi;
        return std::array<double,3>{cx + radius*std::sin(th)*std::cos(ph),
                                    cy + radius*std::sin(th)*std::sin(ph),
                                    cz + radius*std::cos(th)};
    };
    auto facet = [&] (std::array<double,3> const& a, std::array<double,3> const& b,
                      std::array<double,3> const& c)
    {
        ofs << "facet normal 0 0 0\n outer loop\n";
        for (auto const* p : {&a, &b, &c}) {
            ofs << "  vertex " << (*p)[0] << " " << (*p)[1] << " " << (*p)[2] << "\n";
        }
        ofs << " endloop\nendfacet\n";
    };
    for (int it = 0; it < ntheta; ++it) {
        for (int ip = 0; ip < nphi; ++ip) {
            auto p00 = point(it  , ip  );
            auto p01 = point(it  , ip+1);
            auto p10 = point(it+1, ip  );
            auto p11 = point(it+1, ip+1);
            if (it > 0)          { facet(p00, p10, p01); }
            if (it < ntheta - 1) { facet(p01, p10, p11); }
        }
    }
    ofs << "endsolid sphere\n";
}

struct LevelData
{
    IntVect ngrow;
    MultiFab vfrac;
    FabArray<EBCellFlagFab> flag;
};

Vector<LevelData> get_levels (Geometry const& fine_geom, int nlevels, int max_grid_size)
{
    Vector<LevelData> r(nlevels);
    Geometry geom = fine_geom;
    for (int ilev = 0; ilev < nlevels; ++ilev) {
        auto const& level = EB2::IndexSpace::top().getLevel(geom);
        BoxArray ba(geom.Domain());
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        r[ilev].ngrow = level.nGrowVect();
        r[ilev].vfrac.define(ba, dm, 1, r[ilev].ngrow);
        r[ilev].flag.define(ba, dm, 1, r[ilev].ngrow);
        level.fillVolFrac(r[ilev].vfrac, geom);
        level.fillEBCellFlag(r[ilev].flag, geom);
        geom = amrex::coarsen(geom, 2);
    }
    return r;
}

// Count, among the cells outside the lower y face, those that differ from
// the cell just outside the face in the same column, and those that are
// covered.
std::pair<Long,Long> check_columns (LevelData const& ld, Box const& domain)
{
    auto const& va = ld.vfrac.const_arrays();
    auto const& fa = ld.flag.const_arrays();
    auto const dlo = amrex::lbound(domain);
    auto const dhi = amrex::ubound(domain);
    auto const rr = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Long,Long>{},
                              ld.vfrac, ld.ngrow,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept -> GpuTuple<Long,Long>
    {
        amrex::ignore_unused(k);
        bool const inside_others = i >= dlo.x && i <= dhi.x
#if (AMREX_SPACEDIM == 3)
            && k >= dlo.z && k <= dhi.z
#endif
            ;
        if (!inside_others || j >= dlo.y) { return {0, 0}; }
        int const jface = dlo.y - 1;
        auto type = [] (EBCellFlag f) {
            return f.isRegular() ? 0 : (f.isCovered() ? 1 : 2);
        };
        bool const differs = va[b](i,j,k) != va[b](i,jface,k)
            || type(fa[b](i,j,k)) != type(fa[b](i,jface,k));
        return {Long(differs), Long(fa[b](i,j,k).isCovered())};
    });
    Long ndiffer = amrex::get<0>(rr);
    Long ncovered = amrex::get<1>(rr);
    ParallelDescriptor::ReduceLongSum(ndiffer);
    ParallelDescriptor::ReduceLongSum(ncovered);
    return {ndiffer, ncovered};
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int ngrow = 8;
        int required_coarsening_level = 2;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("ngrow", ngrow);
            pp.query("required_coarsening_level", required_coarsening_level);
        }
        int const nlevels = required_coarsening_level + 1;

        Geometry geom(Box(IntVect(0), IntVect(n_cell-1)),
                      RealBox(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.)),
                      CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});

        std::string const stl_file("extend_domain_face_sphere.stl");
        write_sphere_stl(stl_file, 128, 64);
        ParallelDescriptor::Barrier();

        int nerrors = 0;
        auto fail = [&] (std::string const& msg) {
            amrex::Print() << "FAIL: " << msg << "\n";
            ++nerrors;
        };

        // Cells stored outside the domain on each level. With
        // extend_domain_face, every level stores ngrow (rounded up to 16)
        // cells. Otherwise the finest level stores
        // ngrow*2^required_coarsening_level cells and each coarse level holds
        // half of the finer one. No level stores more cells than the length
        // of its domain.
        auto round16 = [] (int n) { return ((n + 15) / 16) * 16; };

        for (std::string const geom_type : {"gshop", "stl"})
        {
            Vector<Vector<LevelData>> levels(2);
            for (int edf = 0; edf < 2; ++edf) {
                bool const extend_domain_face = (edf == 1);
                if (geom_type == "gshop") {
                    EB2::SphereIF sphere(radius, {AMREX_D_DECL(cx,cy,cz)}, false);
                    auto gshop = EB2::makeShop(sphere);
                    EB2::Build(gshop, geom, required_coarsening_level,
                               required_coarsening_level, ngrow, true,
                               extend_domain_face);
                } else {
                    ParmParse pp("eb2");
                    pp.add("geom_type", std::string("stl"));
                    pp.add("stl_file", stl_file);
                    EB2::Build(geom, required_coarsening_level, required_coarsening_level,
                               ngrow, true, extend_domain_face);
                }
                levels[edf] = get_levels(geom, nlevels, max_grid_size);
                EB2::IndexSpace::pop();
            }

            Box domain = geom.Domain();
            for (int ilev = 0; ilev < nlevels; ++ilev) {
                std::string const where = geom_type + " level " + std::to_string(ilev);

                int const ng_true = std::min(round16(ngrow), n_cell >> ilev);
                int const ng_false = std::min(round16(ngrow << required_coarsening_level),
                                              n_cell) >> ilev;
                if (levels[1][ilev].ngrow != IntVect(ng_true)) {
                    fail(where + ": extend_domain_face=true stores "
                         + std::to_string(levels[1][ilev].ngrow[1]) + " cells outside, expected "
                         + std::to_string(ng_true));
                }
                if (levels[0][ilev].ngrow != IntVect(ng_false)) {
                    fail(where + ": extend_domain_face=false stores "
                         + std::to_string(levels[0][ilev].ngrow[1]) + " cells outside, expected "
                         + std::to_string(ng_false));
                }

                auto [ndiff_t, ncov_t] = check_columns(levels[1][ilev], domain);
                amrex::Print() << where << ": extend_domain_face=true:  " << ndiff_t
                               << " cells differ from the face column, " << ncov_t
                               << " covered\n";
                if (ndiff_t != 0) {
                    fail(where + ": cells outside the domain are not an extrusion of the face");
                }
                if (ncov_t == 0) {
                    fail(where + ": no covered cells outside the domain");
                }

                auto [ndiff_f, ncov_f] = check_columns(levels[0][ilev], domain);
                amrex::Print() << where << ": extend_domain_face=false: " << ndiff_f
                               << " cells differ from the face column, " << ncov_f
                               << " covered\n";
                if (ndiff_f == 0) {
                    fail(where + ": extend_domain_face=false gives an extrusion");
                }

                MultiFab diff(levels[0][ilev].vfrac.boxArray(),
                              levels[0][ilev].vfrac.DistributionMap(), 1, 0);
                MultiFab::Copy(diff, levels[1][ilev].vfrac, 0, 0, 1, 0);
                MultiFab::Subtract(diff, levels[0][ilev].vfrac, 0, 0, 1, 0);
                Real const dmax = diff.norm0();
                amrex::Print() << where << ": max volume fraction difference inside the domain: "
                               << dmax << "\n";
                if (dmax != Real(0)) {
                    fail(where + ": volume fractions inside the domain depend on extend_domain_face");
                }

                domain.coarsen(2);
            }
        }

        if (nerrors > 0) {
            amrex::Abort("ExtendDomainFace test failed with " + std::to_string(nerrors)
                         + " errors");
        }
        amrex::Print() << "ExtendDomainFace test passed\n";
    }
    amrex::Finalize();
}
