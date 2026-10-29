/**
 * \file main.cpp
 * \brief Test for the cells EB2 stores outside the domain with extend_domain_face.
 *
 * With EB2::extrude_coarse_levels, the finest level stores ngrow cells outside each
 * non-periodic domain face instead of ngrow * 2^required_coarsening_level, and every coarse
 * level is extended back to ngrow after coarsening by extruding its outermost layer. This
 * test builds the same geometry both ways and checks that
 *
 *  (1) every level fills identical EB data on the domain grown by ngrow-2 (the level set
 *      by sign only, see compareLevels),
 *  (2) every required level stores at least ngrow cells outside the domain,
 *  (3) the finest level stores fewer cells than before,
 *  (4) a checkpoint written the old way still loads, and gives the same data,
 *
 * for an implicit function (a body cut by the domain floor, a box on the floor and a body
 * in a corner), an STL file (3D only), with and without a periodic direction.
 */

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EB2_IndexSpace_STL.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParReduce.H>
#include <AMReX_Print.H>

#if (AMREX_SPACEDIM == 3)
#  include <fstream>
#endif
#include <type_traits>
#include <string>

using namespace amrex;

namespace {

int n_cell = 64;
int required_coarsening_level = 3;
int ngrow = 4;

// The domain grown by ng in the non-periodic directions.
Box grownDomain (Geometry const& geom, int ng)
{
    Box bx = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (!geom.isPeriodic(idim)) { bx.grow(idim, ng); }
    }
    return bx;
}

/**
 * Count the entries that differ between a and b, separately inside inner (cell-centered;
 * converted to the data's index type) and outside it. The old build computes the outermost
 * cells of a coarse level from ghost data, so only the inner region must match exactly.
 */
template <typename MF>
GpuTuple<Long,Long> countDiffs (MF const& a, MF const& b, Box const& inner_cc,
                                bool sign_only = false)
{
    Box const inner = amrex::convert(inner_cc, a.ixType());
    auto const& aa = a.const_arrays();
    auto const& ba = b.const_arrays();
    auto r = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Long,Long>{},
                       a, IntVect(0), a.nComp(),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
            -> GpuTuple<Long,Long>
        {
            bool diff;
            if constexpr (std::is_same_v<MF,MultiFab>) {
                diff = sign_only ? (aa[box_no](i,j,k,n) < 0) != (ba[box_no](i,j,k,n) < 0)
                                 : aa[box_no](i,j,k,n) != ba[box_no](i,j,k,n);
            } else {
                diff = aa[box_no](i,j,k,n) != ba[box_no](i,j,k,n);
            }
            bool const in = inner.contains(i,j,k);
            return { Long(diff && in), Long(diff && !in) };
        });
    Long nin = amrex::get<0>(r);
    Long nout = amrex::get<1>(r);
    ParallelAllReduce::Sum(nin, ParallelContext::CommunicatorSub());
    ParallelAllReduce::Sum(nout, ParallelContext::CommunicatorSub());
    return {nin, nout};
}

int compareLevels (EB2::IndexSpace const& is_ref, EB2::IndexSpace const& is_new,
                   Geometry const& fine_geom, std::string const& label, bool expect_smaller)
{
    int nerrors = 0;
    Geometry geom = fine_geom;
    for (int ilev = 0; ilev <= required_coarsening_level; ++ilev) {
        EB2::Level const& lref = is_ref.getLevel(geom);
        EB2::Level const& lnew = is_new.getLevel(geom);

        BoxArray ba(grownDomain(geom, ngrow));
        ba.maxSize(32);
        DistributionMapping dm(ba);
        Box const inner = grownDomain(geom, ngrow-2);

        Long nin = 0, nout = 0;
        auto check = [&] (char const* name, auto const& a, auto const& b, bool sign_only = false) {
            auto [din, dout] = countDiffs(a, b, inner, sign_only);
            if (din > 0) {
                amrex::Print() << "  ERROR: " << label << " level " << ilev << ": " << name
                               << " differs at " << din << " points\n";
                ++nerrors;
            }
            nin += din;
            nout += dout;
        };

        {
            FabArray<EBCellFlagFab> a(ba, dm, 1, 0), b(ba, dm, 1, 0);
            lref.fillEBCellFlag(a, geom);
            lnew.fillEBCellFlag(b, geom);
            check("cellflag", a, b);
        }
        auto check_cc = [&] (char const* name, int ncomp, auto fill) {
            MultiFab a(ba, dm, ncomp, 0), b(ba, dm, ncomp, 0);
            fill(lref, a);
            fill(lnew, b);
            check(name, a, b);
        };
        check_cc("volfrac", 1, [&] (EB2::Level const& l, MultiFab& mf) { l.fillVolFrac(mf, geom); });
        check_cc("centroid", AMREX_SPACEDIM,
                 [&] (EB2::Level const& l, MultiFab& mf) { l.fillCentroid(mf, geom); });
        check_cc("bndryarea", 1, [&] (EB2::Level const& l, MultiFab& mf) { l.fillBndryArea(mf, geom); });
        check_cc("bndrycent", AMREX_SPACEDIM,
                 [&] (EB2::Level const& l, MultiFab& mf) { l.fillBndryCent(mf, geom); });
        check_cc("bndrynorm", AMREX_SPACEDIM,
                 [&] (EB2::Level const& l, MultiFab& mf) { l.fillBndryNorm(mf, geom); });
        {
            BoxArray const& nba = amrex::convert(ba, IntVect(1));
            MultiFab a(nba, dm, 1, 0), b(nba, dm, 1, 0);
            lref.fillLevelSet(a, geom);
            lnew.fillLevelSet(b, geom);
            // A node shared by a cut box and a covered box is stored twice, as the level set
            // and as the covered value 1, and ParallelCopy takes whichever comes last. That
            // depends on the distribution, which changes with the box sizes, so only the
            // sign is compared.
            check("levelset", a, b, true);
        }
        auto check_face = [&] (char const* name, int ncomp, bool edge, auto fill) {
            Array<MultiFab,AMREX_SPACEDIM> a, b;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                IntVect type = IntVect::TheDimensionVector(idim);
                if (edge) { type = IntVect(1); type[idim] = 0; }
                BoxArray const& fba = amrex::convert(ba, type);
                a[idim].define(fba, dm, ncomp, 0);
                b[idim].define(fba, dm, ncomp, 0);
            }
            fill(lref, amrex::GetArrOfPtrs(a));
            fill(lnew, amrex::GetArrOfPtrs(b));
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                check(name, a[idim], b[idim]);
            }
        };
        using PtrArr = Array<MultiFab*,AMREX_SPACEDIM>;
        check_face("areafrac", 1, false,
                   [&] (EB2::Level const& l, PtrArr const& p) { l.fillAreaFrac(p, geom); });
        check_face("facecent", AMREX_SPACEDIM-1, false,
                   [&] (EB2::Level const& l, PtrArr const& p) { l.fillFaceCent(p, geom); });
        check_face("edgecent", 1, true,
                   [&] (EB2::Level const& l, PtrArr const& p) { l.fillEdgeCent(p, geom); });

        IntVect const ng_new = lnew.nGrowVect();
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (!geom.isPeriodic(idim) && !lnew.isAllRegular() && ng_new[idim] < ngrow) {
                amrex::Print() << "  ERROR: " << label << " level " << ilev
                               << ": nGrowVect " << ng_new << " < " << ngrow << "\n";
                ++nerrors;
            }
        }

        amrex::Print() << "  " << label << " level " << ilev
                       << ": stored cells " << lref.boxArray().numPts()
                       << " -> " << lnew.boxArray().numPts()
                       << ", nGrowVect " << lref.nGrowVect() << " -> " << ng_new
                       << ", differences inside/outside grow(domain," << ngrow-2 << ") "
                       << nin << "/" << nout << "\n";

        geom = amrex::coarsen(geom, IntVect(2));
    }

    // The point of the change: the finest level no longer carries the coarse levels' cells.
    Long const ncells_ref = is_ref.getLevel(fine_geom).boxArray().numPts();
    Long const ncells_new = is_new.getLevel(fine_geom).boxArray().numPts();
    bool const has_wall = !fine_geom.isPeriodic(AMREX_SPACEDIM-1);
    if (expect_smaller && (ncells_new > ncells_ref || (has_wall && ncells_new == ncells_ref))) {
        amrex::Print() << "  ERROR: " << label << ": finest level stores " << ncells_new
                       << " cells, was " << ncells_ref << "\n";
        ++nerrors;
    }

    return nerrors;
}

#if (AMREX_SPACEDIM == 3)
// A pyramid standing on the floor of the unit cube, with its base below the floor, so that
// the floor cuts it along a square and its sides are not grid aligned.
void writePyramidSTL (std::string const& fname)
{
    if (!ParallelDescriptor::IOProcessor()) { return; }
    std::array<std::array<double,3>,5> const v{{{0.3,0.3,-0.2}, {0.7,0.3,-0.2}, {0.7,0.7,-0.2},
                                                 {0.3,0.7,-0.2}, {0.5,0.5,0.45}}};
    std::array<std::array<int,3>,6> const tri{{{0,2,1}, {0,3,2},  // base, facing down
                                               {0,1,4}, {1,2,4}, {2,3,4}, {3,0,4}}};
    std::ofstream os(fname);
    os << "solid pyramid\n";
    for (auto const& t : tri) {
        auto const& a = v[t[0]];
        auto const& b = v[t[1]];
        auto const& c = v[t[2]];
        std::array<double,3> const u{b[0]-a[0], b[1]-a[1], b[2]-a[2]};
        std::array<double,3> const w{c[0]-a[0], c[1]-a[1], c[2]-a[2]};
        os << "facet normal " << u[1]*w[2]-u[2]*w[1] << " " << u[2]*w[0]-u[0]*w[2] << " "
           << u[0]*w[1]-u[1]*w[0] << "\n"
           << "outer loop\n";
        for (int i = 0; i < 3; ++i) {
            auto const& p = v[t[i]];
            os << "vertex " << p[0] << " " << p[1] << " " << p[2] << "\n";
        }
        os << "endloop\n"
           << "endfacet\n";
    }
    os << "endsolid pyramid\n";
}
#endif

// Build twice, the old way and with extrusion, and compare. build() pushes one IndexSpace.
template <typename F>
int buildAndCompare (Geometry const& geom, std::string const& label, F const& build)
{
    bool const extrude = EB2::extrude_coarse_levels;
    amrex::Print() << "  " << label << ": building with ngrow * 2^required_coarsening_level\n";
    EB2::extrude_coarse_levels = false;
    build();
    EB2::IndexSpace const* is_ref = &EB2::IndexSpace::top();
    amrex::Print() << "  " << label << ": building with extruded coarse levels\n";
    EB2::extrude_coarse_levels = true;
    build();
    EB2::IndexSpace const* is_new = &EB2::IndexSpace::top();
    EB2::extrude_coarse_levels = extrude;

    int nerrors = compareLevels(*is_ref, *is_new, geom, label, true);

    // A checkpoint written the old way has more cells outside the domain than the new build
    // asks for. It must still load, and give the same data.
    if (label == "implicit") {
        std::string const fname = "eb_chkpt_" + std::to_string(geom.isPeriodic(0));
        is_ref->getLevel(geom).write_to_chkpt_file(fname, true, EB2::max_grid_size);
        EB2::BuildFromChkptFile(fname, geom, required_coarsening_level,
                                required_coarsening_level, ngrow, true, true);
        nerrors += compareLevels(*is_ref, EB2::IndexSpace::top(), geom, "old checkpoint", false);
        EB2::IndexSpace::pop();
    }

    EB2::IndexSpace::pop();
    EB2::IndexSpace::pop();
    return nerrors;
}

int run (bool periodic_x)
{
    int nerrors = 0;

    RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(int(periodic_x), 0, 0)};
    Geometry geom(Box(IntVect(0), IntVect(n_cell-1)), rb, 0, is_periodic);

    constexpr int floor_dir = AMREX_SPACEDIM-1;

    // Cut by the floor along a circle.
    RealArray center{AMREX_D_DECL(0.32, 0.4, 0.0)};
    center[floor_dir] = -0.1;
    EB2::SphereIF dome(0.28, center, false);

    // A box on the floor whose sides are not grid aligned.
    RealArray lo{AMREX_D_DECL(0.73, 0.71, 0.0)};
    RealArray hi{AMREX_D_DECL(0.91, 0.88, 0.0)};
    lo[floor_dir] = -1.0;
    hi[floor_dir] = 0.2;
    EB2::BoxIF block(lo, hi, false);

    // In the top corner, cut by two domain faces at once. With a periodic x direction it is
    // moved away from the x faces, because it would not be periodic.
    RealArray corner{AMREX_D_DECL(periodic_x ? 0.6 : 1.0, 0.65, 1.0)};
    corner[floor_dir] = 1.0;
    EB2::SphereIF ball(0.23, corner, false);

    auto gshop = EB2::makeShop(EB2::makeUnion(dome, block, ball));

    nerrors += buildAndCompare(geom, "implicit", [&] () {
        EB2::Build(gshop, geom, required_coarsening_level, required_coarsening_level, ngrow);
    });

#if (AMREX_SPACEDIM == 3)
    std::string const stl_file = "pyramid.stl";
    writePyramidSTL(stl_file);
    ParallelDescriptor::Barrier();
    nerrors += buildAndCompare(geom, "stl", [&] () {
        EB2::IndexSpace::push(std::make_unique<EB2::IndexSpaceSTL>
                              (stl_file, 1.0_rt, Array<Real,3>{0.0_rt, 0.0_rt, 0.0_rt}, 0,
                               geom, required_coarsening_level, required_coarsening_level,
                               ngrow, true, true, 0, true, false));
    });
#endif

    return nerrors;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;
        pp.query("n_cell", n_cell);
        pp.query("required_coarsening_level", required_coarsening_level);
        pp.query("ngrow", ngrow);

        int nerrors = 0;
        for (int per = 0; per < 2; ++per) {
            amrex::Print() << "periodic in x = " << per << "\n";
            nerrors += run(per == 1);
        }

        if (nerrors > 0) {
            amrex::Print() << "DomainFaceExtrusion test FAILED with " << nerrors << " errors\n";
            amrex::Abort();
        } else {
            amrex::Print() << "DomainFaceExtrusion test PASSED\n";
        }
    }
    amrex::Finalize();
    return 0;
}
