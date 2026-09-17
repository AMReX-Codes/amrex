// Default gridding must retain the results and validation used before no_chop_dir.
#include <AMReX.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_ParmParse.H>
#include <AMReX_TagBox.H>
#ifdef TEST_AMRLEVEL
#include <AMReX_LevelBld.H>
#endif

#include <functional>
#include <random>
#include <string>

using namespace amrex;

namespace {
class Mesh : public AmrMesh {
public:
    using AmrMesh::AmrMesh;
    using AmrMesh::checkInput;
    using AmrMesh::bfLev;
    Vector<Box> tagged_boxes;
    bool scale_tags = false; // tagged_boxes are level 0 boxes
    Vector<Box> fine_tagged_boxes; // if not empty, used on levels >= 2
    void ErrorEst (int lev, TagBoxArray& tags, Real, int) override {
        if (tagged_boxes.empty()) {
            tags.setVal(tags.boxArray(), TagBox::SET);
        } else {
            Vector<Box> v((lev >= 2 && !fine_tagged_boxes.empty()) ? fine_tagged_boxes : tagged_boxes);
            for (int l = 0; scale_tags && l < lev; ++l) {
                for (auto& b : v) { b.refine(refRatio(l)); }
            }
            tags.setVal(BoxArray(BoxList(std::move(v))), TagBox::SET);
        }
    }
};

// Regrid from level 0 until all levels exist.
void regrid (Mesh& mesh) {
    auto ba = mesh.MakeBaseGrids();
    mesh.SetBoxArray(0, ba);
    mesh.SetDistributionMap(0, DistributionMapping(ba));
    mesh.SetFinestLevel(0);
    Vector<BoxArray> grids;
    for (int iter = 0; iter < mesh.maxLevel(); ++iter) {
        int finest = 0;
        mesh.MakeNewGrids(0, 0., finest, grids);
        for (int lev = 1; lev <= finest; ++lev) {
            mesh.SetBoxArray(lev, grids[lev]);
            mesh.SetDistributionMap(lev, DistributionMapping(grids[lev]));
        }
        mesh.SetFinestLevel(finest);
    }
}

Geometry geometry (IntVect const& size, bool periodic = false) {
    Array<int,AMREX_SPACEDIM> is_per{};
    is_per.fill(periodic ? 1 : 0);
    return Geometry(Box(IntVect(0), size-1),
                    RealBox(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.)), 0, is_per);
}

void rejects (std::function<void()> const& f, std::string const& message) {
    bool caught = false;
    try { f(); }
    catch (std::exception const& e) {
        caught = std::string(e.what()).find(message) != std::string::npos;
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(caught, "Expected rejection: " + message);
}

void test_legacy () {
    AmrInfo info;
    info.max_level = 1;
    info.max_grid_size = {IntVect(32)};
    info.refine_grid_layout = false;
    Mesh mesh(geometry(IntVect(128)), info);
    BoxArray ba(mesh.Geom(0).Domain());
    mesh.ChopGrids(0, ba, 8);
    IntVect chunk(32);
    chunk[AMREX_SPACEDIM-1] = 16;
    BoxArray expected(mesh.Geom(0).Domain());
    expected.maxSize(chunk);
    AMREX_ALWAYS_ASSERT(ba == expected);

    // Unset and any negative no_chop_dir use the same legacy path.
    info.no_chop_dir = -5;
    Mesh negative(geometry(IntVect(128)), info);
    ba = BoxArray(mesh.Geom(0).Domain());
    negative.ChopGrids(0, ba, 8);
    AMREX_ALWAYS_ASSERT(ba == expected);

    // Legacy octree inputs can use BF > MGS and a truncated periodic domain.
    info.max_grid_size = {IntVect(4)};
    Mesh octree(geometry(IntVect(10), true), info);
    AMREX_ALWAYS_ASSERT(octree.bfLev(0) == IntVect(4));
    info.refine_whole_domain_dir = 0;
    rejects([&] { Mesh invalid(geometry(IntVect(10)), info); },
            "Domain size not divisible by blocking_factor/ref_ratio");
    info.refine_whole_domain_dir = -1;

    info.max_grid_size = {IntVect(32)};
    info.blocking_factor = {IntVect(2), IntVect(16)};
    rejects([&] { Mesh invalid(geometry(IntVect(34)), info); },
            "Blocking factors vary too much");
    info.max_grid_size = {IntVect(2), IntVect(32)};
    info.blocking_factor = {IntVect(8)};
    rejects([&] { Mesh invalid(geometry(IntVect(32)), info); },
            "Coarse level blocking factor not a multiple");

    // Integer division here has always been valid, including with assertions on.
    info.max_level = 2;
    info.ref_ratio = {IntVect(6), IntVect(2)};
    info.max_grid_size = {IntVect(48)};
    info.n_error_buf = {IntVect(0)};
    Mesh rr62(geometry(IntVect(16)), info);
    auto base = rr62.MakeBaseGrids();
    rr62.SetBoxArray(0, base);
    rr62.SetDistributionMap(0, DistributionMapping(base));
    rr62.SetFinestLevel(0);
    Vector<BoxArray> grids;
    for (int iter = 0; iter < 3; ++iter) {
        int finest = 0;
        rr62.MakeNewGrids(0, 0., finest, grids);
        for (int lev = 1; lev <= finest; ++lev) {
            rr62.SetBoxArray(lev, grids[lev]);
            rr62.SetDistributionMap(lev, DistributionMapping(grids[lev]));
            expected = BoxArray(rr62.Geom(lev).Domain());
            expected.maxSize(48);
            AMREX_ALWAYS_ASSERT(grids[lev] == expected);
        }
        rr62.SetFinestLevel(finest);
    }
    AMREX_ALWAYS_ASSERT(rr62.finestLevel() == 2);
}

// Without no_chop_dir, odd ratios keep the original rules.
void test_mixed_ratios () {
    AmrInfo info;
    info.max_level = 2;
    info.ref_ratio = {IntVect(3), IntVect(2)};
    info.max_grid_size = {IntVect(32)};
    info.refine_grid_layout = false;
    Mesh mesh(geometry(IntVect(16)), info); // default blocking_factor 8
    regrid(mesh);
    AMREX_ALWAYS_ASSERT(mesh.finestLevel() == 2);
    AMREX_ALWAYS_ASSERT(mesh.boxArray(2).contains(BoxArray(mesh.Geom(2).Domain())));

    // ChopGrids honors max_grid_size in every direction, as before.
    info.max_level = 1;
    info.ref_ratio = {IntVect(3)};
    Mesh rr3(geometry(IntVect(128)), info);
    BoxArray ba(rr3.Geom(0).Domain());
    rr3.ChopGrids(0, ba, 8);
    IntVect chunk(32);
    chunk[AMREX_SPACEDIM-1] = 16;
    BoxArray expected(rr3.Geom(0).Domain());
    expected.maxSize(chunk);
    AMREX_ALWAYS_ASSERT(ba == expected);

#if AMREX_SPACEDIM >= 2
    // With no_chop_dir, ChopGrids still honors max_grid_size in the other
    // directions.
    info.no_chop_dir = AMREX_SPACEDIM-1;
    info.blocking_factor = {IntVect(1), IntVect(24)};
    info.max_grid_size = {IntVect(32), IntVect(96)};
    Mesh nochop(geometry(IntVect(128)), info);
    ba = BoxArray(nochop.Geom(0).Domain());
    nochop.ChopGrids(0, ba, 8);
    for (int i = 0; i < ba.size(); ++i) {
        for (int d = 0; d < AMREX_SPACEDIM-1; ++d) {
            AMREX_ALWAYS_ASSERT(ba[i].length(d) <= 32);
        }
    }
#endif
}

void test_odd_ratio_inputs () {
    AmrInfo info;
    info.no_chop_dir = AMREX_SPACEDIM-1;
    info.max_level = 1;
    info.ref_ratio = {IntVect(3)};
    info.blocking_factor = {IntVect(1), IntVect(16)};
    info.max_grid_size = {IntVect(32), IntVect(64)};
    rejects([&] { Mesh invalid(geometry(IntVect(48)), info); },
            "blocking_factor not power of 2");

    // bfLev(1)*ref_ratio is not divisible by bfLev(2).  Regridding from
    // level 1 must still keep level 3 properly nested in level 2 when the
    // fine tags move to the edge of the level 1 grids.
    info.max_level = 3;
    info.ref_ratio = {IntVect(3), IntVect(3), IntVect(3)};
    IntVect bf3(48);
    bf3[AMREX_SPACEDIM-1] = 24; // max_grid_size does not apply in no_chop_dir
    info.blocking_factor = {IntVect(1), IntVect(24), IntVect(24), bf3};
    info.max_grid_size = {IntVect(32), IntVect(96), IntVect(16), IntVect(96)};
    info.n_error_buf = {IntVect(1)};
    info.refine_grid_layout = false;
    Mesh mesh(geometry(IntVect(64)), info);
    mesh.scale_tags = true;
    mesh.tagged_boxes = {Box(IntVect(16), IntVect(27))};
    regrid(mesh);
    AMREX_ALWAYS_ASSERT(mesh.finestLevel() == 3);
    mesh.fine_tagged_boxes = {Box(IntVect(10), IntVect(13))};
    {
        Vector<BoxArray> grids;
        int finest = 0;
        mesh.MakeNewGrids(1, 0., finest, grids);
        for (int lev = 2; lev <= finest; ++lev) {
            mesh.SetBoxArray(lev, grids[lev]);
            mesh.SetDistributionMap(lev, DistributionMapping(grids[lev]));
        }
        mesh.SetFinestLevel(finest);
    }
    for (int lev = 2; lev <= mesh.finestLevel(); ++lev) {
        BoxArray const& cba = mesh.boxArray(lev-1);
        IntVect const np = mesh.bfLev(lev-1) * mesh.nProper();
        for (int i = 0; i < mesh.boxArray(lev).size(); ++i) {
            Box b = amrex::coarsen(mesh.boxArray(lev)[i], mesh.refRatio(lev-1));
            b.grow(np);
            b &= mesh.Geom(lev-1).Domain();
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(cba.contains(b, true), "Grids not properly nested");
        }
    }
}

void test_tag_shape () {
    auto geom = geometry(IntVect(10));
    BoxArray ba(geom.Domain());
    ba.maxSize(5);
    TagBoxArray tags(ba, DistributionMapping(ba), IntVect(0));
    tags.coarsen(IntVect(2));
    geom.coarsen(IntVect(2));
    tags.mapPeriodicRemoveDuplicates(geom);
    AMREX_ALWAYS_ASSERT(tags.nGrowVect() == IntVect(0));
}

#ifdef TEST_AMRLEVEL
void test_amr_validation () {
    struct Factory : LevelBld {
        void variableSetUp () override {}
        void variableCleanUp () override {}
        AmrLevel* operator() () override { return nullptr; }
        AmrLevel* operator() (Amr&, int, Geometry const&, BoxArray const&,
                             DistributionMapping const&, Real) override { return nullptr; }
    } factory;
    struct Driver : Amr {
        using Amr::Amr;
        using Amr::checkInput;
        void unrefined_direction () { ref_ratio[0][AMREX_SPACEDIM-1] = 1; }
        void odd_ratio () {
            no_chop_dir = 0;
            ref_ratio[0] = IntVect(3);
            ref_ratio[0][AMREX_SPACEDIM-1] = 1;
        }
    };
    RealBox rb(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.));
    Driver driver(&rb, 1, Vector<int>(AMREX_SPACEDIM, 32), 0, &factory);
    driver.SetBlockingFactor(6);
    rejects([&] { driver.checkInput(); }, "blocking_factor not power of 2");
#if AMREX_SPACEDIM >= 2
    driver.SetBlockingFactor(1);
    driver.unrefined_direction();
    IntVect odd_mgs(32);
    odd_mgs[AMREX_SPACEDIM-1] = 33;
    driver.SetMaxGridSize(Vector<IntVect>{IntVect(32), odd_mgs});
    rejects([&] { driver.checkInput(); }, "max_grid_size is not even");
    // With no_chop_dir, an odd ratio elsewhere does not exempt a direction
    // with ratio 1.
    driver.odd_ratio();
    rejects([&] { driver.checkInput(); }, "max_grid_size is not even");
#endif
}
#endif

void test_no_chop () {
    AmrInfo info;
    info.check_input = false;
    info.no_chop_dir = AMREX_SPACEDIM;
    rejects([&] { Mesh invalid(geometry(IntVect(128)), info); }, "no_chop_dir is out of range");
    ParmParse pp("amr");
    pp.add("check_input", false);
    pp.add("no_chop_dir", AMREX_SPACEDIM);
    pp.add("max_level", 0);
    pp.addarr("n_cell", Vector<int>(AMREX_SPACEDIM, 128));
    rejects([] { Mesh invalid; }, "no_chop_dir is out of range");

#if AMREX_SPACEDIM == 3
    info.check_input = true;
    info.no_chop_dir = 2;
    info.blocking_factor = {IntVect(1)};
    info.max_grid_size = {IntVect(128)};
    info.refine_grid_layout_dims = IntVect(1,0,0);
    Mesh columns(geometry(IntVect(128,128,64)), info);
    auto ba = columns.MakeBaseGrids();
    AMREX_ALWAYS_ASSERT(ba.size() >= ParallelDescriptor::NProcs());
    for (int i = 0; i < ba.size(); ++i) {
        AMREX_ALWAYS_ASSERT(ba[i].length(1) == 128 && ba[i].length(2) == 64);
    }
#endif

#if AMREX_SPACEDIM >= 2
    {
        // Level 0 honors max_grid_size in each direction.
        AmrInfo info0;
        info0.no_chop_dir = AMREX_SPACEDIM-1;
        info0.blocking_factor = {IntVect(1)};
        IntVect mgs(64);
        mgs[0] = 16;
        info0.max_grid_size = {mgs};
        info0.refine_grid_layout = false;
        Mesh mesh(geometry(IntVect(64)), info0);
        auto ba0 = mesh.MakeBaseGrids();
        AMREX_ALWAYS_ASSERT(ba0.size() == 4);
        for (int i = 0; i < ba0.size(); ++i) {
            AMREX_ALWAYS_ASSERT(ba0[i].length(0) == 16);
        }
    }
    {
        // An even ratio with a domain not divisible by the blocking factor
        // must not give grids thinner than the blocking factor.
        AmrInfo info1;
        info1.max_level = 1;
        info1.no_chop_dir = AMREX_SPACEDIM-1;
        info1.blocking_factor = {IntVect(1), IntVect(8)};
        info1.max_grid_size = {IntVect(32)};
        info1.n_error_buf = {IntVect(0)};
        info1.refine_grid_layout = false;
        IntVect n_cell(32);
        n_cell[0] = 30;
        n_cell[AMREX_SPACEDIM-1] = 16;
        Mesh mesh(geometry(n_cell), info1);
        IntVect const point(AMREX_D_DECL(29,10,5));
        mesh.tagged_boxes = {Box(point, point)};
        regrid(mesh);
        AMREX_ALWAYS_ASSERT(mesh.finestLevel() == 1);
        BoxArray const& ba1 = mesh.boxArray(1);
        AMREX_ALWAYS_ASSERT(ba1.contains(amrex::refine(Box(point, point), 2)));
        for (int i = 0; i < ba1.size(); ++i) {
            for (int d = 0; d < AMREX_SPACEDIM-1; ++d) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ba1[i].length(d) >= 8, "Thin grid");
            }
        }
    }
#endif
}

#if AMREX_SPACEDIM >= 2
void test_supplied_grids () {
    AmrInfo info;
    info.max_level = 1;
    info.no_chop_dir = AMREX_SPACEDIM-1;
    info.n_error_buf = {IntVect(0)};
    info.refine_grid_layout = false;
    Mesh mesh(geometry(IntVect(32)), info);
    BoxArray ba(mesh.Geom(0).Domain());
    IntVect chunk(7);
    chunk[info.no_chop_dir] = 32;
    ba.maxSize(chunk);
    IntVect const point(AMREX_D_DECL(7,5,5));
    mesh.tagged_boxes = {Box(point, point)};
    mesh.SetBoxArray(0, ba);
    mesh.SetDistributionMap(0, DistributionMapping(ba));
    mesh.SetFinestLevel(0);
    Vector<BoxArray> grids;
    int finest = 0;
    mesh.MakeNewGrids(0, 0., finest, grids);
    AMREX_ALWAYS_ASSERT(finest == 1);
    AMREX_ALWAYS_ASSERT(grids[1].contains(amrex::refine(Box(point, point), 2)));
}

void test_boundary_extension () {
    for (bool new_chop : {false, true}) {
        for (int seed : {1, 19}) {
            AmrInfo info;
            info.max_level = 1;
            info.ref_ratio = {IntVect(3)};
            info.no_chop_dir = AMREX_SPACEDIM-1;
            info.blocking_factor = {IntVect(1), IntVect(24)};
            info.max_grid_size = {IntVect(96)};
            info.n_error_buf = {IntVect(0)};
            info.grid_eff = 1.;
            info.refine_grid_layout = false;
            info.use_new_chop = new_chop;
            Mesh mesh(geometry(IntVect(AMREX_D_DECL(81,83,16))), info);
            // These patterns exercise intersecting extensions at two partial boundaries.
            std::mt19937 gen(seed);
            for (int y = 0; y < 11; ++y) {
                for (int x = 0; x < 11; ++x) {
                    if (gen()%5 < 3) {
                        mesh.tagged_boxes.emplace_back(IntVect(AMREX_D_DECL(8*x,8*y,0)),
                                                       IntVect(AMREX_D_DECL(8*x,8*y,15)));
                    }
                }
            }
            auto ba = mesh.MakeBaseGrids();
            mesh.SetBoxArray(0, ba);
            mesh.SetDistributionMap(0, DistributionMapping(ba));
            mesh.SetFinestLevel(0);
            Vector<BoxArray> grids;
            int finest = 0;
            mesh.MakeNewGrids(0, 0., finest, grids);
            AMREX_ALWAYS_ASSERT(finest == 1 && grids[1].isDisjoint());
            for (int i = 0; i < grids[1].size(); ++i) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(grids[1][i].length().allGE(IntVect(24)),
                                                "Boundary extension was trimmed to a thin grid");
            }
            for (auto const& b : mesh.tagged_boxes) {
                AMREX_ALWAYS_ASSERT(grids[1].contains(amrex::refine(b, 3)));
            }
        }
    }
}
#endif
}

int main (int argc, char* argv[]) {
    amrex::Initialize(argc, argv);
    {
        amrex::system::throw_exception = true;
        test_legacy();
        test_mixed_ratios();
        test_odd_ratio_inputs();
        test_tag_shape();
#if AMREX_SPACEDIM >= 2
        test_supplied_grids();
        test_boundary_extension();
#endif
#ifdef TEST_AMRLEVEL
        test_amr_validation();
#endif
        test_no_chop();
        amrex::Print() << "Default gridding compatibility tests passed\n";
    }
    amrex::Finalize();
}
