#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_BoxList.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

using namespace amrex;

namespace {

// bl must cover the same cells as orig, be disjoint, and have no two
// boxes touching across a face normal to dir.
void check (BoxList const& orig, BoxList const& bl, int dir)
{
    AMREX_ALWAYS_ASSERT(bl.isDisjoint());
    AMREX_ALWAYS_ASSERT(BoxArray(bl).numPts() == BoxArray(orig).numPts());
    AMREX_ALWAYS_ASSERT(bl.contains(orig));
    auto const& v = bl.data();
    for (std::size_t i = 0; i < v.size(); ++i) {
        Box const s = amrex::shift(v[i], dir, 1);
        for (std::size_t j = 0; j < v.size(); ++j) {
            AMREX_ALWAYS_ASSERT(i == j || !s.intersects(v[j]));
        }
    }
}

void test_structured ()
{
    Box const domain(IntVect(0), IntVect(63));
    BoxArray ba(domain);
    ba.maxSize(16);
    BoxList const orig(ba);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        BoxList bl = orig;
        bl.mergeAlongDir(dir);
        check(orig, bl, dir);
        // The pencils are merged back into the domain by the final simplify.
        AMREX_ALWAYS_ASSERT(bl.size() == 1 && bl.front() == domain);
    }
}

#if (AMREX_SPACEDIM >= 2)
// Cuts are staggered across the x-interface, but the union is a single box.
void test_staggered ()
{
    Vector<Box> boxes;
    boxes.emplace_back(IntVect(AMREX_D_DECL( 0, 0,0)), IntVect(AMREX_D_DECL(31,31,7)));
    boxes.emplace_back(IntVect(AMREX_D_DECL(32, 0,0)), IntVect(AMREX_D_DECL(63,15,7)));
    boxes.emplace_back(IntVect(AMREX_D_DECL(32,16,0)), IntVect(AMREX_D_DECL(63,47,7)));
    boxes.emplace_back(IntVect(AMREX_D_DECL( 0,32,0)), IntVect(AMREX_D_DECL(31,47,7)));
    BoxList const orig(std::move(boxes));
    BoxList bl = orig;
    bl.mergeAlongDir(0);
    check(orig, bl, 0);
    AMREX_ALWAYS_ASSERT(bl.size() == 1);
}
#endif

void test_random (int nsplits, ULong seed)
{
    ResetRandomSeed(seed);
    Box const domain(IntVect(0), IntVect(127));
    Vector<Box> boxes{domain};
    for (int isplit = 0; isplit < nsplits; ++isplit) {
        auto const i = Random_int(static_cast<unsigned int>(boxes.size()));
        int const d = static_cast<int>(Random_int(AMREX_SPACEDIM));
        Box& b = boxes[i];
        if (b.length(d) < 2) { continue; }
        int const cut = b.smallEnd(d) + 1
            + static_cast<int>(Random_int(static_cast<unsigned int>(b.length(d)-1)));
        Box const hi = b.chop(d, cut);
        boxes.push_back(hi);
    }
    // Remove some boxes to create holes.
    Vector<Box> kept;
    for (auto const& b : boxes) {
        if (Random_int(5) != 0) { kept.push_back(b); }
    }
    if (kept.empty()) { kept.push_back(boxes[0]); }

    BoxList const orig(std::move(kept));
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        BoxList bl = orig;
        bl.mergeAlongDir(dir);
        check(orig, bl, dir);
        amrex::Print() << "  seed " << seed << " dir " << dir << ": "
                       << orig.size() << " -> " << bl.size() << " boxes\n";
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        test_structured();
#if (AMREX_SPACEDIM >= 2)
        test_staggered();
#endif
        for (ULong seed = 1; seed <= 4; ++seed) {
            test_random(300, seed);
        }
        amrex::Print() << "BoxList::mergeAlongDir tests passed\n";
    }
    amrex::Finalize();
}
