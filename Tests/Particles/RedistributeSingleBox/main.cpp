// On a single process with a single level, box and tile, Redistribute() works in place.
// This test checks that this leaves the container exactly as the full algorithm does
// (same particles, positions and order, in the same tiles), for AoS and pure SoA
// particles, and that it is used when it should be and not otherwise.
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>

#include <string>

using namespace amrex;

namespace
{

template <class PCBase>
class TestParticleContainer
    : public PCBase
{
public:
    using PCBase::PCBase;

    bool redistributeInPlace (IntVect nGrow, bool remove_negative)
    {
        return this->RedistributeInPlace(0, this->finestLevel(), nGrow, remove_negative);
    }
};

enum struct Layout { OneBox, Boxes, Tiling, TwoLevels };
enum struct BoxPart { Domain, LowerHalf, UpperHalf };    //!< of the domain in x
enum struct Move { Inside, AcrossUpperBoundary, IntoGhostCells };
enum struct InPlace { Yes, No, CpuOnly };

struct Scenario
{
    std::string name;
    Layout layout = Layout::OneBox;
    bool periodic = true;
    BoxPart box = BoxPart::Domain;
    int invalidate_every = 3;   //!< 0: none
    Move move = Move::Inside;   //!< of one particle; all others stay inside the box
    int ngrow = 0;
    bool remove_negative = true;
    InPlace in_place = InPlace::Yes;
};

template <class PC>
void modify (PC& pc, Scenario const& s, Real box_lo, Real box_hi)
{
    auto const& geom = pc.Geom(0);
    auto const dx = geom.CellSizeArray();
    auto const rlo = geom.ProbLoArrayInParticleReal();
    auto const rhi = geom.ProbHiArrayInParticleReal();
    int const every = s.invalidate_every;
    Move const move = s.move;

    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
        for (auto& kv : pc.GetParticles(lev)) {
            auto& ptile = kv.second;
            auto ptd = ptile.getParticleTileData();
            ParallelFor(ptile.numParticles(), [=] AMREX_GPU_DEVICE (int i)
            {
                if (every > 0 && i % every == 0) {
                    ptd.id(i).make_invalid();
                    return;
                }
                // shift every particle a bit, keeping it inside the box
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    auto x = ptd.pos(d, i) + ParticleReal(0.25) * ParticleReal(dx[d]);
                    ptd.pos(d, i) = amrex::min(amrex::max(x, rlo[d]), rhi[d]);
                }
                ptd.pos(0, i) = amrex::min(amrex::max(ptd.pos(0, i), ParticleReal(box_lo)),
                                           ParticleReal(box_hi - Real(0.5) * dx[0]));
                if (i == 1 && move == Move::AcrossUpperBoundary) {
                    ptd.pos(0, i) = rhi[0] + ParticleReal(0.5) * ParticleReal(dx[0]);
                } else if (i == 1 && move == Move::IntoGhostCells) {
                    ptd.pos(0, i) = ParticleReal(box_hi + Real(0.5) * dx[0]);
                }
            });
        }
    }
    Gpu::streamSynchronize();
}

//! number of particles that differ between two tiles, comparing all compile-time data
template <class PTile>
Long numDifferent (PTile const& a, PTile const& b)
{
    if (a.numParticles() != b.numParticles()) { return a.numParticles() + 1; }

    auto const pa = a.getConstParticleTileData();
    auto const pb = b.getConstParticleTileData();

    ReduceOps<ReduceOpSum> reduce_op;
    ReduceData<Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(a.numParticles(), reduce_data,
    [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
    {
        auto const sa = pa.getSuperParticle(i);
        auto const sb = pb.getSuperParticle(i);
        bool same = (sa.id() == sb.id()) && (sa.cpu() == sb.cpu());
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { same = same && (sa.pos(d) == sb.pos(d)); }
        using SP = std::decay_t<decltype(sa)>;
        for (int j = 0; j < SP::NReal; ++j) { same = same && (sa.rdata(j) == sb.rdata(j)); }
        for (int j = 0; j < SP::NInt; ++j) { same = same && (sa.idata(j) == sb.idata(j)); }
        return same ? 0 : 1;
    });
    return amrex::get<0>(reduce_data.value(reduce_op));
}

template <class PCBase>
void testScenario (Scenario const& s, std::string const& layout)
{
    using PC = TestParticleContainer<PCBase>;

    int const n_cell = 32;
    int const n_half = n_cell / 2;
    Box const domain(IntVect(0), IntVect(n_cell - 1));
    RealBox real_box;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        real_box.setLo(d, Real(-1.0));
        real_box.setHi(d, Real(1.0));
    }
    Array<int, AMREX_SPACEDIM> is_per{AMREX_D_DECL(int(s.periodic), int(s.periodic), int(s.periodic))};
    Geometry const geom(domain, real_box, CoordSys::cartesian, is_per);

    Box box = domain;
    if (s.box == BoxPart::LowerHalf) { box.setBig(0, n_half - 1); }
    if (s.box == BoxPart::UpperHalf) { box.setSmall(0, n_half); }
    BoxArray ba(box);
    if (s.layout == Layout::Boxes) { ba.maxSize(n_half); }

    Vector<Geometry> geoms{geom};
    Vector<BoxArray> bas{ba};
    if (s.layout == Layout::TwoLevels) {
        geoms.emplace_back(amrex::refine(domain, 2), real_box, CoordSys::cartesian, is_per);
        bas.emplace_back(amrex::refine(Box(IntVect(n_cell / 4), IntVect(3 * n_cell / 4 - 1)), 2));
    }
    Vector<DistributionMapping> dms;
    for (auto const& b : bas) { dms.emplace_back(b); }
    Vector<int> const rr(geoms.size() - 1, 2);

    bool const do_tiling = ParticleContainerBase::do_tiling;
    IntVect const tile_size = ParticleContainerBase::tile_size;
    if (s.layout == Layout::Tiling) {
        ParticleContainerBase::do_tiling = true;
        ParticleContainerBase::tile_size = IntVect(8);
    }

    // initial particles, inside the box
    Real const box_lo = geom.ProbLo(0) + geom.CellSize(0) * Real(box.smallEnd(0));
    Real const box_hi = geom.ProbLo(0) + geom.CellSize(0) * Real(box.bigEnd(0) + 1);
    RealBox init_box = real_box;
    init_box.setLo(0, box_lo);
    init_box.setHi(0, box_hi);
    typename PC::ParticleInitData pdata = {};
    PC a(geoms, dms, bas, rr);
    a.InitRandom(4000, 451, pdata, false, init_box);
    modify(a, s, box_lo, box_hi);

    PC b(geoms, dms, bas, rr);
    b.copyParticles(a, true);
    PC c(geoms, dms, bas, rr);
    c.copyParticles(a, true);

    Long const n_invalid = a.TotalNumberOfParticles(false) - a.TotalNumberOfParticles(true);

    // the decision alone
#ifdef AMREX_USE_GPU
    bool const gpu = true;
#else
    bool const gpu = false;
#endif
    bool const expect = (s.in_place == InPlace::Yes || (s.in_place == InPlace::CpuOnly && !gpu))
                        && ParallelContext::NProcsSub() == 1;
    bool const in_place = c.redistributeInPlace(IntVect(s.ngrow), s.remove_negative);

    // in place vs. the full algorithm, in a deterministic order of moved particles
    a.setStableRedistribute(1);
    b.setStableRedistribute(1);
    bool const early_out = ParticleContainerBase::redistribute_early_out;
    ParticleContainerBase::redistribute_early_out = true;
    a.Redistribute(0, -1, s.ngrow, 0, s.remove_negative);
    ParticleContainerBase::redistribute_early_out = false;
    b.Redistribute(0, -1, s.ngrow, 0, s.remove_negative);
    ParticleContainerBase::redistribute_early_out = early_out;
    ParticleContainerBase::do_tiling = do_tiling;
    ParticleContainerBase::tile_size = tile_size;

    Long n_diff = 0;
    bool same_tiles = a.numLevels() == b.numLevels();
    for (int lev = 0; same_tiles && lev < a.numLevels(); ++lev) {
        same_tiles = same_tiles && a.GetParticles(lev).size() == b.GetParticles(lev).size();
        for (auto const& kv : a.GetParticles(lev)) {
            auto const it = b.GetParticles(lev).find(kv.first);
            if (it == b.GetParticles(lev).end()) { same_tiles = false; continue; }
            n_diff += numDifferent(kv.second, it->second);
        }
    }
    ParallelAllReduce::Sum(n_diff, ParallelContext::CommunicatorSub());
    ParallelAllReduce::And(same_tiles, ParallelContext::CommunicatorSub());

    Long const n_invalid_after = a.TotalNumberOfParticles(false) - a.TotalNumberOfParticles(true);

    amrex::Print() << layout << " " << s.name << ": in place " << in_place
                   << " (expected " << expect << "), particles "
                   << b.TotalNumberOfParticles(false) << ", invalid " << n_invalid
                   << " -> " << n_invalid_after
                   << ", differences to the full Redistribute " << n_diff << "\n";

    AMREX_ALWAYS_ASSERT(in_place == expect);
    AMREX_ALWAYS_ASSERT(same_tiles);
    AMREX_ALWAYS_ASSERT(n_diff == 0);
    AMREX_ALWAYS_ASSERT(n_invalid_after == (s.remove_negative ? 0 : n_invalid));
    AMREX_ALWAYS_ASSERT(a.OK(0, -1, s.ngrow));
}

template <class PCBase>
void testAll (std::string const& layout)
{
    Vector<Scenario> const scenarios{
        {.name="inside"},
        {.name="inside_keep_invalid", .remove_negative=false},
        {.name="all_invalid", .invalidate_every=1},
        {.name="across_periodic_boundary", .move=Move::AcrossUpperBoundary},
        {.name="across_periodic_boundary_ngrow", .move=Move::AcrossUpperBoundary, .ngrow=2},
        {.name="leaves_domain", .periodic=false, .move=Move::AcrossUpperBoundary},
        {.name="ghost_cells_ngrow", .periodic=false, .box=BoxPart::LowerHalf,
         .move=Move::IntoGhostCells, .ngrow=2},
        // GPU: the full algorithm keeps it unshifted outside the domain, CPU: removed
        {.name="kept_unshifted_ngrow", .box=BoxPart::UpperHalf,
         .move=Move::AcrossUpperBoundary, .ngrow=2, .in_place=InPlace::CpuOnly},
        {.name="boxes", .layout=Layout::Boxes, .in_place=InPlace::No},
        {.name="tiling", .layout=Layout::Tiling, .in_place=InPlace::No},
        {.name="two_levels", .layout=Layout::TwoLevels, .in_place=InPlace::No},
    };
    for (auto const& s : scenarios) {
        testScenario<PCBase>(s, layout);
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        testAll<ParticleContainer<1, 1, 2, 1>>("AoS");
        testAll<ParticleContainerPureSoA<AMREX_SPACEDIM + 2, 1>>("SoA");
        amrex::Print() << "PASSED\n";
    }
    amrex::Finalize();
}
