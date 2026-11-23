// A tile created with GetParticles(lev)[index] has all the components of its particle
// container, including the ones added at runtime. This test fills such tiles, without
// DefineAndReturnParticleTile, and checks that the runtime components survive a
// Redistribute that moves particles between boxes (and processes, with MPI). It also
// checks that a container still finds its components by name after it is moved, and
// that a container with a PolymorphicArenaAllocator can create tiles before its arena
// is set.
#include <AMReX.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>

#include <array>
#include <map>
#include <string>
#include <utility>

using namespace amrex;

namespace
{

Geometry makeGeom ()
{
    Box const domain(IntVect(0), IntVect(31));
    RealBox const real_box(std::array<Real, AMREX_SPACEDIM>{AMREX_D_DECL(0., 0., 0.)},
                           std::array<Real, AMREX_SPACEDIM>{AMREX_D_DECL(1., 1., 1.)});
    return Geometry(domain, real_box, CoordSys::cartesian, {AMREX_D_DECL(1, 1, 1)});
}

// the values the runtime components of the particle with this id should have
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE ParticleReal runtimeReal (Long id) { return ParticleReal(id % 4096); }
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int runtimeInt (Long id) { return int(id % 1000) + 1; }

//! Create a tile in pmap, with a key that no MFIter uses, through code that only knows
//! it is a std::map, check that it has all the components of pc, and remove it.
template <class PC, class K, class V, class C, class A>
void checkNewTileThroughMap (PC const& pc, std::map<K, V, C, A>& pmap, int grid)
{
    K const key{grid, -1};
    AMREX_ALWAYS_ASSERT(pmap.count(key) == 0);

    auto [it, inserted] = pmap.try_emplace(key);
    AMREX_ALWAYS_ASSERT(inserted);
    AMREX_ALWAYS_ASSERT(it->second.NumRealComps() == pc.NumRealComps());
    AMREX_ALWAYS_ASSERT(it->second.NumIntComps() == pc.NumIntComps());
    pmap.erase(it);
}

//! Put np particles, spread over the whole domain, into the first tile of this
//! process, created with GetParticles(0)[index].
template <class PC>
void fillFirstTile (PC& pc, int np)
{
    MFIter mfi = pc.MakeMFIter(0);
    if (!mfi.isValid()) { return; }

    auto& ptile = pc.GetParticles(0)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
    AMREX_ALWAYS_ASSERT(ptile.NumRealComps() == pc.NumRealComps());
    AMREX_ALWAYS_ASSERT(ptile.NumIntComps() == pc.NumIntComps());

    ptile.resize(np);
    AMREX_ALWAYS_ASSERT(ptile.GetStructOfArrays().GetRealData(pc.NumRealComps()-1).size() == std::size_t(np));

    // try_emplace, also through a generic std::map reference, sets up new tiles the same way
    checkNewTileThroughMap(pc, pc.GetParticles(0), mfi.index());

    auto const plo = pc.Geom(0).ProbLoArrayInParticleReal();
    auto const phi = pc.Geom(0).ProbHiArrayInParticleReal();
    int const nreal = pc.NumRealComps() - PC::NArrayReal;
    int const nint  = pc.NumIntComps() - PC::NArrayInt;
    int const nreal_ct = PC::NArrayReal;
    Long const id0 = 1 + Long(ParallelDescriptor::MyProc()) * np;
    int const proc = ParallelDescriptor::MyProc();
    auto ptd = ptile.getParticleTileData();
    ParallelForRNG(np, [=] AMREX_GPU_DEVICE (int i, RandomEngine const& engine) noexcept
    {
        Long const id = id0 + i;
        ptd.id(i) = id;
        ptd.cpu(i) = proc;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            ptd.pos(d, i) = plo[d] + (phi[d] - plo[d]) * ParticleReal(Random(engine)) * ParticleReal(0.999);
        }
        for (int j = AMREX_SPACEDIM; j < nreal_ct; ++j) { ptd.rdata(j)[i] = ParticleReal(j); }
        for (int j = 0; j < nreal; ++j) { ptd.m_runtime_rdata[j][i] = runtimeReal(id) + ParticleReal(j); }
        for (int j = 0; j < nint;  ++j) { ptd.m_runtime_idata[j][i] = runtimeInt(id) + j; }
    });
    Gpu::streamSynchronize();
}

//! Check that every runtime component of every particle has its value, and count them.
template <class PC>
Long checkRuntimeComps (PC& pc)
{
    int const nreal = pc.NumRealComps() - PC::NArrayReal;
    int const nint  = pc.NumIntComps() - PC::NArrayInt;

    ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Long, Long> reduce_data(reduce_op);
    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
        for (auto& kv : pc.GetParticles(lev)) {
            auto& ptile = kv.second;
            auto const np = ptile.numParticles();
            auto const& soa = ptile.GetStructOfArrays();
            AMREX_ALWAYS_ASSERT(soa.NumRealComps() == pc.NumRealComps());
            AMREX_ALWAYS_ASSERT(soa.NumIntComps() == pc.NumIntComps());
            for (int j = 0; j < soa.NumRealComps(); ++j) {
                AMREX_ALWAYS_ASSERT(soa.GetRealData(j).size() == std::size_t(np));
            }
            for (int j = 0; j < soa.NumIntComps(); ++j) {
                AMREX_ALWAYS_ASSERT(soa.GetIntData(j).size() == std::size_t(np));
            }
            auto const ptd = ptile.getConstParticleTileData();
            reduce_op.eval(np, reduce_data, [=] AMREX_GPU_DEVICE (int i) -> GpuTuple<Long, Long>
            {
                Long const id = ptd.id(i);
                Long bad = 0;
                for (int j = 0; j < nreal; ++j) {
                    if (ptd.m_runtime_rdata[j][i] != runtimeReal(id) + ParticleReal(j)) { ++bad; }
                }
                for (int j = 0; j < nint; ++j) {
                    if (ptd.m_runtime_idata[j][i] != runtimeInt(id) + j) { ++bad; }
                }
                return {1, bad};
            });
        }
    }
    auto const [count, bad] = reduce_data.value(reduce_op);
    Long result[2] = {count, bad};
    ParallelDescriptor::ReduceLongSum(result, 2);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(result[1] == 0, "runtime components have wrong values");
    return result[0];
}

template <class PC>
void run (std::string const& name, Arena* arena = nullptr)
{
    Geometry const geom = makeGeom();
    BoxArray ba(geom.Domain());
    ba.maxSize(8);
    DistributionMapping const dm(ba);

    int const np = 500;
    Long const total = Long(np) * ParallelDescriptor::NProcs();

    PC pc(geom, dm, ba);
    pc.AddRealComp("runtime_real_0");
    pc.AddRealComp("runtime_real_1");
    pc.AddIntComp("runtime_int_0");

    // a tile created before the arena is set gets its components now and its arena later
    fillFirstTile(pc, np);
    if (arena != nullptr) { pc.SetArena(arena); }

    pc.Redistribute();
    AMREX_ALWAYS_ASSERT(checkRuntimeComps(pc) == total);

    if constexpr (PC::has_polymorphic_allocator) {
        for (auto& kv : pc.GetParticles(0)) {
            AMREX_ALWAYS_ASSERT(kv.second.GetStructOfArrays().GetRealData(0).arena() == arena);
            AMREX_ALWAYS_ASSERT(kv.second.GetStructOfArrays().GetRealData("runtime_real_1").arena() == arena);
        }
    }

    // a moved container, and its tiles, still find the components by name
    PC moved = std::move(pc);
    AMREX_ALWAYS_ASSERT(moved.GetRealCompIndex("runtime_real_1") == PC::NArrayReal + 1);
    for (auto& kv : moved.GetParticles(0)) {
        auto& soa = kv.second.GetStructOfArrays();
        AMREX_ALWAYS_ASSERT(soa.GetRealData("runtime_real_1").size() == std::size_t(kv.second.numParticles()));
        AMREX_ALWAYS_ASSERT(soa.GetIntData("runtime_int_0").size() == std::size_t(kv.second.numParticles()));
    }
    moved.Redistribute();
    AMREX_ALWAYS_ASSERT(checkRuntimeComps(moved) == total);

    Print() << name << ": OK\n";
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        run<ParticleContainer<1, 0, 1, 1>>("AoS");
        run<ParticleContainerPureSoA<AMREX_SPACEDIM + 1, 1>>("pure SoA");
        run<ParticleContainerPureSoA<AMREX_SPACEDIM + 1, 1, PolymorphicArenaAllocator>>(
            "pure SoA, polymorphic arena", The_Arena());
        run<ParticleContainerPureSoA<AMREX_SPACEDIM + 1, 1, PolymorphicArenaAllocator>>(
            "pure SoA, polymorphic pinned arena", The_Pinned_Arena());
    }
    amrex::Finalize();
}
