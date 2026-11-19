// Regression test for particle tiles that exist but hold no real particles.
//
//  1. AddRealComp / AddIntComp must update the runtime component count of
//     every tile, including empty ones (e.g. those left by Redistribute).
//  2. fillNeighbors must define the tiles that receive ghosts but have no
//     real particles, so the ghosts carry the runtime components.
//  3. sumNeighbors must sum back ghosts stored on a tile without real particles.

#include <AMReX.H>
#include <AMReX_MFIter.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace {

Geometry make_geom ()
{
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }
    const Box domain(IntVect(0), IntVect(15));
    Array<int,AMREX_SPACEDIM> is_per{AMREX_D_DECL(0,0,0)};
    return Geometry(domain, real_box, CoordSys::cartesian, is_per);
}

template <typename PC>
void check_runtime_comps (PC& pc, const char* name)
{
    for (int lev = 0; lev < pc.numLevels(); ++lev) {
        for (auto& kv : pc.GetParticles(lev)) {
            auto& ptile = kv.second;
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                ptile.NumRuntimeRealComps() == pc.NumRuntimeRealComps() &&
                ptile.NumRuntimeIntComps()  == pc.NumRuntimeIntComps(),
                std::string(name) + ": tile has a stale number of runtime components");
            auto ptd = ptile.getParticleTileData();
            AMREX_ALWAYS_ASSERT(ptd.m_num_runtime_real == pc.NumRuntimeRealComps());
            AMREX_ALWAYS_ASSERT(ptd.m_num_runtime_int  == pc.NumRuntimeIntComps());
        }
    }
}

template <typename PC>
void test_add_comps_on_empty_tiles (const char* name)
{
    Geometry geom = make_geom();
    BoxArray ba(geom.Domain());
    ba.maxSize(8);
    DistributionMapping dm(ba);

    PC pc(geom, dm, ba);

    // Redistribute defines a (here: empty) tile for every local box.
    pc.Redistribute();
    AMREX_ALWAYS_ASSERT(!pc.GetParticles(0).empty() || ParallelDescriptor::NProcs() > ba.size());

    pc.AddRealComp("w");
    pc.AddIntComp("tag");
    check_runtime_comps(pc, name);

    pc.AddRealComp("w2");
    check_runtime_comps(pc, name);

    amrex::Print() << name << ": AddRealComp/AddIntComp on empty tiles passed.\n";
}

#ifndef AMREX_USE_GPU
// The fixes tested here are in the CPU neighbor code, and sumNeighbors is CPU only.
void test_neighbors_on_empty_tiles ()
{
    using PC = NeighborParticleContainer<1, 0>;
    using ParticleType = PC::ParticleType;

    Geometry geom = make_geom();
    BoxArray ba(geom.Domain());
    ba.maxSize(IntVect(AMREX_D_DECL(8,16,16))); // two grids split in x
    AMREX_ALWAYS_ASSERT(ba.size() == 2);
    DistributionMapping dm(ba);

    PC pc(geom, dm, ba, 1);
    pc.setEnableInverse(true);
    pc.AddRealComp(true);
    const int iw = 0; // SoA index of the runtime real component

    // One particle in grid 0, in the last cell before grid 1, so that it is a
    // neighbor of grid 1. Grid 1 never gets real particles and Redistribute is
    // not called, so no tile is defined for grid 1.
    const auto dx = geom.CellSizeArray();
    if (dm[0] == ParallelDescriptor::MyProc()) {
        auto& ptile = pc.DefineAndReturnParticleTile(0, 0, 0);
        ParticleType p;
        p.id() = ParticleType::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
        p.pos(0) = static_cast<ParticleReal>(7.5*dx[0]);
#if (AMREX_SPACEDIM > 1)
        p.pos(1) = static_cast<ParticleReal>(8.5*dx[1]);
#endif
#if (AMREX_SPACEDIM > 2)
        p.pos(2) = static_cast<ParticleReal>(8.5*dx[2]);
#endif
        p.rdata(0) = 0.0_prt;
        ptile.push_back(p);
        ptile.push_back_real(iw, 3.0_prt);
    }

    pc.fillNeighbors();

    Long total_ghosts = 0;
    for (MFIter mfi = pc.MakeMFIter(0); mfi.isValid(); ++mfi) {
        auto& ptile = pc.ParticlesAt(0, mfi);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ptile.NumRuntimeRealComps() == pc.NumRuntimeRealComps(),
            "fillNeighbors: tile receiving neighbors has no runtime components");
        auto& nbors = pc.GetNeighbors(0, mfi.index(), mfi.LocalTileIndex());
        if (mfi.index() == 1) {
            AMREX_ALWAYS_ASSERT(ptile.numRealParticles() == 0);
            AMREX_ALWAYS_ASSERT(ptile.numNeighborParticles() == 1);
            AMREX_ALWAYS_ASSERT(nbors.numParticles() == 1);
            AMREX_ALWAYS_ASSERT(ptile.GetStructOfArrays().GetRealData(iw)[0] == 3.0_prt);
        } else {
            AMREX_ALWAYS_ASSERT(nbors.numParticles() == 0);
        }
        // accumulate into all ghosts, as an application would
        auto& aos = nbors.GetArrayOfStructs();
        for (int i = 0; i < nbors.numParticles(); ++i) {
            aos[i].rdata(0) += 1.0_prt;
            ++total_ghosts;
        }
    }
    ParallelDescriptor::ReduceLongSum(total_ghosts);
    AMREX_ALWAYS_ASSERT(total_ghosts == 1);

    pc.sumNeighbors(0, 1, 0, 0);

    if (dm[0] == ParallelDescriptor::MyProc()) {
        auto& ptile = pc.ParticlesAt(0, 0, 0);
        AMREX_ALWAYS_ASSERT(ptile.numRealParticles() == 1);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            ptile.GetArrayOfStructs()[0].rdata(0) == 1.0_prt,
            "sumNeighbors: ghost on a tile without real particles was not summed");
    }

    amrex::Print() << "NeighborParticleContainer with empty tiles passed.\n";
}
#endif

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        test_add_comps_on_empty_tiles<ParticleContainer<1, 0, 1, 1>>("AoS");
        test_add_comps_on_empty_tiles<ParticleContainerPureSoA<AMREX_SPACEDIM, 1>>("PureSoA");
#ifndef AMREX_USE_GPU
        test_neighbors_on_empty_tiles();
#endif
    }
    amrex::Finalize();
}
