#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_NeighborParticles.H>

#include <map>
#include <string>
#include <utility>
#include <vector>

using namespace amrex;

// Ghost refresh with and without the mask (setGhostCopyMasked) and with the
// options of the exchange on the CPU (setGhostMpiMasked, setGhostPersistentRecv,
// setGhostDirect, setGhostBatched, setGhostIsend, setGhostCompact), on
// particles that carry per-body contact slots (setGhostSlotGroup).
//
// Every particle has nslots slots. Slot s has a cpu key and an id key (runtime
// ints) and nblocks payload values (runtime reals); the slot is empty when
// both keys are -1. All the values are functions of the particle id and of a
// step number, so a ghost can be checked without looking up its source.
//
// What the refresh has to deliver is the same with and without the mask: the
// components flagged for ghost communication, and the payload of the occupied
// slots. Without the mask the other components are copied too, which is
// checked for the ghosts that stay on the rank. The options change how the
// ghosts are moved, not which ghosts there are or in which order: every
// combination of them is run, after a rebuild and after the updates that
// follow it, and the ghosts of each tile must be the same, in the same
// order, as without any option.

namespace
{
    constexpr int nslots  = 4;
    constexpr int nblocks = 3;

    // runtime reals: payload of block b, slot s at rt_real0 + b*nslots + s,
    // then one real that is not communicated
    constexpr int rt_real0     = 0;
    constexpr int rt_real_free = nblocks*nslots;
    // runtime ints: cpu keys, id keys, then one int that is not communicated
    constexpr int rt_key_cpu = 0;
    constexpr int rt_key_id  = nslots;
    constexpr int rt_int_free = 2*nslots;

    // array components: two reals (the second is not communicated), one int
    constexpr int nar = 2;
    constexpr int nai = 1;

    constexpr ParticleReal sentinel = static_cast<ParticleReal>(-1.0e30);

    ParticleReal rval (Long n, int c, int step)
    {
        return static_cast<ParticleReal>(static_cast<Real>(n % 977) + 0.125*c + 1000.0*step);
    }

    int ival (Long n, int c, int step)
    {
        return static_cast<int>((n*7 + c*13 + step*5) % 10007);
    }

    // 0: empty slot. 1: occupied slot. 2: occupied slot whose id key is -1 but
    // whose cpu key is not, which is not empty.
    int slot_state (Long n, int s, int step)
    {
        if (s == nslots-1 && n % 5 == 0) { return 2; }
        return ((n + 2*s + step) % 3 != 0) ? 1 : 0;
    }

    int key_id_of (Long n, int s, int step)
    {
        int const st = slot_state(n, s, step);
        return (st == 1) ? ival(n, 300+s, 0) + 1 : -1;
    }

    int key_cpu_of (Long n, int s, int step)
    {
        int const st = slot_state(n, s, step);
        if (st == 0) { return -1; }
        return (st == 1) ? ival(n, 400+s, 0) % 4 : 7;
    }

    ParticleReal payload_of (Long n, int b, int s, int step)
    {
        return (slot_state(n, s, step) == 0) ? sentinel : rval(n, 10 + b*nslots + s, step);
    }

    using SlotBase = NeighborParticleContainer<0, 0, nar, nai>;

    class SlotContainer
        : public SlotBase
    {
    public:
        SlotContainer (const Geometry& geom, const DistributionMapping& dm,
                       const BoxArray& ba, int ncells)
            : SlotBase(geom, dm, ba, ncells)
        {
            for (int i = 0; i < nblocks*nslots; ++i) { AddRealComp(true); }
            AddRealComp(false);
            for (int i = 0; i < 2*nslots; ++i) { AddIntComp(true); }
            AddIntComp(false);
            setRealCommComp(RealCommCompStart + 1, false);
            setGhostSlotGroup(nslots, rt_key_cpu, rt_key_id, rt_real0, nblocks);
        }

        // one particle per cell, at its centre
        void initParticles ()
        {
            const int lev = 0;
            const Real* dx  = Geom(lev).CellSize();
            const Real* plo = Geom(lev).ProbLo();

            for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                const Box& tile_box = mfi.tilebox();
                Gpu::HostVector<ParticleType> host_particles;
                for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    AMREX_D_TERM(p.pos(0) = static_cast<ParticleReal>(plo[0] + (Real(iv[0]) + Real(0.5))*dx[0]);,
                                 p.pos(1) = static_cast<ParticleReal>(plo[1] + (Real(iv[1]) + Real(0.5))*dx[1]);,
                                 p.pos(2) = static_cast<ParticleReal>(plo[2] + (Real(iv[2]) + Real(0.5))*dx[2]);)
                    host_particles.push_back(p);
                }

                auto& ptile = DefineAndReturnParticleTile(lev, mfi);
                auto old_size = ptile.GetArrayOfStructs().size();
                ptile.resize(old_size + host_particles.size());
                Gpu::copyAsync(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
                               ptile.GetArrayOfStructs().begin() + old_size);
                Gpu::streamSynchronize();
            }
            setData(0);
        }

        // write the values of step to the real particles
        void setData (int step)
        {
            const int lev = 0;
            for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                auto& ptile = GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
                const int np = ptile.numRealParticles();
                if (np <= 0) { continue; }

                auto const& aos = ptile.GetArrayOfStructs();
                Gpu::HostVector<ParticleType> h_aos(np);
                Gpu::copy(Gpu::deviceToHost, aos().dataPtr(), aos().dataPtr() + np, h_aos.begin());

                auto& soa = ptile.GetStructOfArrays();
                Gpu::HostVector<ParticleReal> hr(np);
                Gpu::HostVector<int> hi(np);

                auto put_real = [&] (int comp, auto&& f)
                {
                    for (int i = 0; i < np; ++i) { hr[i] = f(Long(h_aos[i].id())); }
                    Gpu::copy(Gpu::hostToDevice, hr.begin(), hr.end(), soa.GetRealData(comp).begin());
                };
                auto put_int = [&] (int comp, auto&& f)
                {
                    for (int i = 0; i < np; ++i) { hi[i] = f(Long(h_aos[i].id())); }
                    Gpu::copy(Gpu::hostToDevice, hi.begin(), hi.end(), soa.GetIntData(comp).begin());
                };

                put_real(0, [&] (Long n) { return rval(n, 100, step); });
                put_real(1, [&] (Long n) { return rval(n, 101, step); });
                put_int (0, [&] (Long n) { return ival(n, 200, step); });
                for (int b = 0; b < nblocks; ++b) {
                    for (int s = 0; s < nslots; ++s) {
                        put_real(nar + rt_real0 + b*nslots + s,
                                 [&] (Long n) { return payload_of(n, b, s, step); });
                    }
                }
                put_real(nar + rt_real_free, [&] (Long n) { return rval(n, 102, step); });
                for (int s = 0; s < nslots; ++s) {
                    put_int(nai + rt_key_cpu + s, [&] (Long n) { return key_cpu_of(n, s, step); });
                    put_int(nai + rt_key_id  + s, [&] (Long n) { return key_id_of (n, s, step); });
                }
                put_int(nai + rt_int_free, [&] (Long n) { return ival(n, 201, step); });
            }
        }

        // Check the ghosts against the values of step. Returns their number.
        // full: the components that are not communicated, and the payload of
        // empty slots, must be there too (ghosts that stay on this rank only).
        long checkGhosts (int step, bool full)
        {
            const int lev = 0;
            long nghosts = 0;
            for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                auto& ptile = GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
                const int nreal = ptile.numRealParticles();
                const int nn = ptile.numNeighborParticles();
                if (nn == 0) { continue; }

                auto const& aos = ptile.GetArrayOfStructs();
                Gpu::HostVector<ParticleType> h_aos(nn);
                Gpu::copy(Gpu::deviceToHost, aos().dataPtr() + nreal, aos().dataPtr() + nreal + nn,
                          h_aos.begin());

                auto& soa = ptile.GetStructOfArrays();
                std::map<int, Gpu::HostVector<ParticleReal>> hr;
                std::map<int, Gpu::HostVector<int>> hi;
                for (int c = 0; c < NumRealComps(); ++c) {
                    hr[c].resize(nn);
                    Gpu::copy(Gpu::deviceToHost, soa.GetRealData(c).begin() + nreal,
                              soa.GetRealData(c).begin() + nreal + nn, hr[c].begin());
                }
                for (int c = 0; c < NumIntComps(); ++c) {
                    hi[c].resize(nn);
                    Gpu::copy(Gpu::deviceToHost, soa.GetIntData(c).begin() + nreal,
                              soa.GetIntData(c).begin() + nreal + nn, hi[c].begin());
                }

                for (int i = 0; i < nn; ++i)
                {
                    const Long n = h_aos[i].id();
                    const bool own = full && (h_aos[i].cpu() == ParallelDescriptor::MyProc());

                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hr[0][i] == rval(n, 100, step), "array real 0");
                    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hi[0][i] == ival(n, 200, step), "array int 0");
                    if (own) {
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hr[1][i] == rval(n, 101, step), "array real 1");
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hr[nar + rt_real_free][i] == rval(n, 102, step),
                                                         "runtime real, not communicated");
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hi[nai + rt_int_free][i] == ival(n, 201, step),
                                                         "runtime int, not communicated");
                    }
                    for (int s = 0; s < nslots; ++s) {
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hi[nai + rt_key_cpu + s][i] == key_cpu_of(n, s, step),
                                                         "slot cpu key");
                        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hi[nai + rt_key_id + s][i] == key_id_of(n, s, step),
                                                         "slot id key");
                        if (slot_state(n, s, step) != 0 || own) {
                            for (int b = 0; b < nblocks; ++b) {
                                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                                    hr[nar + rt_real0 + b*nslots + s][i] == payload_of(n, b, s, step),
                                    "slot payload");
                            }
                        }
                    }
                    ++nghosts;
                }
            }
            ParallelDescriptor::ReduceLongSum(nghosts);
            return nghosts;
        }

        // (id, cpu) of the ghosts of every tile in the order they are stored,
        // each tile followed by a (-1, number of ghosts) entry
        std::vector<std::pair<Long,int>> ghostOrder ()
        {
            std::vector<std::pair<Long,int>> order;
            const int lev = 0;
            for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                auto& ptile = GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
                const int nreal = ptile.numRealParticles();
                const int nn = ptile.numNeighborParticles();
                if (nn > 0) {
                    auto const& aos = ptile.GetArrayOfStructs();
                    Gpu::HostVector<ParticleType> h_aos(nn);
                    Gpu::copy(Gpu::deviceToHost, aos().dataPtr() + nreal, aos().dataPtr() + nreal + nn,
                              h_aos.begin());
                    for (int i = 0; i < nn; ++i) {
                        order.emplace_back(Long(h_aos[i].id()), h_aos[i].cpu());
                    }
                }
                order.emplace_back(Long(-1), nn);
            }
            return order;
        }
    };

    // the options of the ghost refresh, as bits
    enum : unsigned {
        OptMasked     = 1u << 0,   // setGhostCopyMasked
        OptMpiMasked  = 1u << 1,   // setGhostMpiMasked
        OptPersistent = 1u << 2,   // setGhostPersistentRecv
        OptDirect     = 1u << 3,   // setGhostDirect
        OptBatched    = 1u << 4,   // setGhostBatched
        OptIsend      = 1u << 5,   // setGhostIsend
        OptCompact    = 1u << 6,   // setGhostCompact
        OptCount      = 1u << 7    // number of combinations
    };

    using GhostOrder = std::vector<std::pair<Long,int>>;

    // Rebuild the ghosts, refresh them three times with new values, rebuild
    // them and refresh them three times more, checking them against the values
    // of the step each time. Returns the order of the ghosts after the second
    // rebuild; with a reference, that order has to be the reference.
    GhostOrder runGhostCopy (unsigned opts, GhostOrder const* reference)
    {
        int const size = 16;
        RealBox real_box;
        for (int n = 0; n < BL_SPACEDIM; n++) {
            real_box.setLo(n, Real(0.0));
            real_box.setHi(n, Real(size));
        }
        const Box domain(IntVect(AMREX_D_DECL(0, 0, 0)),
                         IntVect(AMREX_D_DECL(size-1, size-1, size-1)));
        int is_per[] = {AMREX_D_DECL(1, 1, 1)};
        Geometry geom(domain, &real_box, 0, is_per);

        BoxArray ba(domain);
        ba.maxSize(8);
        DistributionMapping dm(ba);

        SlotContainer pc(geom, dm, ba, 1);
        bool const masked = (opts & OptMasked) != 0;
        pc.setGhostCopyMasked(masked);
        pc.setGhostMpiMasked((opts & OptMpiMasked) != 0);
        pc.setGhostPersistentRecv((opts & OptPersistent) != 0);
        pc.setGhostDirect((opts & OptDirect) != 0);
        pc.setGhostBatched((opts & OptBatched) != 0);
        pc.setGhostIsend((opts & OptIsend) != 0);
        pc.setGhostCompact((opts & OptCompact) != 0);
        AMREX_ALWAYS_ASSERT(pc.ghostCopyMasked() == masked);
        AMREX_ALWAYS_ASSERT(pc.ghostMpiMasked() == ((opts & OptMpiMasked) != 0));
        AMREX_ALWAYS_ASSERT(pc.ghostPersistentRecv() == ((opts & OptPersistent) != 0));
        AMREX_ALWAYS_ASSERT(pc.ghostDirect() == ((opts & OptDirect) != 0));
        AMREX_ALWAYS_ASSERT(pc.ghostBatched() == ((opts & OptBatched) != 0));
        AMREX_ALWAYS_ASSERT(pc.ghostIsend() == ((opts & OptIsend) != 0));
        AMREX_ALWAYS_ASSERT(pc.ghostCompact() == ((opts & OptCompact) != 0));

        // the same ids in every run, so that the order of the ghosts can be compared
        SlotContainer::ParticleType::NextID(1);
        pc.initParticles();
        pc.fillNeighbors();
        long ng = pc.checkGhosts(0, false);
        AMREX_ALWAYS_ASSERT(ng > 0);

        // the slots change state from one step to the next: occupied to
        // empty, empty to occupied, and payload values change
        for (int step = 1; step <= 3; ++step) {
            pc.setData(step);
            pc.updateNeighbors();
            ng = pc.checkGhosts(step, !masked);
            AMREX_ALWAYS_ASSERT(ng > 0);
        }

        // a rebuild makes the layout of the messages and the component lists
        // again, from the particles as they are now
        pc.clearNeighbors();
        pc.fillNeighbors();
        ng = pc.checkGhosts(3, false);
        AMREX_ALWAYS_ASSERT(ng > 0);
        GhostOrder order = pc.ghostOrder();
#ifndef AMREX_USE_GPU
        // on the GPU the order in which the ghosts are stored is not fixed
        if (reference != nullptr) {
            std::string const msg = "ghosts differ from those of the refresh without options, options = "
                                    + std::to_string(opts);
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(order == *reference, msg.c_str());
        }
#else
        amrex::ignore_unused(reference);
#endif
        for (int step = 4; step <= 6; ++step) {
            pc.setData(step);
            pc.updateNeighbors();
            ng = pc.checkGhosts(step, !masked);
            AMREX_ALWAYS_ASSERT(ng > 0);
        }
        return order;
    }
}

void testGhostCopyMasked ()
{
    BL_PROFILE("testGhostCopyMasked");
    amrex::PrintToFile("neighbor_test") << "Running ghost copy test \n";

    // the container reports the options that take effect: not once for each of the runs below
    int const verbose = amrex::Verbose();
    amrex::SetVerbose(0);

    GhostOrder const reference = runGhostCopy(0u, nullptr);
    for (unsigned opts = 1; opts < OptCount; ++opts) {
        runGhostCopy(opts, &reference);
    }
    amrex::SetVerbose(verbose);

    amrex::PrintToFile("neighbor_test") << "  ghost refresh: " << OptCount
                                        << " combinations of the options checked \n";
}
