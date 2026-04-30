#include <AMReX.H>
#include <AMReX_NeighborParticles.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_SPACE.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace amrex;

namespace
{
    constexpr int NumReal = AMREX_SPACEDIM + 2;
    constexpr int NumInt = 2;
    constexpr int MarkerRealComp = AMREX_SPACEDIM;
    constexpr int PayloadRealComp = AMREX_SPACEDIM + 1;
    constexpr int GridIntComp = 0;
    constexpr int MarkerIntComp = 1;
    constexpr ParticleReal Cutoff = 0.2_rt;

    struct TestParams
    {
        IntVect size;
        int max_grid_size = 0;
        int num_ppc = 0;
        int is_periodic = 0;
    };

    struct CheckPair
    {
        template <class P1, class P2>
        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        bool operator() (const P1& p1, const P2& p2) const
        {
            AMREX_D_TERM(Real d0 = p1.pos(0) - p2.pos(0);,
                         Real d1 = p1.pos(1) - p2.pos(1);,
                         Real d2 = p1.pos(2) - p2.pos(2);)
            Real dsquared = AMREX_D_TERM(d0*d0, + d1*d1, + d2*d2);
            return dsquared <= 25.0_rt*Cutoff*Cutoff;
        }
    };

    using ParticleKey = std::pair<Long, int>;

    struct SourceParticleData
    {
        ParticleKey key;
        IntVect cell;
        int grid_id = -1;
        int marker_int = 0;
        int runtime_int = 0;
        std::array<ParticleReal, AMREX_SPACEDIM> pos{};
        ParticleReal marker_real = 0.0_rt;
        ParticleReal payload_real = 0.0_rt;
        ParticleReal runtime_real = 0.0_rt;
    };

    struct PackedSourceParticleData
    {
        Long id = 0;
        int cpu = -1;
        int grid_id = -1;
        int marker_int = 0;
        int runtime_int = 0;
        std::array<int, AMREX_SPACEDIM> cell{};
        std::array<ParticleReal, AMREX_SPACEDIM> pos{};
        ParticleReal marker_real = 0.0_rt;
        ParticleReal payload_real = 0.0_rt;
        ParticleReal runtime_real = 0.0_rt;
    };

    static_assert(std::is_trivially_copyable_v<PackedSourceParticleData>);

    struct InverseContributionData
    {
        ParticleReal marker_real = 0.0_rt;
        ParticleReal payload_real = 0.0_rt;
        ParticleReal runtime_real = 0.0_rt;
        int marker_int = 0;
        int runtime_int = 0;
    };

    struct PackedContributionData
    {
        Long id = 0;
        int cpu = -1;
        ParticleReal marker_real = 0.0_rt;
        ParticleReal payload_real = 0.0_rt;
        ParticleReal runtime_real = 0.0_rt;
        int marker_int = 0;
        int runtime_int = 0;
    };

    static_assert(std::is_trivially_copyable_v<PackedContributionData>);

    struct TileParticleRecord
    {
        ParticleKey key;
        int grid_id = -1;
        IntVect shift = IntVect(0);

        bool operator< (TileParticleRecord const& other) const
        {
            return std::tie(key.first, key.second, grid_id,
                            AMREX_D_DECL(shift[0], shift[1], shift[2]))
                < std::tie(other.key.first, other.key.second, other.grid_id,
                           AMREX_D_DECL(other.shift[0], other.shift[1], other.shift[2]));
        }

        bool operator== (TileParticleRecord const& other) const
        {
            return key == other.key && grid_id == other.grid_id && shift == other.shift;
        }
    };

    void get_position_unit_cell (Real* r, const IntVect& nppc, int i_part)
    {
        int nx = nppc[0];
#if AMREX_SPACEDIM > 1
        int ny = nppc[1];
#else
        int ny = 1;
#endif
#if AMREX_SPACEDIM > 2
        int nz = nppc[2];
#else
        int nz = 1;
#endif

        AMREX_D_TERM(int ix_part = i_part/(ny*nz);,
                     int iy_part = (i_part % (ny*nz)) % ny;,
                     int iz_part = (i_part % (ny*nz)) / ny;)

        AMREX_D_TERM(r[0] = (0.5_rt + ix_part)/nx;,
                     r[1] = (0.5_rt + iy_part)/ny;,
                     r[2] = (0.5_rt + iz_part)/nz;)
    }

    void get_test_params (TestParams& params)
    {
        ParmParse pp("nbor_parts");
        pp.get("size", params.size);
        pp.get("max_grid_size", params.max_grid_size);
        pp.get("num_ppc", params.num_ppc);
        pp.get("is_periodic", params.is_periodic);
    }

    bool almost_equal (ParticleReal lhs, ParticleReal rhs, ParticleReal tol = 1.0e-12_rt)
    {
        return std::abs(lhs-rhs) <= tol;
    }
}

class PureSoANeighborContainer
    : public NeighborParticleContainerPureSoA<NumReal, NumInt>
{
public:
    using Base = NeighborParticleContainerPureSoA<NumReal, NumInt>;
    using ParticleType = typename Base::ParticleType;
    using ParticleTile = typename Base::ParticleTile;

    struct HostTileData
    {
        Gpu::HostVector<uint64_t> idcpu;
        std::array<Gpu::HostVector<ParticleReal>, NumReal> real;
        std::array<Gpu::HostVector<int>, NumInt> idata;
        std::vector<Gpu::HostVector<ParticleReal>> runtime_real;
        std::vector<Gpu::HostVector<int>> runtime_int;
    };

    PureSoANeighborContainer (const Geometry& geom,
                              const DistributionMapping& dmap,
                              const BoxArray& ba,
                              int num_neighbor_cells)
        : Base(geom, dmap, ba, num_neighbor_cells)
    {
        AddRealComp(true);
        AddIntComp(true);
    }

    void InitParticles (const IntVect& nppc, const IntVect& domain_size)
    {
#if AMREX_SPACEDIM == 1
        amrex::ignore_unused(domain_size);
#endif
        const int lev = 0;
        const Real* dx = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();
        const int num_ppc = AMREX_D_TERM(nppc[0], *nppc[1], *nppc[2]);

        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const Box& tile_box = mfi.tilebox();

            Gpu::HostVector<uint64_t> host_idcpu;
            std::array<Gpu::HostVector<ParticleReal>, NumReal> host_real;
            std::array<Gpu::HostVector<int>, NumInt> host_int;
            Gpu::HostVector<ParticleReal> host_runtime_real;
            Gpu::HostVector<int> host_runtime_int;

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
            {
                for (int i_part = 0; i_part < num_ppc; ++i_part) {
                    Real r[AMREX_SPACEDIM];
                    get_position_unit_cell(r, nppc, i_part);

                    int marker = AMREX_D_TERM(iv[0], + domain_size[0]*iv[1],
                                              + domain_size[0]*domain_size[1]*iv[2]);
                    marker = marker*num_ppc + i_part;

                    Long id = ParticleType::NextID();
                    host_idcpu.push_back(0);
                    ParticleIDWrapper(host_idcpu.back()) = id;
                    ParticleCPUWrapper(host_idcpu.back()) = ParallelDescriptor::MyProc();

                    AMREX_D_TERM(host_real[0].push_back(static_cast<ParticleReal>(plo[0] + (iv[0] + r[0])*dx[0]));,
                                 host_real[1].push_back(static_cast<ParticleReal>(plo[1] + (iv[1] + r[1])*dx[1]));,
                                 host_real[2].push_back(static_cast<ParticleReal>(plo[2] + (iv[2] + r[2])*dx[2]));)
                    auto const marker_real = static_cast<ParticleReal>(marker);
                    host_real[MarkerRealComp].push_back(marker_real);
                    host_real[PayloadRealComp].push_back(marker_real + ParticleReal(0.5_rt));

                    host_int[GridIntComp].push_back(mfi.index());
                    host_int[MarkerIntComp].push_back(marker);

                    host_runtime_real.push_back(marker_real + ParticleReal(1.25_rt));
                    host_runtime_int.push_back(marker + 17);
                }
            }

            auto& particle_tile = DefineAndReturnParticleTile(lev, mfi);
            Long old_size = static_cast<Long>(particle_tile.size());
            Long new_size = old_size + static_cast<Long>(host_idcpu.size());
            particle_tile.resize(new_size);

            auto& soa = particle_tile.GetStructOfArrays();
            Gpu::copyAsync(Gpu::hostToDevice,
                           host_idcpu.begin(), host_idcpu.end(),
                           soa.GetIdCPUData().begin() + old_size);

            for (int i = 0; i < NumReal; ++i) {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_real[i].begin(), host_real[i].end(),
                               soa.GetRealData(i).begin() + old_size);
            }
            for (int i = 0; i < NumInt; ++i) {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_int[i].begin(), host_int[i].end(),
                               soa.GetIntData(i).begin() + old_size);
            }
            Gpu::copyAsync(Gpu::hostToDevice,
                           host_runtime_real.begin(), host_runtime_real.end(),
                           soa.GetRealData(NumReal).begin() + old_size);
            Gpu::copyAsync(Gpu::hostToDevice,
                           host_runtime_int.begin(), host_runtime_int.end(),
                           soa.GetIntData(NumInt).begin() + old_size);

            Gpu::streamSynchronize();
        }
        captureOwnedParticles();
    }

    void moveParticles (ParticleReal shift)
    {
        const int lev = 0;
        auto& plev = GetParticles(lev);
        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& ptile = plev[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
            auto ptd = ptile.getParticleTileData();
            const size_t np = ptile.numParticles();
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    ptd.pos(dir, i) += shift;
                }
            });
        }

        for (auto& [key, src] : m_local_source_particles) {
            amrex::ignore_unused(key);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                src.pos[dir] += shift;
                src.cell[dir] = static_cast<int>(
                    amrex::Math::floor((src.pos[dir] - Geom(lev).ProbLoArray()[dir])
                                       * Geom(lev).InvCellSizeArray()[dir]))
                    + Geom(lev).Domain().smallEnd(dir);
            }
        }
    }

    void bumpPayload (int delta)
    {
        const int lev = 0;
        auto& plev = GetParticles(lev);
        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& ptile = plev[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
            auto ptd = ptile.getParticleTileData();
            const size_t np = ptile.numParticles();
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                ptd.rdata(MarkerRealComp)[i] += static_cast<ParticleReal>(delta);
                ptd.rdata(PayloadRealComp)[i] += static_cast<ParticleReal>(delta);
                ptd.idata(MarkerIntComp)[i] += delta;
                ptd.m_runtime_rdata[0][i] += static_cast<ParticleReal>(delta);
                ptd.m_runtime_idata[0][i] += delta;
            });
        }

        for (auto& [key, src] : m_local_source_particles) {
            amrex::ignore_unused(key);
            src.marker_int += delta;
            src.runtime_int += delta;
            src.marker_real += static_cast<ParticleReal>(delta);
            src.payload_real += static_cast<ParticleReal>(delta);
            src.runtime_real += static_cast<ParticleReal>(delta);
        }
    }

    void assertNoNeighbors () const
    {
        for (int lev = 0; lev < numLevels(); ++lev) {
            for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
                auto const& ptile = GetParticles(lev).at(std::make_pair(mfi.index(), mfi.LocalTileIndex()));
                AMREX_ALWAYS_ASSERT(ptile.numNeighborParticles() == 0);
            }
        }
    }

    void checkNeighbors () const
    {
        const int lev = 0;
        auto const& geom = Geom(lev);
        auto const problo = geom.ProbLoArray();
        auto const probhi = geom.ProbHiArray();
        auto const& pshifts = geom.periodicity().shiftIntVect();
        auto const source_particles = gatherSourceParticles();
        auto const& plev = GetParticles(lev);

        auto match_shift = [&] (SourceParticleData const& src,
                                std::array<ParticleReal, AMREX_SPACEDIM> const& pos,
                                Box const& grown_tile_box) -> IntVect
        {
            IntVect matched_shift(std::numeric_limits<int>::max());
            int nmatches = 0;
            for (auto const& shift : pshifts) {
                if (!grown_tile_box.contains(src.cell + shift)) { continue; }

                bool matched = true;
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    ParticleReal expected_pos = src.pos[dir];
                    if (shift[dir] > 0) {
                        expected_pos += static_cast<ParticleReal>(probhi[dir] - problo[dir]);
                    } else if (shift[dir] < 0) {
                        expected_pos -= static_cast<ParticleReal>(probhi[dir] - problo[dir]);
                    }
                    if (!almost_equal(pos[dir], expected_pos)) {
                        matched = false;
                        break;
                    }
                }

                if (matched) {
                    matched_shift = shift;
                    ++nmatches;
                }
            }

            AMREX_ALWAYS_ASSERT(nmatches == 1);
            return matched_shift;
        };

        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto const index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            auto const& ptile = plev.at(index);
            int const np = ptile.numParticles();
            int const np_total = ptile.numTotalParticles();

            HostTileData host;
            copyTileToHost(ptile, np_total, host);

            std::set<TileParticleRecord> expected_records;
            std::set<TileParticleRecord> actual_records;

            Box tile_box = mfi.tilebox();
            Box grown_tile_box = tile_box;
            grown_tile_box.grow(1);

            for (auto const& [key, src] : source_particles) {
                for (auto const& shift : pshifts) {
                    if (!grown_tile_box.contains(src.cell + shift)) { continue; }
                    if (tile_box.contains(src.cell + shift)) { continue; }
                    expected_records.insert(TileParticleRecord{key, src.grid_id, shift});
                }
            }

            for (int i = np; i < np_total; ++i) {
                ParticleKey key{
                    static_cast<Long>(ParticleIDWrapper(host.idcpu[i])),
                    static_cast<int>(ParticleCPUWrapper(host.idcpu[i]))
                };
                auto it = source_particles.find(key);
                AMREX_ALWAYS_ASSERT(it != source_particles.end());
                auto const& src = it->second;

                AMREX_ALWAYS_ASSERT(host.idata[GridIntComp][i] == src.grid_id);
                AMREX_ALWAYS_ASSERT(host.idata[MarkerIntComp][i] == src.marker_int);
                AMREX_ALWAYS_ASSERT(almost_equal(host.real[MarkerRealComp][i], src.marker_real));
                AMREX_ALWAYS_ASSERT(almost_equal(host.real[PayloadRealComp][i], src.payload_real));
                AMREX_ALWAYS_ASSERT(almost_equal(host.runtime_real[0][i], src.runtime_real));
                AMREX_ALWAYS_ASSERT(host.runtime_int[0][i] == src.runtime_int);

                std::array<ParticleReal, AMREX_SPACEDIM> pos{
                    AMREX_D_DECL(host.real[0][i], host.real[1][i], host.real[2][i])
                };
                IntVect shift = match_shift(src, pos, grown_tile_box);
                actual_records.insert(TileParticleRecord{key, src.grid_id, shift});
            }

            AMREX_ALWAYS_ASSERT(actual_records == expected_records);
        }
    }

    void checkNeighborList () const
    {
        const int lev = 0;
        auto const& plev = GetParticles(lev);
        Real const cutoff_sq = 25.0_rt*Cutoff*Cutoff;

        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto const index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            auto const& ptile = plev.at(index);
            int const np = ptile.numParticles();
            int const np_total = ptile.numTotalParticles();

            HostTileData host;
            copyTileToHost(ptile, np_total, host);

            Gpu::HostVector<unsigned int> full_count(np, 0);
            Gpu::HostVector<unsigned int> full_nbors;

            for (int i = 0; i < np; ++i) {
                for (int j = 0; j < np_total; ++j) {
                    if (i == j) { continue; }

                    Real dsquared = 0.0_rt;
                    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                        Real d = host.real[dir][i] - host.real[dir][j];
                        dsquared += d*d;
                    }
                    if (dsquared <= cutoff_sq) {
                        ++full_count[i];
                        full_nbors.push_back(j);
                    }
                }
            }

            auto const& d_counts = m_neighbor_list.at(lev).at(index).GetCounts();
            Gpu::HostVector<unsigned int> h_counts(d_counts.size());
            Gpu::copy(Gpu::deviceToHost, d_counts.begin(), d_counts.end(), h_counts.begin());

            auto const& d_list = m_neighbor_list.at(lev).at(index).GetList();
            Gpu::HostVector<unsigned int> h_list(d_list.size());
            Gpu::copy(Gpu::deviceToHost, d_list.begin(), d_list.end(), h_list.begin());

            for (int i = 0; i < np; ++i) {
                AMREX_ALWAYS_ASSERT(h_counts[i] == full_count[i]);
            }

            unsigned int start = 0;
            for (int i = 0; i < np; ++i) {
                std::sort(full_nbors.begin() + start, full_nbors.begin() + start + full_count[i]);
                std::sort(h_list.begin() + start, h_list.begin() + start + full_count[i]);
                start += full_count[i];
            }

            for (unsigned int i = 0; i < h_list.size(); ++i) {
                AMREX_ALWAYS_ASSERT(h_list[i] == full_nbors[i]);
            }
        }
    }

    void checkInverseSumNeighbors ()
    {
        constexpr ParticleReal marker_real_delta = 1.25_rt;
        constexpr ParticleReal payload_real_delta = 2.5_rt;
        constexpr ParticleReal runtime_real_delta = 3.75_rt;
        constexpr int marker_int_delta = 11;
        constexpr int runtime_int_delta = 13;

        const int lev = 0;
        std::map<ParticleKey, InverseContributionData> local_contributions;

        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& neighb_tile = GetNeighbors(lev, mfi.index(), mfi.LocalTileIndex());
            int const nneighbs = static_cast<int>(neighb_tile.size());
            if (nneighbs == 0) { continue; }

            HostTileData host;
            copyTileToHost(neighb_tile, nneighbs, host);

            for (int i = 0; i < nneighbs; ++i) {
                ParticleKey key{
                    static_cast<Long>(ParticleIDWrapper(host.idcpu[i])),
                    static_cast<int>(ParticleCPUWrapper(host.idcpu[i]))
                };
                auto& contribution = local_contributions[key];
                contribution.marker_real += host.real[MarkerRealComp][i] + marker_real_delta;
                contribution.payload_real += host.real[PayloadRealComp][i] + payload_real_delta;
                contribution.runtime_real += host.runtime_real[0][i] + runtime_real_delta;
                contribution.marker_int += host.idata[MarkerIntComp][i] + marker_int_delta;
                contribution.runtime_int += host.runtime_int[0][i] + runtime_int_delta;
            }

            auto ptd = neighb_tile.getParticleTileData();
            amrex::ParallelFor(nneighbs, [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                ptd.rdata(MarkerRealComp)[i] += marker_real_delta;
                ptd.rdata(PayloadRealComp)[i] += payload_real_delta;
                ptd.m_runtime_rdata[0][i] += runtime_real_delta;
                ptd.idata(MarkerIntComp)[i] += marker_int_delta;
                ptd.m_runtime_idata[0][i] += runtime_int_delta;
            });
        }

        Gpu::streamSynchronize();

        auto const global_contributions = gatherContributionData(local_contributions);

        int total_contributions = 0;
        for (auto const& [key, contribution] : global_contributions) {
            amrex::ignore_unused(key);
            total_contributions += contribution.marker_int;
        }
        AMREX_ALWAYS_ASSERT(total_contributions > 0);

        sumNeighbors(MarkerRealComp, 3, MarkerIntComp, 2);

        auto const& plev = GetParticles(lev);
        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto const index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            auto const& ptile = plev.at(index);
            int const np = ptile.numParticles();
            if (np == 0) { continue; }

            HostTileData host;
            copyTileToHost(ptile, np, host);

            for (int i = 0; i < np; ++i) {
                ParticleKey key{
                    static_cast<Long>(ParticleIDWrapper(host.idcpu[i])),
                    static_cast<int>(ParticleCPUWrapper(host.idcpu[i]))
                };
                auto const src_it = m_local_source_particles.find(key);
                AMREX_ALWAYS_ASSERT(src_it != m_local_source_particles.end());
                auto const& src = src_it->second;

                auto const contribution_it = global_contributions.find(key);
                InverseContributionData contribution;
                if (contribution_it != global_contributions.end()) {
                    contribution = contribution_it->second;
                }

                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    AMREX_ALWAYS_ASSERT(almost_equal(host.real[dir][i], src.pos[dir]));
                }
                AMREX_ALWAYS_ASSERT(host.idata[GridIntComp][i] == src.grid_id);
                AMREX_ALWAYS_ASSERT(host.idata[MarkerIntComp][i] == src.marker_int + contribution.marker_int);
                AMREX_ALWAYS_ASSERT(host.runtime_int[0][i] == src.runtime_int + contribution.runtime_int);
                AMREX_ALWAYS_ASSERT(almost_equal(
                    host.real[MarkerRealComp][i], src.marker_real + contribution.marker_real));
                AMREX_ALWAYS_ASSERT(almost_equal(
                    host.real[PayloadRealComp][i], src.payload_real + contribution.payload_real));
                AMREX_ALWAYS_ASSERT(almost_equal(
                    host.runtime_real[0][i], src.runtime_real + contribution.runtime_real));
            }
        }

        captureOwnedParticles();
    }

private:
    void captureOwnedParticles ()
    {
        const int lev = 0;
        auto const plo = Geom(lev).ProbLoArray();
        auto const dxi = Geom(lev).InvCellSizeArray();
        auto const domain = Geom(lev).Domain();
        auto const& plev = GetParticles(lev);

        m_local_source_particles.clear();

        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto const index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            auto const& ptile = plev.at(index);
            int const np = ptile.numParticles();
            if (np == 0) { continue; }

            HostTileData host;
            copyTileToHost(ptile, np, host);

            for (int i = 0; i < np; ++i) {
                SourceParticleData src;
                src.key = {ParticleIDWrapper(host.idcpu[i]), ParticleCPUWrapper(host.idcpu[i])};
                src.grid_id = host.idata[GridIntComp][i];
                src.marker_int = host.idata[MarkerIntComp][i];
                src.runtime_int = host.runtime_int[0][i];
                src.marker_real = host.real[MarkerRealComp][i];
                src.payload_real = host.real[PayloadRealComp][i];
                src.runtime_real = host.runtime_real[0][i];
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    src.pos[dir] = host.real[dir][i];
                    src.cell[dir] = static_cast<int>(
                        amrex::Math::floor((src.pos[dir] - plo[dir])*dxi[dir])) + domain.smallEnd(dir);
                }

                auto [it, inserted] = m_local_source_particles.emplace(src.key, src);
                amrex::ignore_unused(it);
                AMREX_ALWAYS_ASSERT(inserted);
            }
        }
    }

    void copyTileToHost (ParticleTile const& ptile, int np, HostTileData& host) const
    {
        auto const& soa = ptile.GetStructOfArrays();

        host.idcpu.resize(np);
        Gpu::copy(Gpu::deviceToHost,
                  soa.GetIdCPUData().begin(), soa.GetIdCPUData().begin() + np,
                  host.idcpu.begin());

        for (int i = 0; i < NumReal; ++i) {
            host.real[i].resize(np);
            Gpu::copy(Gpu::deviceToHost,
                      soa.GetRealData(i).begin(), soa.GetRealData(i).begin() + np,
                      host.real[i].begin());
        }

        for (int i = 0; i < NumInt; ++i) {
            host.idata[i].resize(np);
            Gpu::copy(Gpu::deviceToHost,
                      soa.GetIntData(i).begin(), soa.GetIntData(i).begin() + np,
                      host.idata[i].begin());
        }

        host.runtime_real.resize(NumRuntimeRealComps());
        for (int i = 0; i < NumRuntimeRealComps(); ++i) {
            host.runtime_real[i].resize(np);
            Gpu::copy(Gpu::deviceToHost,
                      soa.GetRealData(NumReal + i).begin(), soa.GetRealData(NumReal + i).begin() + np,
                      host.runtime_real[i].begin());
        }

        host.runtime_int.resize(NumRuntimeIntComps());
        for (int i = 0; i < NumRuntimeIntComps(); ++i) {
            host.runtime_int[i].resize(np);
            Gpu::copy(Gpu::deviceToHost,
                      soa.GetIntData(NumInt + i).begin(), soa.GetIntData(NumInt + i).begin() + np,
                      host.runtime_int[i].begin());
        }
    }

    static std::map<ParticleKey, InverseContributionData>
    gatherContributionData (std::map<ParticleKey, InverseContributionData> const& local_contributions)
    {
        std::vector<PackedContributionData> local_packed;
        local_packed.reserve(local_contributions.size());
        for (auto const& [key, contribution] : local_contributions) {
            local_packed.push_back(PackedContributionData{
                key.first, key.second,
                contribution.marker_real, contribution.payload_real, contribution.runtime_real,
                contribution.marker_int, contribution.runtime_int
            });
        }

        auto unpack = [] (std::vector<PackedContributionData> const& packed) {
            std::map<ParticleKey, InverseContributionData> result;
            for (auto const& item : packed) {
                auto& contribution = result[{item.id, item.cpu}];
                contribution.marker_real += item.marker_real;
                contribution.payload_real += item.payload_real;
                contribution.runtime_real += item.runtime_real;
                contribution.marker_int += item.marker_int;
                contribution.runtime_int += item.runtime_int;
            }
            return result;
        };

        if (ParallelDescriptor::NProcs() == 1) {
            return unpack(local_packed);
        }

        int const root = ParallelDescriptor::IOProcessorNumber();
        int const local_nbytes =
            static_cast<int>(local_packed.size() * sizeof(PackedContributionData));
        auto recvcounts = ParallelDescriptor::Gather(local_nbytes, root);

        std::vector<int> displs;
        std::vector<char> recv_bytes;

        if (ParallelDescriptor::MyProc() == root) {
            displs.resize(recvcounts.size(), 0);
            int offset = 0;
            for (int i = 0, n = static_cast<int>(recvcounts.size()); i < n; ++i) {
                displs[i] = offset;
                offset += recvcounts[i];
            }
            recv_bytes.resize(offset);
        }

        auto const* send_ptr = reinterpret_cast<char const*>(local_packed.data());
        ParallelDescriptor::Gatherv(send_ptr, local_nbytes,
                                    recv_bytes.data(), recvcounts, displs, root);

        std::vector<PackedContributionData> global_packed;
        if (ParallelDescriptor::MyProc() == root) {
            AMREX_ALWAYS_ASSERT(recv_bytes.size() % sizeof(PackedContributionData) == 0);
            std::vector<PackedContributionData> gathered(
                recv_bytes.size() / sizeof(PackedContributionData));
            if (!recv_bytes.empty()) {
                std::memcpy(gathered.data(), recv_bytes.data(), recv_bytes.size());
            }

            auto const global_contributions = unpack(gathered);
            global_packed.reserve(global_contributions.size());
            for (auto const& [key, contribution] : global_contributions) {
                global_packed.push_back(PackedContributionData{
                    key.first, key.second,
                    contribution.marker_real, contribution.payload_real, contribution.runtime_real,
                    contribution.marker_int, contribution.runtime_int
                });
            }
        }

        int total_nbytes =
            static_cast<int>(global_packed.size() * sizeof(PackedContributionData));
        ParallelDescriptor::Bcast(&total_nbytes, 1, root);

        std::vector<char> packed_bytes(total_nbytes);
        if (ParallelDescriptor::MyProc() == root && total_nbytes > 0) {
            std::memcpy(packed_bytes.data(), global_packed.data(), total_nbytes);
        }

        if (total_nbytes > 0) {
            ParallelDescriptor::Bcast(packed_bytes.data(), total_nbytes, root);
        }

        AMREX_ALWAYS_ASSERT(total_nbytes % sizeof(PackedContributionData) == 0);
        std::vector<PackedContributionData> unpacked(
            total_nbytes / static_cast<int>(sizeof(PackedContributionData)));
        if (total_nbytes > 0) {
            std::memcpy(unpacked.data(), packed_bytes.data(), total_nbytes);
        }
        return unpack(unpacked);
    }

    std::map<ParticleKey, SourceParticleData> gatherSourceParticles () const
    {
        std::vector<PackedSourceParticleData> local_source_particles;

        local_source_particles.reserve(m_local_source_particles.size());
        for (auto const& [key, src] : m_local_source_particles) {
            PackedSourceParticleData pdata;
            pdata.id = key.first;
            pdata.cpu = key.second;
            pdata.grid_id = src.grid_id;
            pdata.marker_int = src.marker_int;
            pdata.runtime_int = src.runtime_int;
            pdata.marker_real = src.marker_real;
            pdata.payload_real = src.payload_real;
            pdata.runtime_real = src.runtime_real;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                pdata.pos[dir] = src.pos[dir];
                pdata.cell[dir] = src.cell[dir];
            }
            local_source_particles.push_back(pdata);
        }

        auto unpack = [] (std::vector<PackedSourceParticleData> const& packed) {
            std::map<ParticleKey, SourceParticleData> result;
            for (auto const& p : packed) {
                SourceParticleData src;
                src.key = {p.id, p.cpu};
                src.grid_id = p.grid_id;
                src.marker_int = p.marker_int;
                src.runtime_int = p.runtime_int;
                src.marker_real = p.marker_real;
                src.payload_real = p.payload_real;
                src.runtime_real = p.runtime_real;
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    src.cell[dir] = p.cell[dir];
                    src.pos[dir] = p.pos[dir];
                }
                auto [it, inserted] = result.emplace(src.key, src);
                amrex::ignore_unused(it);
                AMREX_ALWAYS_ASSERT(inserted);
            }
            return result;
        };

        if (ParallelDescriptor::NProcs() == 1) {
            return unpack(local_source_particles);
        }

        int const root = ParallelDescriptor::IOProcessorNumber();
        int const local_nbytes =
            static_cast<int>(local_source_particles.size() * sizeof(PackedSourceParticleData));
        auto recvcounts = ParallelDescriptor::Gather(local_nbytes, root);

        std::vector<int> displs;
        std::vector<char> packed_bytes;

        if (ParallelDescriptor::MyProc() == root) {
            displs.resize(recvcounts.size(), 0);
            int offset = 0;
            for (int i = 0, n = static_cast<int>(recvcounts.size()); i < n; ++i) {
                displs[i] = offset;
                offset += recvcounts[i];
            }
            packed_bytes.resize(offset);
        }

        auto const* send_ptr = reinterpret_cast<char const*>(local_source_particles.data());
        ParallelDescriptor::Gatherv(send_ptr, local_nbytes,
                                    packed_bytes.data(), recvcounts, displs, root);

        int total_nbytes = static_cast<int>(packed_bytes.size());
        ParallelDescriptor::Bcast(&total_nbytes, 1, root);
        if (ParallelDescriptor::MyProc() != root) {
            packed_bytes.resize(total_nbytes);
        }
        if (total_nbytes > 0) {
            ParallelDescriptor::Bcast(packed_bytes.data(), total_nbytes, root);
        }

        AMREX_ALWAYS_ASSERT(total_nbytes % sizeof(PackedSourceParticleData) == 0);
        std::vector<PackedSourceParticleData> global_source_particles(
            total_nbytes / static_cast<int>(sizeof(PackedSourceParticleData)));
        if (total_nbytes > 0) {
            std::memcpy(global_source_particles.data(), packed_bytes.data(), total_nbytes);
        }
        return unpack(global_source_particles);
    }

    std::map<ParticleKey, SourceParticleData> m_local_source_particles;
};

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    TestParams params;
    get_test_params(params);
    AMREX_ALWAYS_ASSERT(params.num_ppc == 1);

    RealBox real_box;
    for (int dir = 0; dir < BL_SPACEDIM; ++dir)
    {
        real_box.setLo(dir, 0.0);
        real_box.setHi(dir, params.size[dir]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1, params.size[1]-1, params.size[2]-1));
    Box const domain(domain_lo, domain_hi);

    int is_per[] = {AMREX_D_DECL(params.is_periodic, params.is_periodic, params.is_periodic)};
    Geometry geom(domain, &real_box, 0, is_per);

    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    PureSoANeighborContainer pc(geom, dm, ba, 1);
    pc.InitParticles(IntVect(params.num_ppc), params.size);
    pc.assertNoNeighbors();

    pc.fillNeighbors();
    pc.checkNeighbors();
    pc.buildNeighborList(CheckPair());
    pc.checkNeighborList();

    pc.bumpPayload(1000);
    pc.updateNeighbors();
    pc.checkNeighbors();
    pc.buildNeighborList(CheckPair());
    pc.checkNeighborList();

    pc.moveParticles(0.1_rt);
    pc.bumpPayload(2000);
    pc.updateNeighbors();
    pc.checkNeighbors();
    pc.buildNeighborList(CheckPair());
    pc.checkNeighborList();

#ifndef AMREX_USE_GPU
    // sumNeighbors() is currently only implemented on the CPU path.
    pc.setEnableInverse(true);
    pc.fillNeighbors();
    pc.checkInverseSumNeighbors();
#endif

    amrex::Finalize();
}
