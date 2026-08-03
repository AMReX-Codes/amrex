#include <AMReX_Reduce.H>
#include <AMReX_SPACE.H>

#include "MDParticleContainer.H"
#include "Constants.H"
#include "CheckPair.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <tuple>
#include <type_traits>

using namespace amrex;

namespace
{
    using ParticleKey = std::pair<Long, int>;

    struct SourceParticleData
    {
        ParticleKey key;
        IntVect cell;
        int grid_id = -1;
        std::array<ParticleReal, AMREX_SPACEDIM> pos{};
    };

    struct PackedSourceParticleData
    {
        Long id = 0;
        int cpu = -1;
        int grid_id = -1;
        std::array<int, AMREX_SPACEDIM> cell{};
        std::array<ParticleReal, AMREX_SPACEDIM> pos{};
    };

    static_assert(std::is_trivially_copyable_v<PackedSourceParticleData>);

    struct OwnedParticleData
    {
        ParticleKey key;
        std::array<ParticleReal, AMREX_SPACEDIM> pos{};
        std::array<ParticleReal, PIdx::ncomps> struct_real{};
        int struct_int = 0;
        ParticleReal array_real = 0.0_rt;
        ParticleReal runtime_real = 0.0_rt;
        int array_int = 0;
        int runtime_int = 0;
    };

    struct InverseContributionData
    {
        std::array<ParticleReal, PIdx::ncomps> struct_real{};
        int struct_int = 0;
        ParticleReal array_real = 0.0_rt;
        ParticleReal runtime_real = 0.0_rt;
        int array_int = 0;
        int runtime_int = 0;
    };

    struct PackedContributionData
    {
        Long id = 0;
        int cpu = -1;
        std::array<ParticleReal, PIdx::ncomps> struct_real{};
        int struct_int = 0;
        ParticleReal array_real = 0.0_rt;
        ParticleReal runtime_real = 0.0_rt;
        int array_int = 0;
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

    bool almost_equal (ParticleReal lhs, ParticleReal rhs)
    {
        auto const scale = std::max({ParticleReal(1.0), std::abs(lhs), std::abs(rhs)});
        auto const tol = std::max(ParticleReal(1.0e-12),
                                  ParticleReal(64.0) * std::numeric_limits<ParticleReal>::epsilon());
        return std::abs(lhs - rhs) <= tol * scale;
    }

    void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
    {
        int nx = nppc[0];
#if AMREX_SPACEDIM >= 2
        int ny = nppc[1];
#else
        int ny = 1;
#endif
#if AMREX_SPACEDIM == 3
        int nz = nppc[2];
#else
        int nz = 1;
#endif

        AMREX_D_TERM(int ix_part = i_part/(ny * nz);,
                     int iy_part = (i_part % (ny * nz)) % ny;,
                     int iz_part = (i_part % (ny * nz)) / ny;)

        AMREX_D_TERM(r[0] = (0.5+ix_part)/nx;,
                     r[1] = (0.5+iy_part)/ny;,
                     r[2] = (0.5+iz_part)/nz;)
    }

    void get_gaussian_random_momentum(Real* u, Real u_mean, Real u_std) {
        Real ux_th = amrex::RandomNormal(0.0, u_std);
        Real uy_th = amrex::RandomNormal(0.0, u_std);
        Real uz_th = amrex::RandomNormal(0.0, u_std);

        u[0] = u_mean + ux_th;
        u[1] = u_mean + uy_th;
        u[2] = u_mean + uz_th;
    }

    auto unpack_source_particles (std::vector<PackedSourceParticleData> const& packed_particles)
        -> std::map<ParticleKey, SourceParticleData>
    {
        std::map<ParticleKey, SourceParticleData> source_particles;

        for (auto const& packed : packed_particles) {
            ParticleKey key{packed.id, packed.cpu};
            SourceParticleData pdata;
            pdata.key = key;
            pdata.grid_id = packed.grid_id;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                pdata.cell[idim] = packed.cell[idim];
                pdata.pos[idim] = packed.pos[idim];
            }
            auto [it, inserted] = source_particles.emplace(key, pdata);
            amrex::ignore_unused(it);
            AMREX_ALWAYS_ASSERT(inserted);
        }

        return source_particles;
    }
}

void
MDParticleContainer::
InitParticles(const IntVect& a_num_particles_per_cell,
              Real           a_thermal_momentum_std,
              Real           a_thermal_momentum_mean)
{
    BL_PROFILE("MDParticleContainer::InitParticles");

    amrex::PrintToFile("neighbor_test") << "Generating particles... ";

    const int lev = 0;
    const Real* dx = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                     *a_num_particles_per_cell[1],
                                     *a_num_particles_per_cell[2]);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        const Box& tile_box  = mfi.tilebox();

        Gpu::HostVector<ParticleType> host_particles;
        std::vector<Gpu::HostVector<ParticleReal> > host_real(NumRealComps());
        std::vector<Gpu::HostVector<int> > host_int(NumIntComps());

        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                Real r[AMREX_SPACEDIM];
                Real v[3];

                get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                get_gaussian_random_momentum(v, a_thermal_momentum_mean,
                                             a_thermal_momentum_std);

                AMREX_D_TERM(auto x = static_cast<ParticleReal> (plo[0] + (iv[0] + r[0])*dx[0]);,
                             auto y = static_cast<ParticleReal> (plo[1] + (iv[1] + r[1])*dx[1]);,
                             auto z = static_cast<ParticleReal> (plo[2] + (iv[2] + r[2])*dx[2]);)

                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                AMREX_D_TERM(p.pos(0) = x;,
                             p.pos(1) = y;,
                             p.pos(2) = z;)

                p.rdata(PIdx::vx) = static_cast<ParticleReal> (v[0]);
                p.rdata(PIdx::vy) = static_cast<ParticleReal> (v[1]);
                p.rdata(PIdx::vz) = static_cast<ParticleReal> (v[2]);

                p.rdata(PIdx::ax) = 0.0;
                p.rdata(PIdx::ay) = 0.0;
                p.rdata(PIdx::az) = 0.0;

                p.idata(0) = mfi.index();

                host_particles.push_back(p);
                for (int i = 0; i < NumRealComps(); ++i) {
                    host_real[i].push_back(ParticleReal(mfi.index()));
                }
                for (int i = 0; i < NumIntComps(); ++i) {
                    host_int[i].push_back(mfi.index());
                }
            }
        }

        auto& particle_tile = DefineAndReturnParticleTile(lev, mfi);
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copyAsync(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
                       particle_tile.GetArrayOfStructs().begin() + old_size);

        auto& soa = particle_tile.GetStructOfArrays();
        for (int i = 0; i < NumRealComps(); ++i)
        {
            Gpu::copyAsync(Gpu::hostToDevice,
                           host_real[i].begin(),
                           host_real[i].end(),
                           soa.GetRealData(i).begin() + old_size);
        }

        for (int i = 0; i < NumIntComps(); ++i)
        {
            Gpu::copyAsync(Gpu::hostToDevice,
                           host_int[i].begin(),
                           host_int[i].end(),
                           soa.GetIntData(i).begin() + old_size);
        }

        Gpu::streamSynchronize();
    }

    amrex::PrintToFile("neighbor_test") << " Number of particles is " << this->TotalNumberOfParticles()<< " \n";
    amrex::PrintToFile("neighbor_test") << "done. \n";
}

std::pair<Real, Real> MDParticleContainer::minAndMaxDistance()
{
    BL_PROFILE("MDParticleContainer::minAndMaxDistance");

    const int lev = 0;
    auto& plev  = GetParticles(lev);

    ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op;
    ReduceData<ParticleReal, ParticleReal> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        ParticleReal min_start = std::numeric_limits<ParticleReal>::max();
        ParticleReal max_start = std::numeric_limits<ParticleReal>::lowest();

        reduce_op.eval(aos.numParticles(), reduce_data,
                       [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                       {
                           ParticleType& p1 = pstruct[i];

                           ParticleReal min_d = min_start;
                           ParticleReal max_d = max_start;

                           for (const auto& p2 : nbor_data.getNeighbors(i))
                           {
                               AMREX_D_TERM(ParticleReal dx = p1.pos(0) - p2.pos(0);,
                                            ParticleReal dy = p1.pos(1) - p2.pos(1);,
                                            ParticleReal dz = p1.pos(2) - p2.pos(2);)

                               ParticleReal r2 = AMREX_D_TERM(dx*dx, + dy*dy, + dz*dz);
                               r2 = amrex::max(r2, Params::min_r*Params::min_r);
                               auto r = ParticleReal(std::sqrt(r2));

                               min_d = std::min(min_d, r);
                               max_d = std::max(max_d, r);
                           }
                           return {min_d, max_d};
                       });
    }

    ParticleReal min_d = amrex::get<0>(reduce_data.value(reduce_op));
    ParticleReal max_d = amrex::get<1>(reduce_data.value(reduce_op));
    ParallelDescriptor::ReduceRealMin(min_d, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(max_d, ParallelDescriptor::IOProcessorNumber());

    return std::make_pair(min_d, max_d);
}

void MDParticleContainer::moveParticles(amrex::ParticleReal dx)
{
    BL_PROFILE("MDParticleContainer::moveParticles");

    const int lev = 0;
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();

        auto& ptile = plev[std::make_pair(gid, tid)];
        auto& aos   = ptile.GetArrayOfStructs();
        ParticleType* pstruct = aos.data();

        const size_t np = aos.numParticles();

        // now we move the particles
        AMREX_FOR_1D ( np, i,
        {
            ParticleType& p = pstruct[i];
            AMREX_D_TERM(p.pos(0) += dx;,
                         p.pos(1) += dx;,
                         p.pos(2) += dx;)
        });
    }
}

void MDParticleContainer::writeParticles(int n)
{
    BL_PROFILE("MDParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}

void MDParticleContainer::checkNeighborParticles(bool use_source_grid)
{
    BL_PROFILE("MDParticleContainer::checkNeighborParticles");

    const int lev = 0;
    auto const& geom = Geom(lev);
    auto const dxi = geom.InvCellSizeArray();
    auto const domain = geom.Domain();
    auto const problo = geom.ProbLoArray();
    auto const probhi = geom.ProbHiArray();
    const auto& pshifts = geom.periodicity().shiftIntVect();
    auto& plev  = GetParticles(lev);

    auto ngrids = int(ParticleBoxArray(0).size());
    std::vector<PackedSourceParticleData> local_source_particles;

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& ptile = plev[index];
        auto const& aos = ptile.GetArrayOfStructs();
        int const np = aos.numParticles();

        if (np == 0) { continue; }

        Gpu::HostVector<ParticleType> h_pstruct(np);
        Gpu::copy(Gpu::deviceToHost, aos().dataPtr(), aos().dataPtr() + np, h_pstruct.begin());

        for (int i = 0; i < np; ++i) {
            auto const& p = h_pstruct[i];
            PackedSourceParticleData pdata;
            pdata.id = p.id();
            pdata.cpu = p.cpu();
            pdata.grid_id = p.idata(0);
            auto const cell = getParticleCell(p, problo, dxi, domain);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                pdata.cell[idim] = cell[idim];
                pdata.pos[idim] = p.pos(idim);
            }
            local_source_particles.push_back(pdata);
        }
    }

    std::map<ParticleKey, SourceParticleData> source_particles;
    if (ParallelDescriptor::NProcs() == 1) {
        source_particles = unpack_source_particles(local_source_particles);
    } else {
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
        source_particles = unpack_source_particles(global_source_particles);
    }

    auto match_shift = [&] (ParticleType const& p, SourceParticleData const& src,
                            Box const& grown_tile_box) -> IntVect
    {
        IntVect matched_shift(std::numeric_limits<int>::max());
        int nmatches = 0;
        for (auto const& shift : pshifts) {
            if (!grown_tile_box.contains(src.cell + shift)) { continue; }
            bool matched = true;
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                ParticleReal expected_pos = src.pos[idim];
                if (shift[idim] > 0) {
                    expected_pos += static_cast<ParticleReal>(probhi[idim] - problo[idim]);
                } else if (shift[idim] < 0) {
                    expected_pos -= static_cast<ParticleReal>(probhi[idim] - problo[idim]);
                }
                if (!almost_equal(p.pos(idim), expected_pos)) {
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

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int const gid = mfi.index();
        int const tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);
        auto& ptile = plev[index];
        auto const& aos = ptile.GetArrayOfStructs();
        auto const& soa = ptile.GetStructOfArrays();

        int const np_total = aos.numTotalParticles();
        if (np_total == 0) { continue; }

        Gpu::HostVector<ParticleType> h_pstruct(np_total);
        Gpu::HostVector<ParticleReal> h_rdata(np_total);
        Gpu::HostVector<int> h_idata(np_total);
        Gpu::copy(Gpu::deviceToHost, aos().dataPtr(), aos().dataPtr() + np_total, h_pstruct.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetRealData(0).begin(), soa.GetRealData(0).begin() + np_total,
                  h_rdata.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetIntData(0).begin(), soa.GetIntData(0).begin() + np_total,
                  h_idata.begin());

        Box grown_tile_box = mfi.tilebox();
        if (hasNeighbors()) {
            grown_tile_box.grow(m_num_neighbor_cells);
        }

        std::vector<int> expected_per_grid(ngrids, 0);
        std::vector<int> actual_per_grid(ngrids, 0);
        std::vector<TileParticleRecord> expected_records;
        std::vector<TileParticleRecord> actual_records;

        for (auto const& [key, src] : source_particles) {
            for (auto const& shift : pshifts) {
                if (!grown_tile_box.contains(src.cell + shift)) { continue; }
                int const expected_grid = use_source_grid ? src.grid_id : gid;
                expected_records.push_back({key, expected_grid, shift});
                ++expected_per_grid[expected_grid];
            }
        }

        for (int i = 0; i < np_total; ++i) {
            auto const& p = h_pstruct[i];
            ParticleKey key{p.id(), p.cpu()};
            auto it = source_particles.find(key);
            AMREX_ALWAYS_ASSERT(it != source_particles.end());
            auto const& src = it->second;
            int const expected_grid = use_source_grid ? src.grid_id : gid;

            AMREX_ALWAYS_ASSERT(p.idata(0) == expected_grid);
            AMREX_ALWAYS_ASSERT(static_cast<int>(h_rdata[i]) == expected_grid);
            AMREX_ALWAYS_ASSERT(h_idata[i] == expected_grid);

            IntVect shift = match_shift(p, src, grown_tile_box);
            actual_records.push_back({key, expected_grid, shift});
            ++actual_per_grid[expected_grid];
        }

        for (int igrid = 0; igrid < ngrids; ++igrid) {
            AMREX_ALWAYS_ASSERT(actual_per_grid[igrid] == expected_per_grid[igrid]);
        }

        std::sort(expected_records.begin(), expected_records.end());
        std::sort(actual_records.begin(), actual_records.end());
        AMREX_ALWAYS_ASSERT(actual_records == expected_records);
    }
}

void MDParticleContainer::checkNeighborList()
{
    BL_PROFILE("MDParticleContainer::checkNeighborList");

    const int lev = 0;
    auto& plev  = GetParticles(lev);

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        const int np       = aos.numParticles();
        const int np_total = aos.numTotalParticles();

        amrex::Vector<unsigned int> full_count(np,0);
        amrex::Vector<unsigned int> full_nbors;
        amrex::Gpu::HostVector<ParticleType> h_pstruct(np_total);
        Gpu::copy(Gpu::deviceToHost, aos().dataPtr(), aos().dataPtr() + np_total, h_pstruct.begin());

        // on the host, construct neighbor list using full N^2 search
        for (int i = 0; i < np; i++)
        {
            ParticleType& p1 = h_pstruct[i];

            // Loop over all particles
            for (int j = 0; j < np_total; j++)
            {
                // Don't be your own neighbor.
                if ( i == j ) { continue; }

                ParticleType& p2 = h_pstruct[j];
                AMREX_D_TERM(Real dx = p1.pos(0) - p2.pos(0);,
                             Real dy = p1.pos(1) - p2.pos(1);,
                             Real dz = p1.pos(2) - p2.pos(2);)

                Real r2 = AMREX_D_TERM(dx*dx, + dy*dy, + dz*dz);

                Real cutoff_sq = 25.0*Params::cutoff*Params::cutoff;

                if (r2 <= cutoff_sq)
                {
                    full_count[i] += 1;
                    full_nbors.push_back(j);
                }
            }
        }

        // copy the bin-constructed neighbor list to host
        auto& d_counts = m_neighbor_list[lev][index].GetCounts();
        Gpu::HostVector<unsigned int> h_counts(d_counts.size());
        Gpu::copy(Gpu::deviceToHost, d_counts.begin(), d_counts.end(), h_counts.begin());

        auto& d_list = m_neighbor_list[lev][index].GetList();
        Gpu::HostVector<unsigned int> h_list(d_list.size());
        Gpu::copy(Gpu::deviceToHost, d_list.begin(), d_list.end(), h_list.begin());

        // check answer
        for (int i = 0; i < np; ++i) {
            AMREX_ALWAYS_ASSERT(h_counts[i] == full_count[i]);
        }

        // order not the same, so sort here
        unsigned start = 0;
        for (int i = 0; i < np; ++i) {
            std::sort(full_nbors.data() + start, full_nbors.data() + start + full_count[i]);
            std::sort(h_list.data() + start, h_list.data() + start + full_count[i]);
            start += full_count[i];
        }

        for (unsigned i = 0; i < h_list.size(); ++i) {
            AMREX_ALWAYS_ASSERT(h_list[i] == full_nbors[i]);
        }
    }

    amrex::PrintToFile("neighbor_test") << "All the neighbor list particles match!" << '\n';
}

void MDParticleContainer::checkInverseSumNeighbors ()
{
    BL_PROFILE("MDParticleContainer::checkInverseSumNeighbors");

    constexpr ParticleReal struct_real_delta = 1.5_rt;
    constexpr int struct_int_delta = 9;
    constexpr ParticleReal array_real_delta = 2.5_rt;
    constexpr ParticleReal runtime_real_delta = 3.5_rt;
    constexpr int array_int_delta = 11;
    constexpr int runtime_int_delta = 13;

    const int lev = 0;
    auto& plev = GetParticles(lev);

    std::map<ParticleKey, OwnedParticleData> local_owned_particles;
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        auto const index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& ptile = plev[index];
        auto& aos = ptile.GetArrayOfStructs();
        auto& soa = ptile.GetStructOfArrays();
        int const np = aos.numParticles();
        if (np == 0) { continue; }

        Gpu::HostVector<ParticleType> h_pstruct(np);
        Gpu::HostVector<ParticleReal> h_array_real(np);
        Gpu::HostVector<ParticleReal> h_runtime_real(np);
        Gpu::HostVector<int> h_array_int(np);
        Gpu::HostVector<int> h_runtime_int(np);
        Gpu::copy(Gpu::deviceToHost, aos().dataPtr(), aos().dataPtr() + np, h_pstruct.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetRealData(0).begin(), soa.GetRealData(0).begin() + np,
                  h_array_real.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetRealData(1).begin(), soa.GetRealData(1).begin() + np,
                  h_runtime_real.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetIntData(0).begin(), soa.GetIntData(0).begin() + np,
                  h_array_int.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetIntData(1).begin(), soa.GetIntData(1).begin() + np,
                  h_runtime_int.begin());

        for (int i = 0; i < np; ++i) {
            auto const& p = h_pstruct[i];
            OwnedParticleData pdata;
            pdata.key = {p.id(), p.cpu()};
            AMREX_D_TERM(pdata.pos[0] = p.pos(0);,
                         pdata.pos[1] = p.pos(1);,
                         pdata.pos[2] = p.pos(2);)
            for (int comp = 0; comp < PIdx::ncomps; ++comp) {
                pdata.struct_real[comp] = p.rdata(comp);
            }
            pdata.struct_int = p.idata(0);
            pdata.array_real = h_array_real[i];
            pdata.runtime_real = h_runtime_real[i];
            pdata.array_int = h_array_int[i];
            pdata.runtime_int = h_runtime_int[i];
            auto [it, inserted] = local_owned_particles.emplace(pdata.key, pdata);
            amrex::ignore_unused(it);
            AMREX_ALWAYS_ASSERT(inserted);
        }
    }

    std::map<ParticleKey, InverseContributionData> local_contributions;
    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        auto& neighb_tile = GetNeighbors(lev, mfi.index(), mfi.LocalTileIndex());
        auto& aos = neighb_tile.GetArrayOfStructs();
        auto& soa = neighb_tile.GetStructOfArrays();
        int const nneighbs = static_cast<int>(aos.size());
        if (nneighbs == 0) { continue; }

        Gpu::HostVector<ParticleType> h_pstruct(nneighbs);
        Gpu::HostVector<ParticleReal> h_array_real(nneighbs);
        Gpu::HostVector<ParticleReal> h_runtime_real(nneighbs);
        Gpu::HostVector<int> h_array_int(nneighbs);
        Gpu::HostVector<int> h_runtime_int(nneighbs);
        Gpu::copy(Gpu::deviceToHost, aos().dataPtr(), aos().dataPtr() + nneighbs, h_pstruct.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetRealData(0).begin(), soa.GetRealData(0).begin() + nneighbs,
                  h_array_real.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetRealData(1).begin(), soa.GetRealData(1).begin() + nneighbs,
                  h_runtime_real.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetIntData(0).begin(), soa.GetIntData(0).begin() + nneighbs,
                  h_array_int.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetIntData(1).begin(), soa.GetIntData(1).begin() + nneighbs,
                  h_runtime_int.begin());

        for (int i = 0; i < nneighbs; ++i) {
            auto const& p = h_pstruct[i];
            ParticleKey key{p.id(), p.cpu()};
            auto& contribution = local_contributions[key];
            for (int comp = 0; comp < PIdx::ncomps; ++comp) {
                contribution.struct_real[comp] += p.rdata(comp) + struct_real_delta;
            }
            contribution.struct_int += p.idata(0) + struct_int_delta;
            contribution.array_real += h_array_real[i] + array_real_delta;
            contribution.runtime_real += h_runtime_real[i] + runtime_real_delta;
            contribution.array_int += h_array_int[i] + array_int_delta;
            contribution.runtime_int += h_runtime_int[i] + runtime_int_delta;
        }

        ParticleType* pstruct = aos().dataPtr();
        auto* array_real = soa.GetRealData(0).dataPtr();
        auto* runtime_real = soa.GetRealData(1).dataPtr();
        auto* array_int = soa.GetIntData(0).dataPtr();
        auto* runtime_int = soa.GetIntData(1).dataPtr();
        AMREX_FOR_1D(nneighbs, i,
        {
            ParticleType& p = pstruct[i];
            for (int comp = 0; comp < PIdx::ncomps; ++comp) {
                p.rdata(comp) += struct_real_delta;
            }
            p.idata(0) += struct_int_delta;
            array_real[i] += array_real_delta;
            runtime_real[i] += runtime_real_delta;
            array_int[i] += array_int_delta;
            runtime_int[i] += runtime_int_delta;
        });
    }

    Gpu::streamSynchronize();

    auto gather_contribution_data =
        [&] (std::map<ParticleKey, InverseContributionData> const& local_data)
    {
        std::vector<PackedContributionData> local_packed;
        local_packed.reserve(local_data.size());
        for (auto const& [key, contribution] : local_data) {
            local_packed.push_back(PackedContributionData{
                .id = key.first,
                .cpu = key.second,
                .struct_real = contribution.struct_real,
                .struct_int = contribution.struct_int,
                .array_real = contribution.array_real,
                .runtime_real = contribution.runtime_real,
                .array_int = contribution.array_int,
                .runtime_int = contribution.runtime_int
            });
        }

        auto unpack = [] (std::vector<PackedContributionData> const& packed) {
            std::map<ParticleKey, InverseContributionData> result;
            for (auto const& item : packed) {
                auto& contribution = result[{item.id, item.cpu}];
                for (int comp = 0; comp < PIdx::ncomps; ++comp) {
                    contribution.struct_real[comp] += item.struct_real[comp];
                }
                contribution.struct_int += item.struct_int;
                contribution.array_real += item.array_real;
                contribution.runtime_real += item.runtime_real;
                contribution.array_int += item.array_int;
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

            auto const global_data = unpack(gathered);
            global_packed.reserve(global_data.size());
            for (auto const& [key, contribution] : global_data) {
                global_packed.push_back(PackedContributionData{
                    .id = key.first,
                    .cpu = key.second,
                    .struct_real = contribution.struct_real,
                    .struct_int = contribution.struct_int,
                    .array_real = contribution.array_real,
                    .runtime_real = contribution.runtime_real,
                    .array_int = contribution.array_int,
                    .runtime_int = contribution.runtime_int
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
    };

    auto const global_contributions = gather_contribution_data(local_contributions);
    int total_contributions = 0;
    for (auto const& [key, contribution] : global_contributions) {
        amrex::ignore_unused(key);
        total_contributions += contribution.struct_int;
    }
    AMREX_ALWAYS_ASSERT(total_contributions > 0);

    sumNeighbors(0, PIdx::ncomps + NumRealComps(), 0, 1 + NumIntComps());

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        auto const index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
        auto& ptile = plev[index];
        auto& aos = ptile.GetArrayOfStructs();
        auto& soa = ptile.GetStructOfArrays();
        int const np = aos.numParticles();
        if (np == 0) { continue; }

        Gpu::HostVector<ParticleType> h_pstruct(np);
        Gpu::HostVector<ParticleReal> h_array_real(np);
        Gpu::HostVector<ParticleReal> h_runtime_real(np);
        Gpu::HostVector<int> h_array_int(np);
        Gpu::HostVector<int> h_runtime_int(np);
        Gpu::copy(Gpu::deviceToHost, aos().dataPtr(), aos().dataPtr() + np, h_pstruct.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetRealData(0).begin(), soa.GetRealData(0).begin() + np,
                  h_array_real.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetRealData(1).begin(), soa.GetRealData(1).begin() + np,
                  h_runtime_real.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetIntData(0).begin(), soa.GetIntData(0).begin() + np,
                  h_array_int.begin());
        Gpu::copy(Gpu::deviceToHost, soa.GetIntData(1).begin(), soa.GetIntData(1).begin() + np,
                  h_runtime_int.begin());

        for (int i = 0; i < np; ++i) {
            auto const& p = h_pstruct[i];
            ParticleKey key{p.id(), p.cpu()};
            auto const src_it = local_owned_particles.find(key);
            AMREX_ALWAYS_ASSERT(src_it != local_owned_particles.end());
            auto const& src = src_it->second;

            auto const contribution_it = global_contributions.find(key);
            InverseContributionData contribution;
            if (contribution_it != global_contributions.end()) {
                contribution = contribution_it->second;
            }

            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                AMREX_ALWAYS_ASSERT(almost_equal(p.pos(idim), src.pos[idim]));
            }
            for (int comp = 0; comp < PIdx::ncomps; ++comp) {
                AMREX_ALWAYS_ASSERT(almost_equal(
                    p.rdata(comp), src.struct_real[comp] + contribution.struct_real[comp]));
            }
            AMREX_ALWAYS_ASSERT(p.idata(0) == src.struct_int + contribution.struct_int);
            AMREX_ALWAYS_ASSERT(almost_equal(
                h_array_real[i], src.array_real + contribution.array_real));
            AMREX_ALWAYS_ASSERT(almost_equal(
                h_runtime_real[i], src.runtime_real + contribution.runtime_real));
            AMREX_ALWAYS_ASSERT(h_array_int[i] == src.array_int + contribution.array_int);
            AMREX_ALWAYS_ASSERT(h_runtime_int[i] == src.runtime_int + contribution.runtime_int);
        }
    }
}

void MDParticleContainer::reset_test_id()
{
    BL_PROFILE("MDParticleContainer::reset_test_id");

    const int lev = 0;
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numTotalParticles();

        ParticleType* pstruct = aos().dataPtr();
        auto* rdata = soa.GetRealData(0).dataPtr();
        auto* idata = soa.GetIntData(0).dataPtr();

        AMREX_FOR_1D ( np, i,
        {
            ParticleType& p1 = pstruct[i];
            p1.idata(0) = gid;
            rdata[i] = (ParticleReal) gid;
            idata[i] = gid;
        });
    }
}
