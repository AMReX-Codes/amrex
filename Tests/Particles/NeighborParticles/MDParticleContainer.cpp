#include <AMReX_Reduce.H>
#include <AMReX_SPACE.H>

#include "MDParticleContainer.H"
#include "Constants.H"
#include "CheckPair.H"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <tuple>

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
    auto const plo = geom.ProbLoArray();
    auto const dxi = geom.InvCellSizeArray();
    auto const domain = geom.Domain();
    auto const problo = geom.ProbLoArray();
    auto const probhi = geom.ProbHiArray();
    const auto& pshifts = geom.periodicity().shiftIntVect();
    auto& plev  = GetParticles(lev);

    auto ngrids = int(ParticleBoxArray(0).size());
    std::map<ParticleKey, SourceParticleData> source_particles;

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
            ParticleKey key{p.id(), p.cpu()};
            SourceParticleData pdata;
            pdata.key = key;
            pdata.cell = getParticleCell(p, plo, dxi, domain);
            pdata.grid_id = p.idata(0);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                pdata.pos[idim] = p.pos(idim);
            }
            auto [it, inserted] = source_particles.emplace(key, pdata);
            amrex::ignore_unused(it);
            AMREX_ALWAYS_ASSERT(inserted);
        }
    }

    auto match_shift = [&] (ParticleType const& p, SourceParticleData const& src,
                            Box const& grown_tile_box) -> IntVect
    {
        constexpr auto tol = ParticleReal(1.0e-12);
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
                if (std::abs(p.pos(idim) - expected_pos) > tol) {
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
