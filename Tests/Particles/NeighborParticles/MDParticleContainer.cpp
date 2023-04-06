#include "MDParticleContainer.H"
#include "Constants.H"

#include "CheckPair.H"
#include <AMReX_SPACE.H>

using namespace amrex;

namespace
{
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
                for (int i = 0; i < NumRealComps(); ++i)
                    host_real[i].push_back(ParticleReal(mfi.index()));
                for (int i = 0; i < NumIntComps(); ++i)
                    host_int[i].push_back(mfi.index());
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

    Real min_d = std::numeric_limits<Real>::max();
    Real max_d = std::numeric_limits<Real>::min();

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        Gpu::DeviceScalar<Real> min_d_gpu(min_d);
        Gpu::DeviceScalar<Real> max_d_gpu(max_d);

        Real* pmin_d = min_d_gpu.dataPtr();
        Real* pmax_d = max_d_gpu.dataPtr();

        AMREX_FOR_1D ( np, i,
        {
            ParticleType& p1 = pstruct[i];

            for (const auto& p2 : nbor_data.getNeighbors(i))
            {
                AMREX_D_TERM(Real dx = p1.pos(0) - p2.pos(0);,
                             Real dy = p1.pos(1) - p2.pos(1);,
                             Real dz = p1.pos(2) - p2.pos(2);)

                Real r2 = AMREX_D_TERM(dx*dx, + dy*dy, + dz*dz);
                r2 = amrex::max(r2, Params::min_r*Params::min_r);
                Real r = sqrt(r2);

                Gpu::Atomic::Min(pmin_d, r);
                Gpu::Atomic::Max(pmax_d, r);
            }
        });

        Gpu::Device::streamSynchronize();

        min_d = std::min(min_d, min_d_gpu.dataValue());
        max_d = std::max(max_d, max_d_gpu.dataValue());
    }
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
        ParticleType* pstruct = &(aos[0]);

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

void MDParticleContainer::checkNeighborParticles()
{
    BL_PROFILE("MDParticleContainer::checkNeighborParticles");

    const int lev = 0;
    auto& plev  = GetParticles(lev);

    auto ngrids = int(ParticleBoxArray(0).size());

    amrex::Gpu::DeviceVector<int> d_num_per_grid(ngrids,0);
    amrex::Gpu::HostVector<int> h_num_per_grid(ngrids);

    // CPU version
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        if (gid != 0) continue;
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        auto& soa   = ptile.GetStructOfArrays();
        const int np = aos.numTotalParticles();

        ParticleType* pstruct = aos().dataPtr();
        auto* rdata = soa.GetRealData(0).dataPtr();
        auto* idata = soa.GetIntData(0).dataPtr();

        int* p_num_per_grid = d_num_per_grid.data();
        amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            ParticleType& p1 = pstruct[i];
            Gpu::Atomic::AddNoRet(&(p_num_per_grid[p1.idata(0)]),1);
            AMREX_ALWAYS_ASSERT(p1.idata(0) == (int) rdata[i]);
            AMREX_ALWAYS_ASSERT(p1.idata(0) == idata[i]);
        });

        Gpu::copy(Gpu::deviceToHost,
                  d_num_per_grid.begin(), d_num_per_grid.end(), h_num_per_grid.begin());
        amrex::AllPrintToFile("neighbor_test") << "FOR GRID " << gid << "\n";;

        for (int i = 0; i < ngrids; i++) {
          amrex::AllPrintToFile("neighbor_test") << "   there are " << h_num_per_grid[i] << " with grid id " << i << "\n";
        }

        amrex::AllPrintToFile("neighbor_test") << " \n";
        amrex::AllPrintToFile("neighbor_test") << " \n";
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
                if ( i == j ) continue;

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

    amrex::PrintToFile("neighbor_test") << "All the neighbor list particles match!" << std::endl;
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
