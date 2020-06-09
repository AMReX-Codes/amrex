#include "MDParticleContainer.H"
#include "Constants.H"

#include "CheckPair.H"

using namespace amrex;

namespace
{    
    void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
    {
        int nx = nppc[0];
        int ny = nppc[1];
        int nz = nppc[2];
        
        int ix_part = i_part/(ny * nz);
        int iy_part = (i_part % (ny * nz)) % ny;
        int iz_part = (i_part % (ny * nz)) / ny;
        
        r[0] = (0.5+ix_part)/nx;
        r[1] = (0.5+iy_part)/ny;
        r[2] = (0.5+iz_part)/nz;
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
              const Real     a_thermal_momentum_std,
              const Real     a_thermal_momentum_mean)
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
        
        for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv)) {
            for (int i_part=0; i_part<num_ppc;i_part++) {
                Real r[3];
                Real v[3];
                
                get_position_unit_cell(r, a_num_particles_per_cell, i_part);
                
                get_gaussian_random_momentum(v, a_thermal_momentum_mean,
                                             a_thermal_momentum_std);
                
                Real x = plo[0] + (iv[0] + r[0])*dx[0];
                Real y = plo[1] + (iv[1] + r[1])*dx[1];
                Real z = plo[2] + (iv[2] + r[2])*dx[2];
                
                ParticleType p;
                p.id()  = ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();                
                p.pos(0) = x;
                p.pos(1) = y;
                p.pos(2) = z;
                
                p.rdata(PIdx::vx) = v[0];
                p.rdata(PIdx::vy) = v[1];
                p.rdata(PIdx::vz) = v[2];

                p.rdata(PIdx::ax) = 0.0;
                p.rdata(PIdx::ay) = 0.0;
                p.rdata(PIdx::az) = 0.0;

                p.idata(0) = mfi.index();
                
                host_particles.push_back(p);
            }
        }
        
        auto& particles = GetParticles(lev);
        auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);
        
        Gpu::copy(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
                  particle_tile.GetArrayOfStructs().begin() + old_size);        
    }
    
    amrex::PrintToFile("neighbor_test") << " Number of particles is " << this->TotalNumberOfParticles()<< " \n";
    amrex::PrintToFile("neighbor_test") << "done. \n";
}

std::pair<Real, Real> MDParticleContainer::minAndMaxDistance()
{
    BL_PROFILE("MDParticleContainer::minAndMaxDistance");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
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
                Real dx = p1.pos(0) - p2.pos(0);
                Real dy = p1.pos(1) - p2.pos(1);
                Real dz = p1.pos(2) - p2.pos(2);

                Real r2 = dx*dx + dy*dy + dz*dz;
                r2 = amrex::max(r2, Params::min_r*Params::min_r);
                Real r = sqrt(r2);

		Gpu::Atomic::Min(pmin_d, r);
		Gpu::Atomic::Max(pmax_d, r);
            }
        });

        //	Gpu::Device::streamSynchronize();

	min_d = std::min(min_d, min_d_gpu.dataValue());
	max_d = std::max(max_d, max_d_gpu.dataValue());
    }
    ParallelDescriptor::ReduceRealMin(min_d, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(max_d, ParallelDescriptor::IOProcessorNumber());

    return std::make_pair(min_d, max_d);
}

void MDParticleContainer::moveParticles(amrex::Real dx)
{
    BL_PROFILE("MDParticleContainer::moveParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
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
            p.pos(0) += dx;
            p.pos(1) += dx;
            p.pos(2) += dx;
        });
    }
}

void MDParticleContainer::writeParticles(const int n)
{
    BL_PROFILE("MDParticleContainer::writeParticles");
    const std::string& pltfile = amrex::Concatenate("particles", n, 5);
    WriteAsciiFile(pltfile);
}

void MDParticleContainer::checkNeighborParticles()
{
    BL_PROFILE("MDParticleContainer::checkNeighborParticles");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);

    int ngrids = ParticleBoxArray(0).size();

    amrex::Gpu::ManagedVector<int> d_num_per_grid(ngrids,0);
    int* p_num_per_grid = d_num_per_grid.data();

    // CPU version
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        if (gid != 0) continue;
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numTotalParticles();

        ParticleType* pstruct = aos().dataPtr();

        for (int i = 0; i < np; i++)
        {
            ParticleType& p1 = pstruct[i];
            Gpu::Atomic::Add(&(p_num_per_grid[p1.idata(0)]),1);
        }

        amrex::AllPrintToFile("neighbor_test") << "FOR GRID " << gid << "\n";;

        for (int i = 0; i < ngrids; i++)
          amrex::AllPrintToFile("neighbor_test") << "   there are " << d_num_per_grid[i] << " with grid id " << i << "\n";;

        amrex::AllPrintToFile("neighbor_test") << " \n";
        amrex::AllPrintToFile("neighbor_test") << " \n";
    }

#if 0
    // GPU version
    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        // int mine = 0;

        // amrex::Gpu::DeviceScalar<int> d_mine(mine);
        // int* p_mine = d_mine.dataPtr();

        if (gid != 0) continue;
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numTotalParticles();

        ParticleType* pstruct = aos().dataPtr();

        // now we loop over the neighbor list and compute the forces
        AMREX_FOR_1D ( np, i,
        {
            ParticleType& p1 = pstruct[i];

            // Gpu::Atomic::Add(p_mine,1);

            Gpu::Atomic::Add(&(p_num_per_grid[p1.idata(0)]),1);
        });

        Gpu::Device::synchronize();

        // mine = d_mine.dataValue();

        amrex::AllPrintToFile("neighbor_test") << "FOR GRID " << gid << "\n";;

        for (int i = 0; i < ngrids; i++)
          amrex::AllPrintToFile("neighbor_test") << "   there are " << d_num_per_grid[i] << " with grid id " << i << "\n";;

        amrex::AllPrintToFile("neighbor_test") << " \n";
        amrex::AllPrintToFile("neighbor_test") << " \n";
    }
#endif
}

void MDParticleContainer::checkNeighborList()
{
    BL_PROFILE("MDParticleContainer::checkNeighborList");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);

    int ngrids = ParticleBoxArray(0).size();

    for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        int mine = 0;

        amrex::Gpu::DeviceScalar<int> d_mine(mine);
        int* p_mine = d_mine.dataPtr();

        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        const size_t np       = aos.numParticles();
        const size_t np_total = aos.numTotalParticles();

        amrex::Gpu::ManagedVector<int> d_neighbor_count(np,0);
        int* p_neighbor_count = d_neighbor_count.data();

        amrex::Gpu::ManagedVector<int> d_full_count(np,0);
        int* p_full_count = d_full_count.data();

        auto nbor_data = m_neighbor_list[lev][index].data();
        ParticleType* pstruct = aos().dataPtr();

        // ON DEVIDE: 
        // AMREX_FOR_1D ( np, i,
        // ON HOST:
        // for (int i = 0; i < np; i++)
        for (int i = 0; i < np; i++)
        {
            ParticleType& p1 = pstruct[i];

            amrex::Vector<int> nbor_nbors;
            amrex::Vector<int> full_nbors;

            // Loop over all particles
            for (int j = 0; j < np_total; j++)
            {
                // Don't be your own neighbor.
                if ( i == j ) continue;

                ParticleType& p2 = pstruct[j];
                Real dx = p1.pos(0) - p2.pos(0);
                Real dy = p1.pos(1) - p2.pos(1);
                Real dz = p1.pos(2) - p2.pos(2);
                
                Real r2 = dx*dx + dy*dy + dz*dz;

                Real cutoff_sq = 25.0*Params::cutoff*Params::cutoff; 

		if (r2 <= cutoff_sq)
		{
                   Gpu::Atomic::Add(&(p_full_count[i]),1);
                   full_nbors.push_back(p2.id());
		}
            }

            for (const auto& p2 : nbor_data.getNeighbors(i))
            {               
                Gpu::Atomic::Add(&(p_neighbor_count[i]),1);
                nbor_nbors.push_back(p2.id());
            }

            std::sort(full_nbors.begin(), full_nbors.end());
            std::sort(nbor_nbors.begin(), nbor_nbors.end());

            if (nbor_nbors.size() != full_nbors.size())
            {
               amrex::PrintToFile("neighbor_test") << "Number of neighbors do not match for particle " << i << std::endl;
               amrex::PrintToFile("neighbor_test") << "Neighbor list has " << nbor_nbors.size() << " particles " << std::endl;
               amrex::PrintToFile("neighbor_test") << "Full N^2 list has " << full_nbors.size() << " particles " << std::endl;
               amrex::Abort();
            }

            // amrex::PrintToFile("neighbor_test") << "   there are " << nbor_nbors.size() << " " <<
            //                  full_nbors.size() << " list / full neighbors of particle " << i << std::endl;
            
            // Loop over particles in my neighbor list
            for (int cnt = 0; cnt < nbor_nbors.size(); cnt++)
            {               
                // std::cout << "   NBORS " << nbor_nbors[cnt] << " " << full_nbors[cnt] << std::endl;
                if (nbor_nbors[cnt] != full_nbors[cnt])
                {
                     amrex::PrintToFile("neighbor_test") << "Index of neighbors do not match for particle " << i << std::endl;
                     amrex::PrintToFile("neighbor_test") << "Neighbor list neighbor index: " << nbor_nbors[cnt]  << std::endl;
                     amrex::PrintToFile("neighbor_test") << "Full N^2 list neighbor index: " << full_nbors[cnt]  << std::endl;
                     amrex::Abort();
                }
            }

        } // i
    } // MFIter

    amrex::PrintToFile("neighbor_test") << "All the neighbor list particles match!" << std::endl;
}

void MDParticleContainer::reset_test_id()
{
    BL_PROFILE("MDParticleContainer::reset_test_id");

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev  = GetParticles(lev);

    for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();

        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();
        const size_t np = aos.numTotalParticles();

        ParticleType* pstruct = aos().dataPtr();

        AMREX_FOR_1D ( np, i,
        {
            ParticleType& p1 = pstruct[i];
            p1.idata(0) = gid;
        });
    }
}
