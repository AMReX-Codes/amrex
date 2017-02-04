template <int NR, int NI>
bool    ParticleContainer<NR, NI>::do_tiling = false;

template <int NR, int NI>
IntVect ParticleContainer<NR, NI>::tile_size   { D_DECL(1024000,8,8) };

template <int NR, int NI>
void
ParticleContainer<NR, NI>::allocateData() {
  communicate_comp.resize(m_npartdata, true); 
  int nlevs = m_gdb->finestLevel()+1;
  m_particles.resize(nlevs);
  m_partdata.resize(nlevs);
  m_dummy_mfs.resize(nlevs);
  for (int lev = 0; lev < nlevs; ++ lev)
    {
      auto& partleveldata = m_partdata[lev];
      const BoxArray& ba = m_gdb->ParticleBoxArray(lev);
      const DistributionMapping& dm = m_gdb->ParticleDistributionMap(lev);

      MFInfo Fab_noallocate;
      Fab_noallocate.SetAlloc(false);
      m_dummy_mfs[lev] = std::unique_ptr<MultiFab>(new MultiFab(ba, dm, 1, 0, Fab_noallocate));

      // This is maybe stupid
      // First, allocate for one tile per grid so we can use PIter
      const MultiFab& foo = *(m_dummy_mfs[lev]);
      for (MFIter mfi(foo); mfi.isValid(); ++mfi)
	{
	  int i = mfi.index();
	  partleveldata[i].resize(1);
	  m_particles[lev][i].resize(1);
	}

      // Now resize for the appropriate number of tiles
      PIter it(*this, lev);
      for (MFIter mfi(foo); mfi.isValid(); ++mfi)
	{
	  int i = mfi.index();
	  int num_tiles = it.numTiles();
	  partleveldata[i].resize(num_tiles);
	  m_particles[lev][i].resize(num_tiles);
	}
  
      // Allocate the SoA data on each tile;
      for (PIter it(*this, lev); it.isValid(); ++it) 
	{
	  int gid = it.index();
	  int tid = it.tileIndex();
	  partleveldata[gid][tid].resize(m_npartdata);
	}
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::locateParticle(ParticleType& p, ParticleLocData& pld, 
					  int lev_min, int theEffectiveFinestLevel,
					  int nGrow) {
  if (!ParticleType::Where(p, m_gdb, pld, lev_min, theEffectiveFinestLevel)) {
    if (!ParticleType::PeriodicWhere(p, m_gdb, pld, lev_min, theEffectiveFinestLevel)) {
      if (lev_min != 0) {
	if (!ParticleType::RestrictedWhere(p, m_gdb, pld, nGrow))
	  amrex::Abort("ParticleContainer<NR, NI>::Redistribute(): invalid particle at non-coarse step");
      }
      else {
	//
	// The particle has left the domain; invalidate it.
	// This typically only happens on a coarse timestep.
	//
	p.m_idata.id = -p.m_idata.id;
      }
    }
  }
}

template <int NR, int NI>
long
ParticleContainer<NR, NI>::TotalNumberOfParticles (bool only_valid, bool only_local) const
{
    long nparticles = 0;
    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++) {
        nparticles += NumberOfParticlesAtLevel(lev,only_valid,false);
    }
    if (!only_local) {
	ParallelDescriptor::ReduceLongSum(nparticles);
    }
    return nparticles;
}

template <int NR, int NI>
Array<long>
ParticleContainer<NR, NI>::NumberOfParticlesInGrid (int  lev, 
						    bool only_valid,
						    bool only_local) const
{
    int ngrids = m_gdb->ParticleBoxArray(lev).size();
    Array<long> nparticles(ngrids,0);

    if (lev >= 0 && lev < int(m_particles.size())) {
      for (PIter it(*this, lev); it.isValid(); ++it) {
	int gid = it.index();
	int tid = it.tileIndex();
        const AoS& pbox = m_particles[lev][gid][tid];
	
	if (only_valid) {
	  for (const auto& p : pbox) {
	    if (p.m_idata.id > 0) ++nparticles[gid];
	  }
	}
	else {
	  nparticles[gid] = pbox.size();
	}
      }
      
      if (!only_local)
	ParallelDescriptor::ReduceLongSum(&nparticles[0],ngrids);
    }
    return nparticles;
}

template <int NR, int NI>
long
ParticleContainer<NR, NI>::NumberOfParticlesAtLevel (int  lev,
						     bool only_valid,
						     bool only_local) const
{
    long nparticles = 0;

    if (lev >= 0 && lev < int(m_particles.size())) {
      for (PIter it(*this, lev); it.isValid(); ++it) {
	int gid = it.index();
	int tid = it.tileIndex();
        const AoS& pbox = m_particles[lev][gid][tid];
	
	if (only_valid) {
	  for (const auto& p : pbox) {
	    if (p.m_idata.id > 0) ++nparticles;
	  }
	}
	else {
	  nparticles += pbox.size();
	}
      }
    }
    
    if (!only_local)
      ParallelDescriptor::ReduceLongSum(nparticles);
    
    return nparticles;
}

//
// This includes both valid and invalid particles since invalid particles still take up space.
//

template <int NR, int NI>
void
ParticleContainer<NR, NI>::ByteSpread () const
{
    long cnt = 0;

    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
      const auto& pmap = m_particles[lev];
      for (PIter it(*this, lev); it.isValid(); ++it) {
	int gid = it.index();
	int tid = it.tileIndex();
	cnt += pmap.at(gid)[tid].size();
      }
    }

    long mn = cnt, mx = mn;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    const std::size_t sz = sizeof(ParticleType);

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
    ParallelDescriptor::ReduceLongMin(mn, IOProc);
    ParallelDescriptor::ReduceLongMax(mx, IOProc);
    ParallelDescriptor::ReduceLongSum(cnt,IOProc);

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << "ParticleContainer<NR, NI> byte spread across MPI nodes: ["
                  << mn*sz
		  << " (" << mn << ")"
                  << " ... "
                  << mx*sz
		  << " (" << mx << ")"
                  << "] total particles: (" << cnt << ")\n";
    }
#ifdef BL_LAZY
    });
#endif
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::addOneParticle (int                  id_in,
					   int                  cpu_in, 
					   std::vector<double>& xloc, 
					   std::vector<double>& attributes)
{
    BL_PROFILE("ParticleContainer<NR, NI>::addOneParticle()");
    if (m_particles.size() == 0)
    {
       m_particles.resize(m_gdb->finestLevel()+1);
    }

    ParticleType p;

    p.m_idata.id  = id_in;
    p.m_idata.cpu = cpu_in;

    if (p.m_idata.id <= 0)
        amrex::Abort("Particle ID's must be > 0 in addOneParticle");
 
    if (ParallelDescriptor::MyProc() != p.m_idata.cpu)
        amrex::Abort("cpu_in must equal MyProc() in addOneParticle");

    for (int i = 0; i < BL_SPACEDIM; i++)
       p.m_rdata.pos[i] = xloc[i];

    for (int i = BL_SPACEDIM; i < BL_SPACEDIM + NR; i++)
       p.m_rdata.arr[i] = attributes[i];

    ParticleLocData pld;
    if (!ParticleType::Where(p, m_gdb, pld))
    {
        ParticleType::PeriodicShift(p, m_gdb);

        if (!ParticleType::Where(p, m_gdb, pld))
        {
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                std::cout << "BAD PARTICLE POS(" << d << ") " << p.m_rdata.pos[d] << std::endl;
            }
            amrex::Abort("ParticleContainer<NR, NI>::addOneParticle(): invalid particle");
        }
    }

    m_particles[pld.m_lev][pld.m_grid][pld.m_tile].push_back(p);

    // Note that we will need to call Redistribute once we are done adding particles this way.
    // The Where call above assigns a particle to a grid (pld.m_grid) based on the particle location.
    // However, the particle may not currently live on the processor that owns that grid.
    // The Redistribute routine should ensure that the particle ends up on the right processor.
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::addNParticles(int n_part, Array<double>& x, Array<double>& y, 
#if (BL_SPACEDIM == 3)
                                                      Array<double>& z, 
#endif
                                          int n_attr, Array<double>& attributes)
{
    BL_PROFILE("ParticleContainer<NR, NI>::addNParticles()");

    BL_ASSERT(OK());
    BL_ASSERT(NR >= BL_SPACEDIM+n_attr);

    if (m_particles.size() == 0)
       m_particles.resize(m_gdb->finestLevel()+1);

    ParticleType p;

    for (long j = 0; j < n_part; j++)
    {
       ParticleType p;

       p.m_rdata.pos[j  ] = x[j]; // x
       p.m_rdata.pos[j+1] = y[j]; // y
#if (BL_SPACEDIM == 3)
       p.m_rdata.pos[j+2] = z[j]; // z
#endif

       for (int i = 0; i < n_attr; i++)
          p.m_rdata.arr[BL_SPACEDIM+i] = attributes[n_attr*j + i];

       ParticleLocData pld;
       if (!ParticleType::Where(p, m_gdb, pld))
              amrex::Abort("ParticleContainer<NR, NI>::addNParticles(): invalid particle location");

       p.m_idata.id  = ParticleType::NextID();
       p.m_idata.cpu = ParallelDescriptor::MyProc();

       m_particles[pld.m_lev][pld.m_grid][pld.m_tile].push_back(p);
    }

    Redistribute();
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::MoveRandom ()
{
    //
    // Move particles randomly at all levels
    //
    for (int lev = 0; lev < int(m_particles.size()); lev++)
    {
       MoveRandom(lev);
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::MoveRandom (int lev)
{
    BL_PROFILE("ParticleContainer<NR, NI>::MoveRandom(lev)");
    BL_ASSERT(OK());
    BL_ASSERT(m_gdb != 0);
    // 
    // Move particles up to FRAC*CellSize distance in each coordinate direction.
    //
    const Real FRAC = 0.25;

    static bool first = true;

    static Array<amrex::mt19937> rn;

    if (first)
    {
        first = false;
        //
        // Build and initialize a random number generator per thread.
        //
        int tnum = 1;

#ifdef _OPENMP
        tnum = omp_get_max_threads();
#endif
        rn.resize(tnum);

        for (int i = 0; i < tnum; i++)
        {
            //
            // We want to give each thread across all MPI processes a unique non-zero seed.
            //
            const unsigned long seedbase = 1+tnum*ParallelDescriptor::MyProc();

            rn[i] = amrex::mt19937(seedbase+i);
        }
    }

    AoSMap&       pmap            = m_particles[lev];
    const Real* dx                = m_gdb->Geom(lev).CellSize();
    const Real  dist[BL_SPACEDIM] = { D_DECL(FRAC*dx[0], FRAC*dx[1], FRAC*dx[2]) };

    for (PIter it(*this, lev); it.isValid(); ++it) {
      int gid = it.index();
      int tid = it.tileIndex();
      AoS& pbox = pmap[gid][tid];
      const int n    = pbox.size();
	
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < n; i++)
        {
	  ParticleType& p = pbox[i];
	  
	  if (p.m_idata.id <= 0) continue;
	  
#ifdef _OPENMP
	  int tid = omp_get_thread_num();
#else
	  int tid = 0;
#endif
	  for (int i = 0; i < BL_SPACEDIM; i++)
            {
	      p.m_rdata.pos[i] += dist[i]*(2*rn[tid].d_value()-1);
            }
	  
	  ParticleType::Reset(p, m_gdb, true);
        }
    }
    Redistribute();
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::Increment (MultiFab& mf,
				      int       lev) 
{
  IncrementWithTotal(mf,lev);
}

template <int NR, int NI>
long
ParticleContainer<NR, NI>::IncrementWithTotal (MultiFab& mf,
					       int       lev)
{
  BL_PROFILE("ParticleContainer<NR, NI>::IncrementWithTotal(lev)");
  BL_ASSERT(OK());
  
  if (m_particles.empty()) return 0;
  
  BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));
  
  const AoSMap& pmap = m_particles[lev];
  
  long num_particles_in_domain = 0;
  
  MultiFab* mf_pointer;
  
  if (OnSameGrids(lev, mf))
    {
      // If we are already working with the internal mf defined on the
      // particle_box_array, then we just work with this.
      mf_pointer = &mf;
    }
  else
    {
      // If mf is not defined on the particle_box_array, then we need
      // to make a temporary mf_pointer here and copy it into mf at the end.
      mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),
				m_gdb->ParticleDistributionMap(lev),
				mf.nComp(),mf.nGrow());
    }
  
  for (PIter it(*this, lev); it.isValid(); ++it) {
    int gid = it.index();
    int tid = it.tileIndex();
    const AoS& pbox = pmap[gid][tid];
    FArrayBox&  fab  = (*mf_pointer)[gid];
    for (const auto& p : pbox) {
	ParticleLocData pld;
	Particle<NR, NI>::Where(p, m_gdb, pld);
	if (p.m_idata.id > 0) {
	  BL_ASSERT(pld.m_grid == gid);
	  fab(pld.m_cell) += 1;
	  num_particles_in_domain += 1;
	}
    }
  }
  
  // If mf is not defined on the particle_box_array, then we need
  // to copy here from mf_pointer into mf.   I believe that we don't
  // need any information in ghost cells so we don't copy those.
  if (mf_pointer != &mf) 
    {
      mf.copy(*mf_pointer,0,0,mf.nComp());  
      delete mf_pointer;
    }
  
  ParallelDescriptor::ReduceLongSum(num_particles_in_domain);
  
  return num_particles_in_domain;
}

template <int NR, int NI>
Real
ParticleContainer<NR, NI>::sumParticleMass (int rho_index, int lev) const
{
  BL_PROFILE("ParticleContainer<NR, NI>::sumParticleMass(lev)");
  BL_ASSERT(NR >= 1);
  BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));
  
  Real msum = 0;
  
  const AoSMap& pmap = m_particles[lev];
  for (PIter it(*this, lev); it.isValid(); ++it) {
    int gid = it.index();
    int tid = it.tileIndex();
    const AoS& pbox = pmap.at(gid)[tid];
    for (const auto& p : pbox) {
      if (p.m_idata.id > 0) {
	msum += p.m_rdata.arr[BL_SPACEDIM+rho_index];
      }
    }
  }
  
  ParallelDescriptor::ReduceRealSum(msum);
  
  return msum;
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::GetParticleIDs (Array<int>& part_ids)
{
  BL_PROFILE("ParticleContainer<NR, NI>::GetParticleIDs()");
  //
  // This gives us the starting point into the part_ids array
  // If only one processor (or no MPI), then that's all we need.
  //
  int cnt = 0;
  
#if BL_USE_MPI
  Array<long> cnts(ParallelDescriptor::NProcs());
  
  // This returns the number of particles on this processor
  long lcnt = TotalNumberOfParticles(true,true);
  
  // This accumulates the "lcnt" values into "cnts"
  MPI_Gather(&lcnt,1,              
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     cnts.dataPtr(),
	     1,
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     ParallelDescriptor::IOProcessorNumber(),
	     ParallelDescriptor::Communicator());
  
  ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());
  
  for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++) {
    cnt += cnts[iproc];
  }
  
  std::cout << "PROC CNT " << ParallelDescriptor::MyProc() << " " << cnt << std::endl;
  
  // Each particle takes up 1 int so no need to multiply cnt by anything
#endif

  // This is the total number of particles on *all* processors
  long npart = TotalNumberOfParticles(true,false);

  // Locations
  part_ids.resize(npart, 0);
  
  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
    const auto& pmap = m_particles[lev];
    for (PIter it(*this, lev); it.isValid(); ++it) {
      int gid = it.index();
      int tid = it.tileIndex();
      const AoS& pbx = pmap[gid][tid];
      for (const auto& p : pbx) {
	if (p.m_idata.id > 0) {
	  // Load the ID
	  part_ids[cnt] = p.m_idata.id;
	  
	  // Update counter
	  cnt++;
	}
      }
    }
  }
  
  ParallelDescriptor::ReduceIntSum(part_ids.dataPtr(),part_ids.size()); 
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::GetParticleCPU (Array<int>& part_cpu)
{
  BL_PROFILE("ParticleContainer<NR, NI>::GetParticleCPU()");
  //
  // This gives us the starting point into the part_cpu array
  // If only one processor (or no MPI), then that's all we need.
  //
  int cnt = 0;
  
#if BL_USE_MPI
  Array<long> cnts(ParallelDescriptor::NProcs());
  
  // This returns the number of particles on this processor
  long lcnt = TotalNumberOfParticles(true,true);
  
  // This accumulates the "lcnt" values into "cnts"
  MPI_Gather(&lcnt,1,              
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     cnts.dataPtr(),
	     1,
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     ParallelDescriptor::IOProcessorNumber(),
	     ParallelDescriptor::Communicator());
  
  ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());
  
  for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
    cnt += cnts[iproc];

  // Each particle takes up 1 int so no need to multiply cnt by anything
#endif
  
  // This is the total number of particles on *all* processors
  long npart = TotalNumberOfParticles(true,false);
  
  // Locations
  part_cpu.resize(npart,0);
  
  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
    const auto& pmap = m_particles[lev];
    for (PIter it(*this, lev); it.isValid(); ++it) {
      int gid = it.index();
      int tid = it.tileIndex();
      const AoS& pbx = pmap[gid][tid];
      for (const auto& p : pbx) {
	if (p.m_idata.id > 0) {
	  // Load the ID
	  part_cpu[cnt] = p.m_idata.cpu;
	  
	  // Update counter
	  cnt++;
	}
      }
    }
  }
  
  ParallelDescriptor::ReduceIntSum(part_cpu.dataPtr(),part_cpu.size()); 
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::GetParticleLocations (Array<Real>& part_data)
{
  BL_PROFILE("ParticleContainer<NR, NI>::GetParticleLocations()");
  //
  // This gives us the starting point into the part_data array
  // If only one processor (or no MPI), then that's all we need.
  //
  int cnt = 0;

#if BL_USE_MPI
  Array<long> cnts(ParallelDescriptor::NProcs());
  
  // This returns the number of particles on this processor
  long lcnt = TotalNumberOfParticles(true,true);
  
  // This accumulates the "lcnt" values into "cnts"
  MPI_Gather(&lcnt,1,              
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     cnts.dataPtr(),
	     1,
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     ParallelDescriptor::IOProcessorNumber(),
	     ParallelDescriptor::Communicator());
  
  ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

  for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
    cnt += cnts[iproc];
  
  // Each particle takes up BL_SPACEDIM Reals
  cnt *= (BL_SPACEDIM);
#endif

  // This is the total number of particles on *all* processors
  long npart = TotalNumberOfParticles(true, false);
  
  // Locations
  part_data.resize(BL_SPACEDIM*npart,0);

  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
    const auto& pmap = m_particles[lev];
    for (PIter it(*this, lev); it.isValid(); ++it) {
      int gid = it.index();
      int tid = it.tileIndex();
      const AoS& pbx = pmap[gid][tid];
      for (const auto& p : pbx) {
	if (p.m_idata.id > 0) {
	  // Load positions
	  for (int d=0; d < BL_SPACEDIM; d++)
	    part_data[cnt+d] = p.m_rdata.pos[d];
	  
	  // Update counter
	  cnt += BL_SPACEDIM;
	}
      }
    }
  }
  
  ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}


template <int NR, int NI>
void
ParticleContainer<NR, NI>::GetParticleData (Array<Real>& part_data, int start_comp, int num_comp)
{
  BL_PROFILE("ParticleContainer<NR, NI>::GetParticleData()");
  //
  // This gives us the starting point into the part_data array
  // If only one processor (or no MPI), then that's all we need.
  //
  int cnt = 0;
  
  //
  // Make sure we don't try to get more than we have
  //
  if (start_comp + num_comp > BL_SPACEDIM + NR)
    amrex::Error("Tried to grab too many components in GetParticleData!!");
  
#if BL_USE_MPI
  Array<long> cnts(ParallelDescriptor::NProcs());
  
  // This returns the number of particles on this processor
  long lcnt = TotalNumberOfParticles(true,true);
  
  // This accumulates the "lcnt" values into "cnts"
  MPI_Gather(&lcnt,1,              
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     cnts.dataPtr(),
	     1,
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     ParallelDescriptor::IOProcessorNumber(),
	     ParallelDescriptor::Communicator());
  
  ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());
  
  for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
    cnt += cnts[iproc];
  
  // Each particle takes up num_comp Reals
  cnt*= num_comp;
#endif

  // This is the total number of particles on *all* processors
  long npart = TotalNumberOfParticles(true,false);

  part_data.resize(num_comp*npart,0);

  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
    const auto& pmap = m_particles[lev];
    for (PIter it(*this, lev); it.isValid(); ++it) {
      int gid = it.index();
      int tid = it.tileIndex();
      const AoS& pbx = pmap[gid][tid];
      for (const auto& p : pbx) {
	if (p.m_idata.id > 0) {
	  // Load particle data, whatever it is.
	  for (int d = 0; d < num_comp; d++)
	    part_data[cnt+d] = p.m_rdata.arr[BL_SPACEDIM + start_comp + d];
	  
	  // Update counter
	  cnt += num_comp;
	}
      }
    }
  }
  
  ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::SetAllowParticlesNearBoundary (bool value)
{
  allow_particles_near_boundary = value; 
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::SetParticleLocations (Array<Real>& part_data)
{
  BL_PROFILE("ParticleContainer<NR, NI>::SetParticleLocations()");
  // This gives us the starting point into the part_data array
  // If only one processor (or no MPI), then that's all we need
  int cnt = 0;

#if BL_USE_MPI
  Array<long> cnts(ParallelDescriptor::NProcs());
  
  // This returns the number of particles on this processor
  long lcnt = TotalNumberOfParticles(true,true);
  
  // This accumulates the "lcnt" values into "cnts"
  MPI_Gather(&lcnt,1,              
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     cnts.dataPtr(),
	     1,
	     ParallelDescriptor::Mpi_typemap<long>::type(),
	     ParallelDescriptor::IOProcessorNumber(),
	     ParallelDescriptor::Communicator());

  ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());
  
  for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); iproc++)
    cnt += cnts[iproc];
  
  // Each particle takes up BL_SPACEDIM Reals
  cnt*= BL_SPACEDIM;
#endif
  
  // This is the total number of particles on *all* processors
  long npart = TotalNumberOfParticles(true,false);
  
  // Mass + locations
  if (part_data.size() != npart*BL_SPACEDIM)
    amrex::Abort("Sending in wrong size part_data to SetParticleLocations");

  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
    const auto& pmap = m_particles[lev];
    for (PIter it(*this, lev); it.isValid(); ++it) {
      int gid = it.index();
      int tid = it.tileIndex();
      const AoS& pbx = pmap[gid][tid];
      for (const auto& p : pbx) {
	if (p.m_idata.id > 0) {
	  // Load positions
	  for (int d=0; d < BL_SPACEDIM; d++)
	    p.m_rdata.pos[d] = part_data[cnt+d];
	  
	  // Update counter
	  cnt += BL_SPACEDIM;
	}
      }
    }
  }
}


template <int NR, int NI>
void
ParticleContainer<NR, NI>::RemoveParticlesAtLevel (int level)
{
  BL_PROFILE("ParticleContainer<NR, NI>::RemoveParticlesAtLevel()");
  if (level >= int(this->m_particles.size()))
    return;
  
  if (!this->m_particles[level].empty())
    {
      AoSMap().swap(this->m_particles[level]);
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::RemoveParticlesNotAtFinestLevel ()
{
  BL_PROFILE("ParticleContainer<NR, NI>::RemoveParticlesNotAtFinestLevel()");
  BL_ASSERT(this->m_gdb->finestLevel()+1 == int(this->m_particles.size()));
  
  int cnt = 0;
  
  for (unsigned lev = 0; lev < m_particles.size() - 1; lev++) {
    const auto& pmap = m_particles[lev];
    if (!pmap.empty()) {
      for (PIter it(*this, lev); it.isValid(); ++it) {
	int gid = it.index();
	int tid = it.tileIndex();
	const AoS& pbx = pmap[gid][tid];
	cnt += pbx.size();
      }
      AoSMap().swap(pmap);
    }
  }
  
  //
  // Print how many particles removed on each processor if any were removed.
  //
  if (this->m_verbose > 1)
    {
      int maxcnt = cnt;
      
#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceIntMax(maxcnt);
	  
	  if (maxcnt > 0)
	    {
	      for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
		{
		  if (ParallelDescriptor::MyProc() == i)
		    {
		      if (cnt > 0)
			{
			  std::cout << "Processor "
				    << i
				    << " removed "
				    << cnt
				    << " particles not in finest level" << std::endl;
			}
		    }
		  ParallelDescriptor::Barrier();
		}
	    }
#ifdef BL_LAZY
	});
#endif
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::CreateVirtualParticles (int   level,
						   AoS& virts) const
{
  BL_PROFILE("ParticleContainer<NR, NI>::CreateVirtualParticles()");
  BL_ASSERT(level > 0);
  BL_ASSERT(virts.empty());
  
  if (level >= int(m_particles.size()))
    //
    // This level could exist and simply have no particles.
    //
    return;
  //
  // Read these from the parm file if we haven't done so yet.
  //
  if (aggregation_type == "")
    {
      ParmParse pp("particles");
      aggregation_type = "None";
      pp.query("aggregation_type",aggregation_type);
      aggregation_buffer = 2;
      pp.query("aggregation_buffer",aggregation_buffer);
    }
  //
  // Create a buffer so that particles near the cf border are not aggregated.
  //
  BoxArray buffer = amrex::complementIn(m_gdb->Geom(level).Domain(), m_gdb->ParticleBoxArray(level));
  
  buffer.grow(aggregation_buffer);
  
  const AoSMap& pmap = m_particles[level];

  for (PIter it(*this, level); it.isValid(); ++it) {
    int gid = it.index();
    int tid = it.tileIndex();
    const AoS& pbox = pmap[gid][tid];
  
    //
    // Map for use in Cell aggregation.
    //
    std::map<IntVect,ParticleType> agg_map;
    
    for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
      {
	ParticleLocData pld;
	Particle<NR, NI>::Where(*it, m_gdb, pld);
	if (buffer.contains(pld.m_cell))
	  {
	    //
	    // It's in the no-aggregation buffer.
	    //
	    virts.push_back(*it);
	    //
	    // Set its id to indicate that it's a virt.
	    //
	    virts.back().m_idata.id = VirtualParticleID;
	  }
	else
	  {
	    if (aggregation_type == "None")
	      {
		//
		// No aggregation.  Simply clone the particle.
		//
		virts.push_back(*it);
		//
		// Set its id to indicate that it's a virt.
		//
		virts.back().m_idata.id = VirtualParticleID;
	      }
	    else if (aggregation_type == "Cell")
	      {
		//
		// Note that Cell aggregation assumes that p.m_rdata.arr[BL_SPACEDIM] is mass and
		// that all other components should be combined in a mass-weighted
		// average.
		//
		auto agg_map_it = agg_map.find(pld.m_cell);
		
		if (agg_map_it == agg_map.end())
		  {
		    //
		    // Add the particle.
		    //
		    ParticleType p = *it;
		    //
		    // Set its id to indicate that it's a virt.
		    //
		    p.m_idata.id = VirtualParticleID;
		    agg_map[pld.m_cell] = p;
		  }
		else
		  {
		    BL_ASSERT(agg_map_it != agg_map.end());
		    const ParticleType&  pnew       = *it;
		    ParticleType&        pold       = agg_map_it->second;
		    const Real           old_mass   = pold.m_rdata.arr[BL_SPACEDIM];
		    const Real           new_mass   = pnew.m_rdata.arr[BL_SPACEDIM];
		    const Real           total_mass = old_mass + new_mass;
		    //
		    // Set the position to the center of mass.
		    //
		    for (int i = 0; i < BL_SPACEDIM; i++)
		      {
			pold.m_rdata.pos[i] = (old_mass*pold.m_rdata.pos[i] + new_mass*pnew.m_rdata.pos[i])/total_mass;
		      }
		    BL_ASSERT(ParticleType::Index(pold,m_gdb->Geom(level)) == pld.m_cell);
		    //
		    // Set the metadata (presumably velocity) to the mass-weighted average.
		    //
		    for (int i = BL_SPACEDIM + 1; i < BL_SPACEDIM + NR; i++)
		      {
			pold.m_rdata.arr[i] = (old_mass*pold.m_rdata.arr[i] + new_mass*pnew.m_rdata.arr[i])/total_mass;
		      }
		    pold.m_rdata.arr[BL_SPACEDIM] = total_mass;
		  }
	      }
	    else if (aggregation_type == "Flow")
	      {
		amrex::Abort("Flow aggregation not implemented");
	      }
	    else 
	      {
		amrex::Abort("Unknown Particle Aggregation mode");
                }
	  }
      }
    if (aggregation_type == "Cell")
      {
	//
	// Add the aggregated particles to the virtuals.
	//
	for (const auto& kv : agg_map)
	  {
	    virts.push_back(kv.second);
	  }
      }
  }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::CreateGhostParticles (int   level,
						 int   ngrow,
						 AoS& ghosts) const
{
  BL_PROFILE("ParticleContainer<NR, NI>::CreateGhostParticles()");
  BL_ASSERT(ghosts.empty());
  BL_ASSERT(level < m_gdb->finestLevel());
  
  if (level >= int(m_particles.size()))
    //
    // This level could exist and simply have no particles.
    //
    return;
  
  const BoxArray& fine = m_gdb->ParticleBoxArray(level + 1);
  
  std::vector< std::pair<int,Box> > isects;
  
  const AoSMap& pmap = m_particles[level];
  for (PIter it(*this, level); it.isValid(); ++it) {
    int gid = it.index();
    int tid = it.tileIndex();
    const AoS& pbox = pmap[gid][tid];

    for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
      {
	//
	// Find particle location on the finer level.
	//
	const IntVect& iv = ParticleType::Index(*it,m_gdb->Geom(level+1));
	//
	// Is it in the grown finer level?
	//
	fine.intersections(Box(iv,iv),isects,false,ngrow);
	//
	// Here we add the particle to each potential grid.
	//
	for (const auto& isec : isects)
	  {
	    //
	    // Create a copy.
	    //
	    ParticleType p = *it;
	    ParticleLocData pld;
	    
	    ParticleType::SingleLevelWhere(p, m_gdb, pld, level + 1);
	    
	    //
	    // Set its id to indicate that it's a ghost.
	    //
	    p.m_idata.id = GhostParticleID;
	    
	    //
	    // Store it in the AoS.
	    //
	    ghosts.push_back(p);
	  }
      }
  }
}

//
// This redistributes valid particles and discards invalid ones.
//
template <int NR, int NI>
void
ParticleContainer<NR, NI>::Redistribute (int  lev_min,
					 int  nGrow)
{
  BL_PROFILE("ParticleContainer::Redistribute()");
  const int MyProc    = ParallelDescriptor::MyProc();
  const int  NProcs   = ParallelDescriptor::NProcs();
  Real      strttime  = ParallelDescriptor::second();
  
  //
  // On startup there are cases where Redistribute() could be called
  // with a given finestLevel() where that AmrLevel has yet to be defined.
  //
  int theEffectiveFinestLevel = m_gdb->finestLevel();
  
  while (!m_gdb->LevelDefined(theEffectiveFinestLevel))
    theEffectiveFinestLevel--;
  
  if (int(m_particles.size()) < theEffectiveFinestLevel+1) {
    if (ParallelDescriptor::IOProcessor()) {
      std::cout << "ParticleContainer<NR, NI>::Redistribute() resizing m_particles from "
		<< m_particles.size()
		<< " to "
		<< theEffectiveFinestLevel+1 << '\n';
    }
    m_particles.resize(theEffectiveFinestLevel+1);
  }

  // The valid particles that we don't own.
  std::map<int, AoS> not_ours;

  // struct-of-array data that also needs to be communicated
  std::map<int, SoA> comm_soa_data;
  for (int i = 0; i < NProcs; i++)
    comm_soa_data[i].resize(m_npartdata);

  for (int lev = lev_min, nlevs = m_particles.size(); lev < nlevs; lev++) {
    AoSMap& pmap = m_particles[lev];
    for (PIter it(*this, lev); it.isValid(); ++it) {
      int grid = it.index();
      int tile = it.tileIndex();
      AoS& pbox = pmap[grid][tile];
      unsigned first = 0;
      unsigned npart = pbox.size();
      if (npart != 0) {
	for (unsigned pindex = 0; pindex < npart; ++pindex) {
	  ParticleLocData pld;
	  ParticleType& p = pbox[pindex];
	  if (p.m_idata.id > 0) {
	    locateParticle(p, pld, lev_min, theEffectiveFinestLevel, nGrow);
	    if (p.m_idata.id > 0) {

	      // The owner of the particle is the CPU owning the finest grid
	      // in state data that contains the particle.
	      const int who = m_gdb->ParticleDistributionMap(pld.m_lev)[pld.m_grid];
	      if (who == MyProc) {
		if (pld.m_lev != lev || pld.m_grid != grid || pld.m_tile != tile) {
		  //
		  // We own it but must shift it to another place.
		  //
		  m_particles[pld.m_lev][pld.m_grid][pld.m_tile].push_back(p);
		  for (int comp = 0; comp < m_npartdata; comp++) {
		    Real pdata = m_partdata[lev][grid][tile][comp][pindex];
		    m_partdata[pld.m_lev][pld.m_grid][tile][comp].push_back(pdata);
		  }
		  //
		  // Invalidate the particle so we can reclaim its space.
		  //
		  p.m_idata.id = -p.m_idata.id;
		}
	      }
	      else {
		not_ours[who].push_back(p);
		for (int comp = 0; comp < m_npartdata; comp++) {
		  Real pdata = m_partdata[lev][grid][tile][comp][pindex];
		  if (communicate_comp[comp])
		    comm_soa_data[who][comp].push_back(pdata);
		}
		//
		// Invalidate the particle so we can reclaim its space.
		//
		p.m_idata.id = -p.m_idata.id;
	      }
	    }
	  }
	  // this is a valid particle
	  if (p.m_idata.id > 0) {
	    if (pindex != first) {
	      pbox[first] = p;
	      for (int comp = 0; comp < m_npartdata; comp++) {
		Array<Real>& vec = m_partdata[lev][grid][tile][comp];
		vec[first] = vec[pindex];
	      }
	    }
	    ++first;
	  }
	}
	pbox.erase(pbox.begin() + first, pbox.begin() + npart);
	for (int comp = 0; comp < m_npartdata; comp++) {
	  Array<Real>& vec = m_partdata[lev][grid][tile][comp];
	  vec.erase(vec.begin() + first, vec.begin() + npart);
	}
      }
	
      //
      // Remove any map entries for which the particle container is now empty.
      //
      //if (pmap_it->second.empty()) {
      //	pmap.erase(pmap_it++);
      //      }
      //      else {
      //	++pmap_it;
      //      }
    }
  }
  
  if (int(m_particles.size()) > theEffectiveFinestLevel+1)
    {
      //
      // Looks like we lost an AmrLevel on a regrid.
      //
      if (ParallelDescriptor::IOProcessor())
        {
	  std::cout << "ParticleContainer<NR, NI>::Redistribute() resizing m_particles from "
		    << m_particles.size()
		    << " to "
		    << theEffectiveFinestLevel+1 << '\n';
        }
      BL_ASSERT(int(m_particles.size()) >= 2);
      BL_ASSERT(m_particles[m_particles.size()-1].empty());
      
      m_particles.resize(theEffectiveFinestLevel+1);
    }
  
  if (ParallelDescriptor::NProcs() == 1)
    {
      BL_ASSERT(not_ours.empty());
    }
  else
    {
      RedistributeMPI(not_ours, comm_soa_data);
    }
  
  BL_ASSERT(OK(lev_min, nGrow, theEffectiveFinestLevel));
  
  if (m_verbose > 0)
    {
      Real stoptime = ParallelDescriptor::second() - strttime;
      
      ByteSpread();
      
#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
	  ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
	  if (ParallelDescriptor::IOProcessor())
	    std::cout << "ParticleContainer<NR, NI>::Redistribute() time: " << stoptime << "\n\n";
#ifdef BL_LAZY
	});
#endif
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::RedistributeMPI (std::map<int, AoS>& not_ours,
					    std::map<int, SoA>& comm_soa_data)
{
    BL_PROFILE("ParticleContainer<NR, NI>::RedistributeMPI()");
#if BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    //
    // We may now have particles that are rightfully owned by another CPU.
    //
    Array<int> Snds(NProcs, 0), Rcvs(NProcs, 0);

    int NumSnds = 0, NumRcvs = 0;
    
    for (const auto& kv : not_ours)
      {
        NumSnds       += kv.second.size();
        Snds[kv.first] = kv.second.size();
      }
    
    ParallelDescriptor::ReduceIntMax(NumSnds);

    if (NumSnds == 0)
      //
      // There's no parallel work to do.
      //
      return;

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    typedef std::map<int,int> IntIntMap;

    IntIntMap SndCnts, RcvCnts, rOffset;

    for (int i = 0; i < NProcs; i++)
      if (Snds[i] > 0)
	SndCnts[i] = Snds[i];

    for (int i = 0; i < NProcs; i++)
      {
        if (Rcvs[i] > 0)
	  {
            RcvCnts[i] = Rcvs[i];
            rOffset[i] = NumRcvs;
            NumRcvs   += Rcvs[i];
	  }
      }

    //
    // Don't need these anymore.
    //
    Array<int>().swap(Snds);
    Array<int>().swap(Rcvs);

    //
    // We'll store the particles we're to receive in a map indexed by proc # of receiver.
    //
    std::map<int, AoS> nparticles;

    for (const auto& kv : RcvCnts)
      {
        nparticles[kv.first].resize(kv.second);
      }
    Array<int>         owner(RcvCnts.size());
    Array<int>         index(RcvCnts.size());
    Array<MPI_Status>  stats(RcvCnts.size());
    Array<MPI_Request> rreqs(RcvCnts.size());
    //
    // First send/recv the integer parts of the particles.
    //
    {
      const int SeqNum     = ParallelDescriptor::SeqNum();
      const int iChunkSize = 2 + NI;
      //
      // Allocate data for rcvs as one big chunk.
      //
      Array<int> recvdata(NumRcvs * iChunkSize);
      //
      // Post receives.
      //
      int idx = 0;
      for (auto it = RcvCnts.cbegin(); it != RcvCnts.cend(); ++it, ++idx)
        {
	  const int Who = it->first;
	  const int Cnt = it->second   * iChunkSize;
	  const int Idx = rOffset[Who] * iChunkSize;
	  
	  BL_ASSERT(Cnt > 0);
	  BL_ASSERT(Who >= 0 && Who < NProcs);
	  BL_ASSERT(Cnt < std::numeric_limits<int>::max());
	  
	  owner[idx] = Who;
	  rreqs[idx] = ParallelDescriptor::Arecv(&recvdata[Idx],Cnt,Who,SeqNum).req();
        }
      //
      // Send the integer data.
      //
      Array<int> senddata;
      
      for (const auto& kv : SndCnts)
        {
	  const int Who = kv.first;
	  const int Cnt = kv.second * iChunkSize;
	  
	  BL_ASSERT(Cnt > 0);
	  BL_ASSERT(Who >= 0 && Who < NProcs);
	  BL_ASSERT(Cnt < std::numeric_limits<int>::max());

	  senddata.resize(Cnt);
	  
	  const AoS& pbox = not_ours[Who];
	  
	  int ioff = 0;
	  for (const auto& p : pbox)
            {
	      BL_ASSERT(p.m_idata.id > 0);
	      
	      senddata[ioff+0] = p.m_idata.id;
	      senddata[ioff+1] = p.m_idata.cpu;
	      
	      for (int i = 0; i < NI; ++i) {
		senddata[i+ioff+2] = p.m_idata.arr[i+2];
	      }
	      
	      ioff += iChunkSize;
            }
	  
	  ParallelDescriptor::Send(senddata.dataPtr(),Cnt,Who,SeqNum);
        }
      //
      // Free up this memory ...
      //
      Array<int>().swap(senddata);
      //
      // Now receive and unpack the integer data.
      //
      for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
        {
	  ParallelDescriptor::Waitsome(rreqs, completed, index, stats);
	  
	  for (int k = 0; k < completed; k++)
            {
	      const int  Who  = owner[index[k]];
	      const int  Idx  = rOffset[Who] * iChunkSize;
	      const int* rcvp = &recvdata[Idx];
	      AoS&      pbox = nparticles[Who];
	      
	      BL_ASSERT(int(pbox.size()) == RcvCnts[Who]);
	      
	      for (auto& p : pbox)
                {
		  BL_ASSERT(rcvp != 0);
		  
		  p.m_idata.id   = rcvp[0];
		  p.m_idata.cpu  = rcvp[1];
		  
		  for (int i = 0; i < NI; ++i) {
		    p.m_idata.arr[i+2] = rcvp[i+2];
		  }
		  
		  rcvp += iChunkSize;
                }
            }
        }
    }

    int num_soa_comps = 0;
    for (int i = 0; i < m_npartdata; i++) {
      num_soa_comps += communicate_comp[i];
    }
    
    //
    // Next send/recv the Real parts of the particles.
    //
    {
      const int SeqNum     = ParallelDescriptor::SeqNum();
      const int rChunkSize = BL_SPACEDIM + NR + num_soa_comps;
      //
      // Allocate data for rcvs as one big chunk.
      //
      Array<typename ParticleType::RealType> recvdata(NumRcvs * rChunkSize);
      //
      // Post receives.
      //
      int idx = 0;
      for (auto it = RcvCnts.cbegin(); it != RcvCnts.cend(); ++it, ++idx)
        {
	  const int Who = it->first;
	  const int Cnt = it->second   * rChunkSize;
	  const int Idx = rOffset[Who] * rChunkSize;
	  
	  BL_ASSERT(Cnt > 0);
	  BL_ASSERT(Who >= 0 && Who < NProcs);
	  BL_ASSERT(Cnt < std::numeric_limits<int>::max());
	  
	  rreqs[idx] = ParallelDescriptor::Arecv(&recvdata[Idx],Cnt,Who,SeqNum).req();
        }
      //
      // Send the Real data.
      //
      Array<typename ParticleType::RealType> senddata;
      
      for (const auto& kv : SndCnts)
        {
	  const int Who = kv.first;
	  const int Cnt = kv.second * rChunkSize;
	  
	  BL_ASSERT(Cnt > 0);
	  BL_ASSERT(Who >= 0 && Who < NProcs);
	  BL_ASSERT(Cnt < std::numeric_limits<int>::max());
          
	  senddata.resize(Cnt);
	  
	  AoS& pbox = not_ours[Who];
	  SoA& StructOfArrayData = comm_soa_data[Who];

	  int ioff = 0;
	  int pindex = 0;
	  for (const auto& p : pbox)
            {
	      BL_ASSERT(p.m_idata.id > 0);
	      
	      D_TERM(senddata[ioff+0] = p.m_rdata.pos[0];,
		     senddata[ioff+1] = p.m_rdata.pos[1];,
		     senddata[ioff+2] = p.m_rdata.pos[2];);
	      
	      ioff += BL_SPACEDIM;
	      
	      for (int j = 0; j < NR; j++)
		senddata[ioff+j] = p.m_rdata.arr[BL_SPACEDIM + j];
	      
	      ioff += NR;

	      for (int j = 0; j < m_npartdata; j++) {
		if (communicate_comp[j])
		  senddata[ioff+j] = StructOfArrayData[j][pindex];
	      }

	      ioff += num_soa_comps;
	      pindex++;
            }
	  
	  AoS().swap(pbox);
	  
	  ParallelDescriptor::Send(senddata.dataPtr(), Cnt, Who, SeqNum);
        }
      //
      // Free up this memory ...
      //
      Array<typename ParticleType::RealType>().swap(senddata);
      //
      // Now receive and unpack the Real data.
      //
      for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
        {
	  ParallelDescriptor::Waitsome(rreqs, completed, index, stats);
	  
	  for (int k = 0; k < completed; k++)
            {
	      const int                              Who  = owner[index[k]];
	      const int                              Idx  = rOffset[Who] * rChunkSize;
	      const typename ParticleType::RealType* rcvp = &recvdata[Idx];
	      AoS&                                  pbox = nparticles[Who];
	      
	      BL_ASSERT(int(pbox.size()) == RcvCnts[Who]);
	      
	      ParticleLocData pld;
	      for (auto& p : pbox)
                {
		  BL_ASSERT(rcvp != 0);
		  
		  D_TERM(p.m_rdata.pos[0] = rcvp[0];,
			 p.m_rdata.pos[1] = rcvp[1];,
			 p.m_rdata.pos[2] = rcvp[2];);
		  
		  rcvp += BL_SPACEDIM;
		  
		  for (int j = 0; j < NR; j++)
		    p.m_rdata.arr[BL_SPACEDIM + j] = rcvp[j];

		  rcvp += NR;
		  
		  ParticleType::Where(p, m_gdb, pld);
		  
		  m_particles[pld.m_lev][pld.m_grid][pld.m_tile].push_back(p);
		  
		  for (int j = 0; j < m_npartdata; j++) {
		    if (communicate_comp[j])
		      m_partdata[pld.m_lev][pld.m_grid][pld.m_tile][j].push_back(rcvp[j]);
		    else
		      m_partdata[pld.m_lev][pld.m_grid][pld.m_tile][j].push_back(0.0);		      
		  }
		  
		  rcvp += num_soa_comps;
                }
	      
	      AoS().swap(pbox);
            }
        }
    }
#endif /*BL_USE_MPI*/
}

template <int NR, int NI>
bool
ParticleContainer<NR, NI>::OK (int  lev_min,
			       int  ngrow,
			       int  finest_level) const
{
    BL_PROFILE("ParticleContainer<NR, NI>::OK()");
    if (finest_level == -1)
        finest_level = m_gdb->finestLevel();

    BL_ASSERT(finest_level <= m_gdb->finestLevel());
    //
    // Check that the integer data in each valid particle is what it should be.
    // This includes checking that particles are in the proper place in the particle
    // container based on what Where() says they should be.
    //
    // Particles are copied to avoid accidentally moving them with where.
    //
    for (int lev = lev_min, nlevs=m_particles.size(); lev < nlevs; lev++)
      {
        const AoSMap& pmap = m_particles[lev];
	for (PIter it(*this, lev); it.isValid(); ++it) {
	  int grid = it.index();
	  int tile = it.tileIndex();
	  const AoS& pbox = pmap.at(grid)[tile];
	  for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
            {
	      //
	      // Yes I want to make a copy of the particle.
	      //
	      ParticleType p = *it;
	      
	      if (p.m_idata.id > 0)
                {
		  
		  ParticleLocData pld;
		  Particle<NR, NI>::Where(p, m_gdb, pld);
		  
		  const int     llev  = pld.m_lev;
		  const int     lgrid = pld.m_grid;
		  const int     ltile = pld.m_tile;
		  const IntVect cell  = pld.m_cell;
		  
		  if (!ParticleType::Where(p, m_gdb, pld, lev_min, finest_level))
                    {
		      if (!ParticleType::PeriodicWhere(p, m_gdb, pld, lev_min, finest_level)) 
			{
			  if (!ParticleType::RestrictedWhere(p, m_gdb, pld, ngrow))
			    return false;
			}
                    }
		  if ((lev  != pld.m_lev  || lev  != llev)  ||
		      (grid != pld.m_grid || grid != lgrid) || 
		      cell != pld.m_cell  || ltile != 0)
                    {
		      std::cout << "PARTICLE NUMBER " << p.m_idata.id << '\n';
		      
		      std::cout << "POS IS ";
		      for (int i = 0; i < BL_SPACEDIM; i++)
			std::cout << p.m_rdata.pos[i] << ' ';
		      
		      if (lev != pld.m_lev || lev != llev)
			std::cout << "BAD LEV  " << lev  << " " << pld.m_lev << '\n';
		      
		      if (grid != pld.m_grid || grid != lgrid)
			std::cout << "BAD GRID " << grid << " " << pld.m_grid << '\n';
		      
		      if (cell != pld.m_cell)
			std::cout << "BAD CELL " << cell << " " << pld.m_cell << '\n';
		      
		      return false;
                    }
                }
            }
        }
      }
    
    return true;
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::Checkpoint (const std::string& dir,
				       const std::string& name,
				       bool               is_checkpoint) const
{
    BL_PROFILE("ParticleContainer<NR, NI>::Checkpoint()");
    BL_ASSERT(OK());

    BL_ASSERT(sizeof(typename ParticleType::RealType) == 4 || sizeof(typename ParticleType::RealType) == 8);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();
    //
    // We store the particles in a subdirectory of "dir".
    //
    std::string pdir = dir;

    if (!pdir.empty() && pdir[pdir.size()-1] != '/')
        pdir += '/';

    pdir += name;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!amrex::UtilCreateDirectory(pdir, 0755))
            amrex::CreateDirectoryFailed(pdir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();
    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the I/O processor writes to the header file.
    //
    std::ofstream HdrFile;

    long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
      AoSMap& pmap = m_particles[lev];
      for (PIter it(*this, lev); it.isValid(); ++it) {
	int grid = it.index();
	int tile = it.tileIndex();
	AoS& pbox = pmap[grid][tile];
	for (const auto& p : pbox) {
	  if (p.m_idata.id > 0)
	    //
	    // Only count (and checkpoint) valid particles.
	    //
	    nparticles++;
	}
      }
    }

    ParallelDescriptor::ReduceLongSum(nparticles,IOProc);

    int maxnextid = ParticleType::NextID();

    ParticleType::NextID(maxnextid);

    ParallelDescriptor::ReduceIntMax(maxnextid,IOProc);

    if (ParallelDescriptor::IOProcessor())
      {
        std::string HdrFileName = pdir;
	
        if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
	  HdrFileName += '/';
	
        HdrFileName += "Header";
	
        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);
	
        if (!HdrFile.good())
	  amrex::FileOpenFailed(HdrFileName);
        //
        // First thing written is our Checkpoint/Restart version string.
        // 
        // We append "_single" or "_double" to the version string indicating
        // whether we're using "float" or "double" floating point data in the
        // particles so that we can Restart from the checkpoint files.
        //
        if (sizeof(typename ParticleType::RealType) == 4)
	  {
            HdrFile << ParticleType::Version() << "_single" << '\n';
	  }
        else
	  {
            HdrFile << ParticleType::Version() << "_double" << '\n';
	  }
        //
        // BL_SPACEDIM and N for sanity checking.
        //
        HdrFile << BL_SPACEDIM << '\n';
	
        HdrFile << NR << '\n';
        //
        // The total number of particles.
        //
        HdrFile << nparticles << '\n';
        //
        // The value of nextid that we need to restore on restart.
        //
        HdrFile << maxnextid << '\n';
        //
        // Then the finest level of the AMR hierarchy.
        //
        HdrFile << m_gdb->finestLevel() << '\n';
        //
        // Then the number of grids at each level.
        //
        for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
	  {
            HdrFile << m_gdb->ParticleBoxArray(lev).size() << '\n';
	  }
    }
    //
    // We want to write the data out in parallel.
    //
    // We'll allow up to nOutFiles active writers at a time.
    //
    int nOutFiles(64);
    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) {
      nOutFiles = NProcs;
    }
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
    
    for (int lev = 0; lev <= m_gdb->finestLevel(); lev++)
      {
        const bool gotsome = (NumberOfParticlesAtLevel(lev) > 0);
        //
        // We store the particles at each level in their own subdirectory.
        //
        std::string LevelDir = pdir;
	
        if (gotsome)
	  {
            if (!LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/')
	      LevelDir += '/';
	    
            LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);
	    
            if (ParallelDescriptor::IOProcessor())
	      if (!amrex::UtilCreateDirectory(LevelDir, 0755))
		amrex::CreateDirectoryFailed(LevelDir);
            //
            // Force other processors to wait till directory is built.
            //
            ParallelDescriptor::Barrier();
	  }
	
	MFInfo info;
	info.SetAlloc(false);
	MultiFab state(m_gdb->ParticleBoxArray(lev),
		       m_gdb->ParticleDistributionMap(lev),
		       1,0,info);
        //
        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        //
        Array<int>  which(state.size(),0);
        Array<int > count(state.size(),0);
        Array<long> where(state.size(),0);
	
        if (gotsome)
	  {
            const int   FileNumber   = MyProc % nOutFiles;
            std::string FullFileName = LevelDir;

            FullFileName += '/';
            FullFileName += ParticleType::DataPrefix();
            FullFileName += amrex::Concatenate("", FileNumber, 4);

            std::ofstream ParticleFile;

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            ParticleFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            const int nSets = ((NProcs + (nOutFiles - 1)) / nOutFiles);
            const int mySet = (MyProc / nOutFiles);

            for (int iSet = 0; iSet < nSets; ++iSet)
            {
                if (mySet == iSet)
                {
                    //
                    // Write all the data at this level to the file.
                    //
                    if (iSet == 0)
                        //
                        // First set.
                        //
                        ParticleFile.open(FullFileName.c_str(),
                                          std::ios::out|std::ios::trunc|std::ios::binary);
                    else
                    {
                        ParticleFile.open(FullFileName.c_str(),
                                          std::ios::out|std::ios::app|std::ios::binary);
                        //
                        // Set to the end of the file.
                        //
                        ParticleFile.seekp(0, std::ios::end);
                    }

                    if (!ParticleFile.good())
                        amrex::FileOpenFailed(FullFileName);
                    //
                    // Write out all the valid particles we own at the specified level.
                    // Do it grid block by grid block remembering the seek offset
                    // for the start of writing of each block of data.
                    //
                    WriteParticles(lev, ParticleFile, FileNumber, which, count, where, is_checkpoint);

                    ParticleFile.flush();

                    ParticleFile.close();

                    if (!ParticleFile.good())
                        amrex::Abort("ParticleContainer<NR, NI>::Checkpoint(): problem writing ParticleFile");

                    int iBuff = 0, wakeUpPID = (MyProc + nOutFiles), tag = (MyProc % nOutFiles);

                    if (wakeUpPID < NProcs)
                    {
                        ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
                    }
                }

                if (mySet == (iSet + 1))
                {
                    //
                    // Next set waits.
                    //
                    int iBuff, waitForPID = (MyProc - nOutFiles), tag = (MyProc % nOutFiles);

                    ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
                }
            }

            ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProc);
            ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProc);
            ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProc);
        }

        if (ParallelDescriptor::IOProcessor())
        {
            for (int j = 0; j < state.size(); j++)
            {
                //
                // We now write the which file, the particle count, and the
                // file offset into which the data for each grid was written,
                // to the header file.
                //
                HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
            }

            if (gotsome)
            {
                //
                // Unlink any zero-length data files.
                //
                Array<long> cnt(nOutFiles,0);

                for (int i = 0, N=count.size(); i < N; i++)
                    cnt[which[i]] += count[i];

                for (int i = 0, N=cnt.size(); i < N; i++)
                {
                    if (cnt[i] == 0)
                    {
                        std::string FullFileName = LevelDir;

                        FullFileName += '/';
                        FullFileName += ParticleType::DataPrefix();
                        FullFileName += amrex::Concatenate("", i, 4);

                        amrex::UnlinkFile(FullFileName.c_str());
                    }
                }
            }
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            HdrFile.flush();

            HdrFile.close();

            if (!HdrFile.good())
                amrex::Abort("ParticleContainer<NR, NI>::Checkpoint(): problem writing HdrFile");

            std::cout << "ParticleContainer<NR, NI>::Checkpoint() time: " << stoptime << '\n';
        }
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::WritePlotFile (const std::string& dir,
					 const std::string& name) const
{
    BL_PROFILE("ParticleContainer<NR, NI>::WritePlotFile()");
    BL_ASSERT(OK());
    bool is_checkpoint = false;

    // For yt we need exactly the chk particle format so would need to set is_checkpoint = true
    // Anyway, it's not too bad to have particle ids on disk,
    // think of merger trees or backtracing of particles for nested ics
    // is_checkpoint = true; 
    Checkpoint(dir,name,is_checkpoint);
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::WriteParticles (int            lev,
					  std::ofstream& ofs,
					  int            fnum,
					  Array<int>&    which,
					  Array<int>&    count,
					  Array<long>&   where,
					  bool           is_checkpoint) const
{
    BL_PROFILE("ParticleContainer<NR, NI>::WriteParticles()");
    const AoSMap&     pmap  = m_particles[lev];

    MFInfo info;
    info.SetAlloc(false);
    MultiFab state(m_gdb->ParticleBoxArray(lev),
		   m_gdb->ParticleDistributionMap(lev),
		   1,0,info);

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();
        //
        // Only write out valid particles.
        //
        int cnt = 0;

        auto pmap_it = pmap.find(grid);

        if (pmap_it != pmap.end())
        {
	  for (const auto& p : pmap_it->second[0]) // Tile
            {
                if (p.m_idata.id > 0)
                    cnt++;
            }
        }

        which[grid] = fnum;
        count[grid] = cnt;
        where[grid] = VisMF::FileOffset(ofs);

        if (cnt == 0) continue;

        const AoS& pbox = pmap_it->second[0]; // Tile

        if (is_checkpoint)
        {
            //
            // First write out the integer data in binary.
            // We do not need to write out the m_lev and m_grid
            // info since it's implicit in how the particles
            // are stored.  We can easily recreate them on restart.
            //
            const int iChunkSize = 2;

#ifdef BL_LOWMEMPWRITE
	    int maxItemsToWrite(8192);
	    int cntBufSize(maxItemsToWrite * iChunkSize);
	    int nItems(cnt), nItemsToWrite(0);
            Array<int> istuff(cntBufSize);
            auto it  = pbox.cbegin();
            auto End = pbox.cend();

	    while(nItems > 0) {
              int *iptr = istuff.dataPtr();
	      int itemCount(0);
              for( ; it != End && itemCount < maxItemsToWrite; ++it) {
                if(it->m_idata.id > 0) {

                    iptr[0] = it->m_idata.id;
                    iptr[1] = it->m_idata.cpu;

                    iptr += iChunkSize;
                }
		++itemCount;
              }
	      nItemsToWrite = nItems > maxItemsToWrite ? maxItemsToWrite : nItems;
              ofs.write((char *) istuff.dataPtr(), nItemsToWrite * iChunkSize * sizeof(int));
	      nItems -= nItemsToWrite;
	    }
#else
            Array<int> istuff(cnt*iChunkSize);

            int* iptr = istuff.dataPtr();

	    for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
            {
                if (it->m_idata.id > 0)
                {
                    iptr[0] = it->m_idata.id;
                    iptr[1] = it->m_idata.cpu;

                    iptr += iChunkSize;
                }
            }

            ofs.write((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
#endif
        }


        //
        // Write the Real data in binary.
        //
        const int rChunkSize = BL_SPACEDIM+NR;

#ifdef BL_LOWMEMPWRITE
	int maxItemsToWrite(8192);
	int cntBufSize(maxItemsToWrite * rChunkSize);
	int nItems(cnt), nItemsToWrite(0);
        Array<ParticleType::RealType> rstuff(cntBufSize);
        auto it  = pbox.cbegin();
        auto End = pbox.cend();

	while(nItems > 0) {
          ParticleType::RealType *rptr = rstuff.dataPtr();
	  int itemCount(0);
          for( ; it != End && itemCount < maxItemsToWrite; ++it) {
            if(it->m_idata.id > 0) {
                for (int i = 0; i < BL_SPACEDIM + NR; i++)
                  rptr[i] = it->m.m_rdata.arr[i];
                rptr += rChunkSize;
            }
	    ++itemCount;
          }

	  nItemsToWrite = nItems > maxItemsToWrite ? maxItemsToWrite : nItems;
          ofs.write((char *) rstuff.dataPtr(), nItemsToWrite * rChunkSize * sizeof(typename ParticleType::RealType));
	  nItems -= nItemsToWrite;
	}
#else
        Array<typename ParticleType::RealType> rstuff(cnt*rChunkSize);

        typename ParticleType::RealType* rptr = rstuff.dataPtr();

	for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
        {
            if (it->m_idata.id > 0)
            {
	      for (int i = 0; i < BL_SPACEDIM + NR; i++)
		rptr[i] = it->m.m_rdata.arr[i];
	      rptr += rChunkSize;
            }
        }

        ofs.write((char*)rstuff.dataPtr(),rstuff.size()*sizeof(typename ParticleType::RealType));
#endif
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::Restart (const std::string& dir,
				     const std::string& file,
				     bool is_checkpoint)
{
    BL_PROFILE("ParticleContainer<NR, NI>::Restart()");
    BL_ASSERT(!dir.empty());
    BL_ASSERT(!file.empty());

    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();

    std::string fullname = dir;

    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';

    fullname += file;
    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the IO processor reads the header file.
    //
    // It'll then broadcast() stuff of interest to all CPUs.
    //
    std::ifstream HdrFile;

    std::string HdrFileName = fullname;

    if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
        HdrFileName += '/';

    HdrFileName += "Header";

    HdrFile.open(HdrFileName.c_str(), std::ios::in);

    if (!HdrFile.good())
        amrex::FileOpenFailed(HdrFileName);
    //
    // First value should be the version string.
    //
    Array<char> vbuf(128);

    std::string version;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> version;

        BL_ASSERT(!version.empty());
        BL_ASSERT(vbuf.size() > version.size());

        for (size_t i = 0; i < version.size(); ++i)
            vbuf[i] = version[i];

        vbuf[version.size()] = '\0';
    }

    ParallelDescriptor::Bcast(vbuf.dataPtr(), vbuf.size(), IOProc);
    //
    // What do our version strings mean?
    //
    // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
    // 
    // "Version_One_Dot_One" -- can write out either as either single or double precision.
    //
    // Appended to the latter version string are either "_single" or "_double" to
    // indicate how the particles were written.
    //
    version = vbuf.dataPtr();

    if (version.find("Version_One_Dot_Zero") != std::string::npos)
    {
        Restart_Doit(fullname,HdrFile,"double",is_checkpoint);
    }
    else if (version.find("Version_One_Dot_One") != std::string::npos)
    {
        if (version.find("_single") != std::string::npos)
        {
            Restart_Doit(fullname,HdrFile,"single",is_checkpoint);
        }
        else if (version.find("_double") != std::string::npos)
        {
            Restart_Doit(fullname,HdrFile,"double",is_checkpoint);
        }
        else
        {
            std::string msg("ParticleContainer<NR, NI>::Restart(): bad version string: ");
            msg += version;
            amrex::Error(version.c_str());
        }
    }
    else
    {
        std::string msg("ParticleContainer<NR, NI>::Restart(): unknown version string: ");
        msg += version;
        amrex::Abort(msg.c_str());
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR, NI>::Restart() time: " << stoptime << '\n';
        }
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::Restart_Doit (const std::string& fullname,
					std::ifstream&     HdrFile,
					const std::string& how,
					bool is_checkpoint)
{
    BL_PROFILE("ParticleContainer<NR, NI>::RestartDoit()");
    BL_ASSERT(!fullname.empty());

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    //
    // Next value should be BL_SPACEDIM;
    //
    int dm;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> dm;

        if (dm != BL_SPACEDIM)
            amrex::Abort("ParticleContainer<NR, NI>::Restart(): dm != BL_SPACEDIM");
    }
    ParallelDescriptor::Bcast(&dm, 1, IOProc);
    //
    // Next value should be our "N".
    //
    int n;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> n;

        if (n != NR)
            amrex::Abort("ParticleContainer<NR, NI>::Restart(): n != N");
    }
    ParallelDescriptor::Bcast(&n, 1, IOProc);

    long nparticles;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The total number of particles.
        //
        HdrFile >> nparticles;

        BL_ASSERT(nparticles >= 0);
    }
    ParallelDescriptor::Bcast(&nparticles, 1, IOProc);

    int maxnextid;

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // The value of nextid that we need to restore.
        //
        HdrFile >> maxnextid;

        BL_ASSERT(maxnextid > 0);
    }
    ParallelDescriptor::Bcast(&maxnextid, 1, IOProc);
    //
    // Don't forget to restore it!!!
    //
    ParticleType::NextID(maxnextid);
    //
    // Then the finest level of the AMR hierarchy.
    //
    int finest_level;

    if (ParallelDescriptor::IOProcessor())
    {
        HdrFile >> finest_level;

        BL_ASSERT(finest_level >= 0);
    }
    ParallelDescriptor::Bcast(&finest_level, 1, IOProc);
    //
    // Then the number of grids at each level.
    //
    Array<int> ngrids(finest_level+1);

    BL_ASSERT(finest_level == m_gdb->finestLevel());

    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev = 0; lev <= finest_level; lev++)
        {
            HdrFile >> ngrids[lev];

            BL_ASSERT(ngrids[lev] > 0);
            BL_ASSERT(ngrids[lev] == int(m_gdb->ParticleBoxArray(lev).size()));
        }
    }
    ParallelDescriptor::Bcast(ngrids.dataPtr(), ngrids.size(), IOProc);
    //
    // The rest of HdrFile consists of triples of the form:
    //
    //   which count offset
    //
    // One for each grid at each level from 0 -> finest_level.
    //
    // We rebuild the filename from which and level.
    //
    for (int lev = 0; lev <= finest_level; lev++)
    {
        //
        // Read in the which, count & offset info for this level.
        //
        Array<int>  which(ngrids[lev]);
        Array<int>  count(ngrids[lev]);
        Array<long> where(ngrids[lev]);

        if (ParallelDescriptor::IOProcessor())
        {
            for (int i = 0; i < ngrids[lev]; i++)
            {
                HdrFile >> which[i] >> count[i] >> where[i];
            }
        }
        ParallelDescriptor::Bcast(which.dataPtr(), which.size(), IOProc);
        ParallelDescriptor::Bcast(count.dataPtr(), count.size(), IOProc);
        ParallelDescriptor::Bcast(where.dataPtr(), where.size(), IOProc);

        m_particles.resize(m_gdb->finestLevel()+1);

	MFInfo info;
	info.SetAlloc(false);
	MultiFab state(m_gdb->ParticleBoxArray(lev),
		       m_gdb->ParticleDistributionMap(lev),
		       1,0,info);

        for (MFIter mfi(state); mfi.isValid(); ++mfi)
        {
            const int grid = mfi.index();

            if (count[grid] <= 0) continue;
            //
            // The file names in the header file are relative.
            //
            std::string name = fullname;

            if (!name.empty() && name[name.size()-1] != '/')
                name += '/';

            name += "Level_";
            name += amrex::Concatenate("", lev, 1);
            name += '/';
            name += ParticleType::DataPrefix();
            name += amrex::Concatenate("", which[grid], 4);

            std::ifstream ParticleFile;

            ParticleFile.open(name.c_str(), std::ios::in);

            if (!ParticleFile.good())
                amrex::FileOpenFailed(name);

            ParticleFile.seekg(where[grid], std::ios::beg);

            if (how == "single")
            {
                ReadParticles_SinglePrecision(count[grid],grid,lev,is_checkpoint,ParticleFile);
            }
            else if (how == "double")
            {
                ReadParticles_DoublePrecision(count[grid],grid,lev,is_checkpoint,ParticleFile);
            }
            else
            {
                std::string msg("ParticleContainer<NR, NI>::Restart_Doit(): bad parameter: ");
                msg += how;
                amrex::Error(msg.c_str());
            }
                
            ParticleFile.close();

            if (!ParticleFile.good())
                amrex::Abort("ParticleContainer<NR, NI>::Restart_Doit(): problem reading particles");
        }
    }

    BL_ASSERT(OK());        
}

//
// This one stores real data as doubles.
//

template <int NR, int NI>
void
ParticleContainer<NR, NI>::ReadParticles_DoublePrecision (int            cnt,
							 int            grd,
							 int            lev,
							 bool           is_checkpoint,
							 std::ifstream& ifs)
{
    BL_PROFILE("ParticleContainer<NR, NI>::ReadParticles_DoublePrecision()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < int(m_particles.size()));
    BL_ASSERT(lev >= 0 && lev <= m_gdb->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_gdb->ParticleBoxArray(lev).size());
    //
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    //
    const int iChunkSize = 2;

    Array<int> istuff(cnt*iChunkSize);

    if (is_checkpoint)
        ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the double data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+NR;

    Array<double> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(double));
    //
    // Now reassemble the particles.
    //
    int*            iptr = istuff.dataPtr();
    double*         rptr = rstuff.dataPtr();
    AoS&           pbox = m_particles[lev][grd][0]; // tile
    ParticleType p;

    // If we are restarting from a plotfile instead of a checkpoint file, then we do not
    //    read in the particle id's, so we need to reset the id counter to zero and renumber them
    if (!is_checkpoint)
    {
        int maxnextid = 1;
        ParticleType::NextID(maxnextid);
    }

    for (int i = 0; i < cnt; i++)
    {
        if (is_checkpoint)
        {
            p.m_idata.id   = iptr[0];
            p.m_idata.cpu  = iptr[1];
        }
        else
        {
	  ParticleLocData pld;
	  if (!ParticleType::Where(p, m_gdb, pld))
            {
	      ParticleType::PeriodicShift(p, m_gdb);
	      
	      if (!ParticleType::Where(p, m_gdb, pld))
                {
		  std::cout << "RESTART:BAD PARTICLE ID WOULD BE " << ParticleType::NextID() << '\n';
		  
		  for (int d = 0; d < BL_SPACEDIM; d++)
                    {
		      std::cout << "RESTART:BAD PARTICLE POS(" << d << ") " << p.m_rdata.pos[d] << std::endl;
                    }
		  
		  amrex::Abort("ParticleContainer<NR, NI>::ReadParticles_DoublePrecision(): invalid particle");
                }
            }
	  
	  p.m_idata.id   = ParticleType::NextID();
	  p.m_idata.cpu  = ParallelDescriptor::MyProc();
        }

        BL_ASSERT(p.m_idata.id > 0);

        iptr += iChunkSize;

        D_TERM(p.m_rdata.pos[0] = rptr[0];,
               p.m_rdata.pos[1] = rptr[1];,
               p.m_rdata.pos[2] = rptr[2];);

        for (int i = 0; i < NR; i++)
            p.m_rdata.arr[BL_SPACEDIM+i] = rptr[BL_SPACEDIM+i];

        rptr += rChunkSize;

        pbox.push_back(p);
    }
}

//
// This one stores real data as floats.
//

template <int NR, int NI>
void
ParticleContainer<NR, NI>::ReadParticles_SinglePrecision (int            cnt,
							   int            grd,
							   int            lev,
							   bool           is_checkpoint,
							   std::ifstream& ifs)
{
    BL_PROFILE("ParticleContainer<NR, NI>::ReadParticles_SinglePrecision()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < int(m_particles.size()));
    BL_ASSERT(lev >= 0 && lev <= m_gdb->finestLevel());
    BL_ASSERT(grd >= 0 && grd < m_gdb->ParticleBoxArray(lev).size());
    //
    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    //
    const int iChunkSize = 2;

    Array<int> istuff(cnt*iChunkSize);

    if (is_checkpoint)
        ifs.read((char*)istuff.dataPtr(),istuff.size()*sizeof(int));
    //
    // Then the float data in binary.
    //
    const int rChunkSize = BL_SPACEDIM+NR;

    Array<float> rstuff(cnt*rChunkSize);

    ifs.read((char*)rstuff.dataPtr(),rstuff.size()*sizeof(float));
    //
    // Now reassemble the particles.
    //
    int*   iptr = istuff.dataPtr();
    float* rptr = rstuff.dataPtr();
    AoS&  pbox = m_particles[lev][grd][0]; //Tile

    ParticleLocData pld;
    ParticleType p;

    // If we are restarting from a plotfile instead of a checkpoint file, then we do not
    //    read in the particle id's, so we need to reset the id counter to zero and renumber them
    if (!is_checkpoint)
    {
        int maxnextid = 1;
        ParticleType::NextID(maxnextid);
    }

    for (int i = 0; i < cnt; i++)
    {
        p.m_idata.id   = iptr[0];
        p.m_idata.cpu  = iptr[1];

        if (is_checkpoint)
        {
            p.m_idata.id    = iptr[0];
            p.m_idata.cpu   = iptr[1];
        }
        else
        {
	  if (!ParticleType::Where(p, m_gdb, pld))
            {
	      ParticleType::PeriodicShift(p, m_gdb);

	      if (!ParticleType::Where(p, m_gdb, pld))
                {
                    std::cout << "RESTART:BAD PARTICLE ID WOULD BE " << ParticleType::NextID() << '\n';

                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        std::cout << "RESTART:BAD PARTICLE POS(" << d << ") " << p.m_rdata.pos[d] << std::endl;
                    }

                    amrex::Abort("ParticleContainer<NR, NI>::ReadParticles_SinglePrecision(): invalid particle");
                }
            }

            p.m_idata.id   = ParticleType::NextID();
            p.m_idata.cpu  = ParallelDescriptor::MyProc();
        }

        BL_ASSERT(p.m_idata.id > 0);

        iptr += iChunkSize;

        D_TERM(p.m_rdata.pos[0] = rptr[0];,
               p.m_rdata.pos[1] = rptr[1];,
               p.m_rdata.pos[2] = rptr[2];);

        for (int i = 0; i < NR; i++)
            p.m_rdata.arr[BL_SPACEDIM + i] = rptr[BL_SPACEDIM + i];

        rptr += rChunkSize;

        pbox.push_back(p);
    }
}

template <int NR, int NI>
void
ParticleContainer<NR, NI>::WriteAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer<NR, NI>::WriteAsciiFile()");
    BL_ASSERT(!filename.empty());

    const Real strttime = ParallelDescriptor::second();
    //
    // Count # of valid particles.
    //
    long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
      AoSMap& pmap = m_particles[lev];
      for (PIter it(*this, lev); it.isValid(); ++it) {
	int grid = it.index();
	int tile = it.tileIndex();
	AoS& pbox = pmap[grid][tile];
	for (const auto& p : pbox) {
	  if (p.m_idata.id > 0)
	    //
	    // Only count (and checkpoint) valid particles.
	    //
	    nparticles++;
	}
      }
    }
    
    //
    // And send count to I/O processor.
    //
    ParallelDescriptor::ReduceLongSum(nparticles,ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle count.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles << '\n';
            
        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer<NR, NI>::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
    {
        if (MyProc == i)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            std::ofstream File;

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

	    for (int lev = 0; lev < m_particles.size();  lev++) {
	      AoSMap& pmap = m_particles[lev];
	      for (PIter it(*this, lev); it.isValid(); ++it) {
		int grid = it.index();
		int tile = it.tileIndex();
		AoS& pbox = pmap[grid][tile];
		
		int index = 0;
		for (auto it = pbox.cbegin(); it != pbox.cend(); ++it) {
		    if (it->m_idata.id > 0) {
		      
		      File << it->m_idata.id  << ' ';
		      File << it->m_idata.cpu << ' ';
		      
		      D_TERM(File << it->m_rdata.pos[0] << ' ',
			     << it->m_rdata.pos[1] << ' ',
			     << it->m_rdata.pos[2] << ' ');
		      
		      for (int i = 0; i < m_npartdata; i++)
			File << m_partdata[lev][grid][0][i][index] << ' ';
		      index++;
		      
		      for (int i = BL_SPACEDIM; i < BL_SPACEDIM + NR; i++) {
			char ws = (i == BL_SPACEDIM + NR - 1) ? '\n' : ' ';
			File << it->m_rdata.arr[i] << ws;
		      }
		    }
		}
	      }
            }
	    
            File.flush();
	    
            File.close();

            if (!File.good())
	      amrex::Abort("ParticleContainer<NR, NI>::WriteAsciiFile(): problem writing file");
	    
        }
	
        ParallelDescriptor::Barrier();
    }
    
    if (m_verbose > 1)
      {
        Real stoptime = ParallelDescriptor::second() - strttime;
	
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
	
        if (ParallelDescriptor::IOProcessor())
	  {
            std::cout << "ParticleContainer<NR, NI>::WriteAsciiFile() time: " << stoptime << '\n';
	  }
      }
}

template<int NR, int NI>
inline
void
Particle<NR, NI>::CIC_Fracs (const Real* frac, Real* fracs)
{
    //
    // "frac"  should be dimensioned: Real frac[BL_SPACEDIM]
    //
    // "fracs" should be dimensioned: Real fracs[D_TERM(2,+2,+4)]
    //
#if (BL_SPACEDIM == 1)
    // High
    fracs[0] = frac[0];

    // Low
    fracs[1] = (1-frac[0]);

#elif (BL_SPACEDIM == 2)
    // HH
    fracs[0] = frac[0] * frac[1] ;
    
    // LH
    fracs[1] = (1-frac[0]) * frac[1];
    
    // LL
    fracs[2] = (1-frac[0]) * (1-frac[1]);
    
    // HL
    fracs[3] = frac[0] * (1-frac[1]);

#elif (BL_SPACEDIM == 3)
    // HHH
    fracs[0] = frac[0] * frac[1] * frac[2];

    // LHH
    fracs[1] = (1-frac[0]) * frac[1] * frac[2];

    // LLH
    fracs[2] = (1-frac[0]) * (1-frac[1]) * frac[2];
    
    // HLH
    fracs[3] = frac[0] * (1-frac[1]) * frac[2];

    // HHL
    fracs[4] = frac[0] * frac[1] * (1-frac[2]);
    
    // LHL
    fracs[5] = (1-frac[0]) * frac[1] * (1-frac[2]);

    // LLL
    fracs[6] = (1-frac[0]) * (1-frac[1]) * (1-frac[2]);
    
    // HLL
    fracs[7] = frac[0] * (1-frac[1]) * (1-frac[2]);
#endif
}

template<int NR, int NI>
inline
void
Particle<NR, NI>::CIC_Cells (const IntVect& hicell, IntVect* cells)
{
    //
    // "cells" should be dimensioned: IntVect cells[D_TERM(2,+2,+4)]
    //
    IntVect cell = hicell;

#if (BL_SPACEDIM == 1)
    // High
    cells[0] = cell;

    // Low
    cell[0]  = cell[0] - 1;
    cells[1] = cell;

#elif (BL_SPACEDIM == 2)
    // HH
    cells[0] = cell;
    
    // LH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    
    // LL
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    
    // HL
    cell[0]  = cell[0] + 1;
    cells[3] = cell;

#elif (BL_SPACEDIM == 3)
    // HHH
    cells[0] = cell;

    // LHH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;

    // LLH
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    
    // HLH
    cell[0]  = cell[0] + 1;
    cells[3] = cell;

    cell = hicell;

    // HHL
    cell[2]  = cell[2] - 1;
    cells[4] = cell;
    
    // LHL
    cell[0]  = cell[0] - 1;
    cells[5] = cell;

    // LLL
    cell[1]  = cell[1] - 1;
    cells[6] = cell;
    
    // HLL
    cell[0]  = cell[0] + 1;
    cells[7] = cell;
#endif
}

template<int NR, int NI>
inline
int
Particle<NR, NI>::CIC_Cells_Fracs (const Particle<NR, NI>& p,
				   const Real*         plo,
				   const Real*         dx,
				   Array<Real>&        fracs,
				   Array<IntVect>&     cells)
{
    return Particle<NR, NI>::CIC_Cells_Fracs(p,plo,dx,dx,fracs,cells);
}

//
// This is the multi-level version.
//
// The Array should be empty on input.
//
// There'll be finest_level+1 of them.
//
template <int NR, int NI>
void
ParticleContainer<NR, NI>::AssignDensity (int rho_index, bool sub_cycle,
					   Array<std::unique_ptr<MultiFab> >& mf_to_be_filled, 
					   int lev_min, int ncomp, int finest_level) const
{
    if (rho_index != 0) amrex::Abort("AssignDensity only works if rho_index = 0");

    BL_PROFILE("ParticleContainer<NR, NI>::AssignDensity()");
    BL_ASSERT(NR >= 1);
    BL_ASSERT(NR >= ncomp);
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);

    if (finest_level == -1)
    {
        finest_level = m_gdb->finestLevel();
    }
    while (!m_gdb->LevelDefined(finest_level))
    {
        finest_level--;
    }
    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level, 
    // lev_index is the corresponding index for mf. 
    //

    // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
    mf_to_be_filled.resize(finest_level+1-lev_min);
    for (int lev = lev_min; lev <= finest_level; lev++)
    { 
        const int lev_index = lev - lev_min;
        mf_to_be_filled[lev_index].reset(new MultiFab(m_gdb->boxArray(lev),
						      m_gdb->DistributionMap(lev),
						      ncomp, 1));
	mf_to_be_filled[lev_index]->setVal(0.0);
    }

    // Test whether the grid structure of the boxArray is the same
    //       as the ParticleBoxArray at all levels 
    bool all_grids_the_same = true; 
    for (int lev = lev_min; lev <= finest_level; lev++) {
        if (!OnSameGrids(lev, *mf_to_be_filled[lev-lev_min])) {
	    all_grids_the_same = false;
	    break;
	}
    }

    Array<std::unique_ptr<MultiFab> > mf_part;
    if (!all_grids_the_same)
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1-lev_min);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_part[lev_index].reset(new MultiFab(m_gdb->ParticleBoxArray(lev), 
						  m_gdb->ParticleDistributionMap(lev),
						  ncomp, 1));
	    mf_part[lev_index]->setVal(0.0);
        }
    }

    auto & mf = (all_grids_the_same) ? mf_to_be_filled : mf_part;

    if (finest_level == 0)
    {
        //
        // Just use the far simpler single-level version.
        //
        AssignDensitySingleLevel(rho_index, *mf[0],0,ncomp);
        //
        // I believe that we don't need any information in ghost cells so we don't copy those.
        //
        if ( ! all_grids_the_same) {
            mf_to_be_filled[0]->copy(*mf[0],0,0,ncomp);
	}
        return;
    }
    
    //
    // This is the "data" needed by other MPI procs.
    //
    std::map<int, std::vector<ParticleCommData> > data;

    const Real stime = ParallelDescriptor::second();
    //
    // Minimum M required.
    //
    const int M = D_TERM(2,+2,+4);

    Array<int>     cgrid(M);
    Array<int>    cwhich(M),  fwhich(M);
    Array<Real>    fracs(M),  cfracs(M);
    Array<IntVect> cells(M),  ccells(M), cfshifts(M);

    ParticleCommData pb;
    //
    // I'm going to allocate these badboys here & pass'm into routines that use'm.
    // This should greatly cut down on memory allocation/deallocation.
    //
    Array<IntVect>                    pshifts(27);
    std::vector< std::pair<int,Box> > isects;
    Array<int>                        fgrid(M);
    Array<Real>                       ffracs(M);
    Array<IntVect>                    fcells;
    //
    // "fvalid" contains all the valid region of the MultiFab at this level, together
    // with any ghost cells lying outside the domain, that can be periodically shifted into the
    // valid region.  "compfvalid" is the complement of the "fvalid", while "compfvalid_grown" is 
    // "compfvalid" grown by one.  Using these we can figure out whether or not a cell is in the
    // valid region of our MultiFab as well as whether or not we're at a Fine->Crse boundary.
    //
    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const Geometry& gm        = m_gdb->Geom(lev);
        const Geometry& gm_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1) : gm;
        const Geometry& gm_coarse = (lev > 0) ? m_gdb->Geom(lev-1) : gm;
        const Box&      dm        = gm.Domain();
        const Real*     dx        = gm.CellSize();
        const Real*     plo       = gm.ProbLo();
        const Real*     dx_fine   = (lev < finest_level) ? m_gdb->Geom(lev+1).CellSize() : dx;
        const Real*     dx_coarse = (lev > 0) ? m_gdb->Geom(lev-1).CellSize() : dx;
        const int       lev_index = lev - lev_min;
        const BoxArray& grids     = mf[lev_index]->boxArray();
        const int       dgrow     = (lev == 0) ? 1 : m_gdb->MaxRefRatio(lev-1);

        BoxArray compfvalid, compfvalid_grown, fvalid = mf[lev_index]->boxArray();
        //
        // Do we have Fine->Crse overlap on a periodic boundary?
        // We want to add all ghost cells that can be shifted into valid region.
        //
        BoxList valid;

        for (int i = 0; i < grids.size(); i++)
        {
            if (gm.isAnyPeriodic())
            {
                const Box& dest = amrex::grow(grids[i],dgrow);

                if ( ! dm.contains(dest))
                {
                    for (int j = 0; j < grids.size(); j++)
                    {
                        BL_ASSERT(dm.contains(grids[j]));

                        gm.periodicShift(dest, grids[j], pshifts);

			for (const auto& kiv : pshifts)
                        {
                            const Box& sbx = grids[j] + kiv;
                            const Box& dbx = dest & sbx;

                            BL_ASSERT(dbx.ok());

                            valid.push_back(dbx);
                        }
                    }
                }
            }
        }
        if (valid.isNotEmpty())
        {
            //
            // We've got some Fine->Crse periodic overlap.
            // Don't forget to add the valid boxes too.
            //
            for (int i = 0; i < grids.size(); i++) {
                valid.push_back(grids[i]);
	    }
            fvalid = BoxArray(valid);
            fvalid.removeOverlap();
        }
        //
        // If we're at a lev < finestLevel, this is the coarsened fine BoxArray.
        // We use this for figuring out Crse->Fine issues.
        //
        BoxArray ccba;
        if (lev > 0)
        {
            ccba = m_gdb->boxArray(lev);
            ccba.coarsen(m_gdb->refRatio(lev-1));
        }
        BoxArray cfba;
        if (lev < finest_level)
        {
            cfba = m_gdb->boxArray(lev+1);
            cfba.coarsen(m_gdb->refRatio(lev));

            BL_ASSERT(mf[lev_index]->boxArray().contains(cfba));
        }
        //
        // This is cfba with any shifted ghost cells.
        //
        BoxArray cfvalid = cfba;

        if (lev < finest_level)
        {
            BoxList cvalid;

            const BoxArray& cgrids = mf[lev_index]->boxArray();

            for (int i = 0; i < cfba.size(); i++)
            {
                if (gm.isAnyPeriodic())
                {
                    const Box& dest = amrex::grow(cfba[i],mf[lev_index]->nGrow());

                    if ( ! dm.contains(dest))
                   { 
                        for (int j = 0; j < cgrids.size(); j++)
                        {
                            BL_ASSERT(dm.contains(cgrids[j]));

                            gm.periodicShift(dest, cgrids[j], pshifts);

			    for (const auto& kiv : pshifts)
                            {
                                const Box& sbx = cfba[i] - kiv;

                                cvalid.push_back(sbx);
                            }
                        }
                    }
                }
            }
            if (cvalid.isNotEmpty())
            {
                //
                // We've got some Fine->Crse periodic overlap.
                // Don't forget to add the valid boxes too.
                //
                for (int i = 0; i < cfba.size(); i++) {
                    cvalid.push_back(cfba[i]);
		}
                cfvalid = BoxArray(cvalid);
                cfvalid.removeOverlap();
            }
        }
        //
        // The "+1" is so we enclose the valid region together with any
        //  ghost cells that can be periodically shifted into valid.
        //
        compfvalid = amrex::complementIn(amrex::grow(dm,dgrow+1), fvalid);

        compfvalid_grown = compfvalid;
        compfvalid_grown.grow(1);
        compfvalid_grown.removeOverlap();
            
        if (gm.isAnyPeriodic() && ! gm.isAllPeriodic())
        {
            amrex::Error("AssignDensity: problem must be periodic in no or all directions");
        }
        //
        // If we're at a lev > 0, this is the coarsened BoxArray.
        // We use this for figuring out Fine->Crse issues.
        //
        BoxArray cba;
        if (lev > 0)
        {
            cba = m_gdb->boxArray(lev);
            cba.coarsen(m_gdb->refRatio(lev-1));
        }
        //
        // Do the grids at this level cover the full domain? If they do
        // there can be no Fine->Crse interactions at this level.
        //
        const bool GridsCoverDomain = fvalid.contains(m_gdb->Geom(lev).Domain());

	const AoSMap& pmap = m_particles[lev];
	for (PIter it(*this, lev); it.isValid(); ++it) {
	  int grid = it.index();
	  int tile = it.tileIndex();
	  const AoS& pbx = pmap.at(grid)[tile];
	  FArrayBox&  fab = (*mf[lev_index])[grid];
	  
	  for (const auto& p : pbx)
            {
                if (p.m_idata.id <= 0) {
		  continue;
		}
                //
                // Get "fracs" and "cells" for the particle "p" at this level.
                //
                const int M = ParticleType::CIC_Cells_Fracs(p, plo, dx, fracs, cells);
                //
                // If this is not fully periodic then we have to be careful that no
                // particle's support leaves the domain. We test this by checking the low
                // and high corners respectively.
                //
                if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                    if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                        amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		    }
		}
                //
                // This section differs based on whether we subcycle.
                // Without subcycling we use the "stretchy" support for particles.
                // With subcycling a particles support is strictly defined 
                // by its resident level.
                //
                if (sub_cycle)
                {
                    bool isFiner    = false;
                    bool isBoundary = false;
                    //
                    // First sum the mass in the valid region
                    //
                    for (int i = 0; i < M; i++)
                    {
                        if (cfvalid.contains(cells[i]))
                        {
                            //
                            // Some part of the particle's mass lies in a 
                            // finer region; we'll deal with it shortly.
                            //
                            isFiner    = true;
                            isBoundary = true;
                            continue;
                        }
                        if ( ! fvalid.contains(cells[i]))
                        {
                            //
                            // We're out of the valid region.
                            //
                            isBoundary = true;
                            continue;
                        }
                        //
                        // Sum up mass in first component.
                        //
                        {
                            fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
                        }
                        //
                        // Sum up momenta in next components.
                        //

                        // If the domain is not periodic and we want to let particles
                        //    live near the boundary but "throw away" the contribution that 
                        //    does not fall into the domain ...
                        if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
			     ! gm.Domain().contains(cells[i]))
			{
			  continue;
			}

                        for (int n = 1; n < ncomp; n++) {
                            fab(cells[i],n) += p.m_rdata.arr[BL_SPACEDIM + n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
			}
                    }
                    //
                    // Deal with mass that doesn't belong at this level.
                    // Here we assume proper nesting so that only one special case can
                    // be true for a given particle.
                    //
                    if (isBoundary)
                    {
                        if (isFiner)
                        {
                            BL_ASSERT(lev < finest_level);
                            //
                            // We're at a coarse->fine interface
                            //
                            // get fine cells/fracs
                            //
                            const int MF = ParticleType::CIC_Cells_Fracs(p, plo, dx_fine ,dx, ffracs, fcells);

                            for (int j = 0; j < MF; j++)
                            {
                                //
                                // Make sure this fine cell is valid. Check for periodicity.
                                //
                                const Box bx(fcells[j],fcells[j]);
                                gm_fine.periodicShift(bx, gm_fine.Domain(), pshifts);
                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(int(pshifts.size()) == 1);
                                    fcells[j] = fcells[j] - pshifts[0];
                                }
                                mf[lev_index + 1]->boxArray().intersections(Box(fcells[j],fcells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first; 
                                const int who  = mf[lev_index+1]->DistributionMap()[grid];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
                                    {
                                        (*mf[lev_index+1])[grid](fcells[j],0) += p.m_rdata.arr[BL_SPACEDIM] * ffracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index+1])[grid](fcells[j],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * ffracs[j];
				    }
                                }
                                else
                                {

				  pb.m_lev  = lev+1;
				  pb.m_grid = grid;
				  pb.m_cell = fcells[j];

				  //
				  // Sum up mass in first component.
				  //
				  {
				    pb.m_data[0] = p.m_rdata.arr[BL_SPACEDIM] *  ffracs[j];
				  }
				  
				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * ffracs[j];
				  }
				  
				  data[who].push_back(pb);
                                }
                            }
                        }
                        else if (lev_index > 0)
                        {
                            //
                            // We must be at a fine->coarse interface.
                            //
                            const int MC = ParticleType::CIC_Cells_Fracs(p, plo, dx_coarse, dx, cfracs, ccells);
                            for (int j = 0; j < MC; j++)
                            {
                                //
                                // Make sure this coarse cell isn't in this level's valid region.
                                // This may not matter.
                                //
                                if (cba.contains(ccells[j]))
                                    continue;
                                //
                                // Check for periodicity.
                                //
                                const Box bx(ccells[j],ccells[j]);
                                gm_coarse.periodicShift(bx, gm_coarse.Domain(), pshifts);

                                if ( ! pshifts.empty())
                                {
                                    BL_ASSERT(int(pshifts.size()) == 1);
                                    ccells[j] = ccells[j] - pshifts[0]; 
                                }
                                //
                                // Find its resident grid.
                                //
                                mf[lev_index - 1]->boxArray().intersections(Box(ccells[j],ccells[j]),isects,true,0);
                                if (isects.size() == 0) {
                                    continue;
				}
                                const int grid = isects[0].first;
                                const int who  = mf[lev_index-1]->DistributionMap()[grid];
                                if (who == ParallelDescriptor::MyProc())
                                {
                                    //
                                    // Sum up mass in first component.
                                    //
                                    {
                                        (*mf[lev_index-1])[grid](ccells[j],0) += p.m_rdata.arr[BL_SPACEDIM] * cfracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[grid](ccells[j],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * cfracs[j];
				    }
                                }
                                else
                                {

				  pb.m_lev  = lev-1;
				  pb.m_grid = grid;
				  pb.m_cell = ccells[j];

                                  //
				  // Sum up mass in first component.
				  //
				  {
				    pb.m_data[0] = p.m_rdata.arr[BL_SPACEDIM] * cfracs[j];
				  }
                                  
				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * cfracs[j];
				  }
				  
				  data[who].push_back(pb);
                                }
                            }
                        }
                        else
                        {
                            // The mass is below levels we care about. Ignore it.
                        }
                    }
                }
                else 
                {
                    bool AnyCrseToFine = false;
                    if (lev < finest_level) {
                        AnyCrseToFine = ParticleType::CrseToFine(cfba,cells,cfshifts,gm,cwhich,pshifts);
		    }
                    //
                    // lev_index > 0 means that we don't do F->C for lower levels
                    // This may mean that the mass fraction is off.
                    //
                    bool AnyFineToCrse = false;
                    if (lev_index > 0 && !GridsCoverDomain)
                        AnyFineToCrse = ParticleType::FineToCrse(p,lev,m_gdb,cells,fvalid,compfvalid_grown,ccells,cfracs,fwhich,cgrid,pshifts,isects);

                    BL_ASSERT(!(AnyCrseToFine && AnyFineToCrse));

                    if ( ! AnyCrseToFine && ! AnyFineToCrse)
                    {
                        //
                        // By far the most common case.  Just do it!
                        //
                        for (int i = 0; i < M; i++)
                        {

                            // If the domain is not periodic and we want to let particles
                            //    live near the boundary but "throw away" the contribution that 
                            //    does not fall into the domain ...
                            if (! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i]))
			    {
			      continue;
			    }
                            //
                            // Sum up mass in first component.
                            //
                            {
                                fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
                            }
                            //
                            // Sum up momenta in next components.
                            //
                            for (int n = 1; n < ncomp; n++) {
                                fab(cells[i],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
			    }
                        }
                    }
                    else if (AnyFineToCrse)
                    {
                        Real sum_crse = 0, sum_fine = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // We're at a Fine->Crse boundary.
                                //
                                BL_ASSERT(cgrid[i] >= 0);
                                BL_ASSERT(cgrid[i] < mf[lev_index-1]->size());
                                //
                                // Here we need to update the crse region.  The coarse
                                // region is always going to be updated if we have a
                                // particle in a cell bordering a Fine->Crse boundary.
                                //
                                const int who = mf[lev_index-1]->DistributionMap()[cgrid[i]];

                                if (who == ParallelDescriptor::MyProc())
                                {
                                    if ( ! (*mf[lev_index-1])[cgrid[i]].box().contains(ccells[i])) {
				      continue;
				    }

                                    // If the domain is not periodic and we want to let particles
                                    //    live near the boundary but "throw away" the contribution that 
                                    //    does not fall into the domain ...
                                    if (! gm_coarse.isAllPeriodic() && allow_particles_near_boundary &&
				        ! gm_coarse.Domain().contains(ccells[i]))
				    {
				      continue;
				    }

                                    //
                                    // Sum up mass in first component.
                                    //
                                    {
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],0) += p.m_rdata.arr[BL_SPACEDIM] * cfracs[i];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * cfracs[i];
				    }
                                }
                                else
                                {
				  pb.m_lev  = lev-1;
				  pb.m_grid = cgrid[i];
				  pb.m_cell = ccells[i];

                                  //
				  // Sum up mass in first component.
				  //
				  {
				    pb.m_data[0] = p.m_rdata.arr[BL_SPACEDIM] * cfracs[i];
				  }

				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * cfracs[i];
				  }
				  data[who].push_back(pb);
                                }

                                sum_crse += cfracs[i];
                            }
                        }
                        //
                        // We've updated the Crse cells.  Now we have to update the fine
                        // cells in such a way that the total amount of mass we move
                        // around is precisely p.m_rdata.arr[BL_SPACEDIM]. In other words, the fractions
                        // we use at crse and fine have to sum to zero.  In the fine
                        // case, we have to account for the case where one or more of the
                        // cell indices is not in the valid region of the box containing 
                        // the particle.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {
                            //
                            // Reusing "fwhich" to indicate fine cells that need massaging.
                            //
                            fwhich[i] = true;

                            if ( ! compfvalid_grown.contains(cells[i]))
                            {
                                //
                                // Go ahead and add the full correct amount to these cells.
                                // They can't touch a Fine->Crse boundary.
                                //
                                sum_fine += fracs[i];
                                //
                                // Sum up mass in first component.
                                //
                                {
                                    fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
				}
                                fwhich[i] = false;
                            }
                            else if (compfvalid.contains(cells[i]))
                            {
                                fwhich[i] = false;
                            }
                        }

                        const Real sum_so_far = sum_crse + sum_fine; 

                        BL_ASSERT(sum_so_far > 0);
                        BL_ASSERT(sum_so_far < 1);

                        sum_fine = 0;
                        for (int i = 0; i < M; i++) 
                        {       
                            if (fwhich[i])
                                //
                                // Got to weight cells in this direction differently.
                                //
                                sum_fine += fracs[i];
                        }

                        const Real mult = (1 - sum_so_far) / sum_fine;
                        //
                        // Now add the weighted amount to the fine cells touching the c-f interface.
                        //
                        sum_fine = 0;
                        for (int i = 0; i < M; i++)
                        {
                            if (fwhich[i])
                            {
                                //
                                // Sum up mass in first component.
                                //
                                {
                                    fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM] * fracs[i] * mult;
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i] * mult;
				}

                                sum_fine += fracs[i] * mult;
                            }
                        }

                        BL_ASSERT(std::abs(1-(sum_fine+sum_so_far)) < 1.e-9);
                    }
                    else if (AnyCrseToFine)
                    {
                        Real sum = 0;

                        for (int i = 0; i < M; i++)
                        {
                            if (!cwhich[i])
                            {
                                // If the domain is not periodic and we want to let particles
                                //    live near the boundary but "throw away" the contribution that 
                                //    does not fall into the domain ...
                                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary &&
				     ! gm.Domain().contains(ccells[i]))
				{
				  continue;
				}
                                //
                                // Sum up mass in first component.
                                //
                                {
                                    fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
				}

                                sum += fracs[i];
                            }
                            else
                            {
                                //
                                // We're at a Crse->Fine boundary.
                                //
                                ParticleType::FineCellsToUpdateFromCrse(p,lev,m_gdb,cells[i],cfshifts[i],fgrid,ffracs,fcells,isects);

                                for (int j = 0, nfcells = fcells.size(); j < nfcells; j++)
                                {
                                    const int who = mf[lev_index+1]->DistributionMap()[fgrid[j]];

                                    if (who == ParallelDescriptor::MyProc())
                                    {
                                        //
                                        // Sum up mass in first component.
                                        //
                                        {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],0) += p.m_rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
                                        }
                                        //
                                        // Sum up momenta in next components.
                                        //
                                        for (int n = 1; n < ncomp; n++) {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
					}
                                    }
                                    else
                                    {
				      pb.m_lev  = lev+1;
				      pb.m_grid = fgrid[j];
				      pb.m_cell = fcells[j];

                                      //
				      // Sum up mass in first component.
				      //
				      {
					pb.m_data[0] = p.m_rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
				      }

				      //
				      // Sum up momenta in next components.
				      //
				      for (int n = 1; n < ncomp; n++) {
					pb.m_data[0] = p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i] * ffracs[j];
					}

                                        data[who].push_back(pb);
                                    }

                                    sum += fracs[i] * ffracs[j];
                                }
                            }
                        }

                        BL_ASSERT(std::abs(1-sum) < 1.e-9);
                    }
                }
            }
        }
    }

    //
    // Send any needed data to other MPI processes.
    // This "may" touch ghost cells so we want to do it before
    // the SumBoundary() stuff.
    //
    AssignDensityDoit(rho_index, mf, data, ncomp, lev_min);

    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const int       lev_index = lev - lev_min;
        const Geometry& gm        = m_gdb->Geom(lev);
        const Real*     dx        = gm.CellSize();
        const Real      vol       = D_TERM(dx[0], *dx[1], *dx[2]);

        mf[lev_index]->SumBoundary(gm.periodicity());
        //
        // If ncomp > 1, first divide the momenta (component n) 
        // by the mass (component 0) in order to get velocities.
        // Be careful not to divide by zero.
        //
        for (int n = 1; n < ncomp; n++)
        {
            for (MFIter mfi(*mf[lev_index]); mfi.isValid(); ++mfi)
            {
                (*mf[lev_index])[mfi].protected_divide((*mf[lev_index])[mfi],0,n,1);
            }
        }
        //
        // Only multiply the first component by (1/vol) because this converts mass
        // to density. If there are additional components (like velocity), we don't
        // want to divide those by volume.
        //
        mf[lev_index]->mult(1/vol,0,1);
    }

    //
    // The size of the returned multifab is limited by lev_min and 
    // finest_level. In the following code, lev is the real level,  
    // lev_index is the corresponding index for mf. 
    //
    // I believe that we don't need any information in ghost cells so we don't copy those.
    //
    if ( ! all_grids_the_same)
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_to_be_filled[lev_index]->copy(*mf_part[lev_index],0,0,1);
        }
    
    if (m_verbose > 1)
    {
        Real etime = ParallelDescriptor::second() - stime;

        ParallelDescriptor::ReduceRealMax(etime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR, NI>::AssignDensity(multi-level) time: " << etime << '\n';
        }
    }
}

//
// Used by AssignDensity (Array<std::unique_ptr<MultiFab> >& mf).
//
// Passes data needed by Crse->Fine or Fine->Crse to CPU that needs it.
//
// We store the data that needs to be sent in "data". Note that m_lev is the
// real particle level, while mf may start at a fine level (e.g. lvls 1 and 2).
// Consequently, we must subtract lev_min from m_lev to get the mf lev.
//

template <int NR, int NI>
void
ParticleContainer<NR, NI>::AssignDensityDoit (int               rho_index,
					       Array<std::unique_ptr<MultiFab> >&             mf,
					       std::map<int, std::vector<ParticleCommData> >& data,
					       int               ncomp,
					       int               lev_min) const
{
    if (rho_index != 0) amrex::Abort("AssignDensityDoit only works if rho_index = 0");

    BL_PROFILE("ParticleContainer<NR, NI>::AssignDensityDoit()");
    BL_ASSERT(NR >= ncomp);

    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1)
    {
      BL_ASSERT(data.empty());
      return;
    }

#if BL_USE_MPI
    //
    // We may have data that needs to be sent to another CPU.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    Array<int> Snds(NProcs,0), Rcvs(NProcs,0);

    int NumSnds = 0, NumRcvs = 0;

    for (const auto& kv : data)
    {
        NumSnds       += kv.second.size();
        Snds[kv.first] = kv.second.size();
    }

    ParallelDescriptor::ReduceIntMax(NumSnds);

    if (NumSnds == 0) {
        //
        // There's no parallel work to do.
        //
        return;
    }

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(int),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    typedef std::map<int,int> IntIntMap;

    IntIntMap SndCnts, RcvCnts, rOffset;

    for (int i = 0; i < NProcs; i++) {
        if (Snds[i] > 0) {
            SndCnts[i] = Snds[i];
	}
    }

    for (int i = 0; i < NProcs; i++)
    {
        if (Rcvs[i] > 0)
        {
            RcvCnts[i] = Rcvs[i];
            rOffset[i] = NumRcvs;
            NumRcvs   += Rcvs[i];
        }
    }
    //
    // Don't need these anymore.
    //
    Array<int>().swap(Snds);
    Array<int>().swap(Rcvs);
    //
    // The data we want to receive.
    //
    const int iChunkSize = 2 + BL_SPACEDIM;
    const int rChunkSize = ncomp;

    Array<int>                    irecvdata (NumRcvs*iChunkSize);
    Array<typename ParticleType::RealType> rrecvdata (NumRcvs*rChunkSize);

    Array<int>         index(2*RcvCnts.size());
    Array<MPI_Status>  stats(2*RcvCnts.size());
    Array<MPI_Request> rreqs(2*RcvCnts.size());

    const int SeqNumI = ParallelDescriptor::SeqNum();
    const int SeqNumR = ParallelDescriptor::SeqNum();
    //
    // Post the receives.
    //
    int idx = 0;
    for (auto it = RcvCnts.cbegin(); it != RcvCnts.cend(); ++it, ++idx)
    {
        const int Who  = it->first;
        const int iCnt = it->second   * iChunkSize;
        const int rCnt = it->second   * rChunkSize;
        const int iIdx = rOffset[Who] * iChunkSize;
        const int rIdx = rOffset[Who] * rChunkSize;

        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(iCnt > 0);
        BL_ASSERT(rCnt > 0);
        BL_ASSERT(iCnt < std::numeric_limits<int>::max());
        BL_ASSERT(rCnt < std::numeric_limits<int>::max());

        rreqs[2*idx+0] = ParallelDescriptor::Arecv(&irecvdata[iIdx],iCnt,Who,SeqNumI).req();
        rreqs[2*idx+1] = ParallelDescriptor::Arecv(&rrecvdata[rIdx],rCnt,Who,SeqNumR).req();
    }
    //
    // Send the data.
    //
    Array<int>                             isenddata;
    Array<typename ParticleType::RealType> rsenddata;

    for (const auto& kv : SndCnts)
    {
        const int Who  = kv.first;
        const int iCnt = kv.second * iChunkSize;
        const int rCnt = kv.second * rChunkSize;

        BL_ASSERT(iCnt > 0);
        BL_ASSERT(rCnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(iCnt < std::numeric_limits<int>::max());
        BL_ASSERT(rCnt < std::numeric_limits<int>::max());

        isenddata.resize(iCnt);
        rsenddata.resize(rCnt);

	std::vector<ParticleCommData>& pbox = data[Who];

        int ioff = 0, roff = 0;
	for (const auto& p : pbox)
        {
	  isenddata[ioff+0] = p.m_lev  - lev_min;
	  isenddata[ioff+1] = p.m_grid;

	  D_TERM(isenddata[ioff+2] = p.m_cell[0];,
		 isenddata[ioff+3] = p.m_cell[1];,
		 isenddata[ioff+4] = p.m_cell[2];);

	  ioff += iChunkSize;

	  for (int n = 0; n < ncomp; n++) {
	    rsenddata[roff+n] = p.m_data[n];
	  }

	  roff += ncomp;
        }

	std::vector<ParticleCommData>().swap(pbox);

        ParallelDescriptor::Send(isenddata.dataPtr(),iCnt,Who,SeqNumI);
        ParallelDescriptor::Send(rsenddata.dataPtr(),rCnt,Who,SeqNumR);
    }
    //
    // Receive the data.
    //
    for (int NWaits = rreqs.size(), completed; NWaits > 0; NWaits -= completed)
    {
        ParallelDescriptor::Waitsome(rreqs, completed, index, stats);
    }
    //
    // Now update "mf".
    //
    if (NumRcvs > 0)
    {
        const int*                             idata = irecvdata.dataPtr();
        const typename ParticleType::RealType* rdata = rrecvdata.dataPtr();

        for (int i = 0; i < NumRcvs; i++)
        {
            const int     lev  = idata[0];
            const int     grd  = idata[1];
            const IntVect cell (D_DECL(idata[2],idata[3],idata[4]));

            BL_ASSERT((*mf[lev]).DistributionMap()[grd] == MyProc);
	    BL_ASSERT((*mf[lev])[grd].box().contains(cell));

            for (int n = 0; n < ncomp; n++) {
                (*mf[lev])[grd](cell,n) += rdata[n];
	    }

            idata += iChunkSize;
            rdata += rChunkSize;
        }
    }

#endif /*BL_USE_MPI*/
}

//
// This is the single-level version -- it takes either cell-centered or node-centered MF's
//
template <int NR, int NI>
void
ParticleContainer<NR, NI>::AssignDensitySingleLevel (int rho_index,
						    MultiFab& mf_to_be_filled,
						    int       lev,
						    int       ncomp,
						    int       particle_lvl_offset) const
{
    BL_PROFILE("ParticleContainer<NR, NI>::AssignDensitySingleLevel()");
    BL_ASSERT(NR >= 1);
    BL_ASSERT(ncomp == 1 || ncomp == BL_SPACEDIM+1);

    if (lev >= int(m_particles.size()))
    {
        //
        // Don't do anything if there are no particles at this level.
        //
        return;
    }

    // Keep the same external interface to the applications, but if the
    if (mf_to_be_filled.is_nodal())
    {
        NodalDepositionSingleLevel(rho_index, mf_to_be_filled,lev,ncomp,particle_lvl_offset);
    }
    else if (mf_to_be_filled.boxArray().ixType().cellCentered())
    {
        AssignCellDensitySingleLevel(rho_index, mf_to_be_filled,lev,ncomp,particle_lvl_offset);
    }
    else
    {
	amrex::Abort("AssignCellDensitySingleLevel: mixed type not supported");
    }
}

//
// This is the single-level version for cell-centered density
//
template <int NR, int NI>
void
ParticleContainer<NR, NI>::AssignCellDensitySingleLevel (int rho_index,
							MultiFab& mf_to_be_filled,
							int       lev,
							int       ncomp,
							int       particle_lvl_offset) const
{
    if (rho_index != 0) amrex::Abort("AssignCellDensitySingleLevel only works if rho_index = 0");

    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer = &mf_to_be_filled;
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.
        mf_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev), 
				  m_gdb->ParticleDistributionMap(lev),
				  ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread its effect to an
    //    adjacent grid by first putting the value into ghost cells of its own grid.  The mf->sumBoundary call then
    //    adds the value from one grid's ghost cell to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       amrex::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = m_gdb->Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();
    const AoSMap&   pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
      amrex::Error("AssignDensity: problem must be periodic in no or all directions");
    }
    
    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
      (*mf_pointer)[mfi].setVal(0);
    }

    for (PIter it(*this, lev); it.isValid(); ++it) {
      int grid = it.index();
      int tile = it.tileIndex();
      const AoS& pbx = pmap.at(grid)[tile];
      FArrayBox& fab = (*mf_pointer)[grid];
      auto N = pbx.size();
	
        Array<Real>    fracs;
        Array<IntVect> cells;
	
#ifdef _OPENMP
#pragma omp parallel for default(none) private(fracs,cells) shared(N,plo,dx,dx_particle,gm,fab,ncomp,pbx)
#endif
	for (size_t ip = 0; ip < N; ++ip)
	  {
            const ParticleType& p = pbx[ip];
	    
            if (p.m_idata.id <= 0) {
	      continue;
	    }
	    
            const int M = ParticleType::CIC_Cells_Fracs(p, plo, dx, dx_particle, fracs, cells);
            //
            // If this is not fully periodic then we have to be careful that the
            // particle's support leaves the domain unless we specifically want to ignore
            // any contribution outside the boundary (i.e. if allow_particles_near_boundary = true). 
            // We test this by checking the low and high corners respectively.
            //
            if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
                if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
                    amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
		}
	    }

            for (int i = 0; i < M; i++)
            {
                if ( ! fab.box().contains(cells[i])) {
		  continue;
		}

                // If the domain is not periodic and we want to let particles
                //    live near the boundary but "throw away" the contribution that 
                //    does not fall into the domain ...
                if ( ! gm.isAllPeriodic() && allow_particles_near_boundary && ! gm.Domain().contains(cells[i])) {
		  continue;
		}
                //
                // Sum up mass in first component.
                //
                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                    fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
                }
                // 
                // Sum up momenta in next components.
                //
                for (int n = 1; n < ncomp; n++)
#ifdef _OPENMP
#pragma omp atomic
#endif
		  fab(cells[i],n) += p.m_rdata.arr[BL_SPACEDIM+n] * p.m_rdata.arr[BL_SPACEDIM] * fracs[i];
            }
        }
    }

    mf_pointer->SumBoundary(gm.periodicity());
    //
    // If ncomp > 1, first divide the momenta (component n) 
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    //
    for (int n = 1; n < ncomp; n++)
    {
        for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
        }
    }
    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
	delete mf_pointer;
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR, NI>::AssignDensity(single-level) time: " << stoptime << '\n';
        }
    }
}

//
// This is the single-level version for nodal density
//
template <int NR, int NI>
void
ParticleContainer<NR, NI>::NodalDepositionSingleLevel (int rho_index,
						      MultiFab& mf_to_be_filled,
					              int       lev,
           					      int       ncomp,
	           				      int       particle_lvl_offset) const
{
    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled))
    {
        // If we are already working with the internal mf defined on the 
        // particle_box_array, then we just work with this.
        mf_pointer = &mf_to_be_filled;
    }
    else
    {
        // If mf_to_be_filled is not defined on the particle_box_array, then we need 
        // to make a temporary here and copy into mf_to_be_filled at the end.
        mf_pointer = new MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
						  mf_to_be_filled.boxArray().ixType()),
				  m_gdb->ParticleDistributionMap(lev),
				  ncomp, mf_to_be_filled.nGrow());
    }

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = m_gdb->Geom(lev);
    const Real*     dx          = gm.CellSize();
    const AoSMap&   pmap        = m_particles[lev];
    const int       ngrids      = pmap.size();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) 
        amrex::Error("AssignDensity: problem must be periodic in no or all directions");

    mf_pointer->setVal(0.0);

    Array<IntVect> cells;
    cells.resize(8);

    Array<Real> fracs;
    fracs.resize(8);

#if (BL_SPACEDIM > 1)
    Array<Real> sx;
    sx.resize(2);
    Array<Real> sy;
    sy.resize(2);
#endif
#if (BL_SPACEDIM > 2)
    Array<Real> sz;
    sz.resize(2);
#endif

    for (PIter it(*this, lev); it.isValid(); ++it) {
      int grid = it.index();
      int tile = it.tileIndex();
      const AoS& pbx = pmap.at(grid)[tile];
      FArrayBox& fab = (*mf_pointer)[grid];

      for (const auto& p : pbx)
        {
            if (p.m_idata.id <= 0) {
	      continue;
	    }
	    
	    ParticleLocData pld;
	    Particle<NR, NI>::Where(p, m_gdb, pld);

#if (BL_SPACEDIM == 1)
            cells[0] = pld.m_cell;
            cells[1] = pld.m_cell+IntVect(1);

            Real x = p.m_rdata.pos[0] / dx[0];

            int i = pld.m_cell[0];

            Real xint = x - i;

            for (int i = 0; i < 2; i++)
            {
               fab(cells[0],0) += p.m_rdata.arr[BL_SPACEDIM+rho_index] * (1.0 - xint);
               fab(cells[1],0) += p.m_rdata.arr[BL_SPACEDIM+rho_index] *        xint ;
            }
#elif (BL_SPACEDIM == 2)
            cells[0] = pld.m_cell;
            cells[1] = pld.m_cell+IntVect(1,0);
            cells[2] = pld.m_cell+IntVect(0,1);
            cells[3] = pld.m_cell+IntVect(1,1);

            Real x = p.m_rdata.pos[0] / dx[0];
            Real y = p.m_rdata.pos[1] / dx[1];

            int i = pld.m_cell[0];
            int j = pld.m_cell[1];

            Real xint = x - i;
            Real yint = y - j;

            sx[0] = 1.0-xint;
            sx[1] = xint;
            sy[0] = 1.0-yint;
            sy[1] = yint;

            fracs[0] = sx[0] * sy[0];
            fracs[1] = sx[1] * sy[0];
            fracs[2] = sx[0] * sy[1];
            fracs[3] = sx[1] * sy[1];

            for (int i = 0; i < 4; i++)
            {
               fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM+rho_index] * fracs[i];
            }
#else
            cells[0] = pld.m_cell;
            cells[1] = pld.m_cell+IntVect(1,0,0);
            cells[2] = pld.m_cell+IntVect(0,1,0);
            cells[3] = pld.m_cell+IntVect(1,1,0);
            cells[4] = pld.m_cell+IntVect(0,0,1);
            cells[5] = pld.m_cell+IntVect(1,0,1);
            cells[6] = pld.m_cell+IntVect(0,1,1);
            cells[7] = pld.m_cell+IntVect(1,1,1);

            Real x = p.m_rdata.pos[0] / dx[0];
            Real y = p.m_rdata.pos[1] / dx[1];
            Real z = p.m_rdata.pos[2] / dx[2];

            int i = pld.m_cell[0];
            int j = pld.m_cell[1];
            int k = pld.m_cell[2];

            Real xint = x - i;
            Real yint = y - j;
            Real zint = z - k;

            sx[0] = 1.0-xint;
            sx[1] = xint;
            sy[0] = 1.0-yint;
            sy[1] = yint;
            sz[0] = 1.0-zint;
            sz[1] = zint;

            fracs[0] = sx[0] * sy[0] * sz[0];
            fracs[1] = sx[1] * sy[0] * sz[0];
            fracs[2] = sx[0] * sy[1] * sz[0];
            fracs[3] = sx[1] * sy[1] * sz[0];
            fracs[4] = sx[0] * sy[0] * sz[1];
            fracs[5] = sx[1] * sy[0] * sz[1];
            fracs[6] = sx[0] * sy[1] * sz[1];
            fracs[7] = sx[1] * sy[1] * sz[1];

            for (int i = 0; i < 8; i++)
            {
               fab(cells[i],0) += p.m_rdata.arr[BL_SPACEDIM+rho_index] * fracs[i];
            }
#endif
        }
    }

    mf_pointer->SumBoundary(gm.periodicity());

    //
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    //
    const Real vol = D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1/vol,0,1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
	delete mf_pointer;
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::NodalDepositionSingleLevel time: " << stoptime << '\n';
        }
    }
}

//
// This version takes as input the acceleration vector at cell centers, and has the option of
// returning the acceleration at the particle location in the data array, starting at
// component start_comp_for_accel
//
template <int NR, int NI>
void
ParticleContainer<NR, NI>::moveKick (MultiFab&       acceleration,
				      int             lev,
				      Real            dt,
				      Real            a_new,
				      Real            a_half, 
				      int             start_comp_for_accel)
{
    BL_PROFILE("ParticleContainer::moveKick()");
    BL_ASSERT(NR >= BL_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    const Real strttime  = ParallelDescriptor::second();
    const Real half_dt   = Real(0.5) * dt;
    const Real a_new_inv = 1 / a_new;
    AoSMap&    pmap      = m_particles[lev];

    MultiFab* ac_pointer;
    if (OnSameGrids(lev,acceleration))
    {
        ac_pointer = &acceleration;
    }
    else 
    {
        ac_pointer = new MultiFab(m_gdb->ParticleBoxArray(lev),
				  m_gdb->ParticleDistributionMap(lev),
				  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_pointer); mfi.isValid(); ++mfi)
            ac_pointer->setVal(0.);
        ac_pointer->copy(acceleration,0,0,acceleration.nComp());
        ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (PIter it(*this, lev); it.isValid(); ++it) {
      int grid = it.index();
      int tile = it.tileIndex();
      AoS& pbox = pmap[grid][tile];
      const int n = pbox.size();
      const FArrayBox& gfab = (*ac_pointer)[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < n; i++)
        {
	  ParticleLocData pld;
	  ParticleType& p = pbox[i];

	  if (p.m_idata.id > 0)
            {

	      ParticleType::Where(p, m_gdb, pld);

	      //
	      // Note: rdata.arr[BL_SPACEDIM] is mass, BL_SPACEDIM+1 is v_x, ...
	      //
	      Real grav[BL_SPACEDIM];

	      ParticleType::GetGravity(gfab, m_gdb->Geom(pld.m_lev), p, grav);
	      //
	      // Define (a u)^new = (a u)^half + dt/2 grav^new
	      //
	      D_TERM(p.m_rdata.arr[BL_SPACEDIM+1] *= a_half;,
		     p.m_rdata.arr[BL_SPACEDIM+2] *= a_half;,
		     p.m_rdata.arr[BL_SPACEDIM+3] *= a_half;);

	      D_TERM(p.m_rdata.arr[BL_SPACEDIM+1] += half_dt * grav[0];,
		     p.m_rdata.arr[BL_SPACEDIM+2] += half_dt * grav[1];,
		     p.m_rdata.arr[BL_SPACEDIM+3] += half_dt * grav[2];);

	      D_TERM(p.m_rdata.arr[BL_SPACEDIM+1] *= a_new_inv;,
		     p.m_rdata.arr[BL_SPACEDIM+2] *= a_new_inv;,
		     p.m_rdata.arr[BL_SPACEDIM+3] *= a_new_inv;);

	      if (start_comp_for_accel > BL_SPACEDIM)
                {
		  D_TERM(p.m_rdata.arr[BL_SPACEDIM + start_comp_for_accel  ] = grav[0];,
			 p.m_rdata.arr[BL_SPACEDIM + start_comp_for_accel+1] = grav[1];,
			 p.m_rdata.arr[BL_SPACEDIM + start_comp_for_accel+2] = grav[2];);
                }
            }
        }
    }

    
    if (ac_pointer != &acceleration) delete ac_pointer;

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NR, NI>::moveKick() time: " << stoptime << '\n';
        }
    }
    //
    // No need for Redistribution(), we only change the velocity.
    //
}



