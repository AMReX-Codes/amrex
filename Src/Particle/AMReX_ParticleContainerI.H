
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::do_tiling = false;

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
IntVect
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::tile_size { AMREX_D_DECL(1024000,8,8) };

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> :: SetParticleSize ()
{
    num_real_comm_comps = 0;
    for (int i = 0; i < NumRealComps(); ++i) {
        if (communicate_real_comp[i]) ++num_real_comm_comps;
    }

    num_int_comm_comps = 0;
    for (int i = 0; i < NumIntComps(); ++i) {
        if (communicate_int_comp[i]) ++num_int_comm_comps;
    }

    particle_size = sizeof(ParticleType);
    superparticle_size = particle_size + 
        num_real_comm_comps*sizeof(Real) + num_int_comm_comps*sizeof(int);    
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> :: Initialize ()
{
    levelDirectoriesCreated = false;
    usePrePost = false;
    doUnlink = true;

    SetParticleSize();

    static bool initialized = false;
    if ( ! initialized)
    {
        static_assert(sizeof(ParticleType)%sizeof(RealType) == 0,
                      "sizeof ParticleType is not a multiple of sizeof RealType");
        
        ParmParse pp("particles");
        pp.query("do_tiling", do_tiling);
        Vector<int> tilesize(AMREX_SPACEDIM);
        if (pp.queryarr("tile_size", tilesize, 0, AMREX_SPACEDIM)) {
            for (int i=0; i<AMREX_SPACEDIM; ++i) tile_size[i] = tilesize[i];
        }

        static_assert(std::is_standard_layout<ParticleType>::value
                   && std::is_trivial<ParticleType>::value,
                      "Particle type must be standard layout and trivial.");
        
        pp.query("use_prepost", usePrePost);
        pp.query("do_unlink", doUnlink);

        initialized = true;
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
IntVect
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Index (const ParticleType& p, int lev) const
{
    IntVect iv;
    const Geometry& geom = Geom(lev);

    AMREX_D_TERM(iv[0]=static_cast<int>(floor((p.m_rdata.pos[0]-geom.ProbLo(0))*geom.InvCellSize(0)));,
                 iv[1]=static_cast<int>(floor((p.m_rdata.pos[1]-geom.ProbLo(1))*geom.InvCellSize(1)));,
                 iv[2]=static_cast<int>(floor((p.m_rdata.pos[2]-geom.ProbLo(2))*geom.InvCellSize(2))););

    iv += geom.Domain().smallEnd();

    return iv;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Where (const ParticleType& p,
	 ParticleLocData&    pld,
	 int                 lev_min,
	 int                 lev_max,
	 int                 nGrow,
	 int                 local_grid) const
{

  BL_ASSERT(m_gdb != 0);
  
  if (lev_max == -1)
      lev_max = finestLevel();
  
  BL_ASSERT(lev_max <= finestLevel());

  BL_ASSERT(nGrow == 0 || (nGrow >= 0 && lev_min == lev_max));

  std::vector< std::pair<int, Box> > isects;

  for (int lev = lev_max; lev >= lev_min; lev--) {      
      const IntVect& iv = Index(p, lev);
      if (lev == pld.m_lev) {
          // The fact that we are here means this particle does not belong to any finer grids.
          if (pld.m_grid >= 0) {
              if (pld.m_grown_gridbox.contains(iv)) {
                  pld.m_cell = iv;
                  if (!pld.m_tilebox.contains(iv)) {
		    pld.m_tile = getTileIndex(iv, pld.m_gridbox, do_tiling, tile_size, pld.m_tilebox);
                  }
                  return true;
              }
          }
      }

      int grid;
      const BoxArray& ba = ParticleBoxArray(lev);
      BL_ASSERT(ba.ixType().cellCentered());

      if (local_grid < 0) {
          ba.intersections(Box(iv, iv), isects, true, nGrow);
          grid = isects.empty() ? -1 : isects[0].first;
      } else {
          grid = (*redistribute_mask_ptr)[local_grid](iv, 0);
      }

      if (grid >= 0) {
          const Box& bx = ba.getCellCenteredBox(grid);
	  pld.m_lev  = lev;
	  pld.m_grid = grid;
	  pld.m_tile = getTileIndex(iv, bx, do_tiling, tile_size, pld.m_tilebox);
	  pld.m_cell = iv;
	  pld.m_gridbox = bx;
          pld.m_grown_gridbox = amrex::grow(bx, nGrow);
	  return true;
      }
  }
  
  return false;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::EnforcePeriodicWhere (ParticleType&    p,
			ParticleLocData& pld,
			int              lev_min,
			int              lev_max,
			int              local_grid) const
{

    BL_ASSERT(m_gdb != 0);

    if (!Geom(0).isAnyPeriodic()) return false;

    if (lev_max == -1)
        lev_max = finestLevel();

    BL_ASSERT(lev_max <= finestLevel());

    // Create a copy "dummy" particle to check for periodic outs.
    ParticleType p_prime = p;
    if (PeriodicShift(p_prime)) {
        std::vector< std::pair<int,Box> > isects;
        for (int lev = lev_max; lev >= lev_min; lev--) {

	    int grid;
            IntVect iv;
            const BoxArray& ba = ParticleBoxArray(lev);
            BL_ASSERT(ba.ixType().cellCentered());
            
	    if (local_grid < 0) {
                iv = Index(p_prime, lev);
                ba.intersections(Box(iv, iv), isects, true, 0);
                grid = isects.empty() ? -1 : isects[0].first;
	    } else {
                iv = Index(p_prime, lev);
                if (ba[local_grid].contains(iv))
                {
                    grid = local_grid;
                }
                else
                {
                    grid = (*redistribute_mask_ptr)[local_grid](Index(p, lev), 0);
                }
	    }
            
            if (grid >= 0) {
                AMREX_D_TERM(p.m_rdata.pos[0] = p_prime.m_rdata.pos[0];,
                             p.m_rdata.pos[1] = p_prime.m_rdata.pos[1];,
                             p.m_rdata.pos[2] = p_prime.m_rdata.pos[2];);
                
                const Box& bx = ba.getCellCenteredBox(grid);
                
                pld.m_lev  = lev;
                pld.m_grid = grid;
		pld.m_tile = getTileIndex(iv, bx, do_tiling, tile_size, pld.m_tilebox);
                pld.m_cell = iv;
                pld.m_gridbox = bx;
                pld.m_grown_gridbox = bx;
                return true;
            }
        }
    }
    
    return false;
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::PeriodicShift (ParticleType& p) const
{
    BL_ASSERT(m_gdb != 0);

    const Geometry& geom    = Geom(0);
    const Box&      dmn     = geom.Domain();
    const IntVect&  iv      = Index(p, 0);    
    bool            shifted = false;  
    
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
        if (!geom.isPeriodic(i)) continue;

        if (iv[i] > dmn.bigEnd(i))
        {
            while (p.m_rdata.pos[i] >= geom.ProbHi(i)) p.m_rdata.pos[i] -= geom.ProbLength(i);
            if (p.m_rdata.pos[i] < geom.ProbLo(i)) p.m_rdata.pos[i] = geom.ProbLo(i); // clamp to avoid precision issues;
            shifted = true;
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            while (p.m_rdata.pos[i] <  geom.ProbLo(i)) p.m_rdata.pos[i] += geom.ProbLength(i);

            // clamp to avoid precision issues
            if ( p.m_rdata.pos[i] == geom.ProbHi(i)) p.m_rdata.pos[i] = geom.ProbLo(i);
            if ((p.m_rdata.pos[i] > geom.ProbHi(i))) 
                p.m_rdata.pos[i] = geom.ProbHi(i) - std::numeric_limits<typename ParticleType::RealType>::epsilon();

            shifted = true;
        }
        AMREX_ASSERT( (p.m_rdata.pos[i] >= geom.ProbLo(i) ) and ( p.m_rdata.pos[i] < geom.ProbHi(i) ));
    }

    return shifted;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
ParticleLocData
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
Reset (ParticleType& p,
       bool          update,
       bool          verbose,
       ParticleLocData pld) const
{
    BL_ASSERT(m_gdb != 0);

    bool ok = Where(p, pld);

    if (!ok && Geom(0).isAnyPeriodic())
    {
        // Attempt to shift the particle back into the domain if it
        // crossed a periodic boundary.
      PeriodicShift(p);
      ok = Where(p, pld);
    }
    
    if (!ok) {
        // invalidate the particle.
	if (verbose) {
            amrex::AllPrint()<< "Invalidating out-of-domain particle: " << p << '\n'; 
	}

	BL_ASSERT(p.m_idata.id > 0);

	p.m_idata.id = -p.m_idata.id;
    }

    return pld;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::reserveData ()
{
    int nlevs = maxLevel() + 1;
    m_particles.reserve(nlevs);
    m_dummy_mf.reserve(nlevs);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::resizeData ()
{
    int nlevs = std::max(0, finestLevel()+1);
    m_particles.resize(nlevs);
    m_dummy_mf.resize(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        RedefineDummyMF(lev);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::RedefineDummyMF (int lev) 
{
    if (lev > m_dummy_mf.size()-1) m_dummy_mf.resize(lev+1);
    
    if (m_dummy_mf[lev] == nullptr || 
        ! BoxArray::SameRefs(m_dummy_mf[lev]->boxArray(),
                             ParticleBoxArray(lev))          ||
        ! DistributionMapping::SameRefs(m_dummy_mf[lev]->DistributionMap(), 
                                        ParticleDistributionMap(lev)))
    {
        m_dummy_mf[lev].reset(new MultiFab(ParticleBoxArray(lev),
                                           ParticleDistributionMap(lev),
                                           1,0,MFInfo().SetAlloc(false)));
    };
}  

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::locateParticle (ParticleType& p, ParticleLocData& pld, 
                                                                                   int lev_min, int lev_max, int nGrow, int local_grid) const
{
    bool outside = AMREX_D_TERM(p.m_rdata.pos[0] <  Geom(0).ProbLo(0)
                             || p.m_rdata.pos[0] >= Geom(0).ProbHi(0),
                             || p.m_rdata.pos[1] <  Geom(0).ProbLo(1)
                             || p.m_rdata.pos[1] >= Geom(0).ProbHi(1),
                             || p.m_rdata.pos[2] <  Geom(0).ProbLo(2)
                             || p.m_rdata.pos[2] >= Geom(0).ProbHi(2));

    bool success;
    if (outside)
    {
      // Note that EnforcePeriodicWhere may shift the particle if it is successful.
      success = EnforcePeriodicWhere(p, pld, lev_min, lev_max, local_grid);
      if (!success && lev_min == 0)
      {
          // The particle has left the domain; invalidate it.
          p.m_idata.id = -p.m_idata.id;
          success = true;
      }
    }
    else
    {
        success = Where(p, pld, lev_min, lev_max, 0, local_grid);
    }

    if (!success)
    {
        success = (nGrow > 0) && Where(p, pld, lev_min, lev_min, nGrow);
        pld.m_grown_gridbox = pld.m_gridbox; // reset grown box for subsequent calls.
    }

    if (!success)
    {
        amrex::Abort("ParticleContainer::locateParticle(): invalid particle.");
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
long
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::TotalNumberOfParticles (bool only_valid, bool only_local) const
{
    long nparticles = 0;
    for (int lev = 0; lev <= finestLevel(); lev++) {
        nparticles += NumberOfParticlesAtLevel(lev,only_valid,true);
    }
    if (!only_local) {
	ParallelDescriptor::ReduceLongSum(nparticles);
    }
    return nparticles;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
Vector<long>
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::NumberOfParticlesInGrid (int lev, bool only_valid, bool only_local) const
{
  auto ngrids = ParticleBoxArray(lev).size();
  Vector<long> nparticles(ngrids, 0);

  if (lev >= 0 && lev < int(m_particles.size())) {
    for (const auto& kv : GetParticles(lev)) {
      int gid = kv.first.first;
      const auto& ptile = kv.second;
      
      if (only_valid) {
	for (int k = 0; k < ptile.GetArrayOfStructs().numParticles(); ++k) {
	  const ParticleType& p = ptile.GetArrayOfStructs()[k];
	  if (p.m_idata.id > 0) ++nparticles[gid];
	}
      } else {
	nparticles[gid] += ptile.numParticles();
      }
    }
    
    if (!only_local) ParallelDescriptor::ReduceLongSum(&nparticles[0],ngrids);
  }

  return nparticles;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
long
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::NumberOfParticlesAtLevel (int lev, bool only_valid, bool only_local) const
{
    long nparticles = 0;

    if (lev >= 0 && lev < int(m_particles.size())) {
        for (const auto& kv : GetParticles(lev)) {
            const auto& ptile = kv.second;	
            if (only_valid) {
                for (int k = 0; k < ptile.GetArrayOfStructs().numParticles(); ++k) {
                    const ParticleType& p = ptile.GetArrayOfStructs()[k];
                    if (p.m_idata.id > 0) ++nparticles;
                }
            } else {
                nparticles += ptile.numParticles();
            }
        }
    }

    if (!only_local) ParallelDescriptor::ReduceLongSum(nparticles);
    
    return nparticles;
}

//
// This includes both valid and invalid particles since invalid particles still take up space.
//

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::ByteSpread () const
{
    long cnt = 0;

    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& ptile = kv.second;
            cnt += ptile.numParticles();
        }
    }

    long mn = cnt, mx = mn;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    const std::size_t sz = sizeof(ParticleType)+NumRealComps()*sizeof(Real)+NumIntComps()*sizeof(int);

#ifdef BL_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
    ParallelDescriptor::ReduceLongMin(mn, IOProc);
    ParallelDescriptor::ReduceLongMax(mx, IOProc);
    ParallelDescriptor::ReduceLongSum(cnt,IOProc);

    amrex::Print() << "ParticleContainer byte spread across MPI nodes: ["
                   << mn*sz
                   << " (" << mn << ")"
                   << " ... "
                   << mx*sz
                   << " (" << mx << ")"
                   << "] total particles: (" << cnt << ")\n";
#ifdef BL_LAZY
    });
#endif
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::MoveRandom ()
{
    //
    // Move particles randomly at all levels
    //
    for (int lev = 0; lev < int(m_particles.size()); lev++)
    {
        MoveRandom(lev);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::MoveRandom (int lev)
{
    BL_PROFILE("ParticleContainer::MoveRandom(lev)");
    BL_ASSERT(OK());
    BL_ASSERT(m_gdb != 0);
    // 
    // Move particles up to FRAC*CellSize distance in each coordinate direction.
    //
    const Real FRAC = 0.25;
    auto&       pmap              = m_particles[lev];
    const Real* dx                = Geom(lev).CellSize();
    const Real  dist[AMREX_SPACEDIM] = { AMREX_D_DECL(FRAC*dx[0], FRAC*dx[1], FRAC*dx[2]) };

    for (auto& kv : pmap) {
        auto& aos = kv.second.GetArrayOfStructs();
        const int n = aos.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
	  ParticleType& p = aos[i];
	  
	  if (p.m_idata.id <= 0) continue;
	  
	  for (int j = 0; j < AMREX_SPACEDIM; j++)
              {
                  p.m_rdata.pos[j] += dist[j]*(2*amrex::Random()-1);
              }
	  
	  Reset(p, true);
        }
    }
    Redistribute();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Increment (MultiFab& mf, int lev) 
{
  IncrementWithTotal(mf,lev);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
long
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::IncrementWithTotal (MultiFab& mf, int lev, bool local)
{
  BL_PROFILE("ParticleContainer::IncrementWithTotal(lev)");
  BL_ASSERT(OK());
  
  if (m_particles.empty()) return 0;
  
  BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));
  
  AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

  const auto& pmap = m_particles[lev];
  
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
      mf_pointer = new MultiFab(ParticleBoxArray(lev),
				ParticleDistributionMap(lev),
				mf.nComp(),mf.nGrow());
  }
  
  ParticleLocData pld;
  for (auto& kv : pmap) {
      int gid = kv.first.first;
      const auto& pbox = kv.second.GetArrayOfStructs();
      FArrayBox&  fab  = (*mf_pointer)[gid];
      for (int k = 0; k < pbox.size(); ++ k) {
	const ParticleType& p = pbox[k];
          if (p.m_idata.id > 0) {
              Where(p, pld);
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
  
  if (!local) ParallelDescriptor::ReduceLongSum(num_particles_in_domain);
  
  return num_particles_in_domain;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
Real
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::sumParticleMass (int rho_index, int lev, bool local) const
{
  BL_PROFILE("ParticleContainer::sumParticleMass(lev)");
  BL_ASSERT(NStructReal >= 1);
  BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));
  
  Real msum = 0;
  
  const auto& pmap = m_particles[lev];
  for (const auto& kv : pmap) {
      const auto& pbox = kv.second.GetArrayOfStructs();
      for (int k = 0; k < pbox.size(); ++k) {
	  const ParticleType& p = pbox[k];
	  if (p.m_idata.id > 0) {
	      msum += p.m_rdata.arr[AMREX_SPACEDIM+rho_index];
	  }
      }
  }
  
  if (!local) ParallelDescriptor::ReduceRealSum(msum);
  
  return msum;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::RemoveParticlesAtLevel (int level)
{
    BL_PROFILE("ParticleContainer::RemoveParticlesAtLevel()");
    if (level >= int(this->m_particles.size())) return;
    
    if (!this->m_particles[level].empty())
    {
        ParticleLevel().swap(this->m_particles[level]);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::RemoveParticlesNotAtFinestLevel ()
{
  BL_PROFILE("ParticleContainer::RemoveParticlesNotAtFinestLevel()");
  BL_ASSERT(this->finestLevel()+1 == int(this->m_particles.size()));
  
  long cnt = 0;
  
  for (unsigned lev = 0; lev < m_particles.size() - 1; ++lev) {
      auto& pmap = m_particles[lev];
      if (!pmap.empty()) {
          for (auto& kv : pmap) {
              const auto& pbx = kv.second;
              cnt += pbx.size();
          }
          ParticleLevel().swap(pmap);
      }
  }
  
  //
  // Print how many particles removed on each processor if any were removed.
  //
  if (this->m_verbose > 1 && cnt > 0) {
      amrex::AllPrint() << "Processor " << ParallelDescriptor::MyProc() << " removed " << cnt
                        << " particles not in finest level\n";
  }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CreateVirtualParticles (int level, AoS& virts) const
{
    BL_PROFILE("ParticleContainer::CreateVirtualParticles()");
    BL_ASSERT(level > 0);
    BL_ASSERT(virts.empty());
    
    if (level >= static_cast<int>(m_particles.size()))
        return;

    if (aggregation_type == "")
    {
        ParmParse pp("particles");
        aggregation_type = "None";
        pp.query("aggregation_type", aggregation_type);
        aggregation_buffer = 2;
        pp.query("aggregation_buffer", aggregation_buffer);
    }

    if (aggregation_type == "None");
    else if (aggregation_type == "Cell");
    else if (aggregation_type == "Flow") amrex::Abort("Flow aggregation not implemented");
    else amrex::Abort("Unknown Particle Aggregation mode");
    
    if (aggregation_type == "None")
    { 
        const auto& pmap = m_particles[level];
        for (const auto& kv : pmap)
        {
            const auto& pbox = kv.second.GetArrayOfStructs();
            for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
            {
	        ParticleType p = *it;
	        p.m_idata.id = VirtualParticleID;
                virts.push_back(p);
            }
        }
        return;
    }


    if (aggregation_type == "Cell")
    {
        BoxList bl_buffer;
        bl_buffer.complementIn(Geom(level).Domain(), ParticleBoxArray(level));
        BoxArray buffer(std::move(bl_buffer));
        buffer.grow(aggregation_buffer);
        
        const auto& pmap = m_particles[level];
        for (const auto& kv : pmap)
        {
            const auto& pbox = kv.second.GetArrayOfStructs();
            
            std::map<IntVect,ParticleType> agg_map;
            
            for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
            {
                IntVect cell = Index(*it, level);
                if (buffer.contains(cell))
                {
                    // It's in the no-aggregation buffer.
                    // Set its id to indicate that it's a virt.
		    ParticleType p = *it;
                    p.m_idata.id = VirtualParticleID;
                    virts.push_back(p);
                }
                else
                {
                    //
                    // Note that Cell aggregation assumes that p.m_rdata.arr[AMREX_SPACEDIM] is mass and
                    // that all other components should be combined in a mass-weighted
                    // average.
                    //
                    auto agg_map_it = agg_map.find(cell);
                    
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
                        agg_map[cell] = p;
                    }
                    else
                    {
                        BL_ASSERT(agg_map_it != agg_map.end());
                        const ParticleType&  pnew       = *it;
                        ParticleType&        pold       = agg_map_it->second;
                        const Real           old_mass   = pold.m_rdata.arr[AMREX_SPACEDIM];
                        const Real           new_mass   = pnew.m_rdata.arr[AMREX_SPACEDIM];
                        const Real           total_mass = old_mass + new_mass;
                        //
                        // Set the position to the center of mass.
                        //
                        for (int i = 0; i < AMREX_SPACEDIM; i++)
                        {
                            pold.m_rdata.pos[i] = (old_mass*pold.m_rdata.pos[i] + new_mass*pnew.m_rdata.pos[i])/total_mass;
                        }
                        BL_ASSERT(this->Index(pold, level) == cell);
                        //
                        // Set the metadata (presumably velocity) to the mass-weighted average.
                        //
                        for (int i = AMREX_SPACEDIM + 1; i < AMREX_SPACEDIM + NStructReal; i++)
                        {
                            pold.m_rdata.arr[i] = (old_mass*pold.m_rdata.arr[i] + new_mass*pnew.m_rdata.arr[i])/total_mass;
                        }
                        pold.m_rdata.arr[AMREX_SPACEDIM] = total_mass;
                    }
                }
            }
            
            //
            // Add the aggregated particles to the virtuals.
            //
            for (const auto& agg_particle : agg_map)
            {
                virts.push_back(agg_particle.second);
            }
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CreateGhostParticles (int level, int nGrow, AoS& ghosts) const
{
    BL_PROFILE("ParticleContainer::CreateGhostParticles()");
    BL_ASSERT(ghosts.empty());
    BL_ASSERT(level < finestLevel());
  
    if (level >= static_cast<int>(m_particles.size()))
        return;
  
    const BoxArray& fine = ParticleBoxArray(level + 1);
  
    std::vector< std::pair<int,Box> > isects;
  
    const auto& pmap = m_particles[level];
    for (const auto& kv : pmap)
    {
        const auto& pbox = kv.second.GetArrayOfStructs();
        for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
        {
            const IntVect& iv = Index(*it, level+1);
            fine.intersections(Box(iv,iv),isects,false,nGrow);
            for (const auto& isec : isects)
            {
                amrex::ignore_unused(isec);
                ParticleType p = *it;  // yes, make a copy
                p.m_idata.id = GhostParticleID;	    
                ghosts().push_back(p);
            }
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
clearParticles()
{
    BL_PROFILE("ParticleContainer::clearParticles()");
    
    for (int lev = 0; lev < static_cast<int>(m_particles.size()); ++lev)
    {
        for (auto& kv : m_particles[lev]) 
        {
            kv.second.resize(0);
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
copyParticles(const ParticleContainerType& other, bool local)
{
    BL_PROFILE("ParticleContainer::copyParticles");
    clearParticles();   
    addParticles(other, local);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
addParticles(const ParticleContainerType& other, bool local)
{
    BL_PROFILE("ParticleContainer::addParticles");

    for (int lev = 0; lev < other.numLevels(); ++lev)
    {
        const auto& plevel_other = other.GetParticles(lev);
        auto& plevel = GetParticles(lev);
        for(MFIter mfi = other.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
            
            if(plevel_other.find(index) == plevel_other.end()) continue;
            
            const auto& tile_other = plevel_other.at(index);
            
            if (tile_other.numParticles() == 0) continue;
            
            const auto& aos_other = tile_other.GetArrayOfStructs();
            for (const auto& particle_struct : aos_other)
            {
                plevel[index].push_back(particle_struct);
            }

            const auto& soa_other = tile_other.GetStructOfArrays();
            for (int j = 0; j < NumRealComps(); ++j)
            {
                auto& rdata = soa_other.GetRealData(j);
                for(const auto& real_attrib : rdata)
                {
                    plevel[index].push_back_real(j, real_attrib);
                }                
            }
            for (int j = 0; j < NumIntComps(); ++j)
            {
                auto& idata = soa_other.GetIntData(j);
                for(const auto& int_attrib : idata)
                {
                    plevel[index].push_back_int(j, int_attrib);
                }
            }
        }
    }
    
    if (not local) Redistribute();
}

//
// This redistributes valid particles and discards invalid ones.
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Redistribute (int lev_min, int lev_max, int nGrow, int local)
{
#ifdef AMREX_USE_CUDA
    if ( (nGrow == 0) and Gpu::inLaunchRegion() )
    {
        RedistributeGPU(lev_min, lev_max, nGrow, local);
    }
    else
    {
        RedistributeCPU(lev_min, lev_max, nGrow, local);
    }
#else
    RedistributeCPU(lev_min, lev_max, nGrow, local);
#endif
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::SortParticlesByCell ()
{
#ifdef AMREX_USE_CUDA

    BL_PROFILE("ParticleContainer::SortParticlesByCell()");

    BuildRedistributeMask(0, 1);

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const BoxArray& ba   = ParticleBoxArray(lev);
    const DistributionMapping& dmap = ParticleDistributionMap(lev);
    auto& plev  = m_particles[lev];

    // temporaries
    Gpu::ManagedDeviceVector<amrex::IntVect> cells_tmp;
    Gpu::ManagedDeviceVector<int> index_sequence_tmp;
    ParticleVector aos_r;
    RealVector rdata_r;
    IntVector  idata_r;

    for(MFIter mfi(*redistribute_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();        
        auto& ptile = plev[std::make_pair(gid, tid)];
        auto& aos   = ptile.GetArrayOfStructs();
        auto& soa   = ptile.GetStructOfArrays();
        const size_t np = aos.numParticles();

        cells_tmp.resize(np);
        index_sequence_tmp.resize(np);

        thrust::sequence(thrust::device, index_sequence_tmp.begin(), index_sequence_tmp.end());

        thrust::transform(thrust::device, 
			  aos().begin(), aos().end(), cells_tmp.begin(),
                          functors::assignParticleCell(geom.data()));

        thrust::sort_by_key(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                            cells_tmp.begin(),
                            cells_tmp.end(),
                            index_sequence_tmp.begin());

        //
        // Reorder the particle data 
        //
        {
	    // reorder structs
	    aos_r.resize(np);
            thrust::gather(thrust::device,
                           index_sequence_tmp.begin(), index_sequence_tmp.end(),
                           aos().begin(), aos_r.begin());
            aos().swap(aos_r);

            // reorder real arrays
	    rdata_r.resize(np);
            for (int j = 0; j < NumRealComps(); ++j)
            {
                auto& rdata = ptile.GetStructOfArrays().GetRealData(j);
                thrust::gather(thrust::device,
                               index_sequence_tmp.begin(), index_sequence_tmp.end(),
                               rdata.begin(), rdata_r.begin());
                rdata.swap(rdata_r);
            }

            // reorder int arrays
	    idata_r.resize(np);
            for (int j = 0; j < NumIntComps(); ++j)
            {
                auto& idata = ptile.GetStructOfArrays().GetIntData(j);
                thrust::gather(thrust::device,
                               index_sequence_tmp.begin(), index_sequence_tmp.end(),
                               idata.begin(), idata_r.begin());
                idata.swap(idata_r);
            }
        }
    }
#endif
}

//
// The GPU implementation of Redistribute
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RedistributeGPU (int lev_min, int lev_max, int nGrow, int local)
{
#ifdef AMREX_USE_CUDA

    if (local) AMREX_ASSERT(numParticlesOutOfRange(*this, local) == 0);
    
    // sanity checks
    AMREX_ASSERT(nGrow == 0);
    AMREX_ASSERT(do_tiling == false);

    BL_PROFILE("ParticleContainer::RedistributeGPU()");
    BL_PROFILE_VAR_NS("Redistribute_partition", blp_partition);

    resizeData();

    if (lev_max < 0)
        lev_max = GetParGDB()->finestLevel();
    
    this->defineBufferMap();

    int num_levels = numLevels();
    if (num_levels == 1)
    {
        const auto& ba = ParticleBoxArray(0);
        const auto& geom = Geom(0);
        if ((ba.size() == 1) and geom.isAllPeriodic())
        {
            EnforcePeriodicGPU();
            AMREX_ASSERT(OK());
            return;
        }
    }

    if (not m_particle_locator.isValid(GetParGDB())) m_particle_locator.build(GetParGDB());
    m_particle_locator.setGeometry(GetParGDB());
    auto assign_grid = m_particle_locator.getGridAssignor();

    BL_PROFILE_VAR_START(blp_partition);
    auto& op = redistribute_copy_op;
    redistribute_copy_op.setNumLevels(num_levels);
    Vector<std::map<int, int> > new_sizes(num_levels);    
    int pid = ParallelDescriptor::MyProc();
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        const Geometry& geom = Geom(lev);
        const BoxArray& ba   = ParticleBoxArray(lev);
        const DistributionMapping& dmap = ParticleDistributionMap(lev);

        auto& plev = m_particles[lev];
        const auto plo = Geom(lev).ProbLoArray();
        const auto dxi = Geom(lev).InvCellSizeArray();
        const auto phi = Geom(lev).ProbHiArray();
        const auto is_per = Geom(lev).isPeriodicArray();
        const Box domain = Geom(lev).Domain();

        for (auto& kv : plev)
        {
            int gid = kv.first.first;
            int tid = kv.first.second;
            auto index = std::make_pair(gid, tid);
            
            auto& src_tile = plev[index];
            auto& aos = src_tile.GetArrayOfStructs();
            const size_t np = aos.numParticles();
            
            m_aos_r.resize(np);
            for (int i = 0; i < NArrayReal; ++i) m_rdata_r[i].resize(np);
            for (int i = 0; i < NArrayInt;  ++i) m_idata_r[i].resize(np);
            
            m_grids_tmp.resize(np);
            m_grids_r.resize(np);
            m_levs_r.resize(np);
            m_levs_tmp.resize(np);

            auto grids_ptr = m_grids_tmp.dataPtr();
            auto levs_ptr = m_levs_tmp.dataPtr();
            ParticleType* p_ptr = &(aos[0]);
            amrex::ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
            {
                auto& p = p_ptr[i];
                
                if (p.id() < 0) 
                {    
                    grids_ptr[i] = -1;
                    levs_ptr[i] = -1;
                }
                else
                {
                    enforcePeriodic(p, plo, phi, is_per);
                    auto tup = assign_grid(p);
                    grids_ptr[i] = amrex::get<0>(tup);
                    levs_ptr[i] = amrex::get<1>(tup);
                }
            });

            int num_stay;
            {
                auto p_lev_offsets = BufferMap().levelOffsetsPtr();
                auto p_box_perm = BufferMap().levGridToBucketPtr();
                auto p_pids = BufferMap().bucketToPIDPtr();

                auto particle_stays = [=] AMREX_GPU_DEVICE (int i) -> int
                {
                    return ((grids_ptr[i] == gid) && (levs_ptr[i] == lev)
                            && (p_pids[p_box_perm[p_lev_offsets[lev]+gid]] == pid));
                };
            
                const auto src_data = src_tile.getConstParticleTileData();
                typename ParticleTileType::ParticleTileDataType dst_data;
                dst_data.m_aos = m_aos_r.dataPtr();
                for (int i = 0; i < NArrayReal; ++i)
                    dst_data.m_rdata[i] = m_rdata_r[i].dataPtr();
                for (int i = 0; i < NArrayInt; ++i)
                    dst_data.m_idata[i] = m_idata_r[i].dataPtr();
                auto dst_grids_ptr = m_grids_r.dataPtr();
                auto dst_levs_ptr = m_levs_r.dataPtr();
                
                num_stay = Scan::PrefixSum<int> (np, 
                    [=] AMREX_GPU_DEVICE (int i) -> int 
                    {
                        return particle_stays(i);
                    },
                    [=] AMREX_GPU_DEVICE (int i, int const& s)
                    {
                        if (particle_stays(i)) {
                            dst_data.m_aos[s] = src_data.m_aos[i];
                            for (int j = 0; j < NArrayReal; ++j)
                                dst_data.m_rdata[j][s] = src_data.m_rdata[j][i];
                            for (int j = 0; j < NArrayInt; ++j)
                                dst_data.m_idata[j][s] = src_data.m_idata[j][i];                        
                            dst_grids_ptr[s] = grids_ptr[i];
                            dst_levs_ptr[s] = levs_ptr[i];
                        } else {
                            dst_data.m_aos[np-1-(i-s)] = src_data.m_aos[i];
                            for (int j = 0; j < NArrayReal; ++j)
                                dst_data.m_rdata[j][np-1-(i-s)] = src_data.m_rdata[j][i];
                            for (int j = 0; j < NArrayInt; ++j)
                                dst_data.m_idata[j][np-1-(i-s)] = src_data.m_idata[j][i];
                            dst_grids_ptr[np-1-(i-s)] = grids_ptr[i];
                            dst_levs_ptr[np-1-(i-s)] = levs_ptr[i];
                        }
                    }, 
                    Scan::Type::exclusive);

                aos().swap(m_aos_r);
                for (int i = 0; i < NArrayReal; ++i)
                {
                    auto& rdata = src_tile.GetStructOfArrays().GetRealData(i);
                    rdata.swap(m_rdata_r[i]);
                }
                for (int i = 0; i < NArrayInt; ++i)
                {
                    auto& idata = src_tile.GetStructOfArrays().GetIntData(i);
                    idata.swap(m_idata_r[i]);
                }
            }

            int num_move = np - num_stay;
            new_sizes[lev][gid] = num_stay;
            op.resize(gid, lev, num_move);

            auto p_boxes = op.m_boxes[lev][gid].dataPtr();
            auto p_levs = op.m_levels[lev][gid].dataPtr();
            auto p_src_indices = op.m_src_indices[lev][gid].dataPtr();
            auto p_periodic_shift = op.m_periodic_shift[lev][gid].dataPtr();
            auto dst_grids_ptr = m_grids_r.dataPtr();
            auto dst_levs_ptr  = m_levs_r.dataPtr();
            
	    AMREX_FOR_1D ( num_move, i,
            {
                p_boxes[i] = dst_grids_ptr[i + num_stay];
                p_levs[i] = dst_levs_ptr[i + num_stay];
                p_periodic_shift[i] = IntVect(AMREX_D_DECL(0,0,0));
                p_src_indices[i] = i+num_stay;
            });
        }
    }
    BL_PROFILE_VAR_STOP(blp_partition);

    auto& plan = redistribute_copy_plan;
    
    plan.build(*this, op, local);

    auto& snd_buffer = redistribute_snd_buffer;
    auto& rcv_buffer = redistribute_rcv_buffer;

    packBuffer(*this, op, plan, snd_buffer);

    // clear particles from container
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        auto& plev = m_particles[lev];
        for (auto& kv : plev)
        {
            int gid = kv.first.first;
            int tid = kv.first.second;
            auto index = std::make_pair(gid, tid);
            auto& tile = plev[index];
            tile.resize(new_sizes[lev][gid]);
        }
    }

    // Remove any map entries for which the particle container is now empty.
    for (int lev = lev_min; lev <= lev_max; lev++)
    {
        auto& pmap = m_particles[lev];
        for (auto pmap_it = pmap.begin(); pmap_it != pmap.end(); /* no ++ */)
        {          
            if (pmap_it->second.empty())
            {
                pmap.erase(pmap_it++);
            }
            else
            {
                ++pmap_it;
            }
        }
    }

    if (ParallelDescriptor::UseGpuAwareMpi())
    {
        plan.buildMPIFinish(BufferMap());
        communicateParticlesStart(*this, plan, snd_buffer, rcv_buffer);
        unpackBuffer(*this, plan, snd_buffer, RedistributeUnpackPolicy());
        communicateParticlesFinish(plan);
        unpackRemotes(*this, plan, rcv_buffer, RedistributeUnpackPolicy());        
    }
    else
    {
        pinned_snd_buffer.resize(snd_buffer.size());
        cudaMemcpyAsync(pinned_snd_buffer.dataPtr(), snd_buffer.dataPtr(), snd_buffer.size()*sizeof(SuperParticleType), cudaMemcpyDeviceToHost);    
        plan.buildMPIFinish(BufferMap());
        Gpu::Device::synchronize();
        communicateParticlesStart(*this, plan, pinned_snd_buffer, pinned_rcv_buffer);
        rcv_buffer.resize(pinned_rcv_buffer.size());   
        unpackBuffer(*this, plan, snd_buffer, RedistributeUnpackPolicy());
        communicateParticlesFinish(plan);        
        cudaMemcpyAsync(rcv_buffer.dataPtr(), pinned_rcv_buffer.dataPtr(), pinned_rcv_buffer.size()*sizeof(SuperParticleType), cudaMemcpyHostToDevice);    
        unpackRemotes(*this, plan, rcv_buffer, RedistributeUnpackPolicy());
    }

    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);
#endif    
}

#ifdef AMREX_USE_CUDA
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::EnforcePeriodicGPU ()
{
    BL_PROFILE("ParticleContainer::EnforcePeriodicGPU()");
    const int lev = 0;
    const Geometry& geom = Geom(lev);
    auto& plev = m_particles[lev];
    const auto plo = Geom(lev).ProbLoArray();
    const auto phi = Geom(lev).ProbHiArray();
    const auto is_per = Geom(lev).isPeriodicArray();

    for (auto& kv : plev)
    {
        int gid = kv.first.first;
        int tid = kv.first.second;
        auto index = std::make_pair(gid, tid);
        auto& particles = plev[index];

        const int np = particles.size();
        ParticleType* pstruct = &(particles.GetArrayOfStructs()[0]);
        AMREX_FOR_1D ( np, i,
        {
            enforcePeriodic(pstruct[i], plo, phi, is_per);
        });
    }
}
#endif //AMREX_USE_CUDA

//
// The CPU implementation of Redistribute
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RedistributeCPU (int lev_min, int lev_max, int nGrow, int local)
{
  BL_PROFILE("ParticleContainer::RedistributeCPU()");
    
  const int MyProc    = ParallelDescriptor::MyProc();
  Real      strttime  = amrex::second();
  
  if (local > 0) BuildRedistributeMask(0, local);

  // On startup there are cases where Redistribute() could be called
  // with a given finestLevel() where that AmrLevel has yet to be defined.
  int theEffectiveFinestLevel = m_gdb->finestLevel();

  while (!m_gdb->LevelDefined(theEffectiveFinestLevel))
      theEffectiveFinestLevel--;

  if (int(m_particles.size()) < theEffectiveFinestLevel+1) {
      if (Verbose()) {
          amrex::Print() << "ParticleContainer::Redistribute() resizing containers from "
                         << m_particles.size() << " to " 
                         << theEffectiveFinestLevel + 1 << '\n';
      }
      m_particles.resize(theEffectiveFinestLevel+1);
      m_dummy_mf.resize(theEffectiveFinestLevel+1);
  }
  
  // It is important to do this even if we don't have more levels because we may have changed the 
  // grids at this level in a regrid.
  for (int lev = 0; lev < theEffectiveFinestLevel+1; ++lev)
      RedefineDummyMF(lev);
  
  int nlevs_particles;
  if (lev_max == -1) {
      lev_max = theEffectiveFinestLevel;
      nlevs_particles = m_particles.size() - 1; 
  } else {
      nlevs_particles = lev_max;
  }
  BL_ASSERT(lev_max <= finestLevel());

  // This will hold the valid particles that go to another process
  std::map<int, Vector<char> > not_ours;
  
  int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single
  num_threads = omp_get_num_threads();
#endif
  
  // these are temporary buffers for each thread
  std::map<int, Vector<Vector<char> > > tmp_remote;
  Vector<std::map<std::pair<int, int>, Vector<ParticleVector> > > tmp_local;
  Vector<std::map<std::pair<int, int>, Vector<StructOfArrays<NArrayReal, NArrayInt> > > > soa_local;
  tmp_local.resize(theEffectiveFinestLevel+1);
  soa_local.resize(theEffectiveFinestLevel+1);

  // we resize these buffers outside the parallel region
  for (int lev = lev_min; lev <= lev_max; lev++) {
      for (MFIter mfi(*m_dummy_mf[lev], this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
	   mfi.isValid(); ++mfi) {
          auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
          tmp_local[lev][index].resize(num_threads);
          soa_local[lev][index].resize(num_threads);
          for (int t = 0; t < num_threads; ++t) {
              soa_local[lev][index][t].define(m_num_runtime_real, m_num_runtime_int);
          }
      }
  }
  if (local) {
    for (int i = 0; i < neighbor_procs.size(); ++i)
      tmp_remote[neighbor_procs[i]].resize(num_threads);
  } else {
    for (int i = 0; i < ParallelDescriptor::NProcs(); ++i)
      tmp_remote[i].resize(num_threads);
  }

  // first pass: for each tile in parallel, in each thread copies the particles that
  // need to be moved into it's own, temporary buffer.
  for (int lev = lev_min; lev <= nlevs_particles; lev++) {
      auto& pmap = m_particles[lev];

      Vector<std::pair<int, int> > grid_tile_ids;
      Vector<ParticleTileType*> ptile_ptrs;
      for (auto& kv : pmap)
      {
          grid_tile_ids.push_back(kv.first);
          ptile_ptrs.push_back(&(kv.second));
      }
      
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int pmap_it = 0; pmap_it < static_cast<int>(ptile_ptrs.size()); ++pmap_it)
      {
#ifdef _OPENMP
          int thread_num = omp_get_thread_num();
#else
          int thread_num = 0;
#endif
          int grid = grid_tile_ids[pmap_it].first;
          int tile = grid_tile_ids[pmap_it].second;
          auto& aos = ptile_ptrs[pmap_it]->GetArrayOfStructs();
          auto& soa = ptile_ptrs[pmap_it]->GetStructOfArrays();
          unsigned npart = aos.numParticles();              
          ParticleLocData pld;
          if (npart != 0) {
              long last = npart - 1;
              unsigned pindex = 0;
              while (pindex <= last) {
                  ParticleType& p = aos[pindex];

                  if (p.m_idata.id < 0)
		  {
                      aos[pindex] = aos[last];
                      for (int comp = 0; comp < NumRealComps(); comp++)
                          soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                      for (int comp = 0; comp < NumIntComps(); comp++)
                          soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                      correctCellVectors(last, pindex, grid, aos[pindex]);
                      --last;
                      continue;
                  }
                      
                  locateParticle(p, pld, lev_min, lev_max, nGrow, local ? grid : -1);

                  particlePostLocate(p, pld, lev);
                  
                  if (p.m_idata.id < 0)
                  {
                      aos[pindex] = aos[last];
                      for (int comp = 0; comp < NumRealComps(); comp++)
                          soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                      for (int comp = 0; comp < NumIntComps(); comp++)
                          soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                      correctCellVectors(last, pindex, grid, aos[pindex]);
                      --last;
                      continue;
                  }

                  const int who = ParticleDistributionMap(pld.m_lev)[pld.m_grid];
                  if (who == MyProc) {
                      if (pld.m_lev != lev || pld.m_grid != grid || pld.m_tile != tile) {
                          // We own it but must shift it to another place.
                          auto index = std::make_pair(pld.m_grid, pld.m_tile);
                          BL_ASSERT(tmp_local[pld.m_lev][index].size() == num_threads);
                          tmp_local[pld.m_lev][index][thread_num].push_back(p);
                          for (int comp = 0; comp < NumRealComps(); ++comp) {
                              RealVector& arr = soa_local[pld.m_lev][index][thread_num].GetRealData(comp);
                              arr.push_back(soa.GetRealData(comp)[pindex]);
                          }
                          for (int comp = 0; comp < NumIntComps(); ++comp) {
                              IntVector& arr = soa_local[pld.m_lev][index][thread_num].GetIntData(comp);
                              arr.push_back(soa.GetIntData(comp)[pindex]);
                          }
                          
                          p.m_idata.id = -p.m_idata.id; // Invalidate the particle
                      }
                  }
                  else {
                      auto& particles_to_send = tmp_remote[who][thread_num];
                      auto old_size = particles_to_send.size();
                      auto new_size = old_size + superparticle_size;
                      particles_to_send.resize(new_size);
                      std::memcpy(&particles_to_send[old_size], &p, particle_size);
                      char* dst = &particles_to_send[old_size] + particle_size;
                      for (int comp = 0; comp < NumRealComps(); comp++) {
                          if (communicate_real_comp[comp]) {
                              std::memcpy(dst, &soa.GetRealData(comp)[pindex], sizeof(Real));
                              dst += sizeof(Real);
                          }
                      }
                      for (int comp = 0; comp < NumIntComps(); comp++) {
                          if (communicate_int_comp[comp]) {
			      std::memcpy(dst, &soa.GetIntData(comp)[pindex], sizeof(int));
                              dst += sizeof(int);
                          }
                      }
                      
                      p.m_idata.id = -p.m_idata.id; // Invalidate the particle
                  }
                  
                  if (p.m_idata.id < 0)
		  {
                      aos[pindex] = aos[last];
                      for (int comp = 0; comp < NumRealComps(); comp++)
                          soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                      for (int comp = 0; comp < NumIntComps(); comp++)
                          soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                      correctCellVectors(last, pindex, grid, aos[pindex]);
                      --last;
                      continue;
                  }
                  
                  ++pindex;
              }
              
              aos().erase(aos().begin() + last + 1, aos().begin() + npart);
              for (int comp = 0; comp < NumRealComps(); comp++) {
                  RealVector& rdata = soa.GetRealData(comp);
                  rdata.erase(rdata.begin() + last + 1, rdata.begin() + npart);
              }
              for (int comp = 0; comp < NumIntComps(); comp++) {
                  IntVector& idata = soa.GetIntData(comp);
                  idata.erase(idata.begin() + last + 1, idata.begin() + npart);
              }
          }
      }
  }
  
  for (int lev = lev_min; lev <= lev_max; lev++) {
      auto& pmap = m_particles[lev];
      for (auto pmap_it = pmap.begin(); pmap_it != pmap.end(); /* no ++ */) {
          
          // Remove any map entries for which the particle container is now empty.
          if (pmap_it->second.empty()) {
              pmap.erase(pmap_it++);
          }
          else {
              ++pmap_it;
          }
      }
  }
  
  // Second pass - for each tile in parallel, collect the particles we are owed from all thread's buffers.
  for (int lev = lev_min; lev <= lev_max; lev++) {
      typename std::map<std::pair<int, int>, Vector<ParticleVector > >::iterator pmap_it;
      
      Vector<std::pair<int, int> > grid_tile_ids;
      Vector<Vector<ParticleVector>* > pvec_ptrs;

      // we need to create any missing map entries in serial here
      for (pmap_it=tmp_local[lev].begin(); pmap_it != tmp_local[lev].end(); pmap_it++)
      {
          m_particles[lev][pmap_it->first];
          grid_tile_ids.push_back(pmap_it->first);
          pvec_ptrs.push_back(&(pmap_it->second));
      }

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int pit = 0; pit < static_cast<int>(pvec_ptrs.size()); ++pit)
      {
          auto index = grid_tile_ids[pit];
          auto& ptile = DefineAndReturnParticleTile(lev, index.first, index.second);
          auto& aos = ptile.GetArrayOfStructs();
          auto& soa = ptile.GetStructOfArrays();
          auto& aos_tmp = *(pvec_ptrs[pit]);
          auto& soa_tmp = soa_local[lev][index];
          for (int i = 0; i < num_threads; ++i) {
              aos.insert(aos.end(), aos_tmp[i].begin(), aos_tmp[i].end());
              aos_tmp[i].erase(aos_tmp[i].begin(), aos_tmp[i].end());
              for (int comp = 0; comp < NumRealComps(); ++comp) {
                  RealVector& arr = soa.GetRealData(comp);
                  RealVector& tmp = soa_tmp[i].GetRealData(comp);
                  arr.insert(arr.end(), tmp.begin(), tmp.end());
                  tmp.erase(tmp.begin(), tmp.end());
              }
              for (int comp = 0; comp < NumIntComps(); ++comp) {
                  IntVector& arr = soa.GetIntData(comp);
                  IntVector& tmp = soa_tmp[i].GetIntData(comp);
                  arr.insert(arr.end(), tmp.begin(), tmp.end());
                  tmp.erase(tmp.begin(), tmp.end());
              }
          }
      }
  }

  for (auto& map_it : tmp_remote) {
      int who = map_it.first;
      not_ours[who];
  }

  Vector<int> dest_proc_ids;
  Vector<Vector<Vector<char> >* > pbuff_ptrs;
  for (auto& kv : tmp_remote)
  {
      dest_proc_ids.push_back(kv.first);
      pbuff_ptrs.push_back(&(kv.second));
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (int pmap_it = 0; pmap_it < static_cast<int>(pbuff_ptrs.size()); ++pmap_it)
  {
      int who = dest_proc_ids[pmap_it];
      Vector<Vector<char> >& tmp = *(pbuff_ptrs[pmap_it]);
      for (int i = 0; i < num_threads; ++i) {
          not_ours[who].insert(not_ours[who].end(), tmp[i].begin(), tmp[i].end());
          tmp[i].erase(tmp[i].begin(), tmp[i].end());
      }
  }

  // remove any empty map entries from not_ours
  for (auto pmap_it = not_ours.begin(); pmap_it != not_ours.end(); /* no ++ */) {
      if (pmap_it->second.empty()) {
          not_ours.erase(pmap_it++);
      }
      else {
          ++pmap_it;
      }
  }

  if (int(m_particles.size()) > theEffectiveFinestLevel+1) {
      // Looks like we lost an AmrLevel on a regrid.
      if (m_verbose > 0) {
          amrex::Print() << "ParticleContainer::Redistribute() resizing m_particles from "
                         << m_particles.size() << " to " << theEffectiveFinestLevel+1 << '\n';
      }
      BL_ASSERT(int(m_particles.size()) >= 2);
      
      m_particles.resize(theEffectiveFinestLevel + 1);
      m_dummy_mf.resize(theEffectiveFinestLevel + 1);
  }
  
  if (ParallelDescriptor::NProcs() == 1) {
      BL_ASSERT(not_ours.empty());
  }
  else {
      RedistributeMPI(not_ours, lev_min, lev_max, nGrow, local);
  }
  
  BL_ASSERT(OK(lev_min, lev_max, nGrow));
  
  if (m_verbose > 0) {
      Real stoptime = amrex::second() - strttime;
      
      ByteSpread();
      
#ifdef BL_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
              ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
              amrex::Print() << "ParticleContainer::Redistribute() time: " << stoptime << "\n\n";
#ifdef BL_LAZY
          });
#endif
  }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
defineBufferMap () const
{    
    BL_PROFILE("ParticleContainer::defineBufferMap");
    const int lev = 0;
    const BoxArray& ba = ParticleBoxArray(lev);
    const DistributionMapping& dm = ParticleDistributionMap(lev);

    if (not m_buffer_map.isValid(GetParGDB()))
    {
        m_buffer_map.define(GetParGDB());
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
BuildRedistributeMask (int lev, int nghost) const
{    
    BL_PROFILE("ParticleContainer::BuildRedistributeMask");
    BL_ASSERT(lev == 0);

    if (redistribute_mask_ptr == nullptr ||
        redistribute_mask_nghost < nghost ||
        ! BoxArray::SameRefs(redistribute_mask_ptr->boxArray(), this->ParticleBoxArray(lev)) ||
        ! DistributionMapping::SameRefs(redistribute_mask_ptr->DistributionMap(), this->ParticleDistributionMap(lev)))
    {
        const Geometry& geom = this->Geom(lev);
        const BoxArray& ba = this->ParticleBoxArray(lev);
        const DistributionMapping& dmap = this->ParticleDistributionMap(lev);
        
        redistribute_mask_nghost = nghost;
        redistribute_mask_ptr.reset(new iMultiFab(ba, dmap, 2, nghost));
        redistribute_mask_ptr->setVal(-1, nghost);

        auto tile_size = this->do_tiling ? this->tile_size : IntVect::TheZeroVector();            

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*redistribute_mask_ptr, tile_size); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            (*redistribute_mask_ptr)[mfi].setVal(grid_id, box, 0, 1);
            (*redistribute_mask_ptr)[mfi].setVal(tile_id, box, 1, 1);
        }
        
        redistribute_mask_ptr->FillBoundary(geom.periodicity());

        neighbor_procs.clear();
        for (MFIter mfi(*redistribute_mask_ptr, tile_size); mfi.isValid(); ++mfi) 
        {
            const Box& box = mfi.growntilebox();
            for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv))
            {
                const int grid = (*redistribute_mask_ptr)[mfi](iv, 0);
                if (grid >= 0)
                {
                    const int proc = this->ParticleDistributionMap(lev)[grid];
                    if (proc != ParallelDescriptor::MyProc())  
                        neighbor_procs.push_back(proc);
                }
            }
        }
        RemoveDuplicates(neighbor_procs);
    }    
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
RedistributeMPI (std::map<int, Vector<char> >& not_ours,
                 int lev_min, int lev_max, int nGrow, int local)
{
    BL_PROFILE("ParticleContainer::RedistributeMPI()");
    BL_PROFILE_VAR_NS("RedistributeMPI_locate", blp_locate);
    BL_PROFILE_VAR_NS("RedistributeMPI_copy", blp_copy);

#ifdef BL_USE_MPI

    using buffer_type = unsigned long long;
    
    std::map<int, Vector<buffer_type> > mpi_snd_data;
    for (const auto& kv : not_ours)
    {
        int nbt = (kv.second.size() + sizeof(buffer_type)-1)/sizeof(buffer_type);
        mpi_snd_data[kv.first].resize(nbt);
        std::memcpy((char*) mpi_snd_data[kv.first].data(), kv.second.data(), kv.second.size());
    }
    
    const int NProcs = ParallelDescriptor::NProcs();
    const int NNeighborProcs = neighbor_procs.size();
    
    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;
    if (local > 0)
    {
        AMREX_ALWAYS_ASSERT(lev_min == 0);
        AMREX_ALWAYS_ASSERT(lev_max == 0);
        BuildRedistributeMask(0, local);
        NumSnds = doHandShakeLocal(not_ours, neighbor_procs, Snds, Rcvs);
    }
    else
    {
        NumSnds = doHandShake(not_ours, Snds, Rcvs);
    }

    const int SeqNum = ParallelDescriptor::SeqNum();
    
    if ((not local) and NumSnds == 0)
        return;  // There's no parallel work to do.

    if (local)
    {
        long tot_snds_this_proc = 0;
        long tot_rcvs_this_proc = 0;
        for (int i = 0; i < NNeighborProcs; ++i) {
            tot_snds_this_proc += Snds[neighbor_procs[i]];
            tot_rcvs_this_proc += Rcvs[neighbor_procs[i]];
        }
        if ( (tot_snds_this_proc == 0) and (tot_rcvs_this_proc == 0) ) {
            return; // There's no parallel work to do.
        } 
    }

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    
    std::size_t TotRcvInts = 0;
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvInts);
            TotRcvBytes += Rcvs[i];
            int nbt = (Rcvs[i] + sizeof(buffer_type)-1)/sizeof(buffer_type);
            TotRcvInts += nbt;
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);
    
    // Allocate data for rcvs as one big chunk.
    Vector<unsigned long long> recvdata(TotRcvInts);
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt = (Rcvs[Who] + sizeof(buffer_type)-1)/sizeof(buffer_type);
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (const auto& kv : mpi_snd_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum);
    }
    
    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);
     
	BL_PROFILE_VAR_START(blp_locate);
   
        int npart = TotRcvBytes / superparticle_size;
        
        Vector<int> rcv_levs(npart);
        Vector<int> rcv_grid(npart);
        Vector<int> rcv_tile(npart);

        int ipart = 0;
        ParticleLocData pld;
        for (int j = 0; j < nrcvs; ++j)
        {
            const auto offset = rOffset[j];
            const auto Who    = RcvProc[j];
            const auto Cnt    = Rcvs[Who] / superparticle_size;            
            for (int i = 0; i < Cnt; ++i)
            {
                char* pbuf = ((char*) &recvdata[offset]) + i*superparticle_size;
                ParticleType p;
                std::memcpy(&p, pbuf, sizeof(ParticleType));
                locateParticle(p, pld, lev_min, lev_max, nGrow);
                rcv_levs[ipart] = pld.m_lev;
                rcv_grid[ipart] = pld.m_grid;
                rcv_tile[ipart] = pld.m_tile;
                ++ipart;
            }
        }

	BL_PROFILE_VAR_STOP(blp_locate);

        BL_PROFILE_VAR_START(blp_copy);

#ifndef AMREX_USE_CUDA
        ipart = 0;
        for (int i = 0; i < nrcvs; ++i)
        {
            const auto offset = rOffset[i];
            const auto Who    = RcvProc[i];
            const auto Cnt = Rcvs[Who] / superparticle_size;            
            for (int j = 0; j < Cnt; ++j)
            {                
                auto& ptile = m_particles[rcv_levs[ipart]][std::make_pair(rcv_grid[ipart],
                                                                          rcv_tile[ipart])];
                char* pbuf = ((char*) &recvdata[offset]) + j*superparticle_size;

                ParticleType p;
                std::memcpy(&p, pbuf, sizeof(ParticleType));            
                ptile.push_back(p);
                Real* rdata = (Real*)(pbuf + particle_size);
                for (int comp = 0; comp < NumRealComps(); ++comp) {
                    if (communicate_real_comp[comp]) {
                        ptile.push_back_real(comp, *rdata++);
                    } else {
                        ptile.push_back_real(comp, 0.0);
                    }
                }
            
                int* idata = (int*)(pbuf + particle_size + num_real_comm_comps*sizeof(Real));
                for (int comp = 0; comp < NumIntComps(); ++comp) {
                    if (communicate_int_comp[comp]) {
                        ptile.push_back_int(comp, *idata++);
                    } else {
                        ptile.push_back_int(comp, 0);
                    }
                }
                ++ipart;
            }
        }	    
#else
	Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
	host_particles.reserve(15);
	host_particles.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>,
			std::vector<Cuda::HostVector<Real> > > > host_real_attribs;
	host_real_attribs.reserve(15);
	host_real_attribs.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>,
			std::vector<Cuda::HostVector<int> > > > host_int_attribs;
	host_int_attribs.reserve(15);
	host_int_attribs.resize(finestLevel()+1);

        ipart = 0;
        for (int i = 0; i < nrcvs; ++i)
        {
            const auto offset = rOffset[i];
            const auto Who    = RcvProc[i];
            const auto Cnt = Rcvs[Who] / superparticle_size;            
            for (int j = 0; j < Cnt; ++j)
            {                
                int lev = rcv_levs[ipart];
                std::pair<int, int> ind(std::make_pair(rcv_grid[ipart], rcv_tile[ipart]));

                char* pbuf = ((char*) &recvdata[offset]) + j*superparticle_size;

                ParticleType p;
                std::memcpy(&p, pbuf, sizeof(ParticleType));
	  
                host_real_attribs[lev][ind].resize(NumRealComps());
                host_int_attribs[lev][ind].resize(NumIntComps());
	  
                // add the struct
                host_particles[lev][ind].push_back(p);
	  
                Real* rdata = (Real*)(pbuf + particle_size);
                // add the real...
                for (int comp = 0; comp < NumRealComps(); ++comp) {
                    if (communicate_real_comp[comp]) {
                        host_real_attribs[lev][ind][comp].push_back(*rdata++);
                    } else {
                        host_real_attribs[lev][ind][comp].push_back(0.0);
                    }                    
                }
                
                int* idata = (int*)(pbuf + particle_size + num_real_comm_comps*sizeof(Real));        
                // ... and int array data
                for (int comp = 0; comp < NumIntComps(); ++comp) {
                    if (communicate_int_comp[comp]) {
                        host_int_attribs[lev][ind][comp].push_back(*idata++);
                    } else {
                        host_int_attribs[lev][ind][comp].push_back(0);
                    }
                }
                ++ipart;
            }
        }

	for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
	{
	    for (auto& kv : host_particles[host_lev]) {
	      auto grid = kv.first.first;
	      auto tile = kv.first.second;
	      const auto& src_tile = kv.second;
	      
	      auto& dst_tile = GetParticles(host_lev)[std::make_pair(grid,tile)];
	      auto old_size = dst_tile.GetArrayOfStructs().size();
	      auto new_size = old_size + src_tile.size();
	      dst_tile.resize(new_size);
	      
	      Cuda::thrust_copy(src_tile.begin(),
				src_tile.end(),
				dst_tile.GetArrayOfStructs().begin() + old_size);
	      
	      for (int i = 0; i < NumRealComps(); ++i) {
		Cuda::thrust_copy(host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				  host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				  dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
	      }
	      
	      for (int i = 0; i < NumIntComps(); ++i) {
		Cuda::thrust_copy(host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				  host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				  dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
	      }
	    }
	  }
    
	Gpu::Device::streamSynchronize();
#endif

	BL_PROFILE_VAR_STOP(blp_copy);
    }
#endif /*BL_USE_MPI*/
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::OK (int lev_min, int lev_max, int nGrow) const
{
    BL_PROFILE("ParticleContainer::OK()");

    if (lev_max == -1)
        lev_max = finestLevel();

    return (numParticlesOutOfRange(*this, lev_min, lev_max, nGrow) == 0);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::AddParticlesAtLevel (AoS& particles, int level, int nGrow)
{
    BL_PROFILE("ParticleContainer::AddParticlesAtLevel()");
    if (int(m_particles.size()) < level+1)
        {
            if (Verbose())
            {
                amrex::Print() << "ParticleContainer::AddParticlesAtLevel resizing m_particles from "
                               << m_particles.size()
                               << " to "
                               << level+1 << '\n';
            }
            m_particles.resize(level + 1);
            m_dummy_mf.resize(level+1);
            for (int lev = 0; lev < level+1; ++lev) {
                RedefineDummyMF(lev);
            }
        }    

    ParticleLocData pld;

    for (int i = 0; i < particles.size(); ++i) {
        ParticleType& p = particles[i];

        if (p.id() > 0)
        {
            if (!Where(p, pld, level, level, nGrow))
                amrex::Abort("ParticleContainerAddParticlesAtLevel(): Can't add outside of domain\n");
            m_particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
        }
    }
    Redistribute(level, level, nGrow);
    particles.resize(0);
}

// This is the single-level version for cell-centered density
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
AssignCellDensitySingleLevel (int rho_index,
                              MultiFab& mf_to_be_filled,
                              int       lev,
                              int       ncomp,
                              int       particle_lvl_offset) const
{
    BL_PROFILE("ParticleContainer::AssignCellDensitySingleLevel()");
    
    if (rho_index != 0) amrex::Abort("AssignCellDensitySingleLevel only works if rho_index = 0");
    
    MultiFab* mf_pointer;

    if (OnSameGrids(lev, mf_to_be_filled)) {
      // If we are already working with the internal mf defined on the 
      // particle_box_array, then we just work with this.
      mf_pointer = &mf_to_be_filled;
    }
    else {
      // If mf_to_be_filled is not defined on the particle_box_array, then we need 
      // to make a temporary here and copy into mf_to_be_filled at the end.
      mf_pointer = new MultiFab(ParticleBoxArray(lev), 
				ParticleDistributionMap(lev),
				ncomp, mf_to_be_filled.nGrow());
    }

    // We must have ghost cells for each FAB so that a particle in one grid can spread 
    // its effect to an adjacent grid by first putting the value into ghost cells of its
    // own grid.  The mf->SumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       amrex::Error("Must have at least one ghost cell when in AssignCellDensitySingleLevel");

    const Real      strttime    = amrex::second();

    const auto dxi              = Geom(lev).InvCellSizeArray();
    const auto plo              = Geom(lev).ProbLoArray();
    const auto pdxi             = Geom(lev + particle_lvl_offset).InvCellSizeArray();

    if (Geom(lev).isAnyPeriodic() && ! Geom(lev).isAllPeriodic())
    {
        amrex::Error("AssignCellDensitySingleLevel: problem must be periodic in no or all directions");
    }

    mf_pointer->setVal(0);
    
    using ParConstIter = ParConstIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        FArrayBox local_rho;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            const auto pstruct = particles().data();
            const long np = pti.numParticles();
            FArrayBox& fab = (*mf_pointer)[pti];
            auto rhoarr = fab.array();
#ifdef _OPENMP
            Box tile_box;
            if (Gpu::notInLaunchRegion())
            {
                tile_box = pti.tilebox();
                tile_box.grow(mf_pointer->nGrow());
                local_rho.resize(tile_box,ncomp);
                local_rho = 0.0;
                rhoarr = local_rho.array();
            }
#endif
                        
            if (particle_lvl_offset == 0)
            {
                AMREX_FOR_1D( np, i,
                {
                    amrex_deposit_cic(pstruct[i], ncomp, rhoarr, plo, dxi);
                });
            }
            else
            {
                AMREX_FOR_1D( np, i,
                {
                    amrex_deposit_particle_dx_cic(pstruct[i], ncomp, rhoarr, plo, dxi, pdxi);
                });
            }
                
#ifdef _OPENMP
            if (Gpu::notInLaunchRegion())
            {
                fab.atomicAdd(local_rho, tile_box, tile_box, 0, 0, ncomp);
            }
#endif
        }
    }
    
    mf_pointer->SumBoundary(Geom(lev).periodicity());
    
    // If ncomp > 1, first divide the momenta (component n) 
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    for (int n = 1; n < ncomp; n++)
    {
        for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi)
        {
            (*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
        }
    }
    
    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    const Real* dx = Geom(lev).CellSize();
    const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    
    mf_pointer->mult(1.0/vol, 0, 1, mf_pointer->nGrow());
    
    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled)
    {
        mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
        delete mf_pointer;
    }
    
    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;
        
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        
        amrex::Print() << "ParticleContainer::AssignCellDensitySingleLevel) time: "
                       << stoptime << '\n';
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Interpolate (Vector<std::unique_ptr<MultiFab> >& mesh_data, 
                                                                                int lev_min, int lev_max)
{
    BL_PROFILE("ParticleContainer::Interpolate()");
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        InterpolateSingleLevel(*mesh_data[lev], lev); 
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
InterpolateSingleLevel (MultiFab& mesh_data, int lev)
{
    BL_PROFILE("ParticleContainer::InterpolateSingleLevel()");
    
    if (mesh_data.nGrow() < 1)
        amrex::Error("Must have at least one ghost cell when in InterpolateSingleLevel");
    
    const Geometry& gm = Geom(lev);
    const auto     plo = gm.ProbLoArray();
    const auto     dxi = gm.InvCellSizeArray();

    using ParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (ParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        auto& particles = pti.GetArrayOfStructs();
        auto pstruct = particles().data();
        FArrayBox& fab = mesh_data[pti];
        const auto fabarr = fab.array();
        const Box& box = fab.box();
        const long np = particles.size();
        
        int nComp = fab.nComp();
        AMREX_FOR_1D( np, i,
        {
            amrex_interpolate_cic(pstruct[i], nComp, fabarr, plo, dxi);
        });
    }
}

//
// This version takes as input the acceleration vector at cell centers, and has the option of
// returning the acceleration at the particle location in the data array, starting at
// component start_comp_for_accel
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::moveKick (MultiFab& acceleration, int lev, Real dt, Real a_new, Real a_half, int start_comp_for_accel)
{
    BL_PROFILE("ParticleContainer::moveKick()");
    BL_ASSERT(NStructReal >= AMREX_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    const Real strttime  = amrex::second();
    const Real half_dt   = 0.5_rt * dt;
    const Real a_new_inv = 1_rt / a_new;
    auto&      pmap      = m_particles[lev];

    MultiFab* ac_pointer;
    if (OnSameGrids(lev,acceleration))
    {
        ac_pointer = &acceleration;
    }
    else 
    {
        ac_pointer = new MultiFab(ParticleBoxArray(lev),
				  ParticleDistributionMap(lev),
				  acceleration.nComp(),acceleration.nGrow());
        for (MFIter mfi(*ac_pointer); mfi.isValid(); ++mfi)
            ac_pointer->setVal(0.);
        ac_pointer->copy(acceleration,0,0,acceleration.nComp());
        ac_pointer->FillBoundary(); // DO WE NEED GHOST CELLS FILLED ???
    }

    for (auto& kv : pmap) {
      auto& pbox = kv.second.GetArrayOfStructs();
      const int grid = kv.first.first;
      const int n = pbox.size();
      const FArrayBox& gfab = (*ac_pointer)[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i = 0; i < n; i++)
        {
	  ParticleType& p = pbox[i];

	  if (p.m_idata.id > 0)
            {

	      //
	      // Note: rdata.arr[AMREX_SPACEDIM] is mass, AMREX_SPACEDIM+1 is v_x, ...
	      //
	      Real grav[AMREX_SPACEDIM];

	      ParticleType::GetGravity(gfab, m_gdb->Geom(lev), p, grav);
	      //
	      // Define (a u)^new = (a u)^half + dt/2 grav^new
	      //
	      AMREX_D_TERM(p.m_rdata.arr[AMREX_SPACEDIM+1] *= a_half;,
		     p.m_rdata.arr[AMREX_SPACEDIM+2] *= a_half;,
		     p.m_rdata.arr[AMREX_SPACEDIM+3] *= a_half;);

	      AMREX_D_TERM(p.m_rdata.arr[AMREX_SPACEDIM+1] += half_dt * grav[0];,
		     p.m_rdata.arr[AMREX_SPACEDIM+2] += half_dt * grav[1];,
		     p.m_rdata.arr[AMREX_SPACEDIM+3] += half_dt * grav[2];);

	      AMREX_D_TERM(p.m_rdata.arr[AMREX_SPACEDIM+1] *= a_new_inv;,
		     p.m_rdata.arr[AMREX_SPACEDIM+2] *= a_new_inv;,
		     p.m_rdata.arr[AMREX_SPACEDIM+3] *= a_new_inv;);

	      if (start_comp_for_accel > AMREX_SPACEDIM)
                {
		  AMREX_D_TERM(p.m_rdata.arr[AMREX_SPACEDIM + start_comp_for_accel  ] = grav[0];,
			 p.m_rdata.arr[AMREX_SPACEDIM + start_comp_for_accel+1] = grav[1];,
			 p.m_rdata.arr[AMREX_SPACEDIM + start_comp_for_accel+2] = grav[2];);
                }
            }
        }
    }

    
    if (ac_pointer != &acceleration) delete ac_pointer;

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "ParticleContainer::moveKick() time: " << stoptime << '\n';
    }
}
