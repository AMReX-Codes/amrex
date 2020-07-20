
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
        num_real_comm_comps*sizeof(ParticleReal) + num_int_comm_comps*sizeof(int);
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

    AMREX_D_TERM(iv[0]=static_cast<int>(floor((p.pos(0)-geom.ProbLo(0))*geom.InvCellSize(0)));,
                 iv[1]=static_cast<int>(floor((p.pos(1)-geom.ProbLo(1))*geom.InvCellSize(1)));,
                 iv[2]=static_cast<int>(floor((p.pos(2)-geom.ProbLo(2))*geom.InvCellSize(2))););

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

  AMREX_ASSERT(m_gdb != 0);
  
  if (lev_max == -1)
      lev_max = finestLevel();
  
  AMREX_ASSERT(lev_max <= finestLevel());

  AMREX_ASSERT(nGrow == 0 || (nGrow >= 0 && lev_min == lev_max));

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
      AMREX_ASSERT(ba.ixType().cellCentered());

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

    AMREX_ASSERT(m_gdb != 0);

    if (!Geom(0).isAnyPeriodic()) return false;

    if (lev_max == -1)
        lev_max = finestLevel();

    AMREX_ASSERT(lev_max <= finestLevel());

    // Create a copy "dummy" particle to check for periodic outs.
    ParticleType p_prime = p;
    if (PeriodicShift(p_prime)) {
        std::vector< std::pair<int,Box> > isects;
        for (int lev = lev_max; lev >= lev_min; lev--) {

	    int grid;
            IntVect iv;
            const BoxArray& ba = ParticleBoxArray(lev);
            AMREX_ASSERT(ba.ixType().cellCentered());

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
		    ba.intersections(Box(iv, iv), isects, true, 0);
		    grid = isects.empty() ? -1 : isects[0].first;
		    if(grid == -1)
		    {
		        grid = (*redistribute_mask_ptr)[local_grid](Index(p, lev), 0);
		    }
		}
	    }

            if (grid >= 0) {
                AMREX_D_TERM(p.pos(0) = p_prime.pos(0);,
                             p.pos(1) = p_prime.pos(1);,
                             p.pos(2) = p_prime.pos(2););

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
    const auto& geom = Geom(0);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto is_per = geom.isPeriodicArray();

    return enforcePeriodic(p, plo, phi, is_per);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
ParticleLocData
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
Reset (ParticleType& p,
       bool          update,
       bool          verbose,
       ParticleLocData pld) const
{
    AMREX_ASSERT(m_gdb != 0);

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

	AMREX_ASSERT(p.id() > 0);

	p.id() = -p.id();
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
    bool outside = AMREX_D_TERM(p.pos(0) <  Geom(0).ProbLo(0)
                             || p.pos(0) >= Geom(0).ProbHi(0),
                             || p.pos(1) <  Geom(0).ProbLo(1)
                             || p.pos(1) >= Geom(0).ProbHi(1),
                             || p.pos(2) <  Geom(0).ProbLo(2)
                             || p.pos(2) >= Geom(0).ProbHi(2));

    bool success;
    if (outside)
    {
      // Note that EnforcePeriodicWhere may shift the particle if it is successful.
      success = EnforcePeriodicWhere(p, pld, lev_min, lev_max, local_grid);
      if (!success && lev_min == 0)
      {
          // The particle has left the domain; invalidate it.
          p.id() = -p.id();
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
Long
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::TotalNumberOfParticles (bool only_valid, bool only_local) const
{
    Long nparticles = 0;
    for (int lev = 0; lev <= finestLevel(); lev++) {
        nparticles += NumberOfParticlesAtLevel(lev,only_valid,true);
    }
    if (!only_local) {
        ParallelAllReduce::Sum(nparticles, ParallelContext::CommunicatorSub());
    }
    return nparticles;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
Vector<Long>
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::NumberOfParticlesInGrid (int lev, bool only_valid, bool only_local) const
{
    AMREX_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    LayoutData<Long> np_per_grid_local(ParticleBoxArray(lev),
                                       ParticleDistributionMap(lev));

    for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti)
    {
        int gid = pti.index();
        if (only_valid)
        {
            const auto& ptile = ParticlesAt(lev, pti);
            const auto& aos = ptile.GetArrayOfStructs();
            const auto pstruct = aos().dataPtr();
            const int np = ptile.numParticles();

            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(np, reduce_data,
                           [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                           {
                               return (pstruct[i].id() > 0) ? 1 : 0;
                           });

            int np_valid = amrex::get<0>(reduce_data.value());
            np_per_grid_local[gid] += np_valid;
        } else
        {
            np_per_grid_local[gid] += pti.numParticles();
        }
    }
        
    Vector<Long> nparticles(np_per_grid_local.size(), 0);
    if (only_local)
    {
        for (ParConstIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            nparticles[pti.index()] = np_per_grid_local[pti.index()];
        }
    } 
    else
    {
        ParallelDescriptor::GatherLayoutDataToVector(np_per_grid_local, nparticles,
                                                     ParallelContext::IOProcessorNumberSub());
        ParallelDescriptor::Bcast(&nparticles[0], nparticles.size(), 
                                  ParallelContext::IOProcessorNumberSub());
    }

    return nparticles;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
Long
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::NumberOfParticlesAtLevel (int lev, bool only_valid, bool only_local) const
{
    Long nparticles = 0;

    if (lev >= 0 && lev < int(m_particles.size())) {
        for (const auto& kv : GetParticles(lev)) {
            const auto& ptile = kv.second;
            if (only_valid) {
                for (int k = 0; k < ptile.GetArrayOfStructs().numParticles(); ++k) {
                    const ParticleType& p = ptile.GetArrayOfStructs()[k];
                    if (p.id() > 0) ++nparticles;
                }
            } else {
                nparticles += ptile.numParticles();
            }
        }
    }

    if (!only_local) {
        ParallelAllReduce::Sum(nparticles, ParallelContext::CommunicatorSub());
    }

    return nparticles;
}

//
// This includes both valid and invalid particles since invalid particles still take up space.
//

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::ByteSpread () const
{
    Long cnt = 0;

    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& ptile = kv.second;
            cnt += ptile.numParticles();
        }
    }

    Long mn = cnt, mx = mn;

    const int IOProc = ParallelContext::IOProcessorNumberSub();
    const std::size_t sz = sizeof(ParticleType)+NumRealComps()*sizeof(Real)+NumIntComps()*sizeof(int);

#ifdef AMREX_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelReduce::Min(mn,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Max(mx,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Sum(cnt, IOProc, ParallelContext::CommunicatorSub());

            amrex::Print() << "ParticleContainer byte spread across MPI nodes: ["
                           << mn*sz
                           << " (" << mn << ")"
                           << " ... "
                           << mx*sz
                           << " (" << mx << ")"
                           << "] total particles: (" << cnt << ")\n";
#ifdef AMREX_LAZY
        });
#endif
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::PrintCapacity () const
{
    Long cnt = 0;

    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& ptile = kv.second;
            cnt += ptile.capacity();
        }
    }

    Long mn = cnt, mx = mn;

    const int IOProc = ParallelContext::IOProcessorNumberSub();

#ifdef AMREX_LAZY
    Lazy::QueueReduction( [=] () mutable {
#endif
            ParallelReduce::Min(mn,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Max(mx,  IOProc, ParallelContext::CommunicatorSub());
            ParallelReduce::Sum(cnt, IOProc, ParallelContext::CommunicatorSub());

            amrex::Print() << "ParticleContainer byte spread across MPI nodes: ["
                           << mn
                           << " (" << mn << ")"
                           << " ... "
                           << mx
                           << " (" << mx << ")"
                           << "] total memory: (" << cnt << ")\n";
#ifdef AMREX_LAZY
        });
#endif
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::ShrinkToFit ()
{
    for (unsigned lev = 0; lev < m_particles.size(); lev++) {
        auto& pmap = m_particles[lev];
        for (auto& kv : pmap) {
            auto& ptile = kv.second;
            ptile.shrink_to_fit();
        }
    }
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
    AMREX_ASSERT(OK());
    AMREX_ASSERT(m_gdb != 0);
    // 
    // Move particles up to FRAC*CellSize distance in each coordinate direction.
    //
    const Real FRAC = 0.25;
    auto&       pmap              = m_particles[lev];
    const Real* dx                = Geom(lev).CellSize();
    const Real  dist[AMREX_SPACEDIM] = { AMREX_D_DECL(FRAC*dx[0], FRAC*dx[1], FRAC*dx[2]) };

    for (auto& kv : pmap) {
        auto& aos = kv.second.GetArrayOfStructs();
        const int n = aos.numParticles();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
	  ParticleType& p = aos[i];

	  if (p.id() <= 0) continue;

	  for (int j = 0; j < AMREX_SPACEDIM; j++)
              {
                  p.pos(j) += dist[j]*(2*amrex::Random()-1);
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
Long
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::IncrementWithTotal (MultiFab& mf, int lev, bool local)
{
  BL_PROFILE("ParticleContainer::IncrementWithTotal(lev)");
  AMREX_ASSERT(OK());
  
  if (m_particles.empty()) return 0;
  
  AMREX_ASSERT(lev >= 0 && lev < int(m_particles.size()));
  
  AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

  const auto& pmap = m_particles[lev];
  
  Long num_particles_in_domain = 0;
  
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
      for (int k = 0; k < pbox.numParticles(); ++ k) {
	const ParticleType& p = pbox[k];
        if (p.id() > 0) {
              Where(p, pld);
              AMREX_ASSERT(pld.m_grid == gid);
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

  if (!local) {
      ParallelAllReduce::Sum(num_particles_in_domain, ParallelContext::CommunicatorSub());
  }

  return num_particles_in_domain;
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
  AMREX_ASSERT(this->finestLevel()+1 == int(this->m_particles.size()));

  Long cnt = 0;

  for (unsigned lev = 0; lev < m_particles.size() - 1; ++lev) {
      auto& pmap = m_particles[lev];
      if (!pmap.empty()) {
          for (auto& kv : pmap) {
              const auto& pbx = kv.second;
              cnt += pbx.numParticles();
          }
          ParticleLevel().swap(pmap);
      }
  }
  
  //
  // Print how many particles removed on each processor if any were removed.
  //
  if (this->m_verbose > 1 && cnt > 0) {
      amrex::AllPrint() << "Processor " << ParallelContext::MyProcSub() << " removed " << cnt
                        << " particles not in finest level\n";
  }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CreateVirtualParticles (int level, AoS& virts) const
{
    BL_PROFILE("ParticleContainer::CreateVirtualParticles()");
    AMREX_ASSERT(level > 0);
    AMREX_ASSERT(virts.empty());
    
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
	        p.id() = VirtualParticleID;
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
                    p.id() = VirtualParticleID;
                    virts.push_back(p);
                }
                else
                {
                    //
                    // Note that Cell aggregation assumes that p.rdata(0) is mass and
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
                        p.id() = VirtualParticleID;
                        agg_map[cell] = p;
                    }
                    else
                    {
                        AMREX_ASSERT(agg_map_it != agg_map.end());
                        const ParticleType&  pnew       = *it;
                        ParticleType&        pold       = agg_map_it->second;
                        const Real           old_mass   = pold.rdata(0);
                        const Real           new_mass   = pnew.rdata(0);
                        const Real           total_mass = old_mass + new_mass;
                        //
                        // Set the position to the center of mass.
                        //
                        for (int i = 0; i < AMREX_SPACEDIM; i++)
                        {
                            pold.pos(i) = (old_mass*pold.pos(i) + new_mass*pnew.pos(i))/total_mass;
                        }
                        AMREX_ASSERT(this->Index(pold, level) == cell);
                        //
                        // Set the metadata (presumably velocity) to the mass-weighted average.
                        //
                        for (int i = 1; i < NStructReal; i++)
                        {
                            pold.rdata(i) = (old_mass*pold.rdata(i) + new_mass*pnew.rdata(i))/total_mass;
                        }
                        pold.rdata(0) = total_mass;
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
    AMREX_ASSERT(ghosts.empty());
    AMREX_ASSERT(level < finestLevel());
  
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
                p.id() = GhostParticleID;
                ghosts().push_back(p);
            }
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
clearParticles ()
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
copyParticles (const ParticleContainerType& other, bool local)
{
    using PData = ConstParticleTileData<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    copyParticles(other, [=] AMREX_GPU_HOST_DEVICE (const PData& data, int i) { return 1; }, local);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
addParticles (const ParticleContainerType& other, bool local)
{
    using PData = ConstParticleTileData<NStructReal, NStructInt, NArrayReal, NArrayInt>; 
    addParticles(other, [=] AMREX_GPU_HOST_DEVICE (const PData& data, int i) { return 1; }, local);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F,
          amrex::EnableIf_t<! std::is_integral<F>::value, int> foo>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
copyParticles (const ParticleContainerType& other, F&& f, bool local)
{
    BL_PROFILE("ParticleContainer::copyParticles");
    clearParticles();
    addParticles(other, std::forward<F>(f), local);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class F,
          amrex::EnableIf_t<! std::is_integral<F>::value, int> foo>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
addParticles (const ParticleContainerType& other, F&& f, bool local)
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

            auto& ptile = plevel[index];
            const auto& ptile_other = plevel_other.at(index);
            auto np = ptile_other.numParticles();
            if (np == 0) continue;

            auto dst_index = ptile.numParticles();
            ptile.resize(dst_index + np);

            auto count = amrex::filterParticles(ptile, ptile_other, std::forward<F>(f), 0, dst_index, np);

            ptile.resize(dst_index + count);
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
#ifdef AMREX_USE_GPU
    if ( Gpu::inLaunchRegion() )
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
    BL_PROFILE("ParticleContainer::SortParticlesByCell()");

    for (int lev = 0; lev < numLevels(); ++lev)
    {
        const Geometry& geom = Geom(lev);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();
        const auto domain = geom.Domain();

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& ptile = ParticlesAt(lev, mfi);
            auto& aos   = ptile.GetArrayOfStructs();
            const size_t np = aos.numParticles();
            auto pstruct_ptr = aos().dataPtr();
            
            ParticleTileType ptile_tmp;
            ptile_tmp.define(m_num_runtime_real, m_num_runtime_int);
            ptile_tmp.resize(np);

            m_bins.build(np, pstruct_ptr, mfi.tilebox(),
                       [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> IntVect
                       {
                           return getParticleCell(p, plo, dxi, domain);
                       });
          
            gatherParticles(ptile_tmp, ptile, np, m_bins.permutationPtr());
            ptile.swap(ptile_tmp);
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::SortParticlesByBin (IntVect bin_size)
{
    BL_PROFILE("ParticleContainer::SortParticlesByBin()");

    for (int lev = 0; lev < numLevels(); ++lev)
    {
        const Geometry& geom = Geom(lev);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();
        const auto domain = geom.Domain();

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& ptile = ParticlesAt(lev, mfi);
            auto& aos   = ptile.GetArrayOfStructs();
            const size_t np = aos.numParticles();
            auto pstruct_ptr = aos().dataPtr();

            ParticleTileType ptile_tmp;
            ptile_tmp.define(m_num_runtime_real, m_num_runtime_int);
            ptile_tmp.resize(np);

            const Box& box = mfi.tilebox();
            IntVect lo = box.smallEnd();
			
            m_bins.build(np, pstruct_ptr, mfi.tilebox(),
                       [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> IntVect
                       {
                           return (getParticleCell(p, plo, dxi, domain) - lo) / bin_size;
                       });
          
            gatherParticles(ptile_tmp, ptile, np, m_bins.permutationPtr());
            ptile.swap(ptile_tmp);
        }
    }
}

//
// The GPU implementation of Redistribute
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RedistributeGPU (int lev_min, int lev_max, int nGrow, int local)
{
#ifdef AMREX_USE_GPU

    if (local) AMREX_ASSERT(numParticlesOutOfRange(*this, lev_min, lev_max, local) == 0);

    // sanity check
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
            EnforcePeriodic();
            AMREX_ASSERT(OK());
            return;
        }
    }

    if (not m_particle_locator.isValid(GetParGDB())) m_particle_locator.build(GetParGDB());
    m_particle_locator.setGeometry(GetParGDB());
    auto assign_grid = m_particle_locator.getGridAssignor();

    BL_PROFILE_VAR_START(blp_partition);
    ParticleCopyOp op;
    op.setNumLevels(num_levels);
    Vector<std::map<int, int> > new_sizes(num_levels);
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        const Geometry& geom = Geom(lev);
        const BoxArray& ba   = ParticleBoxArray(lev);

        auto& plev = m_particles[lev];
        for (auto& kv : plev)
        {
            int gid = kv.first.first;
            int tid = kv.first.second;
            auto index = std::make_pair(gid, tid);

            auto& src_tile = plev[index];
            auto& aos = src_tile.GetArrayOfStructs();
            const size_t np = aos.numParticles();

            int num_stay = partitionParticlesByDest(src_tile, assign_grid, BufferMap(),
                                                    geom, lev, gid, tid,
                                                    lev_min, lev_max, nGrow);

            int num_move = np - num_stay;
            new_sizes[lev][gid] = num_stay;
            op.resize(gid, lev, num_move);

            auto p_boxes = op.m_boxes[lev][gid].dataPtr();
            auto p_levs = op.m_levels[lev][gid].dataPtr();
            auto p_src_indices = op.m_src_indices[lev][gid].dataPtr();
            auto p_periodic_shift = op.m_periodic_shift[lev][gid].dataPtr();
            auto p_ptr = &(aos[0]);
            
	    AMREX_FOR_1D ( num_move, i,
            {
                const auto& p = p_ptr[i + num_stay];
                if (p.id() < 0)
                {
                    p_boxes[i] = -1;
                    p_levs[i]  = -1;                    
                }
                else
                {
                    const auto tup = assign_grid(p, lev_min, lev_max, nGrow);
                    p_boxes[i] = amrex::get<0>(tup);
                    p_levs[i]  = amrex::get<1>(tup);
                }
                p_periodic_shift[i] = IntVect(AMREX_D_DECL(0,0,0));
                p_src_indices[i] = i+num_stay;
            });
        }
    }
    BL_PROFILE_VAR_STOP(blp_partition);

    ParticleCopyPlan plan;
    
    plan.build(*this, op, local);

    Gpu::DeviceVector<char> snd_buffer;
    Gpu::DeviceVector<char> rcv_buffer;

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
        Gpu::Device::synchronize();
        Gpu::PinnedVector<char> pinned_snd_buffer;
        Gpu::PinnedVector<char> pinned_rcv_buffer;
        pinned_snd_buffer.resize(snd_buffer.size());
        Gpu::dtoh_memcpy_async(pinned_snd_buffer.dataPtr(), snd_buffer.dataPtr(), snd_buffer.size());
        plan.buildMPIFinish(BufferMap());
        Gpu::Device::synchronize();
        communicateParticlesStart(*this, plan, pinned_snd_buffer, pinned_rcv_buffer);
        rcv_buffer.resize(pinned_rcv_buffer.size());
        unpackBuffer(*this, plan, snd_buffer, RedistributeUnpackPolicy());
        communicateParticlesFinish(plan);
        Gpu::htod_memcpy_async(rcv_buffer.dataPtr(), pinned_rcv_buffer.dataPtr(), pinned_rcv_buffer.size());
        unpackRemotes(*this, plan, rcv_buffer, RedistributeUnpackPolicy());
    }

    Gpu::Device::synchronize();
    AMREX_ASSERT(numParticlesOutOfRange(*this, lev_min, lev_max, nGrow) == 0);
#endif
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::EnforcePeriodic ()
{
    BL_PROFILE("ParticleContainer::EnforcePeriodic()");
    const int lev = 0;
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

        const int np = particles.numParticles();
        ParticleType* pstruct = &(particles.GetArrayOfStructs()[0]);
        AMREX_FOR_1D ( np, i,
        {
            enforcePeriodic(pstruct[i], plo, phi, is_per);
        });
    }
}

//
// The CPU implementation of Redistribute
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RedistributeCPU (int lev_min, int lev_max, int nGrow, int local)
{
  BL_PROFILE("ParticleContainer::RedistributeCPU()");

  const int MyProc    = ParallelContext::MyProcSub();
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
  AMREX_ASSERT(lev_max <= finestLevel());

  // This will hold the valid particles that go to another process
  std::map<int, Vector<char> > not_ours;

  int num_threads = OpenMP::get_max_threads();

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
    for (int i = 0; i < ParallelContext::NProcsSub(); ++i)
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
          int thread_num = OpenMP::get_thread_num();
          int grid = grid_tile_ids[pmap_it].first;
          int tile = grid_tile_ids[pmap_it].second;
          auto& aos = ptile_ptrs[pmap_it]->GetArrayOfStructs();
          auto& soa = ptile_ptrs[pmap_it]->GetStructOfArrays();
          unsigned npart = aos.numParticles();
          ParticleLocData pld;
          if (npart != 0) {
              Long last = npart - 1;
              unsigned pindex = 0;
              while (pindex <= last) {
                  ParticleType& p = aos[pindex];

                  if (p.id() < 0)
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

                  if (p.id() < 0)
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

                  const int who = ParallelContext::global_to_local_rank(ParticleDistributionMap(pld.m_lev)[pld.m_grid]);
                  if (who == MyProc) {
                      if (pld.m_lev != lev || pld.m_grid != grid || pld.m_tile != tile) {
                          // We own it but must shift it to another place.
                          auto index = std::make_pair(pld.m_grid, pld.m_tile);
                          AMREX_ASSERT(tmp_local[pld.m_lev][index].size() == num_threads);
                          tmp_local[pld.m_lev][index][thread_num].push_back(p);
                          for (int comp = 0; comp < NumRealComps(); ++comp) {
                              RealVector& arr = soa_local[pld.m_lev][index][thread_num].GetRealData(comp);
                              arr.push_back(soa.GetRealData(comp)[pindex]);
                          }
                          for (int comp = 0; comp < NumIntComps(); ++comp) {
                              IntVector& arr = soa_local[pld.m_lev][index][thread_num].GetIntData(comp);
                              arr.push_back(soa.GetIntData(comp)[pindex]);
                          }

                          p.id() = -p.id(); // Invalidate the particle
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

                      p.id() = -p.id(); // Invalidate the particle
                  }

                  if (p.id() < 0)
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
      AMREX_ASSERT(int(m_particles.size()) >= 2);

      m_particles.resize(theEffectiveFinestLevel + 1);
      m_dummy_mf.resize(theEffectiveFinestLevel + 1);
  }

  if (ParallelContext::NProcsSub() == 1) {
      AMREX_ASSERT(not_ours.empty());
  }
  else {
      RedistributeMPI(not_ours, lev_min, lev_max, nGrow, local);
  }

  AMREX_ASSERT(OK(lev_min, lev_max, nGrow));

  if (m_verbose > 0) {
      Real stoptime = amrex::second() - strttime;

      ByteSpread();

#ifdef AMREX_LAZY
      Lazy::QueueReduction( [=] () mutable {
#endif
              ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                                  ParallelContext::CommunicatorSub());

              amrex::Print() << "ParticleContainer::Redistribute() time: " << stoptime << "\n\n";
#ifdef AMREX_LAZY
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
    AMREX_ASSERT(lev == 0);

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

        const auto tile_size_do = this->do_tiling ? this->tile_size : IntVect::TheZeroVector();

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*redistribute_mask_ptr, tile_size_do); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.tilebox();
            const int grid_id = mfi.index();
            const int tile_id = mfi.LocalTileIndex();
            (*redistribute_mask_ptr)[mfi].template setVal<RunOn::Host>(grid_id, box, 0, 1);
            (*redistribute_mask_ptr)[mfi].template setVal<RunOn::Host>(tile_id, box, 1, 1);
        }

        redistribute_mask_ptr->FillBoundary(geom.periodicity());

        neighbor_procs.clear();
        for (MFIter mfi(*redistribute_mask_ptr, tile_size_do); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.growntilebox();
            for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv))
            {
                const int grid = (*redistribute_mask_ptr)[mfi](iv, 0);
                if (grid >= 0)
                {
                    const int global_rank = this->ParticleDistributionMap(lev)[grid];
                    const int rank = ParallelContext::global_to_local_rank(global_rank);
                    if (rank != ParallelContext::MyProcSub())
                        neighbor_procs.push_back(rank);
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

#ifdef AMREX_USE_MPI

    using buffer_type = unsigned long long;

    std::map<int, Vector<buffer_type> > mpi_snd_data;
    for (const auto& kv : not_ours)
    {
        int nbt = (kv.second.size() + sizeof(buffer_type)-1)/sizeof(buffer_type);
        mpi_snd_data[kv.first].resize(nbt);
        std::memcpy((char*) mpi_snd_data[kv.first].data(), kv.second.data(), kv.second.size());
    }

    const int NProcs = ParallelContext::NProcsSub();
    const int NNeighborProcs = neighbor_procs.size();

    // We may now have particles that are rightfully owned by another CPU.
    Vector<Long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    Long NumSnds = 0;
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
        Long tot_snds_this_proc = 0;
        Long tot_rcvs_this_proc = 0;
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
        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < size_t(std::numeric_limits<int>::max()));
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum,
                                             ParallelContext::CommunicatorSub()).req();
    }

    // Send.
    for (const auto& kv : mpi_snd_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
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
            for (int i = 0; i < int(Cnt); ++i)
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

#ifndef AMREX_USE_GPU
        ipart = 0;
        for (int i = 0; i < nrcvs; ++i)
        {
            const auto offset = rOffset[i];
            const auto Who    = RcvProc[i];
            const auto Cnt = Rcvs[Who] / superparticle_size;
            for (int j = 0; j < int(Cnt); ++j)
            {
                auto& ptile = m_particles[rcv_levs[ipart]][std::make_pair(rcv_grid[ipart],
                                                                          rcv_tile[ipart])];
                char* pbuf = ((char*) &recvdata[offset]) + j*superparticle_size;

                ParticleType p;
                std::memcpy(&p, pbuf, sizeof(ParticleType));
                pbuf += sizeof(ParticleType);
                ptile.push_back(p);
                for (int comp = 0; comp < NumRealComps(); ++comp) {
                    if (communicate_real_comp[comp]) {
                        Real rdata;
                        std::memcpy(&rdata, pbuf, sizeof(Real));
                        pbuf += sizeof(Real);
                        ptile.push_back_real(comp, rdata);
                    } else {
                        ptile.push_back_real(comp, 0.0);
                    }
                }

                for (int comp = 0; comp < NumIntComps(); ++comp) {
                    if (communicate_int_comp[comp]) {
                        int idata;
                        std::memcpy(&idata, pbuf, sizeof(int));
                        pbuf += sizeof(int);
                        ptile.push_back_int(comp, idata);
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
			std::vector<Gpu::HostVector<Real> > > > host_real_attribs;
	host_real_attribs.reserve(15);
	host_real_attribs.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>,
			std::vector<Gpu::HostVector<int> > > > host_int_attribs;
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
                pbuf += sizeof(ParticleType);

                host_real_attribs[lev][ind].resize(NumRealComps());
                host_int_attribs[lev][ind].resize(NumIntComps());

                // add the struct
                host_particles[lev][ind].push_back(p);
	  
                // add the real...
                for (int comp = 0; comp < NumRealComps(); ++comp) {
                    if (communicate_real_comp[comp]) {
                        Real rdata;
                        std::memcpy(&rdata, pbuf, sizeof(Real));
                        pbuf += sizeof(Real);
                        host_real_attribs[lev][ind][comp].push_back(rdata);
                    } else {
                        host_real_attribs[lev][ind][comp].push_back(0.0);
                    }                    
                }
                
                // ... and int array data
                for (int comp = 0; comp < NumIntComps(); ++comp) {
                    if (communicate_int_comp[comp]) {
                        int idata;
                        std::memcpy(&idata, pbuf, sizeof(int));
                        pbuf += sizeof(int);
                        host_int_attribs[lev][ind][comp].push_back(idata);
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
	      
	      Gpu::copy(Gpu::hostToDevice,
                        src_tile.begin(), src_tile.end(),
                        dst_tile.GetArrayOfStructs().begin() + old_size);
	      
	      for (int i = 0; i < NumRealComps(); ++i) {
                  Gpu::copy(Gpu::hostToDevice,
                            host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                            host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                            dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
	      }
	      
	      for (int i = 0; i < NumIntComps(); ++i) {
                  Gpu::copy(Gpu::hostToDevice,
                            host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                            host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                            dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
	      }
	    }
	  }
    
	Gpu::Device::streamSynchronize();
#endif

	BL_PROFILE_VAR_STOP(blp_copy);
    }
#endif /*AMREX_USE_MPI*/
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

    for (int i = 0; i < particles.numParticles(); ++i) {
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
            const Long np = pti.numParticles();
            FArrayBox& fab = (*mf_pointer)[pti];
            auto rhoarr = fab.array();
#ifdef _OPENMP
            Box tile_box;
            if (Gpu::notInLaunchRegion())
            {
                tile_box = pti.tilebox();
                tile_box.grow(mf_pointer->nGrow());
                local_rho.resize(tile_box,ncomp);
                local_rho.setVal<RunOn::Host>(0.0);
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
                fab.atomicAdd<RunOn::Device>(local_rho, tile_box, tile_box, 0, 0, ncomp);
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
            (*mf_pointer)[mfi].protected_divide<RunOn::Device>((*mf_pointer)[mfi],0,n,1);
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

        ParallelReduce::Max(stoptime, ParallelContext::IOProcessorNumberSub(),
                            ParallelContext::CommunicatorSub());

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
        const Long np = particles.numParticles();

        int nComp = fab.nComp();
        AMREX_FOR_1D( np, i,
        {
            amrex_interpolate_cic(pstruct[i], nComp, fabarr, plo, dxi);
        });
    }
}
