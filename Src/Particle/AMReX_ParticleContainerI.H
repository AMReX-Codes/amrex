
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::do_tiling = false;

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
IntVect
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::tile_size { AMREX_D_DECL(1024000,8,8) };

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> :: SetParticleSize ()
{
    num_real_comm_comps = 0;
    for (int i = 0; i < NArrayReal; ++i) {
        if (communicate_real_comp[i]) ++num_real_comm_comps;
    }

    num_int_comm_comps = 0;
    for (int i = 0; i < NArrayInt; ++i) {
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
        if ( ( not std::is_standard_layout<ParticleType>::value    ) or 
             ( not AMREX_IS_TRIVIALLY_COPYABLE(ParticleType)       )  )
        {
            amrex::Abort("Particle type must be standard layout and trivially copyable.");
        }

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
int
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::getTileIndex(const IntVect& iv, const Box& box, Box& tbx)
{
    if (do_tiling == false) {
        tbx = box;
        return 0;
    } else {
        //
        // This function must be consistent with FabArrayBase::buildTileArray function!!!
        //
        auto tiling_1d = [](int i, int lo, int hi, int tilesize,
                            int& ntile, int& tileidx, int& tlo, int& thi) {
            int ncells = hi-lo+1;
            ntile = std::max(ncells/tilesize, 1);
            int ts_right = ncells/ntile;
            int ts_left  = ts_right+1;
            int nleft = ncells - ntile*ts_right;
	    int ii = i - lo;
            int nbndry = nleft*ts_left;
            if (ii < nbndry) {
                tileidx = ii / ts_left; // tiles on the left of nbndry have size of ts_left
                tlo = lo + tileidx * ts_left;
                thi = tlo + ts_left - 1;
            } else {
                tileidx = nleft + (ii-nbndry) / ts_right;  // tiles on the right: ts_right
                tlo = lo + tileidx * ts_right + nleft;
                thi = tlo + ts_right - 1;
            }
        };
        const IntVect& small = box.smallEnd();
        const IntVect& big   = box.bigEnd();
        IntVect ntiles, ivIndex, tilelo, tilehi;

        AMREX_D_TERM(int iv0 = std::min(std::max(iv[0], small[0]), big[0]);,
               int iv1 = std::min(std::max(iv[1], small[1]), big[1]);,
               int iv2 = std::min(std::max(iv[2], small[2]), big[2]););

        AMREX_D_TERM(tiling_1d(iv0, small[0], big[0], tile_size[0], ntiles[0], ivIndex[0], tilelo[0], tilehi[0]);,
               tiling_1d(iv1, small[1], big[1], tile_size[1], ntiles[1], ivIndex[1], tilelo[1], tilehi[1]);,
               tiling_1d(iv2, small[2], big[2], tile_size[2], ntiles[2], ivIndex[2], tilelo[2], tilehi[2]););

        tbx = Box(tilelo, tilehi);

        return AMREX_D_TERM(ivIndex[0], + ntiles[0]*ivIndex[1], + ntiles[0]*ntiles[1]*ivIndex[2]);
    }
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
                      pld.m_tile = getTileIndex(iv, pld.m_gridbox, pld.m_tilebox);
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
	  pld.m_tile = getTileIndex(iv, bx, pld.m_tilebox);
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
                iv = Index(p, lev);
                grid = (*redistribute_mask_ptr)[local_grid](iv, 0);
                iv = Index(p_prime, lev);
	    }
            
            if (grid >= 0) {
                AMREX_D_TERM(p.m_rdata.pos[0] = p_prime.m_rdata.pos[0];,
                             p.m_rdata.pos[1] = p_prime.m_rdata.pos[1];,
                             p.m_rdata.pos[2] = p_prime.m_rdata.pos[2];);
                
                const Box& bx = ba.getCellCenteredBox(grid);
                
                pld.m_lev  = lev;
                pld.m_grid = grid;
		pld.m_tile = getTileIndex(iv, bx, pld.m_tilebox);
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
            BL_ASSERT(p.m_rdata.pos[i] >= geom.ProbLo(i));            
            shifted = true;
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            while (p.m_rdata.pos[i] <  geom.ProbLo(i)) p.m_rdata.pos[i] += geom.ProbLength(i); 

            // This is necessary when the arithmetic (esp. single precision) is such that the particle 
            //    starts only slightly outside the low boundary (eg slightly negative)
            //    ends up exactly on the high boundary 
            //    needs to get shifted back to the low boundary or it will be deemed out of bounds
            if (p.m_rdata.pos[i] == geom.ProbHi(i)) p.m_rdata.pos[i]  = geom.ProbLo(i);

            BL_ASSERT(p.m_rdata.pos[i] <= geom.ProbHi(i));
            shifted = true;
        }
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
    bool outside = AMREX_D_TERM(   p.m_rdata.pos[0] <  Geometry::ProbLo(0)
                          || p.m_rdata.pos[0] >= Geometry::ProbHi(0),
                          || p.m_rdata.pos[1] <  Geometry::ProbLo(1)
                          || p.m_rdata.pos[1] >= Geometry::ProbHi(1),
                          || p.m_rdata.pos[2] <  Geometry::ProbLo(2)
                          || p.m_rdata.pos[2] >= Geometry::ProbHi(2));

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
                for (const auto& p : ptile.GetArrayOfStructs()) {
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
                for (const auto& p : ptile.GetArrayOfStructs()) {
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
    const std::size_t sz = sizeof(ParticleType) + NArrayReal*sizeof(Real) + NArrayInt*sizeof(int);

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
      for (const auto& p : pbox) {
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
      for (const auto& p : pbox) {
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
                virts.push_back(*it);
                virts.back().m_idata.id = VirtualParticleID;
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
                    virts.push_back(*it);
                    virts.back().m_idata.id = VirtualParticleID;
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

//
// This redistributes valid particles and discards invalid ones.
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Redistribute (int lev_min, int lev_max, int nGrow, int local)
{
#ifdef AMREX_USE_CUDA
    if (local and (lev_min == 0) and (lev_max == 0) and (nGrow == 0))
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

#ifdef AMREX_USE_CUDA

namespace {
    
    struct DeviceBox
    {
        int lo[AMREX_SPACEDIM];
        int hi[AMREX_SPACEDIM];
        
        DeviceBox() 
        {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) 
            {
                lo[i] = 0;
                hi[i] = 0;
            }
        }
        
        DeviceBox(const Box& a_box)
        {
            for (int i = 0; i < AMREX_SPACEDIM; ++i)
            {
                lo[i] = a_box.smallEnd(i);
                hi[i] = a_box.bigEnd(i);
            }
        }
    };

    struct DeviceDomain
    {
        Real   left_edge[AMREX_SPACEDIM];
        Real  right_edge[AMREX_SPACEDIM];
        Real  inverse_dx[AMREX_SPACEDIM];
        Real domain_size[AMREX_SPACEDIM];
        int  is_periodic[AMREX_SPACEDIM];
        int    small_end[AMREX_SPACEDIM];
        int      big_end[AMREX_SPACEDIM];

        DeviceDomain(const Geometry& geom) {
            for (int i = 0; i < AMREX_SPACEDIM; ++i) {
                left_edge[i]   = geom.ProbLo(i);
                right_edge[i]  = geom.ProbHi(i);
                domain_size[i] = right_edge[i] - left_edge[i];
                inverse_dx[i]  = geom.InvCellSize(i);
                is_periodic[i] = geom.isPeriodic(i);
                small_end[i]   = geom.Domain().smallEnd(i);
                big_end[i]     = geom.Domain().bigEnd(i);
            }
        }
    };

    struct applyPeriodicBCs
    {
        DeviceDomain domain;

        __host__ __device__
        applyPeriodicBCs(const DeviceDomain& a_domain)
            : domain(a_domain) {}

        template <typename Particle>
        __host__ __device__
        Particle operator()(const Particle& in) const {
            Particle p = in;
            
            for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                if (not domain.is_periodic[dim]) continue;
                if (p.m_rdata.pos[dim] > domain.right_edge[dim]) { 
                    p.m_rdata.pos[dim] -= domain.domain_size[dim];
                } else if (p.m_rdata.pos[dim] < domain.left_edge[dim]) {
                    p.m_rdata.pos[dim] += domain.domain_size[dim];
                }
            }
            
            return p;
        }
    };

    struct assignParticle
    {        
        DeviceDomain domain;
        DeviceBox box;
        const int* mask_ptr;
        
        __host__ __device__
        assignParticle(const DeviceDomain& a_domain,
                       const DeviceBox&    a_box,
                       const int*          a_mask_ptr) 
            : domain(a_domain), box(a_box), mask_ptr(a_mask_ptr) {}
        
        template <typename Particle>
        __host__ __device__
        int operator()(const Particle& p) const {
            int i = floor((p.m_rdata.pos[0] - domain.left_edge[0]) * domain.inverse_dx[0]);
            i += domain.small_end[0];
#if   AMREX_SPACEDIM > 1
            int j = floor((p.m_rdata.pos[1] - domain.left_edge[1]) * domain.inverse_dx[1]);
            j += domain.small_end[1];
#endif
#if   AMREX_SPACEDIM > 2
            int k = floor((p.m_rdata.pos[2] - domain.left_edge[2]) * domain.inverse_dx[2]);
            k += domain.small_end[2];
#endif
            
            long offset = i - box.lo[0];
#if   AMREX_SPACEDIM==2
            offset += (box.hi[0] - box.lo[0] + 1) * (j - box.lo[1]);
#elif AMREX_SPACEDIM==3
            offset += (box.hi[0] - box.lo[0] + 1) * (j - box.lo[1]
                                                     + (k - box.lo[2]) * (box.hi[1] - box.lo[1] + 1));
#endif
            return mask_ptr[offset];
        }
    };
    
    struct grid_is
    {
        int grid_id;
        
        __host__ __device__
        grid_is(int a_grid_id) : grid_id(a_grid_id) {}
        
        __host__ __device__
        bool operator()(int grid) const
        {
            return grid_id == grid;
        }
    };
    
    struct grid_is_not
    {
        int grid_id;
        
        __host__ __device__
        grid_is_not(int a_grid_id) : grid_id(a_grid_id) {}
    
        __host__ __device__
        bool operator()(int grid) const
        {
            return grid_id != grid;
        }
    };
}

#endif // AMREX_USE_CUDA

//
// The GPU implementation of Redistribute
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RedistributeGPU (int lev_min, int lev_max, int nGrow, int local)
{
#ifdef AMREX_USE_CUDA

    // sanity checks
    AMREX_ASSERT(local);
    AMREX_ASSERT(lev_min == 0);
    AMREX_ASSERT(lev_max == 0);
    AMREX_ASSERT(nGrow == 0);
    AMREX_ASSERT(do_tiling == false);

    BL_PROFILE("ParticleContainer::RedistributeGPU()");
    BL_PROFILE_VAR_NS("compute_grid", blp_transform);
    BL_PROFILE_VAR_NS("compute_partition", blp_compute_partition);
    BL_PROFILE_VAR_NS("copy_partition", blp_copy_partition);
    BL_PROFILE_VAR_NS("copy_buffer", blp_copy_buffer);
    BL_PROFILE_VAR_NS("sort_by_dest", blp_sort);
    BL_PROFILE_VAR_NS("final_copy", blp_copy);
    BL_PROFILE_VAR_NS("resize", blp_resize);

    if (local > 0) BuildRedistributeMask(0, local);

    const int lev = 0;
    const Geometry& geom = Geom(lev);
    const BoxArray& ba   = ParticleBoxArray(lev);
    const DistributionMapping& dmap = ParticleDistributionMap(lev);
    auto& plev  = m_particles[lev];

    m_grids_tmp.resize(0);
    m_index_sequence_tmp.resize(0);

    // temporaries
    m_grids_tmp.resize(0);
    m_index_sequence_tmp.resize(0);
    m_aos_r.resize(0);
    m_grids_r.resize(0);
    m_rdata_r.resize(0);
    m_idata_r.resize(0);

    m_aos_to_redistribute.resize(0);
    for (int k = 0; k < NArrayReal; ++k) m_real_arrays_to_redistribute[k].resize(0);
    for (int k = 0; k < NArrayInt; ++k) m_int_arrays_to_redistribute[k].resize(0);
    m_grids_to_redistribute.resize(0);

    //
    // First pass - figure out what grid each particle should be on, and copy
    // the ones that move into a temporary buffer
    //
    for(MFIter mfi(*redistribute_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        int src_grid = mfi.index();
        const int tid = 0;        
        auto& ptile = plev[std::make_pair(src_grid, tid)];
        auto& aos   = ptile.GetArrayOfStructs();
        auto& soa   = ptile.GetStructOfArrays();
        const size_t old_np = aos.numParticles();

        m_grids_tmp.resize(old_np);
        m_index_sequence_tmp.resize(old_np);

        thrust::sequence(thrust::device, 
			 m_index_sequence_tmp.begin(), 
			 m_index_sequence_tmp.end());        

        BL_PROFILE_VAR_START(blp_transform);
        //
        // Compute the grid each particle belongs to
        //
        thrust::transform(thrust::device, 
                          thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data())),
                          thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data()) + aos.size()),
                          m_grids_tmp.begin(),
                          assignParticle(DeviceDomain(geom),
                                         DeviceBox((*redistribute_mask_ptr)[src_grid].box()),
                                         (*redistribute_mask_ptr)[src_grid].dataPtr()));        

        BL_PROFILE_VAR_STOP(blp_transform);

        BL_PROFILE_VAR_START(blp_compute_partition);
        //
        // Partition the particle data so that particles that stay come first,
        // those that move go to the end.
        //
        auto mid   = thrust::partition(thrust::device,
                                       m_index_sequence_tmp.begin(),
                                       m_index_sequence_tmp.end(), 
                                       m_grids_tmp.begin(), grid_is(src_grid));
        BL_PROFILE_VAR_STOP(blp_compute_partition);

        const size_t num_moved = thrust::distance(mid, m_index_sequence_tmp.end());
        const size_t new_np   = old_np - num_moved;

        if (num_moved == 0) continue;  // early exit if everything is already in right place
        
        const size_t old_size = m_aos_to_redistribute.size();
        const size_t new_size = old_size + num_moved;
        
        m_aos_to_redistribute.resize(new_size);	
        for (int k = 0; k < NArrayReal; ++k) m_real_arrays_to_redistribute[k].resize(new_size);
        for (int k = 0; k < NArrayInt; ++k) m_int_arrays_to_redistribute[k].resize(new_size);
        m_grids_to_redistribute.resize(new_size);

        BL_PROFILE_VAR_START(blp_copy_partition);        
        //
        // Reorder the particle data based on the partition we just computed
        //
        {
	    // reorder structs
	    m_aos_r.resize(old_np);
            thrust::copy_n(thrust::device,
                thrust::make_permutation_iterator(
                thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data())),
                                                 m_index_sequence_tmp.begin()), 
                old_np, m_aos_r.begin());
            thrust::copy(thrust::device, m_aos_r.begin(), m_aos_r.end(), 
                         thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data())));

            // reorder grids
            m_grids_r.resize(old_np);
            thrust::copy_n(thrust::device,
                           thrust::make_permutation_iterator(m_grids_tmp.begin(),
                                                             m_index_sequence_tmp.begin()), 
                           old_np, m_grids_r.begin());
            thrust::copy(thrust::device, m_grids_r.begin(), m_grids_r.end(), m_grids_tmp.begin());

            // reorder real arrays
	    m_rdata_r.resize(old_np);
            for (int j = 0; j < NArrayReal; ++j)
            {
                auto& rdata = ptile.GetStructOfArrays().GetRealData(j);
                thrust::copy_n(thrust::device,
                               thrust::make_permutation_iterator(thrust::device_ptr<Real>(rdata.data()),
                                                       m_index_sequence_tmp.begin()), 
                     old_np, m_rdata_r.begin());
                thrust::copy(thrust::device, m_rdata_r.begin(), m_rdata_r.end(), 
                             thrust::device_ptr<Real>(rdata.data()));
            }

            // reorder int arrays
	    m_idata_r.resize(old_np);
            for (int j = 0; j < NArrayInt; ++j)
            {
                auto& idata = ptile.GetStructOfArrays().GetIntData(j);
                thrust::copy_n(thrust::device,
                               thrust::make_permutation_iterator(thrust::device_ptr<int>(idata.data()),
                                                                 m_index_sequence_tmp.begin()), 
                    old_np, m_idata_r.begin());
                thrust::copy(thrust::device, m_idata_r.begin(), m_idata_r.end(), 
                             thrust::device_ptr<int>(idata.data()));
            }
        }
        BL_PROFILE_VAR_STOP(blp_copy_partition);

        BL_PROFILE_VAR_START(blp_copy_buffer);
        //
        // copy the particle data to be moved into temp buffers
        //
        {
            auto& dst = m_aos_to_redistribute;
            thrust::copy(thrust::device,
                         thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data()) + new_np),
                         thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data()) + old_np),
			 dst.begin() + old_size);
        }
        {
            thrust::copy(thrust::device,
                         m_grids_tmp.begin() + new_np,
                         m_grids_tmp.begin() + old_np,
                         m_grids_to_redistribute.begin() + old_size);
        }
        for (int j = 0; j < NArrayReal; ++j)
        {
  	    auto& src = soa.GetRealData(j);
	    auto& dst = m_real_arrays_to_redistribute[j];
            thrust::copy(thrust::device, 
                         thrust::device_ptr<Real>(src.data() + new_np),
                         thrust::device_ptr<Real>(src.data() + old_np),
			 dst.begin() + old_size);
        }
        for (int j = 0; j < NArrayInt; ++j)
        {
            auto& src = soa.GetIntData(j); 
	    auto& dst = m_int_arrays_to_redistribute[j];
            thrust::copy(thrust::device, 
                         thrust::device_ptr<int>(src.data() + new_np),
                         thrust::device_ptr<int>(src.data() + old_np),
                         dst.begin() + old_size);
        }

        ptile.resize(new_np);

        BL_PROFILE_VAR_STOP(blp_copy_buffer);
    }

    //
    // We now have a temporary buffer that holds all the particles that need to be 
    // copied moved somewhere else. Now we sort those particles by their destination grid
    //
    std::map<int, thrust::device_vector<char> > not_ours;
    const int num_grids = ba.size();
    const int num_to_move = m_grids_to_redistribute.size();

    if (num_to_move > 0)
    {
        BL_PROFILE_VAR_START(blp_sort);

        thrust::device_vector<int> grid_begin(num_grids);
        thrust::device_vector<int> grid_end(num_grids);

        m_index_sequence_tmp.resize(num_to_move);
        thrust::sequence(thrust::device, 
			 m_index_sequence_tmp.begin(),
			 m_index_sequence_tmp.end());
        thrust::sort_by_key(thrust::device,
                            m_grids_to_redistribute.begin(),
                            m_grids_to_redistribute.end(),
                            m_index_sequence_tmp.begin());

        //
        // Reorder particle data using the computed sequence
        //
        {
            // reorder structs
            auto& aos = m_aos_to_redistribute;
            m_aos_r.resize(num_to_move);
            thrust::copy_n(thrust::device,
			   thrust::make_permutation_iterator(aos.begin(),
							     m_index_sequence_tmp.begin()),
			   num_to_move, m_aos_r.begin());
            thrust::copy(thrust::device, m_aos_r.begin(), m_aos_r.end(), aos.begin());

            // reorder real arrays
	    m_rdata_r.resize(num_to_move);
            for (int j = 0; j < NArrayReal; ++j)
            {
                auto& rdata = m_real_arrays_to_redistribute[j];
                thrust::copy_n(thrust::device,
                               thrust::make_permutation_iterator(rdata.begin(),
                                                                 m_index_sequence_tmp.begin()), 
			       num_to_move, m_rdata_r.begin());
                thrust::copy(thrust::device, m_rdata_r.begin(), m_rdata_r.end(), rdata.begin());
            }

            // reorder int arrays
	    m_idata_r.resize(num_to_move);
            for (int j = 0; j < NArrayInt; ++j)
            {
                auto& idata = m_int_arrays_to_redistribute[j];
                thrust::copy_n(thrust::device,
                               thrust::make_permutation_iterator(idata.begin(),
                                                                 m_index_sequence_tmp.begin()), 
			       num_to_move, m_idata_r.begin());
                thrust::copy(thrust::device, m_idata_r.begin(), m_idata_r.end(), idata.begin());
            }
        }

        //
        // Now we compute the start and stop index in the sorted array for each grid
        //
        thrust::counting_iterator<unsigned> search_begin(0);
        thrust::lower_bound(thrust::device,
                            m_grids_to_redistribute.begin(),
                            m_grids_to_redistribute.end(),
                            search_begin,
                            search_begin + num_grids,
                            grid_begin.begin());
        
        thrust::upper_bound(thrust::device,
                            m_grids_to_redistribute.begin(),
                            m_grids_to_redistribute.end(),
                            search_begin,
                            search_begin + num_grids,
                            grid_end.begin());
        
        thrust::host_vector<int> start(grid_begin);
        thrust::host_vector<int> stop(grid_end);
        
        std::map<int, size_t> grid_counts;
        for (int i = 0; i < num_grids; ++i)
        {
            const int dest_proc = dmap[i];
            if (dest_proc != ParallelDescriptor::MyProc())
            {
                grid_counts[dest_proc] += 1;
            }
        }
        BL_PROFILE_VAR_STOP(blp_sort);

        BL_PROFILE_VAR_START(blp_copy);
        
        //
        // Each destination grid, copy the appropriate particle data, passing the non-local data
        // into not_ours
        //
        for (int i = 0; i < num_grids; ++i)
        {
            const int tid = 0;
            auto pair_index = std::make_pair(i, tid);
            
            const size_t num_to_add = stop[i] - start[i];

            const int dest_proc = dmap[i];
            if (dest_proc == ParallelDescriptor::MyProc())  // this is a local copy
            {
                auto& ptile = plev[pair_index];
                const size_t old_size = ptile.numParticles();
                const size_t new_size = old_size + num_to_add;
                ptile.resize(new_size);
                
                // copy structs
                {
		    auto& src = m_aos_to_redistribute;
                    auto& dst = ptile.GetArrayOfStructs();
                    thrust::copy(thrust::device,
                                 src.begin() + start[i],
				 src.begin() + stop[i],
                                 thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(dst().data()) + old_size));
                }
                
                // copy arrays
                for (int j = 0; j < NArrayReal; ++j)
                {
		    auto& src = m_real_arrays_to_redistribute[j];
		    auto& dst = ptile.GetStructOfArrays().GetRealData(j);
                    thrust::copy(thrust::device,
				 src.begin() + start[i],
				 src.begin() + stop[i],
				 thrust::device_ptr<Real>(dst.data() + old_size));
                }
                for (int j = 0; j < NArrayInt; ++j)
                {
		    auto& src = m_int_arrays_to_redistribute[j];
		    auto& dst = ptile.GetStructOfArrays().GetIntData(j);
                    thrust::copy(thrust::device,
				 src.begin() + start[i],
				 src.begin() + stop[i],
				 thrust::device_ptr<int>(dst.data() + old_size));
                }
            }
            else // this is the non-local case
            {
                char* dst;
                const size_t old_size = not_ours[dest_proc].size();
                const size_t new_size
                    = old_size + num_to_add*superparticle_size + sizeof(size_t) + 2*sizeof(int);

                if (old_size == 0)
                {
                    not_ours[dest_proc].resize(new_size + sizeof(size_t));
                    cudaMemcpy(thrust::raw_pointer_cast(not_ours[dest_proc].data()), 
                               &grid_counts[dest_proc], sizeof(size_t), cudaMemcpyHostToDevice);
                    dst = thrust::raw_pointer_cast(
                                                   not_ours[dest_proc].data() + old_size + sizeof(size_t));
                } else
                {
                    not_ours[dest_proc].resize(new_size);
                    dst = thrust::raw_pointer_cast(not_ours[dest_proc].data() + old_size);
                }

                cudaMemcpy(thrust::raw_pointer_cast(dst), 
                           &num_to_add, sizeof(size_t), cudaMemcpyHostToDevice);
                dst += sizeof(size_t);
                
                cudaMemcpy(thrust::raw_pointer_cast(dst), &i, sizeof(int), cudaMemcpyHostToDevice);
                dst += sizeof(int);
                
                cudaMemcpy(thrust::raw_pointer_cast(dst), 
                           &dest_proc, sizeof(int), cudaMemcpyHostToDevice);
                dst += sizeof(int);
                
                // pack structs
                {
		    auto& aos = m_aos_to_redistribute;
                    cudaMemcpy(thrust::raw_pointer_cast(dst), 
                               thrust::raw_pointer_cast(aos.data() + start[i]),
                               num_to_add*sizeof(ParticleType), cudaMemcpyDeviceToDevice);
                    dst += num_to_add*sizeof(ParticleType);
                }

                // pack real arrays
                for (int j = 0; j < NArrayReal; ++j)
                {
                    if (not communicate_real_comp[j]) continue;
                    auto& attrib = m_real_arrays_to_redistribute[j];
                    cudaMemcpy(thrust::raw_pointer_cast(dst), 
                               thrust::raw_pointer_cast(attrib.data() + start[i]),
                               num_to_add*sizeof(Real), cudaMemcpyDeviceToDevice);
                    dst += num_to_add*sizeof(Real);
                }
                
                // pack int arrays           
                for (int j = 0; j < NArrayInt; ++j)
                {
                    if (not communicate_int_comp[j]) continue;
                    auto& attrib = m_int_arrays_to_redistribute[j];
                    cudaMemcpy(thrust::raw_pointer_cast(dst),
                               thrust::raw_pointer_cast(attrib.data() + start[i]),
                               num_to_add*sizeof(int), cudaMemcpyDeviceToDevice);
                    dst += num_to_add*sizeof(int);
                }
            }
        }
        BL_PROFILE_VAR_STOP(blp_copy);
    }

    RedistributeMPIGPU(not_ours);

    EnforcePeriodicGPU();

    AMREX_ASSERT(OK());

#endif // AMREX_USE_CUDA
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

    BuildRedistributeMask(0, 1);

    for (MFIter mfi(*redistribute_mask_ptr, false); mfi.isValid(); ++mfi)
    {
        const int grid_id = mfi.index();
        const int tile_id = mfi.LocalTileIndex();
        auto& particles = ParticlesAt(lev, grid_id, tile_id);
        auto& aos = particles.GetArrayOfStructs();

        thrust::transform(thrust::device, 
                          thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data())),
                          thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data()) + aos.size()),
                          thrust::device_ptr<ParticleType>(static_cast<ParticleType*>(aos().data())),
                          applyPeriodicBCs(geom));
    }
}
#endif //AMREX_USE_CUDA

#ifdef AMREX_USE_CUDA
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::RedistributeMPIGPU (std::map<int, thrust::device_vector<char> >& not_ours)
{
    BL_PROFILE("ParticleContainer::RedistributeMPIGPU()");

#if BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    const int lev = 0;

    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;

    for (const auto& kv : not_ours)
    {
        const int np = kv.second.size(); 
        Snds[kv.first] = np;
        NumSnds += np;
    }
    
    ParallelDescriptor::ReduceLongMax(NumSnds);

    if (NumSnds == 0) return;

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());
    
    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 ParallelDescriptor::Communicator()) );
    
    BL_ASSERT(Rcvs[ParallelDescriptor::MyProc()] == 0);
    
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += Rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();
    
    // Allocate data for rcvs as one big chunk.
    thrust::device_vector<char> recvdata(TotRcvBytes);
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = Rcvs[Who];
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv(thrust::raw_pointer_cast(recvdata.data() + offset),
                                             Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (const auto& kv : not_ours) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send(thrust::raw_pointer_cast(kv.second.data()),
                                 Cnt, Who, SeqNum);
    }

    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);

        for (int i = 0; i < nrcvs; ++i) {
            const int offset = rOffset[i];
            char* buffer = thrust::raw_pointer_cast(recvdata.data() + offset);
            size_t num_grids, num_particles;
            int gid, pid;
            cudaMemcpy(&num_grids, buffer, sizeof(size_t), cudaMemcpyDeviceToHost); 
            buffer += sizeof(size_t);

            for (int g = 0; g < num_grids; ++g) {
                cudaMemcpy(&num_particles, buffer, sizeof(size_t), cudaMemcpyDeviceToHost);
                buffer += sizeof(size_t);
                cudaMemcpy(&gid, buffer, sizeof(int), cudaMemcpyDeviceToHost);
                buffer += sizeof(int);
                cudaMemcpy(&pid, buffer, sizeof(int), cudaMemcpyDeviceToHost);
                buffer += sizeof(int);

                if (num_particles == 0) continue;

                AMREX_ALWAYS_ASSERT(pid == ParallelDescriptor::MyProc());
                {
                    const int tid = 0;
                    auto pair_index = std::make_pair(gid, tid);
                    auto& ptile = m_particles[lev][pair_index];
                    auto& aos = ptile.GetArrayOfStructs();
                    auto& soa = ptile.GetStructOfArrays();
                    const size_t old_size = ptile.numParticles();
                    const size_t new_size = old_size + num_particles;
                    ptile.resize(new_size);
                   
                    //copy structs
                    cudaMemcpy(static_cast<ParticleType*>(aos().data()) + old_size,
                               buffer, num_particles*sizeof(ParticleType),
                               cudaMemcpyDeviceToDevice);
                    buffer += num_particles*sizeof(ParticleType);

                    // copy real arrays 
                    for (int j = 0; j < NArrayReal; ++j) {
                        if (not communicate_real_comp[j]) continue;
                        auto& attrib = soa.GetRealData(j);
                        cudaMemcpy(attrib.data() + old_size, buffer, num_particles*sizeof(Real),
                                   cudaMemcpyDeviceToDevice);
                        buffer += num_particles*sizeof(Real);
                    }

                    // copy int arrays 
                    for (int j = 0; j < NArrayInt; ++j) {
                        if (not communicate_int_comp[j]) continue;
                        auto& attrib = soa.GetIntData(j);
                        cudaMemcpy(attrib.data() + old_size, buffer, num_particles*sizeof(int),
                                   cudaMemcpyDeviceToDevice);
                        buffer += num_particles*sizeof(int);
                    }
                }            
            }
        }
    }        
#endif // MPI    
}
#endif // AMREX_USE_CUDA

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
                      for (int comp = 0; comp < NArrayReal; comp++)
                          soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                      for (int comp = 0; comp < NArrayInt; comp++)
                          soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                      correctCellVectors(last, pindex, grid, aos[pindex]);
                      --last;
                      continue;
                  }
                      
                  locateParticle(p, pld, lev_min, lev_max, nGrow, local ? grid : -1);
                  
                  if (p.m_idata.id < 0)
                  {
                      aos[pindex] = aos[last];
                      for (int comp = 0; comp < NArrayReal; comp++)
                          soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                      for (int comp = 0; comp < NArrayInt; comp++)
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
                          for (int comp = 0; comp < NArrayReal; ++comp) {
                              RealVector& arr = soa_local[pld.m_lev][index][thread_num].GetRealData(comp);
                              arr.push_back(soa.GetRealData(comp)[pindex]);
                          }
                          for (int comp = 0; comp < NArrayInt; ++comp) {
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
                      for (int comp = 0; comp < NArrayReal; comp++) {
                          if (communicate_real_comp[comp]) {
                              std::memcpy(dst, &soa.GetRealData(comp)[pindex], sizeof(Real));
                              dst += sizeof(Real);
                          }
                      }
                      for (int comp = 0; comp < NArrayInt; comp++) {
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
                      for (int comp = 0; comp < NArrayReal; comp++)
                          soa.GetRealData(comp)[pindex] = soa.GetRealData(comp)[last];
                      for (int comp = 0; comp < NArrayInt; comp++)
                          soa.GetIntData(comp)[pindex] = soa.GetIntData(comp)[last];
                      correctCellVectors(last, pindex, grid, aos[pindex]);
                      --last;
                      continue;
                  }
                  
                  ++pindex;
              }
              
              aos().erase(aos().begin() + last + 1, aos().begin() + npart);
              for (int comp = 0; comp < NArrayReal; comp++) {
                  RealVector& rdata = soa.GetRealData(comp);
                  rdata.erase(rdata.begin() + last + 1, rdata.begin() + npart);
              }
              for (int comp = 0; comp < NArrayInt; comp++) {
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
          auto& aos = m_particles[lev][index].GetArrayOfStructs();
          auto& soa = m_particles[lev][index].GetStructOfArrays();
          auto& aos_tmp = *(pvec_ptrs[pit]);
          auto& soa_tmp = soa_local[lev][index];
          for (int i = 0; i < num_threads; ++i) {
              aos.insert(aos.end(), aos_tmp[i].begin(), aos_tmp[i].end());
              aos_tmp[i].erase(aos_tmp[i].begin(), aos_tmp[i].end());
              for (int comp = 0; comp < NArrayReal; ++comp) {
                  RealVector& arr = soa.GetRealData(comp);
                  RealVector& tmp = soa_tmp[i].GetRealData(comp);
                  arr.insert(arr.end(), tmp.begin(), tmp.end());
                  tmp.erase(tmp.begin(), tmp.end());
              }
              for (int comp = 0; comp < NArrayInt; ++comp) {
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
BuildRedistributeMask (int lev, int nghost) const {
    
    BL_PROFILE("ParticleContainer::BuildRedistributeMask");
    BL_ASSERT(lev == 0);

    if (redistribute_mask_ptr == nullptr ||
        ! BoxArray::SameRefs(redistribute_mask_ptr->boxArray(), this->ParticleBoxArray(lev)) ||
        ! DistributionMapping::SameRefs(redistribute_mask_ptr->DistributionMap(), this->ParticleDistributionMap(lev)))
        {
            const Geometry& geom = this->Geom(lev);
            const BoxArray& ba = this->ParticleBoxArray(lev);
            const DistributionMapping& dmap = this->ParticleDistributionMap(lev);
            
            redistribute_mask_ptr.reset(new iMultiFab(ba, dmap, 2, nghost));
            redistribute_mask_ptr->setVal(-1, nghost);
            
#ifdef _OPENMP
#pragma omp parallel
#endif
            for (MFIter mfi(*redistribute_mask_ptr, this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
                 mfi.isValid(); ++mfi) {
                const Box& box = mfi.tilebox();
                const int grid_id = mfi.index();
                const int tile_id = mfi.LocalTileIndex();
                redistribute_mask_ptr->setVal(grid_id, box, 0, 1);
                redistribute_mask_ptr->setVal(tile_id, box, 1, 1);
            }
            
            redistribute_mask_ptr->FillBoundary(geom.periodicity());
            
            neighbor_procs.clear();
            for (MFIter mfi(*redistribute_mask_ptr, this->do_tiling ? this->tile_size : IntVect::TheZeroVector());
                 mfi.isValid(); ++mfi) {
                const Box& box = mfi.growntilebox();
                for (IntVect iv = box.smallEnd(); iv <= box.bigEnd(); box.next(iv)) {
                    const int grid = (*redistribute_mask_ptr)[mfi](iv, 0);
                    if (grid >= 0) {
                        const int proc = this->ParticleDistributionMap(lev)[grid];
                        if (proc != ParallelDescriptor::MyProc())  
                            neighbor_procs.push_back(proc);
                    }
                }
            }
            
            RemoveDuplicates(neighbor_procs);

#ifdef AMREX_USE_CUDA
	    Box bx = ba[0];	    
	    m_aos_to_redistribute.reserve(bx.numPts());
	    for (int k = 0; k < NArrayReal; ++k) 
	        m_real_arrays_to_redistribute[k].reserve(bx.numPts());
	    for (int k = 0; k < NArrayInt; ++k)  
	        m_int_arrays_to_redistribute[k].reserve(bx.numPts());
	    m_grids_to_redistribute.reserve(bx.numPts());
#endif
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

#if BL_USE_MPI

    const int NProcs = ParallelDescriptor::NProcs();
    const int NNeighborProcs = neighbor_procs.size();
    
    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;
    if (local > 0) {
        AMREX_ALWAYS_ASSERT(lev_min == 0);
        AMREX_ALWAYS_ASSERT(lev_max == 0);
        BuildRedistributeMask(0, local);
        NumSnds = doHandShakeLocal(not_ours, neighbor_procs, Snds, Rcvs);
    }
    else {
        NumSnds = doHandShake(not_ours, Snds, Rcvs);
    }

    const int SeqNum = ParallelDescriptor::SeqNum();
    
    if ((not local) and NumSnds == 0)
        return;  // There's no parallel work to do.

    if (local) {
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
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += Rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);
    
    // Allocate data for rcvs as one big chunk.
    Vector<char> recvdata(TotRcvBytes);
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = Rcvs[Who];
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (const auto& kv : not_ours) {
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
   
        if (recvdata.size() % superparticle_size != 0) {
            if (m_verbose) {
                amrex::AllPrint() << "ParticleContainer::RedistributeMPI: sizes = "
                                  << recvdata.size() << ", " << superparticle_size << "\n";
            }
            amrex::Abort("ParticleContainer::RedistributeMPI: How did this happen?");
        }

        int npart = recvdata.size() / superparticle_size;
        
        Vector<int> rcv_levs(npart);
        Vector<int> rcv_grid(npart);
        Vector<int> rcv_tile(npart);

        ParticleLocData pld;
#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif                
        for (int i = 0; i < npart; ++i)
        {
            char* pbuf = recvdata.data() + i*superparticle_size;
            ParticleType p;
            std::memcpy(&p, pbuf, sizeof(ParticleType));
            locateParticle(p, pld, lev_min, lev_max, nGrow);
            rcv_levs[i] = pld.m_lev;
            rcv_grid[i] = pld.m_grid;
            rcv_tile[i] = pld.m_tile;
        }

	BL_PROFILE_VAR_STOP(blp_locate);

        BL_PROFILE_VAR_START(blp_copy);

        for (int j = 0; j < npart; ++j)
        {
            auto& ptile = m_particles[rcv_levs[j]][std::make_pair(rcv_grid[j], rcv_tile[j])];
            
            char* pbuf = recvdata.data() + j*superparticle_size;
            ParticleType p;
            std::memcpy(&p, pbuf, sizeof(ParticleType));            
            ptile.push_back(p);

            Real* rdata = (Real*)(pbuf + particle_size);
            for (int comp = 0; comp < NArrayReal; ++comp) {
                if (communicate_real_comp[comp]) {
                    ptile.push_back_real(comp, *rdata++);
                } else {
                    ptile.push_back_real(comp, 0.0);
                }
            }
            
            int* idata = (int*)(pbuf + particle_size + num_real_comm_comps*sizeof(Real));
            for (int comp = 0; comp < NArrayInt; ++comp) {
                if (communicate_int_comp[comp]) {
                    ptile.push_back_int(comp, *idata++);
                } else {
                    ptile.push_back_int(comp, 0);
                }
            }            
        }

	BL_PROFILE_VAR_STOP(blp_copy);
    }
#endif /*BL_USE_MPI*/
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::OK (int lev_min, int lev_max, int nGrow) const
{
#ifdef AMREX_USE_CUDA
    if ( (finestLevel() == 0) and (nGrow == 0) )
    {
        return OKGPU(lev_min, lev_max, nGrow);
    }
    else
    {
        return OKCPU(lev_min, lev_max, nGrow);
    }
#else
    return OKCPU(lev_min, lev_max, nGrow);
#endif    
}

#ifdef AMREX_USE_CUDA
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::OKGPU (int lev_min, int lev_max, int nGrow) const
{
    BL_PROFILE("ParticleContainer::OKGPU()");

    if (lev_max == -1)
        lev_max = finestLevel();

    AMREX_ASSERT(lev_max <= finestLevel());

    AMREX_ASSERT(lev_min == 0);
    AMREX_ASSERT(lev_max == 0);
    AMREX_ASSERT(nGrow == 0);
    AMREX_ASSERT(do_tiling == false);

    const int lev = 0;
    const Geometry& geom = Geom(lev);

    BuildRedistributeMask(0, 1);

    thrust::device_vector<int> grid_indices;

    long total_np = 0;
    long total_wrong = 0;
    for(MFIter mfi(*redistribute_mask_ptr, false); mfi.isValid(); ++mfi) {
        int i = mfi.index();
        const int tid = 0;

        const auto& ptile = ParticlesAt(lev, i, tid);
        auto& aos   = ptile.GetArrayOfStructs();
        const int np = ptile.numParticles();

        total_np += np;

        if (np == 0) continue;
        
        grid_indices.resize(np);
        
        thrust::transform(thrust::device, 
                          thrust::device_ptr<const ParticleType>(static_cast<const ParticleType*>(aos().data())),
                          thrust::device_ptr<const ParticleType>(static_cast<const ParticleType*>(aos().data()) 
                                                   + aos.size()),
                          grid_indices.begin(),
                          assignParticle(DeviceDomain(geom),
                                         DeviceBox((*redistribute_mask_ptr)[i].box()),
                                         (*redistribute_mask_ptr)[i].dataPtr()));

        int count = thrust::count_if(grid_indices.begin(),
                                     grid_indices.end(),
                                     grid_is_not(i));

        total_wrong += count;
    }

    ParallelDescriptor::ReduceLongMax(total_np);
    ParallelDescriptor::ReduceLongMax(total_wrong);

    amrex::Print() << "I have " << total_np << " particles in OK(). \n"; 
    amrex::Print() << "I have " << total_wrong << " particles in the wrong place. \n"; 

    return (total_wrong == 0);
}
#endif

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::OKCPU (int lev_min, int lev_max, int nGrow) const
{
    BL_PROFILE("ParticleContainer::OKCPU()");
    if (lev_max == -1)
        lev_max = finestLevel();

    AMREX_ASSERT(lev_max <= finestLevel());

    ParticleLocData pld;
    for (int lev = lev_min; lev <= lev_max; lev++)
    {
        const auto& pmap = m_particles[lev];
	for (const auto& kv : pmap)
        {
            const int grid = kv.first.first;
            const int tile = kv.first.second;
            const auto& aos = kv.second.GetArrayOfStructs();
            const auto& soa = kv.second.GetStructOfArrays();
            
            int np = aos.numParticles();
            for (int i = 0; i < NArrayReal; i++) {
                BL_ASSERT(np == soa.GetRealData(i).size());
            }
            for (int i = 0; i < NArrayInt; i++) {
                BL_ASSERT(np == soa.GetIntData(i).size());
            }

            const BoxArray& ba = ParticleBoxArray(lev);
            BL_ASSERT(ba.ixType().cellCentered());
            for (const auto& p : aos)
            {
                if (p.m_idata.id > 0)
                {
                    if (grid < 0 || grid >= ba.size()) return false;

                    //
                    // First, make sure the particle COULD be in this container
                    // 
                    const IntVect& iv = Index(p, lev);
                    const Box& gridbox = ba.getCellCenteredBox(grid);
                    if (!amrex::grow(gridbox, nGrow).contains(iv))
                    {
                        return false;
                    }
                    Box tbx;
                    if (getTileIndex(iv, gridbox, tbx) != tile)
                    {
                        return false;
                    }

                    //
                    // Then, we need to make sure it cannot be stored in finer level 
                    // or valid box of current level.
                    //
                    if (Where(p, pld, lev_min, lev_max))
                    {
                        if (lev != pld.m_lev  || grid != pld.m_grid || tile != pld.m_tile)
                        {
                            if (m_verbose) {
                                amrex::AllPrint() << "PARTICLE NUMBER " << p.m_idata.id << '\n'
                                                  << "POS  " << AMREX_D_TERM(p.m_rdata.pos[0], << p.m_rdata.pos[1], << p.m_rdata.pos[2]) << "\n"
                                                  << "LEV  " << lev  << " " << pld.m_lev << '\n'
                                                  << "GRID " << grid << " " << pld.m_grid << '\n';
                            }
                            return false;
                        }
                    }
                }
            }
        }
    }
    
    return true;
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

    while (!particles.empty()) {
        ParticleType& p = particles.back();

        if (p.id() > 0)
        {
            if (!Where(p, pld, level, level, nGrow))
                amrex::Abort("ParticleContainerAddParticlesAtLevel(): Can't add outside of domain\n");
            m_particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
        }

        particles.pop_back();
    }
    Redistribute(level, level, nGrow);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteParticleRealData (void* data, size_t size,
                         std::ostream& os, const RealDescriptor& rd) const
{
    if (sizeof(typename ParticleType::RealType) == 4) {
        writeFloatData((float*) data, size, os, ParticleRealDescriptor);
    } 
    else if (sizeof(typename ParticleType::RealType) == 8) {
        writeDoubleData((double*) data, size, os, ParticleRealDescriptor);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::ReadParticleRealData (void* data, size_t size,
                         std::istream& is, const RealDescriptor& rd)
{
    if (sizeof(typename ParticleType::RealType) == 4) {
        readFloatData((float*) data, size, is, ParticleRealDescriptor);
    } 
    else if (sizeof(typename ParticleType::RealType) == 8) {
        readDoubleData((double*) data, size, is, ParticleRealDescriptor);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Checkpoint (const std::string& dir, const std::string& name, bool is_checkpoint,
              const Vector<std::string>& real_comp_names,
              const Vector<std::string>& int_comp_names) const
{
    BL_PROFILE("ParticleContainer::Checkpoint()");
    BL_ASSERT(OK());

    BL_ASSERT(sizeof(typename ParticleType::RealType) == 4 || sizeof(typename ParticleType::RealType) == 8);

    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProcNumber   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = amrex::second();

    std::string pdir = dir;
    if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
    pdir += name;
    
    if ( ! levelDirectoriesCreated) {
        if (ParallelDescriptor::IOProcessor()) 
            if ( ! amrex::UtilCreateDirectory(pdir, 0755)) 
                amrex::CreateDirectoryFailed(pdir);
        ParallelDescriptor::Barrier();
    }
    
    std::ofstream HdrFile;

    long nparticles = 0;
    int maxnextid;

    if(usePrePost) {
      nparticles = nparticlesPrePost;
      maxnextid  = maxnextidPrePost;
    } else {
      nparticles = 0;
      maxnextid  = ParticleType::NextID();

      for (int lev = 0; lev < m_particles.size();  lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
            for (const auto& p : aos) {
                // Only count (and checkpoint) valid particles.
                if (p.m_idata.id > 0) nparticles++;
            }
        }
      }
      ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);

      ParticleType::NextID(maxnextid);
      ParallelDescriptor::ReduceIntMax(maxnextid, IOProcNumber);
    }


    if (ParallelDescriptor::IOProcessor())
    {
        std::string HdrFileName = pdir;
	
        if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/') {
	  HdrFileName += '/';
	}
	
        HdrFileName += "Header";
        HdrFileNamePrePost = HdrFileName;
	
        HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);
	
        if ( ! HdrFile.good()) {
	  amrex::FileOpenFailed(HdrFileName);
	}
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
        // AMREX_SPACEDIM and N for sanity checking.
        //
        HdrFile << AMREX_SPACEDIM << '\n';
	
	// The number of extra real parameters
        HdrFile << NStructReal + NArrayReal << '\n';

        // Real component names
        if (real_comp_names.size() == 0) {
            for (int i = 0; i < NStructReal + NArrayReal; ++i ) {
                HdrFile << "real_comp" << i << '\n';
            }
        } else {
            BL_ASSERT(real_comp_names.size() == NStructReal + NArrayReal);
            for (int i = 0; i < NStructReal + NArrayReal; ++i ) {
                HdrFile << real_comp_names[i] << '\n';
            }
        }

	// The number of extra int parameters
        HdrFile << NStructInt + NArrayInt << '\n';

        // int component names
        if (int_comp_names.size() == 0) {
            for (int i = 0; i < NStructInt + NArrayInt; ++i ) {
                HdrFile << "int_comp" << i << '\n';
            }
        } else {
            BL_ASSERT(int_comp_names.size() == NStructInt + NArrayInt);
            for (int i = 0; i < NStructInt + NArrayInt; ++i ) {
                HdrFile << int_comp_names[i] << '\n';
            }
        }

        HdrFile << is_checkpoint << '\n';

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
        HdrFile << finestLevel() << '\n';
        //
        // Then the number of grids at each level.
        //
        for (int lev = 0; lev <= finestLevel(); lev++)
	  {
            HdrFile << ParticleBoxArray(lev).size() << '\n';
	  }
    }
    //
    // We want to write the data out in parallel.
    //
    // We'll allow up to nOutFiles active writers at a time.
    //
    int nOutFiles(256);

    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) {
      nOutFiles = NProcs;
    }
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
    nOutFilesPrePost = nOutFiles;

    for (int lev = 0; lev <= finestLevel(); lev++)
      {
        bool gotsome;
	if(usePrePost) {
          gotsome = (nParticlesAtLevelPrePost[lev] > 0);
	} else {
          gotsome = (NumberOfParticlesAtLevel(lev) > 0);
	}
        //
        // We store the particles at each level in their own subdirectory.
        //
        std::string LevelDir = pdir;
	
        if (gotsome)
	  {
            if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') {
	      LevelDir += '/';
	    }
	    
            LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);
	    
	    if ( ! levelDirectoriesCreated) {
              if (ParallelDescriptor::IOProcessor()) 
	        if ( ! amrex::UtilCreateDirectory(LevelDir, 0755)) 
		  amrex::CreateDirectoryFailed(LevelDir);
              //
              // Force other processors to wait until directory is built.
              //
              ParallelDescriptor::Barrier();
	    }
	  }
	
        // Write out the header for each particle
        if (gotsome and ParallelDescriptor::IOProcessor()) {
            std::string HeaderFileName = LevelDir;
            HeaderFileName += "/Particle_H";
            std::ofstream ParticleHeader(HeaderFileName);
            
            ParticleBoxArray(lev).writeOn(ParticleHeader);
            ParticleHeader << '\n';

            ParticleHeader.flush();
            ParticleHeader.close();
        }
        
	MFInfo info;
	info.SetAlloc(false);
	MultiFab state(ParticleBoxArray(lev),
		       ParticleDistributionMap(lev),
		       1,0,info);
        //
        // We eventually want to write out the file name and the offset
        // into that file into which each grid of particles is written.
        //
        Vector<int>  which(state.size(),0);
        Vector<int > count(state.size(),0);
        Vector<long> where(state.size(),0);
	
	std::string filePrefix(LevelDir);
	filePrefix += '/';
	filePrefix += ParticleType::DataPrefix();
	if(usePrePost) {
	  filePrefixPrePost[lev] = filePrefix;
	}
	bool groupSets(false), setBuf(true);

        if (gotsome)
	{
	    for(NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi)
	    {
	      std::ofstream& myStream = (std::ofstream&) nfi.Stream();
	      //
	      // Write out all the valid particles we own at the specified level.
	      // Do it grid block by grid block remembering the seek offset
	      // for the start of writing of each block of data.
	      //
	      WriteParticles(lev, myStream, nfi.FileNumber(), which, count, where, is_checkpoint);
	    }

	    if(usePrePost) {
              whichPrePost[lev] = which;
              countPrePost[lev] = count;
              wherePrePost[lev] = where;
	    } else {
              ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProcNumber);
              ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProcNumber);
              ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProcNumber);
	    }
        }

        if (ParallelDescriptor::IOProcessor())
        {
	  if(usePrePost) {
	    // ---- write to the header and unlink in CheckpointPost
	  } else {
            for (int j = 0; j < state.size(); j++)
            {
                //
                // We now write the which file, the particle count, and the
                // file offset into which the data for each grid was written,
                // to the header file.
                //
                HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
            }

            if (gotsome && doUnlink)
            {
		BL_PROFILE_VAR("PC<NNNN>::Checkpoint:unlink", unlink);
                //
                // Unlink any zero-length data files.
                //
                Vector<long> cnt(nOutFiles,0);

                for (int i = 0, N=count.size(); i < N; i++) {
                    cnt[which[i]] += count[i];
		}

                for (int i = 0, N=cnt.size(); i < N; i++)
                {
                    if (cnt[i] == 0)
                    {
                        std::string FullFileName = NFilesIter::FileName(i, filePrefix);
                        amrex::UnlinkFile(FullFileName.c_str());
                    }
                }
            }

	  }

        }
    }            // ---- end for(lev...)

    if (ParallelDescriptor::IOProcessor())
    {
      HdrFile.flush();
      HdrFile.close();
      if ( ! HdrFile.good())
      {
        amrex::Abort("ParticleContainer::Checkpoint(): problem writing HdrFile");
      }
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime, IOProcNumber);

        amrex::Print() << "ParticleContainer::Checkpoint() time: " << stoptime << '\n';
    }
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CheckpointPre ()
{
    if( ! usePrePost) {
        return;
    }

    BL_PROFILE("ParticleContainer::CheckpointPre()");
    
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    long nparticles = 0;
    int  maxnextid  = ParticleType::NextID();
    
    for (int lev = 0; lev < m_particles.size();  lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
            for (const auto& p : aos) {
                if (p.m_idata.id > 0) {
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
		}
            }
        }
    }
    ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);

    ParticleType::NextID(maxnextid);
    ParallelDescriptor::ReduceIntMax(maxnextid, IOProcNumber);
    
    nparticlesPrePost = nparticles;
    maxnextidPrePost  = maxnextid;
    
    nParticlesAtLevelPrePost.clear();
    nParticlesAtLevelPrePost.resize(finestLevel() + 1, 0);
    for(int lev(0); lev <= finestLevel(); ++lev) {
        nParticlesAtLevelPrePost[lev] = NumberOfParticlesAtLevel(lev);
    }
    
    whichPrePost.clear();
    whichPrePost.resize(finestLevel() + 1);
    countPrePost.clear();
    countPrePost.resize(finestLevel() + 1);
    wherePrePost.clear();
    wherePrePost.resize(finestLevel() + 1);
    
    filePrefixPrePost.clear();
    filePrefixPrePost.resize(finestLevel() + 1);
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::CheckpointPost ()
{
    if( ! usePrePost) {
      return;
    }

    BL_PROFILE("ParticleContainer::CheckpointPost()");

    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    std::ofstream HdrFile;
    HdrFile.open(HdrFileNamePrePost.c_str(), std::ios::out | std::ios::app);

    for(int lev(0); lev <= finestLevel(); ++lev) {
      ParallelDescriptor::ReduceIntSum (whichPrePost[lev].dataPtr(), whichPrePost[lev].size(), IOProcNumber);
      ParallelDescriptor::ReduceIntSum (countPrePost[lev].dataPtr(), countPrePost[lev].size(), IOProcNumber);
      ParallelDescriptor::ReduceLongSum(wherePrePost[lev].dataPtr(), wherePrePost[lev].size(), IOProcNumber);


      if(ParallelDescriptor::IOProcessor()) {
        for(int j(0); j < whichPrePost[lev].size(); ++j) {
          HdrFile << whichPrePost[lev][j] << ' ' << countPrePost[lev][j] << ' ' << wherePrePost[lev][j] << '\n';
        }

        const bool gotsome = (nParticlesAtLevelPrePost[lev] > 0);
        if(gotsome && doUnlink) {
            BL_PROFILE_VAR("PC<NNNN>::Checkpoint:unlink", unlink_post);
          // Unlink any zero-length data files.
          Vector<long> cnt(nOutFilesPrePost,0);

          for(int i(0), N = countPrePost[lev].size(); i < N; ++i) {
            cnt[whichPrePost[lev][i]] += countPrePost[lev][i];
	  }

          for(int i(0), N = cnt.size(); i < N; ++i) {
            if(cnt[i] == 0) {
              std::string FullFileName = NFilesIter::FileName(i, filePrefixPrePost[lev]);
              amrex::UnlinkFile(FullFileName.c_str());
            }
          }
        }
      }
    }

    if(ParallelDescriptor::IOProcessor()) {
      HdrFile.flush();
      HdrFile.close();
      if( ! HdrFile.good()) {
        amrex::Abort("ParticleContainer::CheckpointPost(): problem writing HdrFile");
      }
    }
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
WritePlotFile ( const std::string&        dir,
                const std::string&        name,
                const Vector<std::string>& real_comp_names,
                const Vector<std::string>& int_comp_names) const
{
    BL_PROFILE("ParticleContainer::WritePlotFile()");
    BL_ASSERT(OK());
    bool is_checkpoint = false;

    // For yt we need exactly the chk particle format so would need to set is_checkpoint = true
    // Anyway, it's not too bad to have particle ids on disk,
    // think of merger trees or backtracing of particles for nested ics
    // is_checkpoint = true; 
    Checkpoint(dir,name,is_checkpoint,real_comp_names,int_comp_names);
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFilePre ()
{
  CheckpointPre();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WritePlotFilePost ()
{
  CheckpointPost();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteParticles (int            lev,
                                                                                   std::ofstream& ofs,
                                                                                   int            fnum,
                                                                                   Vector<int>&    which,
                                                                                   Vector<int>&    count,
                                                                                   Vector<long>&   where,
                                                                                   bool           is_checkpoint) const
{
    BL_PROFILE("ParticleContainer::WriteParticles()");

    // For a each grid, the tiles it contains
    std::map<int, Vector<int> > tile_map;

    for (const auto& kv : m_particles[lev])
    {
        const int grid = kv.first.first;
        const int tile = kv.first.second;
        tile_map[grid].push_back(tile);

        // Only write out valid particles.
        int cnt = 0;	
        for (const auto& p : kv.second.GetArrayOfStructs()) {
            if (p.m_idata.id > 0) {
                cnt++;
	    }
	}

        count[grid] += cnt;
    }
	
    MFInfo info;
    info.SetAlloc(false);
    MultiFab state(ParticleBoxArray(lev),
		   ParticleDistributionMap(lev),
		   1,0,info);

    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
      const int grid = mfi.index();
      
      which[grid] = fnum;
      where[grid] = VisMF::FileOffset(ofs);
      
      if (count[grid] == 0) {
        continue;
      }
      
      if (is_checkpoint) {
	// First write out the integer data in binary.
	const int iChunkSize = 2 + NStructInt + NArrayInt;
	Vector<int> istuff(count[grid]*iChunkSize);
	int* iptr = istuff.dataPtr();

	for (unsigned i = 0; i < tile_map[grid].size(); i++) {
            const auto& pbox = m_particles[lev].at(std::make_pair(grid, tile_map[grid][i]));
            int pindex = 0;
            for (const auto& p : pbox.GetArrayOfStructs()) {
                if (p.m_idata.id > 0) {
                    for (int j = 0; j < 2 + NStructInt; j++) {
                        iptr[j] = p.m_idata.arr[j];
		    }
                    iptr += 2 + NStructInt;
                    const auto& soa  = pbox.GetStructOfArrays();
                    for (int j = 0; j < NArrayInt; j++) {
                        iptr[j] = soa.GetIntData(j)[pindex];
                    }
                    iptr += NArrayInt;
                }
                ++pindex;
            }
	}
        
        writeIntData(istuff.dataPtr(), istuff.size(), ofs);
        ofs.flush();  // Some systems require this flush() (probably due to a bug)
      }
      
      // Write the Real data in binary.
      const int rChunkSize = AMREX_SPACEDIM + NStructReal + NArrayReal;
      Vector<typename ParticleType::RealType> rstuff(count[grid]*rChunkSize);
      typename ParticleType::RealType* rptr = rstuff.dataPtr();
      
      for (unsigned i = 0; i < tile_map[grid].size(); i++) {
          const auto& pbox = m_particles[lev].at(std::make_pair(grid, tile_map[grid][i]));
          int pindex = 0;
          for (const auto& p : pbox.GetArrayOfStructs()) {
              if (p.m_idata.id > 0) {
                  for (int j = 0; j < AMREX_SPACEDIM + NStructReal; j++) {
                      rptr[j] = p.m_rdata.arr[j];
                  }
                  rptr += AMREX_SPACEDIM + NStructReal;
                  const auto& soa  = pbox.GetStructOfArrays();
                  for (int j = 0; j < NArrayReal; j++) {
                      rptr[j] = (typename ParticleType::RealType) soa.GetRealData(j)[pindex];
                  }
                  rptr += NArrayReal;
              }
              ++pindex;
          }
      }
      WriteParticleRealData(rstuff.dataPtr(), rstuff.size(), ofs, ParticleRealDescriptor);
      ofs.flush();  // Some systems require this flush() (probably due to a bug)              
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::Restart (const std::string& dir, const std::string& file, bool is_checkpoint)
{
  BL_PROFILE("ParticleContainer::Restart()");
  BL_ASSERT(!dir.empty());
  BL_ASSERT(!file.empty());
  
  const Real strttime = amrex::second();
  
  int DATA_Digits_Read(5);
  ParmParse pp("particles");
  pp.query("datadigits_read",DATA_Digits_Read);

  std::string fullname = dir;
  if (!fullname.empty() && fullname[fullname.size()-1] != '/')
    fullname += '/';
  fullname += file;
  std::string HdrFileName = fullname;
  if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
    HdrFileName += '/';
  HdrFileName += "Header";
  
  Vector<char> fileCharPtr;
  ParallelDescriptor::ReadAndBcastFile(HdrFileName, fileCharPtr);
  std::string fileCharPtrString(fileCharPtr.dataPtr());
  std::istringstream HdrFile(fileCharPtrString, std::istringstream::in);
  
  std::string version;
  HdrFile >> version;
  BL_ASSERT(!version.empty());
  
  // What do our version strings mean?
  // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
  // "Version_One_Dot_One" -- can write out either as either single or double precision.
  // Appended to the latter version string are either "_single" or "_double" to
  // indicate how the particles were written.
  // "Version_Two_Dot_Zero" -- this is the AMReX particle file format
  std::string how;
  if (version.find("Version_One_Dot_Zero") != std::string::npos) {
    how = "double";
  }
  else if (version.find("Version_One_Dot_One")  != std::string::npos or
           version.find("Version_Two_Dot_Zero") != std::string::npos) {
    if (version.find("_single") != std::string::npos) {
      how = "single";
    }
    else if (version.find("_double") != std::string::npos) {
      how = "double";
    }
    else {
      std::string msg("ParticleContainer::Restart(): bad version string: ");
      msg += version;
      amrex::Error(version.c_str());
    }
  }
  else {
    std::string msg("ParticleContainer::Restart(): unknown version string: ");
    msg += version;
    amrex::Abort(msg.c_str());
  }
  
  int dm;
  HdrFile >> dm;
  if (dm != AMREX_SPACEDIM)
    amrex::Abort("ParticleContainer::Restart(): dm != AMREX_SPACEDIM");
  
  int nr;
  HdrFile >> nr;
  if (nr != NStructReal + NArrayReal)
    amrex::Abort("ParticleContainer::Restart(): nr != NStructReal + NArrayReal");

  std::string comp_name;
  for (int i = 0; i < nr; ++i)
      HdrFile >> comp_name;

  int ni;
  HdrFile >> ni;
  if (ni != NStructInt + NArrayInt)
    amrex::Abort("ParticleContainer::Restart(): ni != NStructInt");
  
  for (int i = 0; i < ni; ++i)
      HdrFile >> comp_name;

  bool checkpoint;
  HdrFile >> checkpoint;

  long nparticles;
  HdrFile >> nparticles;
  BL_ASSERT(nparticles >= 0);
  
  int maxnextid;
  HdrFile >> maxnextid;
  BL_ASSERT(maxnextid > 0);
  ParticleType::NextID(maxnextid);
  
  int finest_level_in_file;
  HdrFile >> finest_level_in_file;
  BL_ASSERT(finest_level_in_file >= 0);

  // Determine whether this is a dual-grid restart or not.
  Vector<BoxArray> particle_box_arrays(finest_level_in_file + 1);
  bool dual_grid = false;
  for (int lev = 0; lev <= finest_level_in_file; lev++) {
      std::string phdr_name = fullname;
      phdr_name = amrex::Concatenate(phdr_name + "/Level_", lev, 1);
      phdr_name += "/Particle_H";

      if (not amrex::FileExists(phdr_name)) {
          dual_grid = false;
          break;
      }

      Vector<char> phdr_chars;
      ParallelDescriptor::ReadAndBcastFile(phdr_name, phdr_chars);
      std::string phdr_string(phdr_chars.dataPtr());
      std::istringstream phdr_file(phdr_string, std::istringstream::in);

      particle_box_arrays[lev].readFrom(phdr_file);
      if (not particle_box_arrays[lev].CellEqual(ParticleBoxArray(lev))) dual_grid = true;
  }
  
  if (dual_grid) {
      for (int lev = 0; lev <= finestLevel(); lev++) {
          SetParticleBoxArray(lev, particle_box_arrays[lev]);
          DistributionMapping pdm(particle_box_arrays[lev]);
          SetParticleDistributionMap(lev, pdm);
      }
  }
  
  Vector<int> ngrids(finest_level_in_file+1);
  for (int lev = 0; lev <= finest_level_in_file; lev++) {
      HdrFile >> ngrids[lev];
      BL_ASSERT(ngrids[lev] > 0);
      if (lev <= finestLevel()) {
          BL_ASSERT(ngrids[lev] == int(ParticleBoxArray(lev).size()));
      }
  }

  resizeData();
  
  if (finest_level_in_file > finestLevel()) {
      m_particles.resize(finest_level_in_file+1);
  }

  for (int lev = 0; lev <= finest_level_in_file; lev++) {
      Vector<int>  which(ngrids[lev]);
      Vector<int>  count(ngrids[lev]);
      Vector<long> where(ngrids[lev]);
      for (int i = 0; i < ngrids[lev]; i++) {
          HdrFile >> which[i] >> count[i] >> where[i];
      }
    
      Vector<int> grids_to_read;
      if (lev <= finestLevel()) {
          for (MFIter mfi(*m_dummy_mf[lev]); mfi.isValid(); ++mfi) {
              grids_to_read.push_back(mfi.index());
          }
      } else {
          
          // we lost a level on restart. we still need to read in particles
          // on finer levels, and put them in the right place via Redistribute()
          
          const int rank = ParallelDescriptor::MyProc();
          const int NReaders = ParticleType::MaxReaders();
          if (rank >= NReaders) return;
          
          const int Navg = ngrids[lev] / NReaders;
          const int Nleft = ngrids[lev] - Navg * NReaders;
          
          int lo, hi;
          if (rank < Nleft) {
              lo = rank*(Navg + 1);
              hi = lo + Navg + 1;
          } 
          else {
              lo = rank * Navg + Nleft;
              hi = lo + Navg;
          }
          
          for (int i = lo; i < hi; ++i) {
              grids_to_read.push_back(i);
          }
      }

      for(int igrid = 0; igrid < static_cast<int>(grids_to_read.size()); ++igrid) {
          const int grid = grids_to_read[igrid];
          
          if (count[grid] <= 0) continue;
          
          // The file names in the header file are relative.
          std::string name = fullname;
          
          if (!name.empty() && name[name.size()-1] != '/')
              name += '/';
          
          name += "Level_";
          name += amrex::Concatenate("", lev, 1);
          name += '/';
          name += ParticleType::DataPrefix();
          name += amrex::Concatenate("", which[grid], DATA_Digits_Read);
          
          std::ifstream ParticleFile;
          
          ParticleFile.open(name.c_str(), std::ios::in);
          
          if (!ParticleFile.good())
              amrex::FileOpenFailed(name);
          
          ParticleFile.seekg(where[grid], std::ios::beg);
          
          if (how == "single") {
              ReadParticles<float>(count[grid], grid, lev, is_checkpoint, ParticleFile);
          }
          else if (how == "double") {
              ReadParticles<double>(count[grid], grid, lev, is_checkpoint, ParticleFile);
          }
          else {
              std::string msg("ParticleContainer::Restart(): bad parameter: ");
              msg += how;
              amrex::Error(msg.c_str());
          }
          
          ParticleFile.close();
          
          if (!ParticleFile.good())
              amrex::Abort("ParticleContainer::Restart(): problem reading particles");
      }
  }
  
  Redistribute();
  
  BL_ASSERT(OK());
  
  if (m_verbose > 1) {
      Real stoptime = amrex::second() - strttime;	
      ParallelDescriptor::ReduceRealMax(stoptime, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "ParticleContainer::Restart() time: " << stoptime << '\n';
  }
}

// Read a batch of particles from the checkpoint file
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
template <class RTYPE>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::ReadParticles (int            cnt,
                                                                                  int            grd,
                                                                                  int            lev,
                                                                                  bool           is_checkpoint,
                                                                                  std::ifstream& ifs) 
{
    BL_PROFILE("ParticleContainer::ReadParticles()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < int(m_particles.size()));

    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    const int iChunkSize = 2 + NStructInt + NArrayInt;
    Vector<int> istuff(cnt*iChunkSize);
    if (is_checkpoint)
        readIntData(istuff.dataPtr(), istuff.size(), ifs, FPC::NativeIntDescriptor());

    // Then the real data in binary.
    const int rChunkSize = AMREX_SPACEDIM + NStructReal + NArrayReal;
    Vector<RTYPE> rstuff(cnt*rChunkSize);
    ReadParticleRealData(rstuff.dataPtr(), rstuff.size(), ifs, ParticleRealDescriptor);

    // Now reassemble the particles.
    int*   iptr = istuff.dataPtr();
    RTYPE* rptr = rstuff.dataPtr();

    // If we are restarting from a plotfile instead of a checkpoint file, then we do not
    //    read in the particle id's, so we need to reset the id counter to zero and renumber them
    if (!is_checkpoint) {
      int maxnextid = 1;
      ParticleType::NextID(maxnextid);
    }

    ParticleType p;
    ParticleLocData pld;
    for (int i = 0; i < cnt; i++) {
      if (is_checkpoint) {
	p.m_idata.id   = iptr[0];
	p.m_idata.cpu  = iptr[1];
      }
      else {
	p.m_idata.id   = ParticleType::NextID();
	p.m_idata.cpu  = ParallelDescriptor::MyProc();
      }

      BL_ASSERT(p.m_idata.id > 0);

      for (int j = 0; j < NStructInt; j++)
          p.m_idata.arr[2+j] = iptr[2+j];

      iptr += 2 + NStructInt;

      AMREX_D_TERM(p.m_rdata.pos[0] = rptr[0];,
	     p.m_rdata.pos[1] = rptr[1];,
	     p.m_rdata.pos[2] = rptr[2];);

      for (int j = 0; j < NStructReal; j++)
	p.m_rdata.arr[AMREX_SPACEDIM+j] = rptr[AMREX_SPACEDIM+j];

      rptr += AMREX_SPACEDIM + NStructReal;
      
      locateParticle(p, pld, 0, finestLevel(), 0);

      auto& ptile = m_particles[lev][std::make_pair(grd, pld.m_tile)];

      ptile.push_back(p);

      for (int j = 0; j < NArrayReal; j++) {
          ptile.push_back_real(j, rptr[j]);
      }

      rptr += NArrayReal;

      for (int j = 0; j < NArrayInt; j++) {
          ptile.push_back_int(j, iptr[j]);
      }

      iptr += NArrayInt;

    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer::WriteAsciiFile()");
    BL_ASSERT(!filename.empty());

    const Real strttime = amrex::second();
    //
    // Count # of valid particles.
    //
    long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
        auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
            for (const auto& p : aos) {
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
        // Have I/O processor open file and write out particle metadata.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles  << '\n';
        File << NStructReal << '\n';
        File << NStructInt  << '\n';
        File << NArrayReal  << '\n';
        File << NArrayInt   << '\n';
            
        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
    {
        if (MyProc == proc)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            std::ofstream File;

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

	    for (int lev = 0; lev < m_particles.size();  lev++) {
	      auto& pmap = m_particles[lev];
	      for (const auto& kv : pmap) {
                const auto& aos = kv.second.GetArrayOfStructs();
                const auto& soa = kv.second.GetStructOfArrays();

		int index = 0;
		for (auto it = aos.cbegin(); it != aos.cend(); ++it) {
		    if (it->m_idata.id > 0) {

                        // write out the particle struct first... 
                        AMREX_D_TERM(File << it->m_rdata.pos[0] << ' ',
                               << it->m_rdata.pos[1] << ' ',
                               << it->m_rdata.pos[2] << ' ');

                        for (int i = AMREX_SPACEDIM; i < AMREX_SPACEDIM + NStructReal; i++)
                            File << it->m_rdata.arr[i] << ' ';

                        File << it->m_idata.id  << ' ';
                        File << it->m_idata.cpu << ' ';
                        
                        for (int i = 2; i < 2 + NStructInt; i++)
                            File << it->m_idata.arr[i] << ' ';
		      
                        // then the particle attributes.
                        for (int i = 0; i < NArrayReal; i++)
                            File << soa.GetRealData(i)[index] << ' ';

                        for (int i = 0; i < NArrayInt; i++)
                            File << soa.GetIntData(i)[index] << ' ';

                        File << '\n';                            
                        
                        index++;		      
                    }
                }
              }
            }
	    
            File.flush();
	    
            File.close();
            
            if (!File.good())
                amrex::Abort("ParticleContainer::WriteAsciiFile(): problem writing file");
	    
        }
	
        ParallelDescriptor::Barrier();
    }
    
    if (m_verbose > 1)
        {
            Real stoptime = amrex::second() - strttime;
            
            ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
            
            amrex::Print() << "ParticleContainer::WriteAsciiFile() time: " << stoptime << '\n';
        }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::WriteCoarsenedAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer::WriteCoarsenedAsciiFile()");
    BL_ASSERT(!filename.empty());

    const Real strttime = amrex::second();
 
    //
    // Count # of valid particles.
    //
    long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
        auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
            for (const auto& p : aos) {
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
            amrex::Abort("ParticleContainer::WriteCoarsenedAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
    {
        if (MyProc == proc)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            std::ofstream File;

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

	    for (int lev = 0; lev < m_particles.size();  lev++) {
	      auto& pmap = m_particles[lev];
	      for (auto& kv : pmap) {
                  auto& aos = kv.second.GetArrayOfStructs();
                  auto& soa = kv.second.GetStructOfArrays();
                  
                  int index = 0;
                  ParticleLocData pld;
                  for (auto it = aos.begin(); it != aos.end(); ++it) {
                      locateParticle(*it, pld, 0, finestLevel(), 0);
                      // Only keep particles in even cells
                      if (it->id() > 0 &&
                          (pld.m_cell[0])%2 == 0 && (pld.m_cell[1])%2 == 0 && (pld.m_cell[2])%2 == 0)
                      {
                          
                          // Only keep particles in even cells               
                          if (it->m_idata.id > 0) {
                              
                              File << it->m_idata.id  << ' ';
                              File << it->m_idata.cpu << ' ';
                              
                              AMREX_D_TERM(File << it->m_rdata.pos[0] << ' ',
                                     << it->m_rdata.pos[1] << ' ',
                                     << it->m_rdata.pos[2] << ' ');
                              
                              for (int i = 0; i < NArrayReal; i++) {
                                  File << soa.GetRealData(i)[index] << ' ';
                              }
                              index++;
                              
                              for (int i = AMREX_SPACEDIM; i < AMREX_SPACEDIM + NStructReal; i++) {
                                  char ws = (i == AMREX_SPACEDIM + NStructReal - 1) ? '\n' : ' ';
                                  if (i == AMREX_SPACEDIM) {
                                      // Multiply mass by 8 since we are only taking 1/8 of the 
                                      // total particles and want to keep the mass in the domain the same.
                                      File << 8.0* it->m_rdata.arr[i] << ws;
                                  }
                                  else {
                                      File << it->m_rdata.arr[i] << ws;
                                  }
                              }
                          }
                      }
                  }
              }
            }

            File.flush();

            File.close();

            if (!File.good())
                amrex::Abort("ParticleContainer::WriteCoarsenedAsciiFile(): problem writing file");

        }

        ParallelDescriptor::Barrier();
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "ParticleContainer::WriteCoarsenedAsciiFile() time: " << stoptime << '\n';
    }
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
    // own grid.  The mf->sumBoundary call then adds the value from one grid's ghost cell
    // to another grid's valid region.
    if (mf_pointer->nGrow() < 1) 
       amrex::Error("Must have at least one ghost cell when in AssignDensitySingleLevel");

#ifdef _OPENMP
    const int       ng          = mf_pointer->nGrow();
#endif
    const Real      strttime    = amrex::second();
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
      amrex::Error("AssignDensity: problem must be periodic in no or all directions");
    }
    
    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }

    using ParConstIter = ParConstIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        FArrayBox local_rho;
        for (ParConstIter pti(*this, lev); pti.isValid(); ++pti) {
            const auto& particles = pti.GetArrayOfStructs();
            int nstride = particles.dataShape().first;
            const long np = pti.numParticles();
            FArrayBox& fab = (*mf_pointer)[pti];
            Real* data_ptr;
            const int *lo, *hi;
#ifdef _OPENMP
            Box tile_box = pti.tilebox();
            tile_box.grow(ng);
            local_rho.resize(tile_box,ncomp);
            local_rho = 0.0;
            data_ptr = local_rho.dataPtr();
            lo = tile_box.loVect();
            hi = tile_box.hiVect();
#else
            const Box& box = fab.box();
            data_ptr = fab.dataPtr();
            lo = box.loVect();
            hi = box.hiVect();
#endif

            if (dx == dx_particle) {
                amrex_deposit_cic(particles.data(), nstride, np, ncomp, 
                                  data_ptr, lo, hi, plo, dx);
            } else {
                amrex_deposit_particle_dx_cic(particles.data(), nstride, np, ncomp,
                                              data_ptr, lo, hi, plo, dx, dx_particle);
            }
                

#ifdef _OPENMP
            amrex_atomic_accumulate_fab(BL_TO_FORTRAN_3D(local_rho), 
                                        BL_TO_FORTRAN_3D(fab), ncomp);
#endif

        }
    }

    mf_pointer->SumBoundary(gm.periodicity());
    
    // If ncomp > 1, first divide the momenta (component n) 
    // by the mass (component 0) in order to get velocities.
    // Be careful not to divide by zero.
    for (int n = 1; n < ncomp; n++){
      for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
	(*mf_pointer)[mfi].protected_divide((*mf_pointer)[mfi],0,n,1);
      }
    }

    // Only multiply the first component by (1/vol) because this converts mass
    // to density. If there are additional components (like velocity), we don't
    // want to divide those by volume.
    const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    mf_pointer->mult(1.0/vol, 0, 1, mf_pointer->nGrow());

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled. I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled) {
      mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
      delete mf_pointer;
    }
    
    if (m_verbose > 1) {
      Real stoptime = amrex::second() - strttime;
      
      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
      
      amrex::Print() << "ParticleContainer::AssignDensity(single-level) time: " << stoptime << '\n';
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
    
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx          = gm.CellSize();

    using ParIter = ParIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (ParIter pti(*this, lev); pti.isValid(); ++pti) {
        auto& particles = pti.GetArrayOfStructs();
        FArrayBox& fab = mesh_data[pti];
        const Box& box = fab.box();
        const long N = particles.size();
        int nstride = particles.dataShape().first;
        int nComp = fab.nComp();
        amrex_interpolate_cic(particles.data(), nstride, N, 
                              fab.dataPtr(), box.loVect(), box.hiVect(), nComp, plo, dx);
    }
}

//
// This version takes as input the acceleration vector at cell centers, and has the option of
// returning the acceleration at the particle location in the data array, starting at
// component start_comp_for_accel
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::moveKick (MultiFab&       acceleration,
                                                                             int             lev,
                                                                             Real            dt,
                                                                             Real            a_new,
                                                                             Real            a_half, 
                                                                             int             start_comp_for_accel)
{
    BL_PROFILE("ParticleContainer::moveKick()");
    BL_ASSERT(NStructReal >= AMREX_SPACEDIM+1);
    BL_ASSERT(lev >= 0 && lev < int(m_particles.size()));

    const Real strttime  = amrex::second();
    const Real half_dt   = Real(0.5) * dt;
    const Real a_new_inv = 1 / a_new;
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
