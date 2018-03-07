
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::do_tiling = false;

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
IntVect
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::tile_size   { AMREX_D_DECL(1024000,8,8) };

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> :: Initialize ()
{
    levelDirectoriesCreated = false;
    usePrePost = false;
    doUnlink = true;

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
        if ( ! std::is_pod<ParticleType>::value) {
            amrex::Abort("Particle is not POD");
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

    AMREX_D_TERM(iv[0]=floor((p.m_rdata.pos[0]-geom.ProbLo(0))*geom.InvCellSize(0));,
                 iv[1]=floor((p.m_rdata.pos[1]-geom.ProbLo(1))*geom.InvCellSize(1));,
                 iv[2]=floor((p.m_rdata.pos[2]-geom.ProbLo(2))*geom.InvCellSize(2)););

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
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Where (const ParticleType& p,
                                                                          ParticleLocData&    pld,
                                                                          int                 lev_min,
                                                                          int                 lev_max,
                                                                          int                 nGrow) const
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

      const BoxArray& ba = ParticleBoxArray(lev);
      BL_ASSERT(ba.ixType().cellCentered());      
      ba.intersections(Box(iv, iv), isects, true, nGrow);
      
      if (!isects.empty()) {
          const Box& bx = ba.getCellCenteredBox(isects[0].first);
	  pld.m_lev  = lev;
	  pld.m_grid = isects[0].first;
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
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::EnforcePeriodicWhere (ParticleType&    p,
                                                                                         ParticleLocData& pld,
                                                                                         int              lev_min,
                                                                                         int              lev_max) const
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
            const IntVect& iv = Index(p_prime, lev);
            const BoxArray& ba = ParticleBoxArray(lev);
            BL_ASSERT(ba.ixType().cellCentered());
            
            ba.intersections(Box(iv,iv),isects,true,0);
            
            if (!isects.empty()) {
                AMREX_D_TERM(p.m_rdata.pos[0] = p_prime.m_rdata.pos[0];,
                             p.m_rdata.pos[1] = p_prime.m_rdata.pos[1];,
                             p.m_rdata.pos[2] = p_prime.m_rdata.pos[2];);
                
                const Box& bx = ba.getCellCenteredBox(isects[0].first);
                
                pld.m_lev  = lev;
                pld.m_grid = isects[0].first;
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
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::PeriodicShift (ParticleType& p) const
{
    //    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::PeriodicShift()");

    // This routine should only be called when Where() returns false.
    BL_ASSERT(m_gdb != 0);

    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    const Geometry& geom    = Geom(0);
    const Box&      dmn     = geom.Domain();
    const IntVect&  iv      = Index(p, 0);
    bool            shifted = false;  

    for (int i = 0; i < AMREX_SPACEDIM; i++) {
        if (!geom.isPeriodic(i)) continue;

        if (iv[i] > dmn.bigEnd(i)) {
            if (p.m_rdata.pos[i] == geom.ProbHi(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.m_rdata.pos[i] += .125*geom.CellSize(i);
            
            p.m_rdata.pos[i] -= geom.ProbLength(i);
            
            if (p.m_rdata.pos[i] <= geom.ProbLo(i))
                //
                // This can happen due to precision issues.
                //
                p.m_rdata.pos[i] += .125*geom.CellSize(i);
            
            BL_ASSERT(p.m_rdata.pos[i] >= geom.ProbLo(i));
            
            shifted = true;
        }
        else if (iv[i] < dmn.smallEnd(i)) {
            if (p.m_rdata.pos[i] == geom.ProbLo(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.m_rdata.pos[i] -= .125*geom.CellSize(i);
            
            p.m_rdata.pos[i] += geom.ProbLength(i);
            
            if (p.m_rdata.pos[i] >= geom.ProbHi(i))
                //
                // This can happen due to precision issues.
                //
                p.m_rdata.pos[i] -= .125*geom.CellSize(i);
            
            BL_ASSERT(p.m_rdata.pos[i] <= geom.ProbHi(i));
            
            shifted = true;
        }
    }
    //
    // The particle may still be outside the domain in the case
    // where we aren't periodic on the face out which it travelled.
    //
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
    //    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::Reset()");
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
                                                                                   int lev_min, int lev_max, int nGrow) const
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
        success = EnforcePeriodicWhere(p, pld, lev_min, lev_max);
        if (!success && lev_min == 0)
        {
            // The particle has left the domain; invalidate it.
            p.m_idata.id = -p.m_idata.id;
            success = true;
        }
    }
    else
    {
        success = Where(p, pld, lev_min, lev_max);
    }

    if (!success)
    {
        success = (nGrow > 0) && Where(p, pld, lev_min, lev_min, nGrow);
        pld.m_grown_gridbox = pld.m_gridbox; // reset grown box for subsequent calls.
    }

    if (!success)
    {
        amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::locateParticle(): invalid particle.");
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

    amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> byte spread across MPI nodes: ["
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
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::MoveRandom(lev)");
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
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::IncrementWithTotal(lev)");
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
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::sumParticleMass(lev)");
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
std::pair<long,long>
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::StartIndexInGlobalArray () const
{
    // This returns the number of particles on this processor
    long lcnt = TotalNumberOfParticles(true,true);

#if BL_USE_MPI
    Vector<long> cnts(ParallelDescriptor::NProcs());
        
    // This accumulates the "lcnt" values into "cnts"
    MPI_Gather(&lcnt,1,              
               ParallelDescriptor::Mpi_typemap<long>::type(),
               cnts.dataPtr(),
               1,
               ParallelDescriptor::Mpi_typemap<long>::type(),
               ParallelDescriptor::IOProcessorNumber(),
               ParallelDescriptor::Communicator());
    
    ParallelDescriptor::Bcast(cnts.dataPtr(), cnts.size(), ParallelDescriptor::IOProcessorNumber());

    long start = 0L;
    for (int iproc = 0; iproc < ParallelDescriptor::MyProc(); ++iproc) {
        start += cnts[iproc];
    }
    long tot = std::accumulate(cnts.begin(), cnts.end(), 0L);
#else
    long start = 0L;
    long tot = lcnt;
#endif
    return std::make_pair(start,tot);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetParticleIDs (Vector<int>& part_ids)
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::GetParticleIDs()");

  long start, npart;
  std::tie(start,npart) = StartIndexInGlobalArray();

  // Locations
  part_ids.resize(npart, 0);
  
  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
      const auto& pmap = m_particles[lev];
      for (auto& kv : pmap) {
          const auto& pbx = kv.second.GetArrayOfStructs();
          for (const auto& p : pbx) {
              if (p.m_idata.id > 0) {
                  part_ids[start++] = p.m_idata.id;
              }
          }
      }
  }
  
  ParallelDescriptor::ReduceIntSum(part_ids.dataPtr(),part_ids.size()); 
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetParticleCPU (Vector<int>& part_cpu)
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::GetParticleCPU()");

  long start, npart;
  std::tie(start,npart) = StartIndexInGlobalArray();

  // Locations
  part_cpu.resize(npart,0);
  
  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
      const auto& pmap = m_particles[lev];
      for (auto& kv : pmap) {
          const auto& pbx = kv.second.GetArrayOfStructs();
          for (const auto& p : pbx) {
              if (p.m_idata.id > 0) {
                  part_cpu[start++] = p.m_idata.cpu;
              }
          }
      }
  }
  
  ParallelDescriptor::ReduceIntSum(part_cpu.dataPtr(),part_cpu.size()); 
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetParticleLocations (Vector<Real>& part_data)
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::GetParticleLocations()");
 
  long start, npart;
  std::tie(start,npart) = StartIndexInGlobalArray();
  
  // Each particle takes up AMREX_SPACEDIM Reals
  start *= (AMREX_SPACEDIM);
  
  // Locations
  part_data.resize(AMREX_SPACEDIM*npart,0);

  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
      const auto& pmap = m_particles[lev];
      for (auto& kv : pmap) {
          const auto& pbx = kv.second.GetArrayOfStructs();
          for (const auto& p : pbx) {
              if (p.m_idata.id > 0) {
                  // Load positions
                  for (int d=0; d < AMREX_SPACEDIM; d++)
                      part_data[start++] = p.m_rdata.pos[d];
              }
          }
      }
  }
  
  ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetParticleData (Vector<Real>& part_data, int start_comp, int num_comp)
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::GetParticleData()");

  long start, npart;
  std::tie(start,npart) = StartIndexInGlobalArray();
  
  // Each particle takes up num_comp Reals
  start*= num_comp;

  part_data.resize(num_comp*npart,0);

  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
      const auto& pmap = m_particles[lev];
      for (auto& kv : pmap) {
          const auto& pbx = kv.second.GetArrayOfStructs();
          for (const auto& p : pbx) {
              if (p.m_idata.id > 0) {
                  // Load particle data, whatever it is.
                  for (int d = 0; d < num_comp; d++)
                      part_data[start++] = p.m_rdata.arr[AMREX_SPACEDIM + start_comp + d];
              }
          }
      }
  }
  
  ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::GetArrayData (Vector<Real>& part_data, int start_comp, int num_comp)
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::GetArrayData()");

  long start, npart;
  std::tie(start,npart) = StartIndexInGlobalArray();  
  part_data.resize(num_comp*npart,0);

  // Each particle takes up num_comp Reals
  start*= num_comp;
  
  for (int comp = 0; comp < num_comp; ++comp) {
      for (unsigned lev = 0; lev < m_particles.size(); lev++) {
          const auto& pmap = m_particles[lev];
          for (auto& kv : pmap) {
              const auto& pbx = kv.second.GetArrayOfStructs();
              const auto& soa = kv.second.GetStructOfArrays();
              const Vector<Real>& arr = soa[start_comp + comp];
              for (unsigned i = 0; i < arr.size(); ++i) {
                  if (pbx[i].m_idata.id > 0) {
                      part_data[start++] = arr[i];
                  }
              }
          }
      }
  }
  
  ParallelDescriptor::ReduceRealSum(part_data.dataPtr(),part_data.size()); 
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetAllowParticlesNearBoundary (bool value)
{
  allow_particles_near_boundary = value; 
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::SetParticleLocations (Vector<Real>& part_data)
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::SetParticleLocations()");

  long start, npart;
  std::tie(start,npart) = StartIndexInGlobalArray();
  
  // Each particle takes up AMREX_SPACEDIM Reals
  start*= AMREX_SPACEDIM;
  
  // Mass + locations
  if (part_data.size() != npart*AMREX_SPACEDIM)
    amrex::Abort("Sending in wrong size part_data to SetParticleLocations");

  for (unsigned lev = 0; lev < m_particles.size(); lev++) {
      auto& pmap = m_particles[lev];
      for (auto& kv : pmap) {
          auto& pbx = kv.second.GetArrayOfStructs();
          for (auto& p : pbx) {
              if (p.m_idata.id > 0) {
                  // Load positions
                  for (int d=0; d < AMREX_SPACEDIM; d++)
                      p.m_rdata.pos[d] = part_data[start++];
              }
          }
      }
  }
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::RemoveParticlesAtLevel (int level)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::RemoveParticlesAtLevel()");
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
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::RemoveParticlesNotAtFinestLevel()");
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
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::CreateVirtualParticles (int   level,
                                                       AoS& virts) const
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::CreateVirtualParticles()");
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
  BoxArray buffer = amrex::complementIn(Geom(level).Domain(), ParticleBoxArray(level));
  
  buffer.grow(aggregation_buffer);
  
  const auto& pmap = m_particles[level];
  for (const auto& kv : pmap) {
    const auto& pbox = kv.second.GetArrayOfStructs();
    
    //
    // Map for use in Cell aggregation.
    //
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
	    if (aggregation_type == "None")
	      {
		//
		// No aggregation.  Simply clone the particle.
		// Set its id to indicate that it's a virt.
                  virts.push_back(*it);
                  virts.back().m_idata.id = VirtualParticleID;
	      }
	    else if (aggregation_type == "Cell")
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
	for (const auto& agg_particle : agg_map)
	  {
	    virts.push_back(agg_particle.second);
	  }
      }
  }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::CreateGhostParticles (int   level,
                                                     int   nGrow,
                                                     AoS& ghosts) const
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::CreateGhostParticles()");
  BL_ASSERT(ghosts.empty());
  BL_ASSERT(level < finestLevel());
  
  if (level >= int(m_particles.size()))
    //
    // This level could exist and simply have no particles.
    //
    return;
  
  const BoxArray& fine = ParticleBoxArray(level + 1);
  
  std::vector< std::pair<int,Box> > isects;
  
  const auto& pmap = m_particles[level];
  for (const auto& kv : pmap) {
    const auto& pbox = kv.second.GetArrayOfStructs();
    for (auto it = pbox.cbegin(); it != pbox.cend(); ++it)
      {
	//
	// Find particle location on the finer level.
	//
	const IntVect& iv = Index(*it, level+1);
	//
	// Is it in the grown finer level?
	//
	fine.intersections(Box(iv,iv),isects,false,nGrow);
	//
	// Here we add the particle to each potential grid.
	//
	for (const auto& isec : isects)
	  {
	    //
	    // Create a copy.
	    //
	    ParticleType p = *it;
	    
	    //
	    // Set its id to indicate that it's a ghost.
	    //
	    p.m_idata.id = GhostParticleID;
	    
	    //
	    // Store it in the AoS.
	    //
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
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Redistribute (int lev_min, int lev_max, int nGrow, bool local)
{

  BL_PROFILE("ParticleContainer::Redistribute()");
  const int MyProc    = ParallelDescriptor::MyProc();
  Real      strttime  = ParallelDescriptor::second();
  
  // On startup there are cases where Redistribute() could be called
  // with a given finestLevel() where that AmrLevel has yet to be defined.
  int theEffectiveFinestLevel = m_gdb->finestLevel();

  while (!m_gdb->LevelDefined(theEffectiveFinestLevel))
      theEffectiveFinestLevel--;

  if (int(m_particles.size()) < theEffectiveFinestLevel+1) {
      if (Verbose()) {
          amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Redistribute() resizing containers from "
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
  Vector<std::map<std::pair<int, int>, Vector<Vector<ParticleType> > > > tmp_local;
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
  for (int i = 0; i < ParallelDescriptor::NProcs(); ++i) {
      tmp_remote[i].resize(num_threads);
  }

  // first pass: for each tile in parallel, in each thread copies the particles that
  // need to be moved into it's own, temporary buffer.
  for (int lev = lev_min; lev <= nlevs_particles; lev++) {
      auto& pmap = m_particles[lev];
      typename std::map<std::pair<int, int>, ParticleTileType >::iterator pmap_it;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single nowait
#endif
      for (pmap_it = pmap.begin(); pmap_it != pmap.end(); pmap_it++) {
#ifdef _OPENMP
#pragma omp task firstprivate(pmap_it)
#endif
          {
#ifdef _OPENMP
              int thread_num = omp_get_thread_num();
#else
              int thread_num = 0;
#endif
              int grid = pmap_it->first.first;
              int tile = pmap_it->first.second;
              auto& aos = pmap_it->second.GetArrayOfStructs();
              auto& soa = pmap_it->second.GetStructOfArrays();
              unsigned first = 0;
              unsigned npart = aos.numParticles();              
              ParticleLocData pld;
              if (npart != 0) {
                  for (unsigned pindex = 0; pindex < npart; ++pindex) {
                      ParticleType& p = aos[pindex];
                      
                      if (p.m_idata.id < 0) continue;                      
                      locateParticle(p, pld, lev_min, lev_max, nGrow);
                      if (p.m_idata.id < 0) continue;                      
                      
                      // The owner of the particle is the CPU owning the finest grid
                      // in state data that contains the particle.
                      const int who = ParticleDistributionMap(pld.m_lev)[pld.m_grid];
                      if (who == MyProc) {
                          if (pld.m_lev != lev || pld.m_grid != grid || pld.m_tile != tile) {
                              // We own it but must shift it to another place.
                              auto index = std::make_pair(pld.m_grid, pld.m_tile);
                              BL_ASSERT(tmp_local[pld.m_lev][index].size() == num_threads);
                              tmp_local[pld.m_lev][index][thread_num].push_back(p);
                              for (int comp = 0; comp < NArrayReal; ++comp) {
                                  Vector<Real>& arr = soa_local[pld.m_lev][index][thread_num].GetRealData(comp);
                                  arr.push_back(soa.GetRealData(comp)[pindex]);
                              }
                              for (int comp = 0; comp < NArrayInt; ++comp) {
                                  Vector<int>& arr = soa_local[pld.m_lev][index][thread_num].GetIntData(comp);
                                  arr.push_back(soa.GetIntData(comp)[pindex]);
                              }
                              
                              // Invalidate the particle so we can reclaim its space.
                              p.m_idata.id = -p.m_idata.id;
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
                          // Invalidate the particle so we can reclaim its space.
                          p.m_idata.id = -p.m_idata.id;
                      }
                      
                      // this is a valid particle
                      if (p.m_idata.id > 0) {
                          if (pindex != first) {
                              aos[first] = p;
                              for (int comp = 0; comp < NArrayReal; comp++)
                                  soa.GetRealData(comp)[first] = soa.GetRealData(comp)[pindex];
                              for (int comp = 0; comp < NArrayInt; comp++)
                                  soa.GetIntData(comp)[first] = soa.GetIntData(comp)[pindex];
                          }
                          ++first;
                      }
                  }
              
                  aos().erase(aos().begin() + first, aos().begin() + npart);
                  for (int comp = 0; comp < NArrayReal; comp++) {
                      Vector<Real>& rdata = soa.GetRealData(comp);
                      rdata.erase(rdata.begin() + first, rdata.begin() + npart);
                  }
                  for (int comp = 0; comp < NArrayInt; comp++) {
                      Vector<int>& idata = soa.GetIntData(comp);
                      idata.erase(idata.begin() + first, idata.begin() + npart);
                  }
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
        typename std::map<std::pair<int, int>, Vector<Vector<ParticleType> > >::iterator pmap_it;

        // we need to create any missing map entries in serial here
        for (pmap_it=tmp_local[lev].begin(); pmap_it != tmp_local[lev].end(); pmap_it++)
            m_particles[lev][pmap_it->first];

#ifdef _OPENMP
#pragma omp parallel
#pragma omp single nowait
#endif
      for (pmap_it=tmp_local[lev].begin(); pmap_it != tmp_local[lev].end(); pmap_it++) {
#ifdef _OPENMP
#pragma omp task firstprivate(pmap_it)
#endif
          {
              auto index = pmap_it->first;
              auto& aos = m_particles[lev][index].GetArrayOfStructs();
              auto& soa = m_particles[lev][index].GetStructOfArrays();
              auto& aos_tmp = pmap_it->second;
              auto& soa_tmp = soa_local[lev][index];
              for (int i = 0; i < num_threads; ++i) {
                  aos.insert(aos.end(), aos_tmp[i].begin(), aos_tmp[i].end());
                  aos_tmp[i].erase(aos_tmp[i].begin(), aos_tmp[i].end());
                  for (int comp = 0; comp < NArrayReal; ++comp) {
                      Vector<Real>& arr = soa.GetRealData(comp);
                      Vector<Real>& tmp = soa_tmp[i].GetRealData(comp);
                      arr.insert(arr.end(), tmp.begin(), tmp.end());
                      tmp.erase(tmp.begin(), tmp.end());
                  }
                  for (int comp = 0; comp < NArrayInt; ++comp) {
                      Vector<int>& arr = soa.GetIntData(comp);
                      Vector<int>& tmp = soa_tmp[i].GetIntData(comp);
                      arr.insert(arr.end(), tmp.begin(), tmp.end());
                      tmp.erase(tmp.begin(), tmp.end());
                  }
              }
          }
      }
  }

  for (auto& map_it : tmp_remote) {
      int who = map_it.first;
      not_ours[who];
  }  

  std::map<int, Vector<Vector<char> > >::iterator it;
#ifdef _OPENMP
#pragma omp parallel
#pragma omp single nowait
#endif
  for (it=tmp_remote.begin(); it != tmp_remote.end(); it++) {
#ifdef _OPENMP
#pragma omp task firstprivate(it)
#endif
      {
          int who = it->first;
          Vector<Vector<char> >& tmp = it->second;
          for (int i = 0; i < num_threads; ++i) {
              not_ours[who].insert(not_ours[who].end(), tmp[i].begin(), tmp[i].end());
              tmp[i].erase(tmp[i].begin(), tmp[i].end());
          }
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
      amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Redistribute() resizing m_particles from "
                     << m_particles.size() << " to " << theEffectiveFinestLevel+1 << '\n';
      BL_ASSERT(int(m_particles.size()) >= 2);
      BL_ASSERT(m_particles[m_particles.size()-1].empty());
      
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
      Real stoptime = ParallelDescriptor::second() - strttime;
      
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

namespace {

#ifdef BL_USE_MPI    
    long CountSnds(const std::map<int, Vector<char> >& not_ours, Vector<long>& Snds)
    {
        long NumSnds = 0;        
        for (const auto& kv : not_ours)
        {
            NumSnds       += kv.second.size();
            Snds[kv.first] = kv.second.size();
        }
        
        ParallelDescriptor::ReduceLongMax(NumSnds);

        return NumSnds;
    }

    long doHandShake(const std::map<int, Vector<char> >& not_ours,
                     Vector<long>& Snds, Vector<long>& Rcvs)
    {
        long NumSnds = CountSnds(not_ours, Snds);
        if (NumSnds == 0) return NumSnds;

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

        return NumSnds;
    }

    long doHandShakeLocal(const std::map<int, Vector<char> >& not_ours,
                          const Vector<int>& neighbor_procs, Vector<long>& Snds, Vector<long>& Rcvs)
    {
        long NumSnds = CountSnds(not_ours, Snds);
        if (NumSnds == 0) return NumSnds;
        
        const int num_rcvs = neighbor_procs.size();
        Vector<MPI_Status>  stats(num_rcvs);
        Vector<MPI_Request> rreqs(num_rcvs);
        
        const int SeqNum = ParallelDescriptor::SeqNum();
        
        // Post receives
        for (int i = 0; i < num_rcvs; ++i) {
            const int Who = neighbor_procs[i];
            const long Cnt = 1;
            
            BL_ASSERT(Who >= 0 && Who < ParallelDescriptor::NProcs());
            
            rreqs[i] = ParallelDescriptor::Arecv(&Rcvs[Who], Cnt, Who, SeqNum).req();
        }
        
        // Send.
        for (int i = 0; i < num_rcvs; ++i) {
            const int Who = neighbor_procs[i];
            const long Cnt = 1;
            
            BL_ASSERT(Who >= 0 && Who < ParallelDescriptor::NProcs());
            
            ParallelDescriptor::Send(&Snds[Who], Cnt, Who, SeqNum);        
        }
        
        if (num_rcvs > 0) {
            BL_MPI_REQUIRE( MPI_Waitall(num_rcvs, rreqs.data(), stats.data()) );
        }
        
        return NumSnds;
    }
#endif  // BL_USE_MPI
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
BuildRedistributeMask(int lev, int nghost) {
    
    BL_PROFILE("NeighborParticleContainer::BuildRedistributeMask");
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
                        const int tile = (*redistribute_mask_ptr)[mfi](iv, 1);
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
                 int lev_min, int lev_max, int nGrow, bool local)
{
    BL_PROFILE("ParticleContainer::RedistributeMPI()");
#if BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    const int NProcs = ParallelDescriptor::NProcs();
    
    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;
    if (local) {
        AMREX_ALWAYS_ASSERT(lev_min == 0);
        AMREX_ALWAYS_ASSERT(lev_max == 0);        
        BuildRedistributeMask(0);
        NumSnds = doHandShakeLocal(not_ours, neighbor_procs, Snds, Rcvs);
    }
    else {
        NumSnds = doHandShake(not_ours, Snds, Rcvs);
    }
    if (NumSnds == 0)
        return;  // There's no parallel work to do.

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
        BL_MPI_REQUIRE( MPI_Waitall(nrcvs, rreqs.data(), stats.data()) );
        
        ParticleLocData pld;
        
        if (recvdata.size() % superparticle_size != 0) {
            amrex::AllPrint() << "ParticleContainer::RedistributeMPI: sizes = "
                              << recvdata.size() << ", " << superparticle_size << "\n";
            amrex::Abort("ParticleContainer::RedistributeMPI: How did this happen?");
        }
        
        char* pbuf = recvdata.data();
        int npart = recvdata.size() / superparticle_size;
        for (int i = 0; i < npart; ++i) {
            ParticleType p;
            std::memcpy(&p, pbuf, particle_size);
            
            locateParticle(p, pld, lev_min, lev_max, nGrow);
            
            auto& ptile = m_particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)];
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
            
            pbuf += superparticle_size;
        }
    }
#endif /*BL_USE_MPI*/
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::OK (int lev_min, int lev_max, int nGrow) const
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::OK()");
    if (lev_max == -1)
        lev_max = finestLevel();

    BL_ASSERT(lev_max <= finestLevel());
    //
    // Check that each valid particle is in the proper container.
    //
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
                            amrex::AllPrint() << "PARTICLE NUMBER " << p.m_idata.id << '\n'
                                              << "POS  " << AMREX_D_TERM(p.m_rdata.pos[0], << p.m_rdata.pos[1], << p.m_rdata.pos[2]) << "\n"
                                              << "LEV  " << lev  << " " << pld.m_lev << '\n'
                                              << "GRID " << grid << " " << pld.m_grid << '\n';
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
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::AddParticlesAtLevel()");
    if (int(m_particles.size()) < level+1)
        {
            if (ParallelDescriptor::IOProcessor())
                {
                    std::cout << "ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::AddParticlesAtLevel resizing m_particles from "
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
                amrex::Abort("ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::AddParticlesAtLevel(): Can't add outside of domain\n");
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
::Checkpoint (const std::string&        dir,
              const std::string&        name,
              bool                      is_checkpoint,
              const Vector<std::string>& real_comp_names,
              const Vector<std::string>& int_comp_names) const
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::Checkpoint()");
    BL_ASSERT(OK());

    BL_ASSERT(sizeof(typename ParticleType::RealType) == 4 || sizeof(typename ParticleType::RealType) == 8);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProcNumber   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();
    //
    // We store the particles in a subdirectory of "dir".
    //
    std::string pdir = dir;

    if ( ! pdir.empty() && pdir[pdir.size()-1] != '/') {
        pdir += '/';
    }

    pdir += name;
    //
    // Make the particle directories if they don't already exist.
    //
    amrex::Print() << "Particle::Checkpoint: pdir levelDirectoriesCreated = "
                   << pdir << "  " << levelDirectoriesCreated << '\n';
    if ( ! levelDirectoriesCreated) {
      if (ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "IOIOIOIO:CD  Particle::Checkpoint:0:  creating directory:  "
	               << pdir << '\n';
        if ( ! amrex::UtilCreateDirectory(pdir, 0755)) {
            amrex::CreateDirectoryFailed(pdir);
	}
      }
      // Force other processors to wait until directory is built.
      ParallelDescriptor::Barrier();
    }

    //
    // The header contains the info we need to read back in the particles.
    //
    // Only the I/O processor writes to the header file.
    //
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
            std::cout << "PC " << NStructReal << " " << NArrayReal << " " << real_comp_names.size() 
            //BL_ASSERT(real_comp_names.size() == NStructReal + NArrayReal);
            for (int i = 0; i < real_comp_names.size(); ++i ) {
            //for (int i = 0; i < NStructReal + NArrayReal; ++i ) {
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
              if (ParallelDescriptor::IOProcessor()) {
                amrex::Print() << "IOIOIOIO:CD  Particle::Checkpoint:1:  creating directory:  "
		               << LevelDir << '\n';
	        if ( ! amrex::UtilCreateDirectory(LevelDir, 0755)) {
		  amrex::CreateDirectoryFailed(LevelDir);
	        }
	      }
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
			amrex::Print() << "::IOIOIOIO:  Unlinking file:  "<< FullFileName << '\n';

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
        amrex::Abort("ParticleContainer<NSR, NSI, NAR, NAI>::Checkpoint(): problem writing HdrFile");
      }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime, IOProcNumber);

        amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, "
	               << "NArrayInt>::Checkpoint() time: " << stoptime << '\n';
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

    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::CheckpointPre()");

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

    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::CheckpointPost()");

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
	      std::cout << ParallelDescriptor::MyProc() << "::IOIOIOIO:  Unlinking file:  "
	                << FullFileName << '\n';

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
        amrex::Abort("ParticleContainer<NSR, NSI, NAR, NAI>::CheckpointPost(): problem writing HdrFile");
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
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::WritePlotFile()");
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
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::WriteParticles()");

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
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart (const std::string& dir,
                                                                            const std::string& file,
                                                                            bool is_checkpoint)
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::Restart()");
  BL_ASSERT(!dir.empty());
  BL_ASSERT(!file.empty());
  
  const Real strttime = ParallelDescriptor::second();
  
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
      std::string msg("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart(): bad version string: ");
      msg += version;
      amrex::Error(version.c_str());
    }
  }
  else {
    std::string msg("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart(): unknown version string: ");
    msg += version;
    amrex::Abort(msg.c_str());
  }
  
  int dm;
  HdrFile >> dm;
  if (dm != AMREX_SPACEDIM)
    amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart(): dm != AMREX_SPACEDIM");
  
  int nr;
  HdrFile >> nr;
  if (nr != NStructReal + NArrayReal)
    amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart(): nr != NStructReal + NArrayReal");

  std::string comp_name;
  for (int i = 0; i < nr; ++i)
      HdrFile >> comp_name;

  int ni;
  HdrFile >> ni;
  if (ni != NStructInt + NArrayInt)
    amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart(): ni != NStructInt");
  
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
  
  int finest_level;
  HdrFile >> finest_level;
  BL_ASSERT(finest_level >= 0);
  
  Vector<int> ngrids(finest_level+1);
  BL_ASSERT(finest_level == finestLevel());
  for (int lev = 0; lev <= finest_level; lev++) {
    HdrFile >> ngrids[lev];
    BL_ASSERT(ngrids[lev] > 0);
    BL_ASSERT(ngrids[lev] == int(ParticleBoxArray(lev).size()));
  }

  resizeData();
  
  for (int lev = 0; lev <= finest_level; lev++) {
    Vector<int>  which(ngrids[lev]);
    Vector<int>  count(ngrids[lev]);
    Vector<long> where(ngrids[lev]);
    for (int i = 0; i < ngrids[lev]; i++) {
      HdrFile >> which[i] >> count[i] >> where[i];
    }
    
    for (MFIter mfi(*m_dummy_mf[lev]); mfi.isValid(); ++mfi) {
      const int grid = mfi.index();
      
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
	std::string msg("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart(): bad parameter: ");
	msg += how;
	amrex::Error(msg.c_str());
      }
      
      ParticleFile.close();
      
      if (!ParticleFile.good())
	amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart(): problem reading particles");
    }
  }

  BL_ASSERT(OK());
  
  if (m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;	
    ParallelDescriptor::ReduceRealMax(stoptime, ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::Restart() time: " << stoptime << '\n';
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
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::ReadParticles()");
    BL_ASSERT(cnt > 0);
    BL_ASSERT(lev < int(m_particles.size()));
    BL_ASSERT(lev >= 0 && lev <= finestLevel());
    BL_ASSERT(grd >= 0 && grd < ParticleBoxArray(lev).size());

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
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::WriteAsciiFile()");
    BL_ASSERT(!filename.empty());

    const Real strttime = ParallelDescriptor::second();
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

//      File << nparticles  << '\n';
//      File << NStructReal << '\n';
//      File << NStructInt  << '\n';
//      File << NArrayReal  << '\n';
//      File << NArrayInt   << '\n';
            
        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteAsciiFile(): problem writing file");
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
            std::ofstream File;

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

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
                amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteAsciiFile(): problem writing file");
	    
        }
	
        ParallelDescriptor::Barrier();
    }
    
    if (m_verbose > 1)
        {
            Real stoptime = ParallelDescriptor::second() - strttime;
            
            ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
            
            amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::WriteAsciiFile() time: " << stoptime << '\n';
        }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::WriteCoarsenedAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::WriteCoarsenedAsciiFile()");
    BL_ASSERT(!filename.empty());

    const Real strttime = ParallelDescriptor::second();
 
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
            amrex::Abort("ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::WriteCoarsenedAsciiFile(): problem writing file");
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
            std::ofstream File;

            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

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
                amrex::Abort("ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::WriteCoarsenedAsciiFile(): problem writing file");

        }

        ParallelDescriptor::Barrier();
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<NStructReal,NStructInt,NArrayReal, NArrayInt>::WriteCoarsenedAsciiFile() time: " << stoptime << '\n';
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::CrseToFine (const BoxArray&       cfba,
                                                                               const Vector<IntVect>& cells,
                                                                               Vector<IntVect>&       cfshifts,
                                                                               const Geometry&       gm,
                                                                               Vector<int>&           which,
                                                                               Vector<IntVect>&       pshifts) const
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::CrseToFine()");
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle, will we cross a  crse->fine boundary of the level
    // with coarsened fine BoxArray "cfba".  "cells" are as calculated from
    // CIC_Cells_Fracs().
    //
    const int M = cells.size();

    which.resize(M);
    cfshifts.resize(M);

    for (int i = 0; i < M; i++)
        which[i] =  0;

    bool result = false;

    for (int i = 0; i < M; i++)
    {
        if (cfba.contains(cells[i]))
        {
            result      = true;
            which[i]    = 1;
            cfshifts[i] = IntVect::TheZeroVector();
        }
        else if (!gm.Domain().contains(cells[i]))
        {
            BL_ASSERT(gm.isAnyPeriodic());
            //
            // Can the cell be shifted into cfba?
            //
            const Box bx(cells[i],cells[i]);

            gm.periodicShift(bx, gm.Domain(), pshifts);

            if (!pshifts.empty())
            {
                BL_ASSERT(pshifts.size() == 1);

                const Box& dbx = bx - pshifts[0];

                BL_ASSERT(dbx.ok());

                if (cfba.contains(dbx))
                {
                    //
                    // Note that pshifts[0] is from the coarse perspective.
                    // We'll later need to multiply it by ref ratio to use
                    // at the fine level.
                    //
                    result      = true;
                    which[i]    = 1;
                    cfshifts[i] = pshifts[0];
                }
            }
        }
    }

    return result;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
bool
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::FineToCrse (const ParticleType&                p,
					   int                                flev,
					   const Vector<IntVect>&              fcells,
					   const BoxArray&                    fvalid,
					   const BoxArray&                    compfvalid_grown,
					   Vector<IntVect>&                    ccells,
					   Vector<Real>&                       cfracs,
					   Vector<int>&                        which,
					   Vector<int>&                        cgrid,
					   Vector<IntVect>&                    pshifts,
					   std::vector< std::pair<int,Box> >& isects) const
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::FineToCrse()");
    BL_ASSERT(m_gdb != 0);
    BL_ASSERT(flev > 0);
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle we'll cross a fine->crse boundary.  Note that crossing
    // a periodic boundary, where the periodic shift lies in our valid region,
    // is not considered a Fine->Crse crossing.
    //
    const int M = fcells.size();

    which.resize(M);
    cgrid.resize(M);
    ccells.resize(M);
    cfracs.resize(M);

    for (int i = 0; i < M; i++)
    {
        cgrid[i] = -1;
        which[i] =  0;
    }

    ParticleLocData pld;
    Where(p, pld);

    const Box& ibx = amrex::grow(ParticleBoxArray(flev)[pld.m_grid],-1);

    BL_ASSERT(ibx.ok());

    if (ibx.contains(pld.m_cell))
        //
        // We're strictly contained in our valid box.
        // We can't cross a fine->crse boundary.
        //
        return false;

    if (!compfvalid_grown.contains(pld.m_cell))
        //
        // We're strictly contained in our "valid" region. Note that the valid
        // region contains any periodically shifted ghost cells that intersect
        // valid region.
        //
        return false;
    //
    // Otherwise ...
    //
    const Geometry& cgm = Geom(flev-1);
    const IntVect&  rr  = m_gdb->refRatio(flev-1);
    const BoxArray& cba = ParticleBoxArray(flev-1);

    Particle<NStructReal, NStructInt>::CIC_Cells_Fracs(p, cgm.ProbLo(), cgm.CellSize(), cfracs, ccells);

    bool result = false;

    for (int i = 0; i < M; i++)
    {
        IntVect ccell_refined = ccells[i]*rr;
        //
        // We've got to protect against the case when we're at the low
        // end of the domain because coarsening & refining don't work right
        // when indices go negative.
        //
        for (int dm = 0; dm < AMREX_SPACEDIM; dm++)
            ccell_refined[dm] = std::max(ccell_refined[dm], -1);

        if (!fvalid.contains(ccell_refined))
        {
            result   = true;
            which[i] = 1;

            Box cbx(ccells[i],ccells[i]);
    
            if (!cgm.Domain().contains(ccells[i]))
            {
                //
                // We must be at a periodic boundary.
                // Find valid box into which we can be periodically shifted.
                //
                BL_ASSERT(cgm.isAnyPeriodic());

                cgm.periodicShift(cbx, cgm.Domain(), pshifts);

                BL_ASSERT(pshifts.size() == 1);

                cbx -= pshifts[0];

                ccells[i] -= pshifts[0];
                BL_ASSERT(cbx.ok());
                BL_ASSERT(cgm.Domain().contains(cbx));
            }
            //
            // Which grid at the crse level do we need to update?
            //
            cba.intersections(cbx,isects,true,0);

            BL_ASSERT(!isects.empty());

            cgrid[i] = isects[0].first;  // The grid ID at crse level that we hit.
        }
    }

    return result;
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::FineCellsToUpdateFromCrse (const ParticleType&                p,
						      int                                lev,
						      const IntVect&                     ccell,
						      const IntVect&                     cshift,
						      Vector<int>&                        fgrid,
						      Vector<Real>&                       ffrac,
						      Vector<IntVect>&                    fcells,
						      std::vector< std::pair<int,Box> >& isects) const
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::FineCellsToUpdateFromCrse()");
    BL_ASSERT(lev >= 0);
    BL_ASSERT(lev < finestLevel());

    const Box&      fbx = amrex::refine(Box(ccell,ccell),m_gdb->refRatio(lev));
    const BoxArray& fba = ParticleBoxArray(lev+1);
    const Real*     plo = Geom(lev).ProbLo();
    const Real*     dx  = Geom(lev).CellSize();
    const Real*     fdx = Geom(lev+1).CellSize();

    if (cshift == IntVect::TheZeroVector())
    {
        BL_ASSERT(fba.contains(fbx));
    }
    //
    // Instead of clear()ing these we'll do a resize(0).
    // This'll preserve their capacity so that we'll only need
    // to do any memory allocation when their capacity needs to increase.
    //
    fgrid.resize(0);
    ffrac.resize(0);
    fcells.resize(0);
    //
    // Which fine cells does particle "p" (that wants to update "ccell") do we
    // touch at the finer level?
    //
    for (IntVect iv = fbx.smallEnd(); iv <= fbx.bigEnd(); fbx.next(iv))
    {
        bool touches = true;

        for (int k = 0; k < AMREX_SPACEDIM; k++)
        {
            const Real celllo = iv[k]  * fdx[k] + plo[k];
            const Real cellhi = celllo + fdx[k];

            if ((p.m_rdata.pos[k] < celllo) && (celllo > (p.m_rdata.pos[k] + dx[k]/2)))
                touches = false;

            if ((p.m_rdata.pos[k] > cellhi) && (cellhi < (p.m_rdata.pos[k] - dx[k]/2)))
                touches = false;
        }

        if (touches)
        {
            fcells.push_back(iv);
        }
    }

    Real sum_fine = 0;
    //
    // We need to figure out the fine fractions and the fine grid needed updating.
    //
    for (int j = 0; j < fcells.size(); j++)
    {
        IntVect& iv = fcells[j];

        Real the_frac = 1;

        for (int k = 0; k < AMREX_SPACEDIM; k++)
        {
            const Real celllo = (iv[k] * fdx[k] + plo[k]);

            if (p.m_rdata.pos[k] <= celllo)
            {
                const Real isecthi = p.m_rdata.pos[k] + dx[k]/2;

                the_frac *= std::min((isecthi - celllo),fdx[k]);
            }
            else
            {
                const Real cellhi  = (iv[k]+1) * fdx[k] + plo[k];
                const Real isectlo = p.m_rdata.pos[k] - dx[k]/2;

                the_frac *= std::min((cellhi - isectlo),fdx[k]);
            }
        }

        ffrac.push_back(the_frac);

        sum_fine += the_frac;

        if (cshift != IntVect::TheZeroVector())
        {
            //
            // Update to the correct fine cell needing updating.
            // Note that "cshift" is from the coarse perspective.
            //
            const IntVect& fshift = cshift * m_gdb->refRatio(lev);
            //
            // Update fcells[j] to indicate a shifted fine cell needing updating.
            //
            iv -= fshift;
        }

        fba.intersections(Box(iv,iv),isects,true,0);

        BL_ASSERT(!isects.empty());

        fgrid.push_back(isects[0].first);
    }

    BL_ASSERT(ffrac.size() == fcells.size());
    BL_ASSERT(fgrid.size() == fcells.size());
    //
    // Now adjust the fine fractions so they sum to one.
    //
    for (int j = 0; j < ffrac.size(); j++)
        ffrac[j] /= sum_fine;
}


//
// This is the multi-level version.
// The Array should be empty on input.
// There'll be finest_level+1 of them.
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensity (int rho_index, bool sub_cycle,
					      Vector<std::unique_ptr<MultiFab> >& mf_to_be_filled, 
					      int lev_min, int ncomp, int finest_level) const
{
    if (rho_index != 0) amrex::Abort("AssignDensity only works if rho_index = 0");

    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::AssignDensity()");
    BL_ASSERT(NStructReal >= 1);
    BL_ASSERT(NStructReal >= ncomp);
    BL_ASSERT(ncomp == 1 || ncomp == AMREX_SPACEDIM+1);

    if (finest_level == -1) {
      finest_level = finestLevel();
    }
    while (!m_gdb->LevelDefined(finest_level)) {
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

    Vector<std::unique_ptr<MultiFab> > mf_part;
    if (!all_grids_the_same)
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1-lev_min);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            const int lev_index = lev - lev_min;
            mf_part[lev_index].reset(new MultiFab(ParticleBoxArray(lev), 
						  ParticleDistributionMap(lev),
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
    std::map<int, Vector<ParticleCommData> > data;

    const Real stime = ParallelDescriptor::second();
    //
    // Minimum NUMCELLS required.
    //
    const int NUM_CELLS = AMREX_D_TERM(2,+2,+4);

    Vector<int>     cgrid(NUM_CELLS);
    Vector<int>    cwhich(NUM_CELLS),  fwhich(NUM_CELLS);
    Vector<Real>    fracs(NUM_CELLS),  cfracs(NUM_CELLS);
    Vector<IntVect> cells(NUM_CELLS),  ccells(NUM_CELLS), cfshifts(NUM_CELLS);

    ParticleCommData pb;
    //
    // I'm going to allocate these badboys here & pass'm into routines that use'm.
    // This should greatly cut down on memory allocation/deallocation.
    //
    Vector<IntVect>                    pshifts(27);
    std::vector< std::pair<int,Box> > isects;
    Vector<int>                        fgrid(NUM_CELLS);
    Vector<Real>                       ffracs(NUM_CELLS);
    Vector<IntVect>                    fcells;
    //
    // "fvalid" contains all the valid region of the MultiFab at this level, together
    // with any ghost cells lying outside the domain, that can be periodically shifted into the
    // valid region.  "compfvalid" is the complement of the "fvalid", while "compfvalid_grown" is 
    // "compfvalid" grown by one.  Using these we can figure out whether or not a cell is in the
    // valid region of our MultiFab as well as whether or not we're at a Fine->Crse boundary.
    //
    for (int lev = lev_min; lev <= finest_level; lev++)
    {
        const Geometry& gm        = Geom(lev);
        const Geometry& gm_fine   = (lev < finest_level) ? Geom(lev+1) : gm;
        const Geometry& gm_coarse = (lev > 0) ? Geom(lev-1) : gm;
        const Box&      dm        = gm.Domain();
        const Real*     dx        = gm.CellSize();
        const Real*     plo       = gm.ProbLo();
        const Real*     dx_fine   = (lev < finest_level) ? Geom(lev+1).CellSize() : dx;
        const Real*     dx_coarse = (lev > 0) ? Geom(lev-1).CellSize() : dx;
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
        const bool GridsCoverDomain = fvalid.numPts() == Geom(lev).Domain().numPts();

	const auto& pmap = m_particles[lev];
	for (const auto& kv : pmap) {
	  const int gid = kv.first.first;
	  const auto& aos = kv.second.GetArrayOfStructs();
	  FArrayBox&  fab = (*mf[lev_index])[gid];
	  for (const auto& p : aos)
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
                            fab(cells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
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
                            fab(cells[i],n) += p.m_rdata.arr[AMREX_SPACEDIM + n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
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
                                        (*mf[lev_index+1])[grid](fcells[j],0) += p.m_rdata.arr[AMREX_SPACEDIM] * ffracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index+1])[grid](fcells[j],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * ffracs[j];
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
				    pb.m_data[0] = p.m_rdata.arr[AMREX_SPACEDIM] *  ffracs[j];
				  }
				  
				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * ffracs[j];
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
                                        (*mf[lev_index-1])[grid](ccells[j],0) += p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[j];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[grid](ccells[j],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[j];
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
				    pb.m_data[0] = p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[j];
				  }
                                  
				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[j];
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
                        AnyCrseToFine = CrseToFine(cfba,cells,cfshifts,gm,cwhich,pshifts);
		    }
                    //
                    // lev_index > 0 means that we don't do F->C for lower levels
                    // This may mean that the mass fraction is off.
                    //
                    bool AnyFineToCrse = false;
                    if (lev_index > 0 && !GridsCoverDomain)
                        AnyFineToCrse = FineToCrse(p,lev,cells,fvalid,compfvalid_grown,ccells,cfracs,fwhich,cgrid,pshifts,isects);

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
                                fab(cells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
                            }
                            //
                            // Sum up momenta in next components.
                            //
                            for (int n = 1; n < ncomp; n++) {
                                fab(cells[i],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
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
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[i];
                                    }
                                    //
                                    // Sum up momenta in next components.
                                    //
                                    for (int n = 1; n < ncomp; n++) {
                                        (*mf[lev_index-1])[cgrid[i]](ccells[i],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[i];
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
				    pb.m_data[0] = p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[i];
				  }

				  //
				  // Sum up momenta in next components.
				  //
				  for (int n = 1; n < ncomp; n++) {
				    pb.m_data[n] = p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * cfracs[i];
				  }
				  data[who].push_back(pb);
                                }

                                sum_crse += cfracs[i];
                            }
                        }
                        //
                        // We've updated the Crse cells.  Now we have to update the fine
                        // cells in such a way that the total amount of mass we move
                        // around is precisely p.m_rdata.arr[AMREX_SPACEDIM]. In other words, the fractions
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
                                    fab(cells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
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
                                    fab(cells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i] * mult;
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i] * mult;
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
                                    fab(cells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
                                }
                                //
                                // Sum up momenta in next components.
                                //
                                for (int n = 1; n < ncomp; n++) {
                                    fab(cells[i],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
				}

                                sum += fracs[i];
                            }
                            else
                            {
                                //
                                // We're at a Crse->Fine boundary.
                                //
                                FineCellsToUpdateFromCrse(p,lev,cells[i],cfshifts[i],fgrid,ffracs,fcells,isects);

                                for (int j = 0, nfcells = fcells.size(); j < nfcells; j++)
                                {
                                    const int who = mf[lev_index+1]->DistributionMap()[fgrid[j]];

                                    if (who == ParallelDescriptor::MyProc())
                                    {
                                        //
                                        // Sum up mass in first component.
                                        //
                                        {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],0) += p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i] * ffracs[j];
                                        }
                                        //
                                        // Sum up momenta in next components.
                                        //
                                        for (int n = 1; n < ncomp; n++) {
                                            (*mf[lev_index+1])[fgrid[j]](fcells[j],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i] * ffracs[j];
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
					pb.m_data[0] = p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i] * ffracs[j];
				      }

				      //
				      // Sum up momenta in next components.
				      //
				      for (int n = 1; n < ncomp; n++) {
					pb.m_data[0] = p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i] * ffracs[j];
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
        const Geometry& gm        = Geom(lev);
        const Real*     dx        = gm.CellSize();
        const Real      vol       = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

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

        amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensity(multi-level) time: " << etime << '\n';
    }
}

//
// Used by AssignDensity (Vector<std::unique_ptr<MultiFab> >& mf).
//
// Passes data needed by Crse->Fine or Fine->Crse to CPU that needs it.
//
// We store the data that needs to be sent in "data". Note that m_lev is the
// real particle level, while mf may start at a fine level (e.g. lvls 1 and 2).
// Consequently, we must subtract lev_min from m_lev to get the mf lev.
//

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensityDoit (int               rho_index,
					       Vector<std::unique_ptr<MultiFab> >&             mf,
					       std::map<int, Vector<ParticleCommData> >& data,
					       int               ncomp,
					       int               lev_min) const
{
    if (rho_index != 0) amrex::Abort("AssignDensityDoit only works if rho_index = 0");

    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::AssignDensityDoit()");
    BL_ASSERT(NStructReal >= ncomp);

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

    Vector<int> Snds(NProcs,0), Rcvs(NProcs,0);

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
    Vector<int>().swap(Snds);
    Vector<int>().swap(Rcvs);
    //
    // The data we want to receive.
    //
    const int iChunkSize = 2 + AMREX_SPACEDIM;
    const int rChunkSize = ncomp;

    Vector<int>                    irecvdata (NumRcvs*iChunkSize);
    Vector<typename ParticleType::RealType> rrecvdata (NumRcvs*rChunkSize);

    Vector<int>         index(2*RcvCnts.size());
    Vector<MPI_Status>  stats(2*RcvCnts.size());
    Vector<MPI_Request> rreqs(2*RcvCnts.size());

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
    Vector<int>                             isenddata;
    Vector<typename ParticleType::RealType> rsenddata;

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

	auto& pbox = data[Who];

        int ioff = 0, roff = 0;
	for (const auto& p : pbox)
        {
	  isenddata[ioff+0] = p.m_lev  - lev_min;
	  isenddata[ioff+1] = p.m_grid;

	  AMREX_D_TERM(isenddata[ioff+2] = p.m_cell[0];,
		 isenddata[ioff+3] = p.m_cell[1];,
		 isenddata[ioff+4] = p.m_cell[2];);

	  ioff += iChunkSize;

	  for (int n = 0; n < ncomp; n++) {
	    rsenddata[roff+n] = p.m_data[n];
	  }

	  roff += ncomp;
        }

	Vector<ParticleCommData>().swap(pbox);

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
            const IntVect cell (AMREX_D_DECL(idata[2],idata[3],idata[4]));

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

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensityFort (int rho_index,
                                                                                      Vector<std::unique_ptr<MultiFab> >& mf_to_be_filled, 
                                                                                      int lev_min, int ncomp, int finest_level) const
{

    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::AssignDensityFort()");
    
    if (rho_index != 0) amrex::Abort("AssignDensity only works if rho_index = 0");
    
    BL_ASSERT(NStructReal >= 1);
    BL_ASSERT(NStructReal >= ncomp);
    BL_ASSERT(ncomp == 1 || ncomp == AMREX_SPACEDIM+1);

    if (finest_level == -1) {
        finest_level = finestLevel();
    }
    while (!m_gdb->LevelDefined(finest_level)) {
        finest_level--;
    }

    // Create the space for mf_to_be_filled, regardless of whether we'll need a temporary mf
    mf_to_be_filled.resize(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; lev++)
    { 
        mf_to_be_filled[lev].reset(new MultiFab(m_gdb->boxArray(lev),
				                m_gdb->DistributionMap(lev),
					        ncomp, 2));
	mf_to_be_filled[lev]->setVal(0.0);
    }

    // Test whether the grid structure of the boxArray is the same
    //       as the ParticleBoxArray at all levels 
    bool all_grids_the_same = true; 
    for (int lev = lev_min; lev <= finest_level; lev++) {
        if (!OnSameGrids(lev, *mf_to_be_filled[lev])) {
	    all_grids_the_same = false;
	    break;
	}
    }

    Vector<std::unique_ptr<MultiFab> > mf_part;
    if (!all_grids_the_same)
    { 
        // Create the space for the temporary, mf_part
        mf_part.resize(finest_level+1);
        for (int lev = lev_min; lev <= finest_level; lev++)
        {
            mf_part[lev].reset(new MultiFab(ParticleBoxArray(lev), 
					    ParticleDistributionMap(lev),
					    ncomp, 2));
	    mf_part[lev]->setVal(0.0);
        }
    }

    auto & mf = (all_grids_the_same) ? mf_to_be_filled : mf_part;

    if (finest_level == 0)
    {
        //
        // Just use the far simpler single-level version.
        //
        AssignCellDensitySingleLevelFort(rho_index, *mf[0], 0, ncomp);
        //
        // I believe that we don't need any information in ghost cells so we don't copy those.
        //
        if ( ! all_grids_the_same) {
            mf_to_be_filled[0]->copy(*mf[0],0,0,ncomp);
	}
        return;
    }

    PhysBCFunct cphysbc, fphysbc;
    int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaries
    int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
    Vector<BCRec> bcs(1, BCRec(lo_bc, hi_bc));
    PCInterp mapper;
    
    Vector<std::unique_ptr<MultiFab> > tmp(finest_level+1);
    for (int lev = lev_min; lev <= finest_level; ++lev) {
        const BoxArray& ba = mf[lev]->boxArray();
        const DistributionMapping& dm = mf[lev]->DistributionMap();
        tmp[lev].reset(new MultiFab(ba, dm, 1, 0));
        tmp[lev]->setVal(0.0);
    }

    for (int lev = lev_min; lev <= finest_level; ++lev) {

        AssignCellDensitySingleLevelFort(rho_index, *mf[lev], lev, 1, 0);

        if (lev < finest_level) {
            amrex::InterpFromCoarseLevel(*tmp[lev+1], 0.0, *mf[lev],
                                         rho_index, rho_index, ncomp,
                                         m_gdb->Geom(lev), m_gdb->Geom(lev+1),
                                         cphysbc, fphysbc,
                                         m_gdb->refRatio(lev), &mapper, bcs);
        }
        
        if (lev > lev_min) {
            // Note - this will double count the mass on the coarse level in 
            // regions covered by the fine level, but this will be corrected
            // below in the call to average_down.
            amrex::sum_fine_to_coarse(*mf[lev], *mf[lev-1],
                                      rho_index, 1, m_gdb->refRatio(lev-1),
                                      m_gdb->Geom(lev-1), m_gdb->Geom(lev));
        }
        
        mf[lev]->plus(*tmp[lev], rho_index, ncomp, 0);
    }
    
    for (int lev = finest_level - 1; lev >= lev_min; --lev) {
        amrex::average_down(*mf[lev+1], *mf[lev], rho_index, ncomp, m_gdb->refRatio(lev));
    }
    
    if (!all_grids_the_same) {
        for (int lev = lev_min; lev <= finest_level; lev++) {
            mf_to_be_filled[lev]->copy(*mf_part[lev],0,0,1);
        }
    }
    if (lev_min > 0) {
        int nlevels = finest_level - lev_min + 1;
        for (int i = 0; i < nlevels; i++)
            {
                mf_to_be_filled[i] = std::move(mf_to_be_filled[i+lev_min]);
            }
        mf_to_be_filled.resize(nlevels);
    }
}

//
// This is the single-level version -- it takes either cell-centered or node-centered MF's
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensitySingleLevel (int rho_index,
                                                                                             MultiFab& mf_to_be_filled,
                                                                                             int       lev,
                                                                                             int       ncomp,
                                                                                             int       particle_lvl_offset) const
{
  BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::AssignDensitySingleLevel()");
  BL_ASSERT(NStructReal >= 1);
  BL_ASSERT(ncomp == 1 || ncomp == AMREX_SPACEDIM+1);
  
  if (lev >= int(m_particles.size()))
    {
      // Don't do anything if there are no particles at this level.
      return;
    }

  // Keep the same external interface to the applications, but if the
  if (mf_to_be_filled.is_nodal()) {
    NodalDepositionSingleLevel(rho_index, mf_to_be_filled,lev,ncomp,particle_lvl_offset);
  }
  else if (mf_to_be_filled.boxArray().ixType().cellCentered()) {
    AssignCellDensitySingleLevel(rho_index, mf_to_be_filled,lev,ncomp,particle_lvl_offset);
  }
  else {
    amrex::Abort("AssignCellDensitySingleLevel: mixed type not supported");
  }
}

// This is the single-level version for cell-centered density
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
AssignCellDensitySingleLevelFort (int rho_index,
                                  MultiFab& mf_to_be_filled,
                                  int       lev,
                                  int       ncomp,
                                  int       particle_lvl_offset) const
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::AssignCellDensitySingleLevelFort()");
    
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
    const Real      strttime    = ParallelDescriptor::second();
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
            const Box& box = fab.box();
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
      Real stoptime = ParallelDescriptor::second() - strttime;
      
      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
      
      amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensity(single-level) time: " << stoptime << '\n';
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InterpolateFort (Vector<std::unique_ptr<MultiFab> >& mesh_data, 
                                                                                    int lev_min, int lev_max)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InterpolateFort()");
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        InterpolateSingleLevelFort(*mesh_data[lev], lev); 
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
InterpolateSingleLevelFort (MultiFab& mesh_data, int lev)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InterpolateSingleLevelFort()");
    
    if (mesh_data.nGrow() < 1)
        amrex::Error("Must have at least one ghost cell when in InterpolateSingleLevelFort");
    
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
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::AssignCellDensitySingleLevel()");
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

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = Geom(lev);
    const Real*     plo         = gm.ProbLo();
    const Real*     dx_particle = Geom(lev + particle_lvl_offset).CellSize();
    const Real*     dx          = gm.CellSize();
    const auto&     pmap        = m_particles[lev];

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) {
      amrex::Error("AssignDensity: problem must be periodic in no or all directions");
    }
    
    for (MFIter mfi(*mf_pointer); mfi.isValid(); ++mfi) {
        (*mf_pointer)[mfi].setVal(0);
    }

    for (const auto& kv : pmap) {
      const int grid = kv.first.first;
      const auto& pbx = kv.second.GetArrayOfStructs();
      FArrayBox& fab = (*mf_pointer)[grid];
      auto N = pbx.size();
	
        Vector<Real>    fracs;
        Vector<IntVect> cells;
	
#ifdef _OPENMP
#pragma omp parallel for default(none) private(fracs,cells) shared(N,plo,dx,dx_particle,gm,fab,ncomp,pbx)
#endif
	for (size_t ip = 0; ip < N; ++ip) {
	  const ParticleType& p = pbx[ip];
	  
	  if (p.m_idata.id <= 0)
	    continue;
	    
	  const int M = ParticleType::CIC_Cells_Fracs(p, plo, dx, dx_particle, fracs, cells);

	  // If this is not fully periodic then we have to be careful that the
	  // particle's support leaves the domain unless we specifically want to ignore
	  // any contribution outside the boundary (i.e. if allow_particles_near_boundary = true). 
	  // We test this by checking the low and high corners respectively.
	  if ( ! gm.isAllPeriodic() && ! allow_particles_near_boundary) {
	    if ( ! gm.Domain().contains(cells[0]) || ! gm.Domain().contains(cells[M-1])) {
	      amrex::Error("AssignDensity: if not periodic, all particles must stay away from the domain boundary");
	    }
	  }
	  
	  for (int i = 0; i < M; i++) {
	    if ( !fab.box().contains(cells[i]) )
	      continue;

	    // If the domain is not periodic and we want to let particles
	    // live near the boundary but "throw away" the contribution that 
	    // does not fall into the domain ...
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
	      fab(cells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
	    }
	    // 
	    // Sum up momenta in next components.
	    //
	    for (int n = 1; n < ncomp; n++)
#ifdef _OPENMP
#pragma omp atomic
#endif
	      fab(cells[i],n) += p.m_rdata.arr[AMREX_SPACEDIM+n] * p.m_rdata.arr[AMREX_SPACEDIM] * fracs[i];
	  }
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

    mf_pointer->mult(1/vol, 0, 1);

    // If mf_to_be_filled is not defined on the particle_box_array, then we need
    // to copy here from mf_pointer into mf_to_be_filled.   I believe that we don't
    // need any information in ghost cells so we don't copy those.
    if (mf_pointer != &mf_to_be_filled) {
      mf_to_be_filled.copy(*mf_pointer,0,0,ncomp);
      delete mf_pointer;
    }
    
    if (m_verbose > 1) {
      Real stoptime = ParallelDescriptor::second() - strttime;
      
      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
      
      amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::AssignDensity(single-level) time: " << stoptime << '\n';
    }
}

//
// This is the single-level version for nodal density
//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::NodalDepositionSingleLevel (int rho_index,
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
        mf_pointer = new MultiFab(amrex::convert(ParticleBoxArray(lev),
						  mf_to_be_filled.boxArray().ixType()),
				  ParticleDistributionMap(lev),
				  ncomp, mf_to_be_filled.nGrow());
    }

    const Real      strttime    = ParallelDescriptor::second();
    const Geometry& gm          = Geom(lev);
    const Real*     dx          = gm.CellSize();
    const auto&     pmap        = m_particles[lev];

    if (gm.isAnyPeriodic() && ! gm.isAllPeriodic()) 
        amrex::Error("AssignDensity: problem must be periodic in no or all directions");

    mf_pointer->setVal(0.0);

    Vector<IntVect> cells;
    cells.resize(8);

    Vector<Real> fracs;
    fracs.resize(8);

#if (AMREX_SPACEDIM > 1)
    Vector<Real> sx;
    sx.resize(2);
    Vector<Real> sy;
    sy.resize(2);
#endif
#if (AMREX_SPACEDIM > 2)
    Vector<Real> sz;
    sz.resize(2);
#endif

    ParticleLocData pld;

    for (const auto& kv : pmap) {
      const auto& pbx = kv.second.GetArrayOfStructs();
      const int grid = kv.first.first;
      FArrayBox& fab = (*mf_pointer)[grid];

      for (const auto& p : pbx)
        {
            if (p.m_idata.id <= 0) {
	      continue;
	    }
	    
	    Where(p, pld);

#if (AMREX_SPACEDIM == 1)
            cells[0] = pld.m_cell;
            cells[1] = pld.m_cell+IntVect(1);

            Real x = p.m_rdata.pos[0] / dx[0];

            int i = pld.m_cell[0];

            Real xint = x - i;

            for (int i = 0; i < 2; i++)
            {
               fab(cells[0],0) += p.m_rdata.arr[AMREX_SPACEDIM+rho_index] * (1.0 - xint);
               fab(cells[1],0) += p.m_rdata.arr[AMREX_SPACEDIM+rho_index] *        xint ;
            }
#elif (AMREX_SPACEDIM == 2)
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

            for (int ind = 0; ind < 4; ind++)
            {
               fab(cells[ind],0) += p.m_rdata.arr[AMREX_SPACEDIM+rho_index] * fracs[ind];
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
               fab(cells[i],0) += p.m_rdata.arr[AMREX_SPACEDIM+rho_index] * fracs[i];
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
    const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

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

        amrex::Print() << "ParticleContainer<N>::NodalDepositionSingleLevel time: " << stoptime << '\n';
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

    const Real strttime  = ParallelDescriptor::second();
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
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::moveKick() time: " << stoptime << '\n';
    }
    //
    // No need for Redistribute(), we only change the velocity.
    //
}
