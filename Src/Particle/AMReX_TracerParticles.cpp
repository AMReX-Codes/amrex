
#include <AMReX_TracerParticles.H>

namespace amrex {

//
// Uses midpoint method to advance particles using umac.
//
void
TracerParticleContainer::AdvectWithUmac (MultiFab* umac, int lev, Real dt)
{
    BL_PROFILE("TracerParticleContainer::AdvectWithUmac()");
    BL_ASSERT(OK(lev, lev, umac[0].nGrow()-1));
    BL_ASSERT(lev >= 0 && lev < GetParticles().size());

    AMREX_D_TERM(BL_ASSERT(umac[0].nGrow() >= 1);,
           BL_ASSERT(umac[1].nGrow() >= 1);,
           BL_ASSERT(umac[2].nGrow() >= 1););

    AMREX_D_TERM(BL_ASSERT(!umac[0].contains_nan());,
           BL_ASSERT(!umac[1].contains_nan());,
           BL_ASSERT(!umac[2].contains_nan()););

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();
    const Real*     plo      = geom.ProbLo();

    Vector<std::unique_ptr<MultiFab> > raii_umac(AMREX_SPACEDIM);
    Vector<MultiFab*> umac_pointer(AMREX_SPACEDIM);
    // We assume that if umac[0]'s boxArray matches then the others will too...
    if (OnSameGrids(lev, umac[0]))
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
	    umac_pointer[i] = &umac[i];
	}
    }
    else
    {
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        {
	    int ng = umac[i].nGrow();
	    raii_umac[i].reset(new MultiFab(amrex::convert(m_gdb->ParticleBoxArray(lev),
                                                           IntVect::TheDimensionVector(i)),
					    m_gdb->ParticleDistributionMap(lev),
					    umac[i].nComp(), ng));
					    
	    umac_pointer[i] = raii_umac[i].get();
	    umac_pointer[i]->copy(umac[i],0,0,umac[i].nComp(),ng,ng);
        }
    }

    for (int ipass = 0; ipass < 2; ipass++)
    {
        auto& pmap = GetParticles(lev);
	for (auto& kv : pmap) {
	  int grid = kv.first.first;
	  auto& pbox = kv.second.GetArrayOfStructs();
	  const int n = pbox.size();

	  FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
						 &((*umac_pointer[1])[grid]),
						 &((*umac_pointer[2])[grid])) };

#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];

                if (p.m_idata.id <= 0) continue;

                const IntVect& cc_cell = Index(p, lev);

                const Real len[AMREX_SPACEDIM] = { AMREX_D_DECL((p.m_rdata.pos[0]-plo[0])/dx[0] + Real(0.5),
                                                       (p.m_rdata.pos[1]-plo[1])/dx[1] + Real(0.5),
                                                       (p.m_rdata.pos[2]-plo[2])/dx[2] + Real(0.5)) };

                const IntVect cell(AMREX_D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

                const Real frac[AMREX_SPACEDIM] = { AMREX_D_DECL(len[0]-cell[0], len[1]-cell[1], len[2]-cell[2]) };

                for (int d = 0; d < AMREX_SPACEDIM; d++)
                {
                    IntVect ecell = cell;

                    ecell[d] = cc_cell[d] + 1;

                    Real efrac[AMREX_SPACEDIM] = { AMREX_D_DECL(frac[0], frac[1], frac[2]) };

                    efrac[d] = (p.m_rdata.pos[d]-plo[d])/dx[d] - cc_cell[d];

                    for (int j = 0; j < AMREX_SPACEDIM; j++)
                    {
                        if (efrac[j] > 1) efrac[j] = 1;
                        if (efrac[j] < 0) efrac[j] = 0;
                    }

                    Real vel = ParticleType::InterpDoit(*fab[d], ecell, efrac, 0);

                    if (ipass == 0)
                    {
                        //
                        // Save old position and the vel & predict location at dt/2.
                        //
                        p.m_rdata.arr[AMREX_SPACEDIM+d] = p.m_rdata.pos[d];
                        p.m_rdata.pos[d] += 0.5*dt*vel;
                    }
                    else
                    {
		        //
                        // Update to final time using the orig position and the vel at dt/2.
                        //
		        p.m_rdata.pos[d]  = p.m_rdata.arr[AMREX_SPACEDIM+d] + dt*vel;
                        // Save the velocity for use in Timestamp().
			p.m_rdata.arr[AMREX_SPACEDIM+d] = vel;
                    }
                }
            }
        }
    }
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "TracerParticleContainer::AdvectWithUmac() time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
	});
#endif
    }
}

//
// Uses midpoint method to advance particles using cell-centered velocity
//
void
TracerParticleContainer::AdvectWithUcc (const MultiFab& Ucc, int lev, Real dt)
{
    BL_ASSERT(Ucc.nGrow() > 0);
    BL_ASSERT(OK(lev, lev, Ucc.nGrow()-1));
    BL_ASSERT(lev >= 0 && lev < GetParticles().size());

    BL_ASSERT(!Ucc.contains_nan());

    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(lev);

    BL_ASSERT(OnSameGrids(lev, Ucc));

    int idx[AMREX_SPACEDIM] = {AMREX_D_DECL(0,1,2)};

    for (int ipass = 0; ipass < 2; ipass++)
    {
        auto& pmap = GetParticles(lev);
	for (auto& kv : pmap) {
	  int grid = kv.first.first;
	  auto& pbox = kv.second.GetArrayOfStructs();
	  const int n    = pbox.size();
	  const FArrayBox& fab = Ucc[grid];
	    
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < n; i++)
            {
                ParticleType& p = pbox[i];

                if (p.m_idata.id <= 0) continue;

		Real v[AMREX_SPACEDIM];

		ParticleType::Interp(p, geom, fab, idx, v, AMREX_SPACEDIM);

		if (ipass == 0) {
		    //
		    // Save old position and the vel & predict location at dt/2.
		    //
		    for (int d = 0; d < AMREX_SPACEDIM; d++)
		    {
			p.m_rdata.arr[AMREX_SPACEDIM+d] = p.m_rdata.pos[d];
                        p.m_rdata.pos[d] += 0.5*dt*v[d];
                    }
		} else {
		    //
		    // Update to final time using the orig position and the vel at dt/2.
		    //
		    for (int d = 0; d < AMREX_SPACEDIM; d++)
		    {
                        p.m_rdata.pos[d]  = p.m_rdata.arr[AMREX_SPACEDIM+d] + dt*v[d];
                        // Save the velocity for use in Timestamp().
			p.m_rdata.arr[AMREX_SPACEDIM+d] = v[d];
                    }
                }
            }
        }
    }
    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "TracerParticleContainer::AdvectWithUcc() time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
	});
#endif
    }
}

void
TracerParticleContainer::Timestamp (const std::string&      basename,
				    const MultiFab&         mf,
				    int                     lev,
				    Real                    time,
				    const std::vector<int>& indices)
{
    BL_PROFILE("TracerParticleContainer::Timestamp()");
    //
    // basename -> base filename for the output file
    // mf       -> the multifab
    // lev      -> level to check for particles
    // time     -> simulation time (will be recorded in Timestamp file)
    // indices  -> indices into mf that we output
    //
    BL_ASSERT(lev >= 0);
    BL_ASSERT(time >= 0);
    BL_ASSERT(!basename.empty());
    BL_ASSERT(lev <= m_gdb->finestLevel());

    const Real strttime = ParallelDescriptor::second();

    const int   MyProc    = ParallelDescriptor::MyProc();
    const int   NProcs    = ParallelDescriptor::NProcs();
    // We'll spread the output over this many files.
    int nOutFiles(64);
    ParmParse pp("particles");
    pp.query("particles_nfiles",nOutFiles);
    if(nOutFiles == -1) {
      nOutFiles = NProcs;
    }
    nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
    const int   nSets     = ((NProcs + (nOutFiles - 1)) / nOutFiles);
    const int   mySet     = (MyProc / nOutFiles);

    for (int iSet = 0; iSet < nSets; ++iSet)
      {
        if (mySet == iSet)
	  {
            //
            // Do we have any particles at this level that need writing?
            //
            bool gotwork = false;
	    
            const auto& pmap = GetParticles(lev);
	    for (auto& kv : pmap) {
              const auto& pbox = kv.second.GetArrayOfStructs();
	      for (const auto& p : pbox) {
		if (p.m_idata.id > 0) {
		  gotwork = true;
		  break;
		}
	      }
	      if (gotwork) break;
	    }

            if (gotwork)
	      {
                std::string FileName = amrex::Concatenate(basename + '_', MyProc % nOutFiles, 2);
		
                std::ofstream TimeStampFile;
		
                VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

                TimeStampFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

                TimeStampFile.open(FileName.c_str(), std::ios::out|std::ios::app|std::ios::binary);

                TimeStampFile.setf(std::ios_base::scientific,std::ios_base::floatfield);

                TimeStampFile.precision(10);

                TimeStampFile.seekp(0, std::ios::end);

                if (!TimeStampFile.good())
                    amrex::FileOpenFailed(FileName);

                const int       M  = indices.size();
                const BoxArray& ba = mf.boxArray();

                std::vector<Real> vals(M);

		for (auto& kv : pmap) {
		  int grid = kv.first.first;
		  const auto& pbox = kv.second.GetArrayOfStructs();
		  const Box&       bx   = ba[grid];
		  const FArrayBox& fab  = mf[grid];

		  for (const auto& p : pbox)
                    {
		      if (p.m_idata.id <= 0) continue;
		      
		      const IntVect& iv = Index(p,lev);
		      
		      if (!bx.contains(iv) && !ba.contains(iv)) continue;
		      
		      TimeStampFile << p.m_idata.id  << ' ' << p.m_idata.cpu << ' ';
		      
		      AMREX_D_TERM(TimeStampFile << p.m_rdata.pos[0] << ' ';,
			     TimeStampFile << p.m_rdata.pos[1] << ' ';,
			     TimeStampFile << p.m_rdata.pos[2] << ' ';);
		      
		      TimeStampFile << time;
		      //
		      // AdvectWithUmac stores the velocity in rdata ...
		      //
		      AMREX_D_TERM(TimeStampFile << ' ' << p.m_rdata.arr[AMREX_SPACEDIM+0];,
			     TimeStampFile << ' ' << p.m_rdata.arr[AMREX_SPACEDIM+1];,
			     TimeStampFile << ' ' << p.m_rdata.arr[AMREX_SPACEDIM+2];);
		      
		      if (M > 0)
                        {
			  ParticleType::Interp(p,m_gdb->Geom(lev),fab,&indices[0],&vals[0],M);
			  
			  for (int i = 0; i < M; i++)
                            {
			      TimeStampFile << ' ' << vals[i];
                            }
                        }
		      
		      TimeStampFile << '\n';
                    }
                }
		
                TimeStampFile.flush();
                TimeStampFile.close();
            }

            const int iBuff     = 0;
            const int wakeUpPID = (MyProc + nOutFiles);
            const int tag       = (MyProc % nOutFiles);

            if (wakeUpPID < NProcs)
                ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
        if (mySet == (iSet + 1))
        {
            //
            // Next set waits.
            //
            int       iBuff;
            const int waitForPID = (MyProc - nOutFiles);
            const int tag        = (MyProc % nOutFiles);

            ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
        }
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "TracerParticleContainer::Timestamp: lev: " << lev << " time: " << stoptime << '\n';
        }
#ifdef BL_LAZY
        });
#endif
    }
}

}
