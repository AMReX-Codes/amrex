#include <AMReX_TracerParticle_mod_K.H>
#include <AMReX_TracerParticles.H>
#include "AMReX_TracerParticles.H"
#include "AMReX_TracerParticle_mod_K.H"
#include <AMReX_Print.H>
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

    const Real      strttime = amrex::second();
    const Geometry& geom     = m_gdb->Geom(lev);
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();
    
    Vector<std::unique_ptr<MultiFab> > raii_umac(AMREX_SPACEDIM);
    Vector<MultiFab*> umac_pointer(AMREX_SPACEDIM);
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
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n = aos.size();
            auto p_pbox = aos().data();
            const FArrayBox* fab[AMREX_SPACEDIM] = { AMREX_D_DECL(&((*umac_pointer[0])[grid]),
                                                                  &((*umac_pointer[1])[grid]),
                                                                  &((*umac_pointer[2])[grid])) };
	  
            //array of these pointers to pass to the GPU
            amrex::GpuArray<amrex::Array4<const Real>, AMREX_SPACEDIM>
            const umacarr {AMREX_D_DECL((*fab[0]).array(),
                                        (*fab[1]).array(),
                                        (*fab[2]).array() )};

            amrex::ParallelFor(n,   
                               [=] AMREX_GPU_DEVICE (int i)      
            {					
                ParticleType& p = p_pbox[i];
                if (p.m_idata.id <= 0) return;
                Real v[AMREX_SPACEDIM];
                mac_interpolate(p, plo, dxi, umacarr, v);
                if (ipass == 0)
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.pos[dim];     
                        p.m_rdata.pos[dim] += 0.5*dt*v[dim];   		 
                    }    		  
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.m_rdata.pos[dim]  = p.m_rdata.arr[AMREX_SPACEDIM+dim] + dt*v[dim];
                        p.m_rdata.arr[AMREX_SPACEDIM+dim] = v[dim];
                    }
                }
            });
        }
    }
    
    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "TracerParticleContainer::AdvectWithUmac() time: " << stoptime << '\n';
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
    BL_PROFILE("TracerParticleContainer::AdvectWithUcc()");
    BL_ASSERT(Ucc.nGrow() > 0);
    BL_ASSERT(OK(lev, lev, Ucc.nGrow()-1));
    BL_ASSERT(lev >= 0 && lev < GetParticles().size());
    BL_ASSERT(!Ucc.contains_nan());

    const Real          strttime = amrex::second();
    const Geometry&     geom     = m_gdb->Geom(lev);
    const auto          plo      = geom.ProbLoArray();
    const auto          dxi      = geom.InvCellSizeArray();

    BL_ASSERT(OnSameGrids(lev, Ucc));

    for (int ipass = 0; ipass < 2; ipass++)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (ParIterType pti(*this, lev); pti.isValid(); ++pti)
        {
            int grid    = pti.index();
            auto& ptile = ParticlesAt(lev, pti);
            auto& aos  = ptile.GetArrayOfStructs();
            const int n          = aos.size();
            const FArrayBox& fab = Ucc[grid];
            const auto uccarr = fab.array();
            auto  p_pbox = aos().data();

            amrex::ParallelFor(n,
                               [=] AMREX_GPU_DEVICE (int i)
            {
                ParticleType& p  = p_pbox[i];
                if (p.m_idata.id <= 0) return;
                Real v[AMREX_SPACEDIM];
                
                cic_interpolate(p, plo, dxi, uccarr, v);
                
                if (ipass == 0)
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.m_rdata.arr[AMREX_SPACEDIM+dim] = p.m_rdata.pos[dim];
                        p.m_rdata.pos[dim] += 0.5*dt*v[dim];
                    }                  
                }
                else
                {
                    for (int dim=0; dim < AMREX_SPACEDIM; dim++)
                    {
                        p.m_rdata.pos[dim]  = p.m_rdata.arr[AMREX_SPACEDIM+dim] + dt*v[dim];
                        p.m_rdata.arr[AMREX_SPACEDIM+dim] = v[dim];
                    }
                }
            });
        }
    }
    
    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "TracerParticleContainer::AdvectWithUcc() time: " << stoptime << '\n';
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

    const Real strttime = amrex::second();

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
	      for (int k = 0; k < pbox.size(); ++k)
	      {
		const ParticleType& p = pbox[k];
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

                VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

                std::ofstream TimeStampFile;

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

		  for (int k = 0; k < pbox.size(); ++k)
		    {
		      const ParticleType& p = pbox[k];

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
        Real stoptime = amrex::second() - strttime;

#ifdef BL_LAZY
        Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "TracerParticleContainer::Timestamp: lev: " << lev << " time: " << stoptime << '\n';
#ifdef BL_LAZY
        });
#endif
    }
}
}
