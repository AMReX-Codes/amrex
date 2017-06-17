#include "MyParticleContainer.H"

#include "deposit_F.H"

#ifdef CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif // CUDA

using namespace amrex;

MyParticleContainer::MyParticleContainer(const Geometry            & geom,
					 const DistributionMapping & dmap,
					 const BoxArray            & ba)
  : ParticleContainer<1+2*BL_SPACEDIM> (geom, dmap, ba)
{}

void MyParticleContainer::InitParticles(int num_particles, Real mass) {
  bool serialize = true;
  int iseed = 451;
  ParticleInitData pdata = {mass, 0, 0, 0, 0, 0, 0};
  InitRandom(num_particles, iseed, pdata, serialize); 
}

void MyParticleContainer::Deposit(MultiFab& partMF) {

  BL_PROFILE("Particle Deposit.");

  const int lev = 0;
  const Geometry& gm  = Geom(lev);
  const Real*     plo = gm.ProbLo();
  const Real*     dx  = gm.CellSize();  

  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {    
    const auto& particles = pti.GetArrayOfStructs();
    int nstride = particles.dataShape().first;
    const long np  = pti.numParticles();    
    FArrayBox& rhofab = partMF[pti];
    const Box& box    = rhofab.box();        

#if CUDA
    cuda_deposit_cic(particles.data(), nstride, np, pti.index(),
                     rhofab.dataPtr(), box.loVect(), box.hiVect(), 
                     plo, dx);
#else 
    deposit_cic(particles.data(), nstride, np,
  		rhofab.dataPtr(), box.loVect(), box.hiVect(), 
  		plo, dx);
#endif // CUDA    
  }
  partMF.SumBoundary(gm.periodicity());
}

void MyParticleContainer::Interpolate(MultiFab& acc) {
  
  BL_PROFILE("Particle Interpolate.");
  
  const int lev = 0;
  const Geometry& gm  = Geom(lev);
  const Real*     plo = gm.ProbLo();
  const Real*     dx  = gm.CellSize();  

  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {    
    const auto& particles = pti.GetArrayOfStructs();
    int nstride = particles.dataShape().first;
    const long np  = pti.numParticles();    
    FArrayBox& accfab = acc[pti];
    const Box& box    = accfab.box();

#if CUDA    
    cuda_interpolate_cic(particles.data(), nstride, np, pti.index(),
			 accfab.dataPtr(), box.loVect(), box.hiVect(), 
			 plo, dx);
#else
    interpolate_cic(particles.data(), nstride, np,
		    accfab.dataPtr(), box.loVect(), box.hiVect(), 
		    plo, dx);
#endif // CUDA    
  }
}

void MyParticleContainer::Push() {
  
  BL_PROFILE("Particle Push.");

  const int lev = 0;
  const Geometry& gm  = Geom(lev);
  const Real*     plo = gm.ProbLo();
  const Real*     dx  = gm.CellSize();
  
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      const auto& particles = pti.GetArrayOfStructs();
      int nstride = particles.dataShape().first;
      const long np  = pti.numParticles();
#if CUDA
      push_particles(particles.data(), nstride, np);
#else
      cuda_push_particles(particles.data(), nstride, np);
#endif // CUDA
  }
}
