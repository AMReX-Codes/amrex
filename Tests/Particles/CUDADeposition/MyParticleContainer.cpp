#include <cuda_runtime_api.h>
#include <cuda.h>

#include "MyParticleContainer.H"

#include "deposit_F.H"

using namespace amrex;

MyParticleContainer::MyParticleContainer(const Geometry            & geom,
					 const DistributionMapping & dmap,
					 const BoxArray            & ba)
  : ParticleContainer<1> (geom, dmap, ba)
{
  m_np = 0;
  psize = sizeof(ParticleType);
}

void MyParticleContainer::InitParticles(int num_particles, Real mass) {
  bool serialize = true;
  int iseed = 451;
  ParticleInitData pdata = {mass};
  InitRandom(num_particles, iseed, pdata, serialize);
  m_np = num_particles;
  
  cudaMalloc((void**) &device_particles, m_np*psize); 
}

void MyParticleContainer::CopyParticlesToDevice() {
  const int lev = 0;
  int offset = 0;
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
    const auto& particles = pti.GetArrayOfStructs();
    const long np  = pti.numParticles();
    cudaMemcpy(device_particles + offset,
	       particles.data(), np*psize, cudaMemcpyHostToDevice);
    offset += np;
  }
}

void MyParticleContainer::CopyParticlesFromDevice() {
  const int lev = 0;
  int offset = 0;
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
    auto& particles = pti.GetArrayOfStructs();
    const long np  = pti.numParticles();
    cudaMemcpy(particles.data(), device_particles + offset,
	       np*psize, cudaMemcpyDeviceToHost);
    offset += np;
  }
}

void MyParticleContainer::Deposit(MultiFab& partMF) {

  CopyParticlesToDevice();

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
    
    deposit_cic((Real*) device_particles, nstride, np,
  		rhofab.dataPtr(), box.loVect(), box.hiVect(), 
  		plo, dx);
  }
  
  partMF.SumBoundary(gm.periodicity());

  CopyParticlesFromDevice();
}
