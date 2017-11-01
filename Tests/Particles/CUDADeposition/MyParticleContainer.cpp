#include "MyParticleContainer.H"

#include "deposit_F.H"

#ifdef AMREX_USE_CUDA
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif // CUDA

using namespace amrex;

MyParticleContainer::MyParticleContainer(const Geometry            & geom,
					 const DistributionMapping & dmap,
					 const BoxArray            & ba)
  : ParticleContainer<1+2*BL_SPACEDIM> (geom, dmap, ba)
{
  m_np = 0;
  m_ngrids = 0;
  psize = sizeof(ParticleType);
}

void MyParticleContainer::InitParticles(int num_particles, Real mass) {
  bool serialize = true;
  int iseed = 451;
  ParticleInitData pdata = {mass, 0, 0, 0, 0, 0, 0};
  InitRandom(num_particles, iseed, pdata, serialize);
  m_np = num_particles;

#ifdef AMREX_USE_CUDA  

  cudaError_t err_code;
  err_code = cudaMalloc((void**) &device_particles, m_np*psize); 
  if (err_code != cudaSuccess) {
    amrex::Print() << "Could not allocate device memory in InitParticles. \n";
    amrex::Abort();
  }

  const int lev = 0;
  int offset = 0;
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
    const auto& particles = pti.GetArrayOfStructs();
    const long np  = pti.numParticles();
    particle_counts.push_back(np);
    particle_offsets.push_back(offset);
    offset += np;
  }

  m_ngrids = particle_counts.size();
  err_code = cudaMalloc((void**) &device_particle_offsets, m_ngrids*sizeof(int)); 
  if (err_code != cudaSuccess) {
    amrex::Print() << "Could not allocate device memory in InitParticles. \n";
    amrex::Abort();
  }
  
  err_code = cudaMalloc((void**) &device_particle_counts,  m_ngrids*sizeof(int)); 
  if (err_code != cudaSuccess) {
    amrex::Print() << "Could not allocate device memory in InitParticles. \n";
    amrex::Abort();
  }

#endif
}

void MyParticleContainer::CopyParticlesToDevice() {
#ifdef AMREX_USE_CUDA
  const int lev = 0;
  int offset = 0;
  particle_counts.clear();
  particle_offsets.clear();
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
    const auto& particles = pti.GetArrayOfStructs();
    const long np  = pti.numParticles();
    cudaMemcpy(device_particles + offset,
	       particles.data(), np*psize, cudaMemcpyHostToDevice);
    particle_counts.push_back(np);
    particle_offsets.push_back(offset);
    offset += np;
  }

  cudaMemcpy(device_particle_counts, particle_counts.data(), 
	     m_ngrids*sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(device_particle_offsets, particle_offsets.data(),
	     m_ngrids*sizeof(int), cudaMemcpyHostToDevice);

#endif
}

void MyParticleContainer::CopyParticlesFromDevice() {
#ifdef AMREX_USE_CUDA
  const int lev = 0;
  int offset = 0;
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
    auto& particles = pti.GetArrayOfStructs();
    const long np  = pti.numParticles();
    cudaMemcpy(particles.data(), device_particles + offset,
	       np*psize, cudaMemcpyDeviceToHost);
    offset += np;
  }
#endif
}

void MyParticleContainer::Deposit(MultiFab& partMF, MultiFab& acc) {

  BL_PROFILE("Particle Process.");

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
    FArrayBox& accfab = acc[pti];
    const Box& box    = rhofab.box();        

#if CUDA
    cuda_deposit_cic((Real*) device_particles, nstride, m_np, np, 
		device_particle_counts, device_particle_offsets, 
		m_ngrids, pti.uniqueIndex(),
  		rhofab.dataPtr(), box.loVect(), box.hiVect(), 
  		plo, dx);
#else 
    deposit_cic(particles.data(), nstride, np,
  		rhofab.dataPtr(), box.loVect(), box.hiVect(), 
  		plo, dx);
#endif // CUDA    
  }

  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {    
    const auto& particles = pti.GetArrayOfStructs();
    int nstride = particles.dataShape().first;
    const long np  = pti.numParticles();    
    FArrayBox& rhofab = partMF[pti];
    FArrayBox& accfab = acc[pti];
    const Box& box    = rhofab.box();        

#if CUDA    
    cuda_interpolate_cic((Real*) device_particles, nstride, m_np, np,
			 device_particle_counts, device_particle_offsets, 
			 m_ngrids, pti.uniqueIndex(),
			 accfab.dataPtr(), box.loVect(), box.hiVect(), 
			 plo, dx);
#else
    interpolate_cic(particles.data(), nstride, m_np,
		    accfab.dataPtr(), box.loVect(), box.hiVect(), 
		    plo, dx);
#endif // CUDA    
  }

#if CUDA
    cuda_push_particles((Real*) device_particles, 11, m_np);
#else
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {    
    const auto& particles = pti.GetArrayOfStructs();
    int nstride = particles.dataShape().first;
    const long np  = pti.numParticles();    
    FArrayBox& rhofab = partMF[pti];
    FArrayBox& accfab = acc[pti];
    const Box& box    = rhofab.box();        
    
    push_particles(particles.data(), nstride, np);    
  }
#endif // CUDA    
    
  partMF.SumBoundary(gm.periodicity());
  
  //  CopyParticlesFromDevice();

}
