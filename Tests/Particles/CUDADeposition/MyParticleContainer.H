#ifndef _MYPARTICLECONTAINER_H_
#define _MYPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

namespace amrex {

class MyParticleContainer
    : public ParticleContainer<1>
{
public:
  
  using MyParIter = amrex::ParIter<1>;
  
  MyParticleContainer(const Geometry            & geom, 
		      const DistributionMapping & dmap,
		      const BoxArray            & ba);
  
  void InitParticles(int num_particles, Real mass);
  
  void Deposit(MultiFab& partMF);
};

}

#endif // _MYPARTICLECONTAINER_H_
