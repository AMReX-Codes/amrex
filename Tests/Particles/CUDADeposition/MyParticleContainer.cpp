#include "MyParticleContainer.H"

#include "deposit_F.H"

using namespace amrex;

MyParticleContainer::MyParticleContainer(const Geometry            & geom,
					 const DistributionMapping & dmap,
					 const BoxArray            & ba)
  : ParticleContainer<1> (geom, dmap, ba)
{}

void MyParticleContainer::InitParticles(int num_particles, Real mass) {
  bool serialize = true;
  int iseed = 451;
  ParticleInitData pdata = {mass};
  InitRandom(num_particles, iseed, pdata, serialize);
}

void MyParticleContainer::Deposit(MultiFab& partMF) {

  const int lev = 0;
  const Geometry& gm          = Geom(lev);
  const Real*     plo         = gm.ProbLo();
  const Real*     dx          = gm.CellSize();
  
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {    
    const auto& particles = pti.GetArrayOfStructs();
    int nstride = particles.dataShape().first;
    const long np  = pti.numParticles();    
    FArrayBox& rhofab = partMF[pti];
    const Box& box    = rhofab.box();        
    deposit_cic(particles.data(), nstride, np,
		rhofab.dataPtr(), box.loVect(), box.hiVect(), 
		plo, dx);
  }
}
