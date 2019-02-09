#ifndef MD_K_H_
#define MD_K_H_

#include "Constants.H"
#include "MDParticleContainer.H"

AMREX_GPU_DEVICE AMREX_INLINE
bool check_pair(const MDParticleContainer::ParticleType& p1, 
                const MDParticleContainer::ParticleType& p2)
{
    amrex::Real d0 = (p1.pos(0) - p2.pos(0));
    amrex::Real d1 = (p1.pos(1) - p2.pos(1));
    amrex::Real d2 = (p1.pos(2) - p2.pos(2));    
    amrex::Real dsquared = d0*d0 + d1*d1 + d2*d2;   
    return (dsquared <= Params::cutoff*Params::cutoff);
}

#endif
