#ifndef _EBPARTICLECONTAINER_H_
#define _EBPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"

namespace amrex {
    
    ///
    /// This class demonstrates how to have particles interact with 
    /// embedded boundaries. Particles are initialized inside a 
    /// circle (or sphere in 3D) with randomly oriented velocities.
    /// Particles do not interact with each other, but they do
    /// specularly reflect off of the embedded boundaries.
    ///
    class EBParticleContainer
        : public ParticleContainer<2*BL_SPACEDIM>
    {
    public:

        using MyParIter      = ParIter<2*BL_SPACEDIM>;
        using MyConstParIter = ParConstIter<2*BL_SPACEDIM>;

        EBParticleContainer(const Geometry            & geom, 
                            const DistributionMapping & dmap,
                            const BoxArray            & ba);
        
        void InitParticles(const Real& radius, const RealVect& center);
        
        void moveParticles(const Real dt);
        
        void bounceWalls(const MultiFab& dummy,
                         const MultiFab* volfrac,
                         const MultiCutFab* bndrycent,
                         std::array<const MultiCutFab*, AMREX_SPACEDIM>& areafrac,
                         const Real dt);
        
        void writeParticles(int n);
    };
    
}

#endif // _EBPARTICLECONTAINER_H_
