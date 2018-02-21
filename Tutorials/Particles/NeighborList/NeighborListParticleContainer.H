#ifndef _NEIGHBORLISTPARTICLECONTAINER_H_
#define _NEIGHBORLISTPARTICLECONTAINER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_NeighborParticles.H"

namespace amrex {

class NeighborListParticleContainer
    : public NeighborParticleContainer<2*BL_SPACEDIM, 0>
{
public:

    ///
    /// This particle container fills a mask for quickly computing
    /// neighbor grids / tiles for a given particle
    ///
    NeighborListParticleContainer(const Geometry            & geom, 
                                  const DistributionMapping & dmap,
                                  const BoxArray            & ba,
                                  int                         nneighbor);

    ///
    /// Init one particle per cell with random velocities
    ///
    void InitParticles();

    ///
    /// Compute the short range forces on a tile's worth of particles.
    /// fillNeighbors must have already been called.
    ///
    void computeForces();

    ///
    /// Compute the short range forces on a tile's worth of particles using
    /// the neighbor list instead of the N^2 approach.
    /// fillNeighbors must have already been called.
    ///
    void computeForcesNL();

    ///
    /// Move the particles according to their forces, reflecting at domain boundaries
    ///
    void moveParticles(const Real dt);

    ///
    /// Save the particle data in an ASCII format
    ///
    void writeParticles(int n);

private:
    
    inline virtual bool check_pair(const ParticleType& p1, const ParticleType& p2) final {
        return AMREX_D_TERM(   (p1.pos(0) - p2.pos(0))*(p1.pos(0) - p2.pos(0)) ,
                             + (p1.pos(1) - p2.pos(1))*(p1.pos(1) - p2.pos(1)) ,
                             + (p1.pos(2) - p2.pos(2))*(p1.pos(2) - p2.pos(2)) ) <= cutoff*cutoff;
    }
    
    static constexpr Real cutoff = 1.e-2;
    static constexpr Real min_r  = 1.e-4;    
};
    
}

#endif // _NEIGHBORLISTPARTICLECONTAINER_H_
