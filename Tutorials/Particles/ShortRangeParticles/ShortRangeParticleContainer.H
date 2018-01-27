#ifndef _SHORTRANGEPARTICLECONTAINER_H_
#define _SHORTRANGEPARTICLECONTAINER_H_

#include "AMReX_NeighborParticles.H"

namespace amrex {

///
/// The particles in this container carry velocity and acceleration in addition to the 
/// usual components.
///
/// In 2D, the struct is: x y vx vy ax ay id cpu
/// In 3D: x y z vx vy vz ax ay az id cpu
///
class ShortRangeParticleContainer
    : public NeighborParticleContainer<2*BL_SPACEDIM, 0>
{
public:

    ///
    /// This particle container fills a mask for quickly computing
    /// neighbor grids / tiles for a given particle
    ///
    ShortRangeParticleContainer(const Geometry            & geom, 
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
    /// Move the particles according to their forces, reflecting at domain boundaries
    ///
    void moveParticles(const Real dt);

    ///
    /// Save the particle data in an ASCII format
    ///
    void writeParticles(int n);
};

}

#endif // _SHORTRANGEPARTICLECONTAINER_H_
