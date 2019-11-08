#ifndef WARPX_PARTICLES_COLLISION_COULOMBCOLLISIONS_H_
#define WARPX_PARTICLES_COLLISION_COULOMBCOLLISIONS_H_

#include "WarpXParticleContainer.H"

void
doCoulombCollisionsWithinTile ( int lev, amrex::MFIter const& mfi,
      std::unique_ptr< WarpXParticleContainer>& species1,
      std::unique_ptr< WarpXParticleContainer>& species2 );

#endif // WARPX_PARTICLES_COLLISION_COULOMBCOLLISIONS_H_
