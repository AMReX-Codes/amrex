#ifndef WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_

#include "WarpXParticleContainer.H"
#include <AMReX_DenseBins.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

class CollisionType
{
public:
    int m_species1;
    int m_species2;

    CollisionType(
        const std::vector<std::string>& species_names,
        std::string collision_name);

    static void doCoulombCollisionsWithinTile ( int lev, amrex::MFIter const& mfi,
        std::unique_ptr<WarpXParticleContainer>& species1,
        std::unique_ptr<WarpXParticleContainer>& species2 );

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
