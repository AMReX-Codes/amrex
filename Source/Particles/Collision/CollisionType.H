#ifndef WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_

#include "WarpXParticleContainer.H"
#include <AMReX_DenseBins.H>
#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

class CollisionType
{
public:
    int  m_species1_index;
    int  m_species2_index;
    bool m_isSameSpecies;
    amrex::Real m_CoulombLog;

    CollisionType(
        const std::vector<std::string>& species_names,
        std::string const collision_name);

    static void doCoulombCollisionsWithinTile (
        int const lev, amrex::MFIter const& mfi,
        std::unique_ptr<WarpXParticleContainer>& species1,
        std::unique_ptr<WarpXParticleContainer>& species2,
        bool const isSameSpecies, amrex::Real const CoulombLog );

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
