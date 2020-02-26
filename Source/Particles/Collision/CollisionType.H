/* Copyright 2019 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
#define WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_

#include "Particles/WarpXParticleContainer.H"
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

    /** Perform all binary collisions within a tile
     *
     * @param lev AMR level of the tile
     * @param mfi iterator for multifab
     * @param species1/2 pointer to species container
     * @param isSameSpecies true if collision is between same species
     * @param CoulombLog user input Coulomb logrithm
     *
     */

    static void doCoulombCollisionsWithinTile (
        int const lev, amrex::MFIter const& mfi,
        std::unique_ptr<WarpXParticleContainer>& species1,
        std::unique_ptr<WarpXParticleContainer>& species2,
        bool const isSameSpecies, amrex::Real const CoulombLog );

};

#endif // WARPX_PARTICLES_COLLISION_COLLISIONTYPE_H_
