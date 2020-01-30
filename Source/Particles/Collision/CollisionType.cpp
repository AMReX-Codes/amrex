/* Copyright 2019 Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "CollisionType.H"
#include "ShuffleFisherYates.H"
#include "ElasticCollisionPerez.H"
#include <WarpX.H>

CollisionType::CollisionType(
    const std::vector<std::string>& species_names,
    std::string const collision_name)
{

    #if defined WARPX_DIM_RZ
    amrex::Abort("Collisions only work in Cartesian geometry for now.");
    #endif

    // read collision species
    std::vector<std::string> collision_species;
    amrex::ParmParse pp(collision_name);
    pp.getarr("species", collision_species);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(collision_species.size() == 2,
    "Collision species must name exactly two species.");

    // default Coulomb log, if < 0, will be computed automatically
    m_CoulombLog = -1.0;
    pp.query("CoulombLog", m_CoulombLog);

    for (int i=0; i<species_names.size(); i++)
    {
        if (species_names[i] == collision_species[0])
        { m_species1_index = i; }
        if (species_names[i] == collision_species[1])
        { m_species2_index = i; }
    }

    if (collision_species[0] == collision_species[1])
        m_isSameSpecies = true;
    else
        m_isSameSpecies = false;

}

using namespace amrex;
// Define shortcuts for frequently-used type names
using ParticleType = WarpXParticleContainer::ParticleType;
using ParticleTileType = WarpXParticleContainer::ParticleTileType;
using ParticleBins = DenseBins<ParticleType>;
using index_type = ParticleBins::index_type;

namespace {

    /* Find the particles and count the particles that are in each cell.
       Note that this does *not* rearrange particle arrays */
    ParticleBins
    findParticlesInEachCell( int const lev, MFIter const& mfi,
                             ParticleTileType const& ptile) {

        // Extract particle structures for this tile
        int const np = ptile.numParticles();
        ParticleType const* particle_ptr = ptile.GetArrayOfStructs()().data();

        // Extract box properties
        Geometry const& geom = WarpX::GetInstance().Geom(lev);
        Box const& cbx = mfi.tilebox(IntVect::TheZeroVector()); //Cell-centered box
        const auto lo = lbound(cbx);
        const auto dxi = geom.InvCellSizeArray();
        const auto plo = geom.ProbLoArray();

        // Find particles that are in each cell;
        // results are stored in the object `bins`.
        ParticleBins bins;
        bins.build(np, particle_ptr, cbx,
            // Pass lambda function that returns the cell index
            [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p) noexcept -> IntVect
            {
                return IntVect(AMREX_D_DECL((p.pos(0)-plo[0])*dxi[0] - lo.x,
                                            (p.pos(1)-plo[1])*dxi[1] - lo.y,
                                            (p.pos(2)-plo[2])*dxi[2] - lo.z));
            });

        return bins;
    }

}

/** Perform all binary collisions within a tile
 *
 * @param lev AMR level of the tile
 * @param mfi iterator for multifab
 * @param species1/2 pointer to species container
 * @param isSameSpecies true if collision is between same species
 * @param CoulombLog user input Coulomb logrithm
 *
 */
void CollisionType::doCoulombCollisionsWithinTile
    ( int const lev, MFIter const& mfi,
    std::unique_ptr<WarpXParticleContainer>& species_1,
    std::unique_ptr<WarpXParticleContainer>& species_2,
    bool const isSameSpecies, Real const CoulombLog )
{

    if ( isSameSpecies ) // species_1 == species_2
    {
        // Extract particles in the tile that `mfi` points to
        ParticleTileType& ptile_1 = species_1->ParticlesAt(lev, mfi);

        // Find the particles that are in each cell of this tile
        ParticleBins bins_1 = findParticlesInEachCell( lev, mfi, ptile_1 );

        // Loop over cells, and collide the particles in each cell

        // Extract low-level data
        int const n_cells = bins_1.numBins();
        // - Species 1
        auto& soa_1 = ptile_1.GetStructOfArrays();
        ParticleReal * const AMREX_RESTRICT ux_1 =
            soa_1.GetRealData(PIdx::ux).data();
        ParticleReal * const AMREX_RESTRICT uy_1 =
            soa_1.GetRealData(PIdx::uy).data();
        ParticleReal * const AMREX_RESTRICT uz_1  =
            soa_1.GetRealData(PIdx::uz).data();
        ParticleReal const * const AMREX_RESTRICT w_1 =
            soa_1.GetRealData(PIdx::w).data();
        index_type* indices_1 = bins_1.permutationPtr();
        index_type const* cell_offsets_1 = bins_1.offsetsPtr();
        Real q1 = species_1->getCharge();
        Real m1 = species_1->getMass();

        const Real dt = WarpX::GetInstance().getdt(lev);
        Geometry const& geom = WarpX::GetInstance().Geom(lev);
        #if (AMREX_SPACEDIM == 2)
        auto dV = geom.CellSize(0) * geom.CellSize(1);
        #elif (AMREX_SPACEDIM == 3)
        auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
        #endif

        // Loop over cells
        amrex::ParallelFor( n_cells,
            [=] AMREX_GPU_DEVICE (int i_cell) noexcept
            {
                // The particles from species1 that are in the cell `i_cell` are
                // given by the `indices_1[cell_start_1:cell_stop_1]`
                index_type const cell_start_1 = cell_offsets_1[i_cell];
                index_type const cell_stop_1  = cell_offsets_1[i_cell+1];
                index_type const cell_half_1 = (cell_start_1+cell_stop_1)/2;

                // Do not collide if there is only one particle in the cell
                if ( cell_stop_1 - cell_start_1 >= 2 )
                {
                    // shuffle
                    ShuffleFisherYates(
                        indices_1, cell_start_1, cell_half_1 );

                    // Call the function in order to perform collisions
                    ElasticCollisionPerez(
                        cell_start_1, cell_half_1,
                        cell_half_1, cell_stop_1,
                        indices_1, indices_1,
                        ux_1, uy_1, uz_1, ux_1, uy_1, uz_1, w_1, w_1,
                        q1, q1, m1, m1, Real(-1.0), Real(-1.0),
                        dt, CoulombLog, dV );
                }
            }
        );
    }
    else // species_1 != species_2
    {
        // Extract particles in the tile that `mfi` points to
        ParticleTileType& ptile_1 = species_1->ParticlesAt(lev, mfi);
        ParticleTileType& ptile_2 = species_2->ParticlesAt(lev, mfi);

        // Find the particles that are in each cell of this tile
        ParticleBins bins_1 = findParticlesInEachCell( lev, mfi, ptile_1 );
        ParticleBins bins_2 = findParticlesInEachCell( lev, mfi, ptile_2 );

        // Loop over cells, and collide the particles in each cell

        // Extract low-level data
        int const n_cells = bins_1.numBins();
        // - Species 1
        auto& soa_1 = ptile_1.GetStructOfArrays();
        ParticleReal * const AMREX_RESTRICT ux_1 =
            soa_1.GetRealData(PIdx::ux).data();
        ParticleReal * const AMREX_RESTRICT uy_1 =
            soa_1.GetRealData(PIdx::uy).data();
        ParticleReal * const AMREX_RESTRICT uz_1 =
            soa_1.GetRealData(PIdx::uz).data();
        ParticleReal const * const AMREX_RESTRICT w_1 =
            soa_1.GetRealData(PIdx::w).data();
        index_type* indices_1 = bins_1.permutationPtr();
        index_type const* cell_offsets_1 = bins_1.offsetsPtr();
        Real q1 = species_1->getCharge();
        Real m1 = species_1->getMass();
        // - Species 2
        auto& soa_2 = ptile_2.GetStructOfArrays();
        Real* ux_2  = soa_2.GetRealData(PIdx::ux).data();
        Real* uy_2  = soa_2.GetRealData(PIdx::uy).data();
        Real* uz_2  = soa_2.GetRealData(PIdx::uz).data();
        Real* w_2   = soa_2.GetRealData(PIdx::w).data();
        index_type* indices_2 = bins_2.permutationPtr();
        index_type const* cell_offsets_2 = bins_2.offsetsPtr();
        Real q2 = species_2->getCharge();
        Real m2 = species_2->getMass();

        const Real dt = WarpX::GetInstance().getdt(lev);
        Geometry const& geom = WarpX::GetInstance().Geom(lev);
        #if (AMREX_SPACEDIM == 2)
        auto dV = geom.CellSize(0) * geom.CellSize(1);
        #elif (AMREX_SPACEDIM == 3)
        auto dV = geom.CellSize(0) * geom.CellSize(1) * geom.CellSize(2);
        #endif

        // Loop over cells
        amrex::ParallelFor( n_cells,
            [=] AMREX_GPU_DEVICE (int i_cell) noexcept
            {
                // The particles from species1 that are in the cell `i_cell` are
                // given by the `indices_1[cell_start_1:cell_stop_1]`
                index_type const cell_start_1 = cell_offsets_1[i_cell];
                index_type const cell_stop_1  = cell_offsets_1[i_cell+1];
                // Same for species 2
                index_type const cell_start_2 = cell_offsets_2[i_cell];
                index_type const cell_stop_2  = cell_offsets_2[i_cell+1];

                // ux from species1 can be accessed like this:
                // ux_1[ indices_1[i] ], where i is between
                // cell_start_1 (inclusive) and cell_start_2 (exclusive)

                // Do not collide if one species is missing in the cell
                if ( cell_stop_1 - cell_start_1 >= 1 &&
                     cell_stop_2 - cell_start_2 >= 1 )
                {
                    // shuffle
                    ShuffleFisherYates(indices_1, cell_start_1, cell_stop_1);
                    ShuffleFisherYates(indices_2, cell_start_2, cell_stop_2);

                    // Call the function in order to perform collisions
                    ElasticCollisionPerez(
                        cell_start_1, cell_stop_1, cell_start_2, cell_stop_2,
                        indices_1, indices_2,
                        ux_1, uy_1, uz_1, ux_2, uy_2, uz_2, w_1, w_2,
                        q1, q2, m1, m2, Real(-1.0), Real(-1.0),
                        dt, CoulombLog, dV );
                }
            }
        );
    } // end if ( isSameSpecies)

}
