#include <WarpX.H>
#include "CollisionType.H"
#include "ElasticCollisionPerez.H"

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
    findParticlesInEachCell( int lev, MFIter const& mfi,
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

        // Find particles that are in each cell ;
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


/* \brief Perform Coulomb collisions within one particle tile */
void CollisionType::doCoulombCollisionsWithinTile
    ( int lev, MFIter const& mfi,
    std::unique_ptr<WarpXParticleContainer>& species_1,
    std::unique_ptr<WarpXParticleContainer>& species_2 )
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
    auto& soa_1= ptile_1.GetStructOfArrays();
    Real* ux_1 = soa_1.GetRealData(PIdx::ux).data();
    Real* uy_1 = soa_1.GetRealData(PIdx::ux).data();
    Real* uz_1 = soa_1.GetRealData(PIdx::ux).data();
    Real* w_1 = soa_1.GetRealData(PIdx::w).data();
    index_type* indices_1 = bins_1.permutationPtr();
    index_type const* cell_offsets_1 = bins_1.offsetsPtr();
    Real q1 = species_1->getCharge();
    Real m1 = species_1->getMass();
    // - Species 2
    auto& soa_2= ptile_2.GetStructOfArrays();
    Real* ux_2 = soa_2.GetRealData(PIdx::ux).data();
    Real* uy_2 = soa_2.GetRealData(PIdx::ux).data();
    Real* uz_2 = soa_2.GetRealData(PIdx::ux).data();
    Real* w_2 = soa_2.GetRealData(PIdx::w).data();
    index_type* indices_2 = bins_2.permutationPtr();
    index_type const* cell_offsets_2 = bins_2.offsetsPtr();
    Real q2 = species_2->getCharge();
    Real m2 = species_2->getMass();

    const Real dt = WarpX::GetInstance().getdt(lev);
    Geometry const& geom = WarpX::GetInstance().Geom(lev);
    const Real dV = geom.CellSize(0)*geom.CellSize(1)*geom.CellSize(2);

    // Loop over cells
    amrex::ParallelFor( n_cells,
        [=] AMREX_GPU_DEVICE (int i_cell) noexcept
        {
            // The particles from species1 that are in the cell `i_cell` are
            // given by the `indices_1[cell_start_1:cell_stop_1]`
            index_type const cell_start_1 = cell_offsets_1[i_cell];
            index_type const cell_stop_1 = cell_offsets_1[i_cell+1];
            // Same for species 2
            index_type const cell_start_2 = cell_offsets_2[i_cell];
            index_type const cell_stop_2 = cell_offsets_2[i_cell+1];

            // ux from species1 can be accessed like this:
            // ux_1[ indices_1[i] ], where i is between
            // cell_start_1 (inclusive) and cell_start_2 (exclusive)

            // Call the function in order to perform collisions

            ElasticCollisionPerez(
                cell_start_1, cell_stop_1, cell_start_2, cell_stop_2,
                indices_1, indices_2,
                ux_1, uy_1, uz_1, ux_2, uy_2, uz_2, w_1, w_2,
                q1, q2, m1, m2, -1.0, -1.0, dt, -1.0, dV);

        }
    );

}

CollisionType::CollisionType(
    const std::vector<std::string>& species_names,
    std::string collision_name)
{

    std::vector<std::string> collision_species;

    amrex::ParmParse pp(collision_name);
    pp.getarr("species", collision_species);

    for (int i=0; i<species_names.size(); i++)
    {
        if (species_names[i] == collision_species[0])
	{ m_species1 = i; }
        if (species_names[i] == collision_species[1])
	{ m_species2 = i; }
    }

}


