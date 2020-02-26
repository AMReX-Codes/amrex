/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "SortingUtils.H"
#include "Particles/PhysicalParticleContainer.H"
#include "WarpX.H"

#include <AMReX_Particles.H>


using namespace amrex;

/* \brief Determine which particles deposit/gather in the buffer, and
 *        and reorder the particle arrays accordingly
 *
 *  More specifically:
 *  - Modify `nfine_current` and `nfine_gather` (in place)
 *     so that they correspond to the number of particles
 *     that deposit/gather in the fine patch respectively.
 *  - Reorder the particle arrays,
 *     so that the `nfine_current`/`nfine_gather` first particles
 *     deposit/gather in the fine patch
 *     and (thus) the `np-nfine_current`/`np-nfine_gather` last particles
 *     deposit/gather in the buffer
 *
 * \param nfine_current number of particles that deposit to the fine patch
 *         (modified by this function)
 * \param nfine_gather number of particles that gather into the fine patch
 *         (modified by this function)
 * \param np total number of particles in this tile
 * \param pti object that holds the particle information for this tile
 * \param lev current refinement level
 * \param current_masks indicates, for each cell, whether that cell is
 *       in the deposition buffers or in the interior of the fine patch
 * \param gather_masks indicates, for each cell, whether that cell is
 *       in the gather buffers or in the interior of the fine patch
 * \param uxp, uyp, uzp, wp references to the particle momenta and weight
 *         (modified by this function)
 */
void
PhysicalParticleContainer::PartitionParticlesInBuffers(
    long& nfine_current, long& nfine_gather, long const np,
    WarpXParIter& pti, int const lev,
    iMultiFab const* current_masks,
    iMultiFab const* gather_masks,
    RealVector& uxp, RealVector& uyp, RealVector& uzp, RealVector& wp)
{
    WARPX_PROFILE("PPC::Evolve::partition");

    // Initialize temporary arrays
    Gpu::DeviceVector<int> inexflag;
    inexflag.resize(np);
    Gpu::DeviceVector<long> pid;
    pid.resize(np);

    // First, partition particles into the larger buffer

    // - Select the larger buffer
    iMultiFab const* bmasks =
        (WarpX::n_field_gather_buffer >= WarpX::n_current_deposition_buffer) ?
        gather_masks : current_masks;
    // - For each particle, find whether it is in the larger buffer,
    //   by looking up the mask. Store the answer in `inexflag`.
    amrex::ParallelFor( np, fillBufferFlag(pti, bmasks, inexflag, Geom(lev)) );
    // - Find the indices that reorder particles so that the last particles
    //   are in the larger buffer
    fillWithConsecutiveIntegers( pid );
    auto const sep = stablePartition( pid.begin(), pid.end(), inexflag );
    // At the end of this step, `pid` contains the indices that should be used to
    // reorder the particles, and `sep` is the position in the array that
    // separates the particles that deposit/gather on the fine patch (first part)
    // and the particles that deposit/gather in the buffers (last part)
    long const n_fine = iteratorDistance(pid.begin(), sep);
    // Number of particles on fine patch, i.e. outside of the larger buffer

    // Second, among particles that are in the larger buffer, partition
    // particles into the smaller buffer

    if (WarpX::n_current_deposition_buffer == WarpX::n_field_gather_buffer) {
        // No need to do anything if the buffers have the same size
        nfine_current = nfine_gather = iteratorDistance(pid.begin(), sep);
    } else if (sep == pid.end()) {
        // No need to do anything if there are no particles in the larger buffer
        nfine_current = nfine_gather = np;
    } else {
        int n_buf;
        if (bmasks == gather_masks) {
            nfine_gather = n_fine;
            bmasks = current_masks;
            n_buf = WarpX::n_current_deposition_buffer;
        } else {
            nfine_current = n_fine;
            bmasks = gather_masks;
            n_buf = WarpX::n_field_gather_buffer;
        }
        if (n_buf > 0)
        {
            // - For each particle in the large buffer, find whether it is in
            // the smaller buffer, by looking up the mask. Store the answer in `inexflag`.
            amrex::ParallelFor( np - n_fine,
               fillBufferFlagRemainingParticles(pti, bmasks, inexflag, Geom(lev), pid, n_fine) );
            auto const sep2 = stablePartition( sep, pid.end(), inexflag );

            if (bmasks == gather_masks) {
                nfine_gather = iteratorDistance(pid.begin(), sep2);
            } else {
                nfine_current = iteratorDistance(pid.begin(), sep2);
            }
        }
    }

    // only deposit / gather to coarsest grid
    if (m_deposit_on_main_grid && lev > 0) {
        nfine_current = 0;
    }
    if (m_gather_from_main_grid && lev > 0) {
        nfine_gather = 0;
    }

    // Reorder the actual particle array, using the `pid` indices
    if (nfine_current != np || nfine_gather != np)
    {
        // Temporary array for particle AoS
        ParticleVector particle_tmp;
        particle_tmp.resize(np);

        // Copy particle AoS
        auto& aos = pti.GetArrayOfStructs();
        amrex::ParallelFor( np,
            copyAndReorder<ParticleType>( aos(), particle_tmp, pid ) );
        std::swap(aos(), particle_tmp);

        // Temporary array for particle individual attributes
        RealVector tmp;
        tmp.resize(np);

        // Copy individual attributes
        amrex::ParallelFor( np, copyAndReorder<Real>( wp, tmp, pid ) );
        std::swap(wp, tmp);
        amrex::ParallelFor( np, copyAndReorder<Real>( uxp, tmp, pid ) );
        std::swap(uxp, tmp);
        amrex::ParallelFor( np, copyAndReorder<Real>( uyp, tmp, pid ) );
        std::swap(uyp, tmp);
        amrex::ParallelFor( np, copyAndReorder<Real>( uzp, tmp, pid ) );
        std::swap(uzp, tmp);

        // Make sure that the temporary arrays are not destroyed before
        // the GPU kernels finish running
        Gpu::streamSynchronize();
    }
    // Make sure that the temporary arrays are not destroyed before
    // the GPU kernels finish running
    Gpu::streamSynchronize();
}
