#include <SortingUtils.H>
#include <PhysicalParticleContainer.H>
#include <WarpX.H>
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
    BL_PROFILE("PPC::Evolve::partition");
    Gpu::ManagedDeviceVector<int> inexflag;
    inexflag.resize(np);

    // Select the largest buffer first
    iMultiFab const* bmasks =
        (WarpX::n_field_gather_buffer >= WarpX::n_current_deposition_buffer) ?
        gather_masks : current_masks;

    // For each particle, find whether it is in the large buffer, by looking up the mask
    const Array4<const int>& msk = (*bmasks)[pti].array();
    auto& aos = pti.GetArrayOfStructs();
    ParticleType * AMREX_RESTRICT pstructs = &(aos[0]);
    int * inexflag_ptr = inexflag.dataPtr();
    const Geometry& geom = Geom(lev);
    const Real prob_lo_x = geom.ProbLo(0);
    const Real prob_lo_y = geom.ProbLo(1);
    const Real prob_lo_z = geom.ProbLo(2);
    const Real inv_dx = geom.InvCellSize(0);
    const Real inv_dy = geom.InvCellSize(1);
    const Real inv_dz = geom.InvCellSize(2);
    const IntVect domain_small_end = geom.Domain().smallEnd();
    amrex::ParallelFor( np,
         [=] AMREX_GPU_DEVICE (long i) {

             const IntVect& iv = findParticleIndex(pstructs[i],
                 prob_lo_x, prob_lo_y, prob_lo_z,
                 inv_dx, inv_dy, inv_dz, domain_small_end );
             inexflag_ptr[i] = msk(iv);
         }
     );

    // Partition the particles according whether they are in the large buffer or not
    Gpu::ManagedDeviceVector<long> pid;
    pid.resize(np);
    fillWithConsecutiveIntegers( pid );
    auto sep = std::stable_partition(pid.begin(), pid.end(),
                                     [&inexflag](long id) { return inexflag[id]; });

    if (WarpX::n_current_deposition_buffer == WarpX::n_field_gather_buffer) {
        nfine_current = nfine_gather = std::distance(pid.begin(), sep);
    } else if (sep != pid.end()) {
        int n_buf;
        if (bmasks == gather_masks) {
            nfine_gather = std::distance(pid.begin(), sep);
            bmasks = current_masks;
            n_buf = WarpX::n_current_deposition_buffer;
        } else {
            nfine_current = std::distance(pid.begin(), sep);
            bmasks = gather_masks;
            n_buf = WarpX::n_field_gather_buffer;
        }
        if (n_buf > 0)
        {
            const auto& msk2 = (*bmasks)[pti];
            for (auto it = sep; it != pid.end(); ++it) {
                const long id = *it;
                const IntVect& iv = Index(aos[id], lev);
                inexflag[id] = msk2(iv);
            }

            auto sep2 = std::stable_partition(sep, pid.end(),
                                              [&inexflag](long id) {return inexflag[id];});
            if (bmasks == gather_masks) {
                nfine_gather = std::distance(pid.begin(), sep2);
            } else {
                nfine_current = std::distance(pid.begin(), sep2);
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

    if (nfine_current != np || nfine_gather != np)
    {

        ParticleVector particle_tmp;
        particle_tmp.resize(np);

        for (long ip = 0; ip < np; ++ip) {
            particle_tmp[ip] = aos[pid[ip]];
        }
        std::swap(aos(), particle_tmp);

        RealVector tmp;
        tmp.resize(np);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = wp[pid[ip]];
        }
        std::swap(wp, tmp);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = uxp[pid[ip]];
        }
        std::swap(uxp, tmp);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = uyp[pid[ip]];
        }
        std::swap(uyp, tmp);

        for (long ip = 0; ip < np; ++ip) {
            tmp[ip] = uzp[pid[ip]];
        }
        std::swap(uzp, tmp);
    }

}
