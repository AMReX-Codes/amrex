/**
 * \file AMReX_EB_StateRedistUtils.cpp
 */

#include <AMReX_EB_Redistribution.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_GpuMemory.H>

// Neighbor index lookup tables in constant memory for efficient GPU access
// These map neighbor indices to (i,j,k) offsets
#if (AMREX_SPACEDIM == 2)
namespace {
    AMREX_GPU_CONSTANT int const imap[9] = {0,-1, 0, 1,-1, 1,-1, 0, 1};
    AMREX_GPU_CONSTANT int const jmap[9] = {0,-1,-1,-1, 0, 0, 1, 1, 1};
    AMREX_GPU_CONSTANT int const kmap[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
}
#elif (AMREX_SPACEDIM == 3)
namespace {
    AMREX_GPU_CONSTANT int const imap[27] = {0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    AMREX_GPU_CONSTANT int const jmap[27] = {0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    AMREX_GPU_CONSTANT int const kmap[27] = {0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
}
#endif

namespace amrex {

void
MakeStateRedistUtils ( Box const& bx,
                       Array4<EBCellFlag const> const& flag,
                       Array4<Real const> const& vfrac,
                       Array4<Real const> const& ccent,
                       Array4<int  const> const& itracker,
                       Array4<Real      > const& nrs,
                       Array4<Real      > const& alpha,
                       Array4<Real      > const& nbhd_vol,
                       Array4<Real      > const& cent_hat,
                       Geometry const& lev_geom,
                       Real target_vol)
{
    // Note that itracker has {4 in 2D, 8 in 3D} components and all are initialized to zero
    // We will add to the first component every time this cell is included in a merged neighborhood,
    //    either by merging or being merged
    //
    // In 2D, we identify the cells in the remaining three components with the following ordering
    //
    // ^  6 7 8
    // |  4   5
    // j  1 2 3
    //   i --->
    //
    // In 3D, We identify the cells in the remaining three components with the following ordering
    //
    //    at k-1   |   at k  |   at k+1
    //
    // ^  15 16 17 |  6 7 8  |  24 25 26
    // |  12 13 14 |  4   5  |  21 22 23
    // j  9  10 11 |  1 2 3  |  18 19 20
    //   i --->
    //
    // Note the first component of each of these arrays should never be used
    // imap, jmap, kmap are now in constant memory (see file scope declarations above)
    //

    AMREX_D_TERM(const auto& is_periodic_x = lev_geom.isPeriodic(0);,
                 const auto& is_periodic_y = lev_geom.isPeriodic(1);,
                 const auto& is_periodic_z = lev_geom.isPeriodic(2););

    Box const& bxg3 = amrex::grow(bx,3);
    Box const& bxg4 = amrex::grow(bx,4);
    Box const& bxg5 = amrex::grow(bx,5);

    const Box domain = lev_geom.Domain();

    Box domain_per_grown = domain;
    if (is_periodic_x) { domain_per_grown.grow(0,5); }
    if (is_periodic_y) { domain_per_grown.grow(1,5); }
#if (AMREX_SPACEDIM == 3)
    if (is_periodic_z) { domain_per_grown.grow(2,5); }
#endif

    // Cells with a nonempty tracker are sparse for the usual EB workload. Append
    // their coordinates to a compact list while initializing nrs and alpha.
    Long const num_candidates = bxg5.numPts();
    Gpu::AsyncVector<Dim3> active_cells_v(num_candidates);
    Dim3* const active_cells = active_cells_v.data();
    Gpu::DeviceScalar<int> active_count_v(0);
    int* const active_count = active_count_v.dataPtr();

    amrex::ParallelFor(bxg5,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        nrs(i,j,k) = 1.0_rt;

        if (bxg4.contains(i,j,k)) {
            auto const cell_flag = flag(i,j,k);
            Real const value = cell_flag.isCovered() ? 0.0_rt : 1.0_rt;
            alpha(i,j,k,0) = value;

            if (bxg3.contains(i,j,k)) {
                if (cell_flag.isRegular()) {
                    AMREX_D_TERM(cent_hat(i,j,k,0) = 0.0_rt;,
                                 cent_hat(i,j,k,1) = 0.0_rt;,
                                 cent_hat(i,j,k,2) = 0.0_rt;);
                } else if (cell_flag.isCovered()) {
                    AMREX_D_TERM(cent_hat(i,j,k,0) = eb_covered_val;,
                                 cent_hat(i,j,k,1) = eb_covered_val;,
                                 cent_hat(i,j,k,2) = eb_covered_val;);
                } else {
                    AMREX_D_TERM(cent_hat(i,j,k,0) = ccent(i,j,k,0);,
                                 cent_hat(i,j,k,1) = ccent(i,j,k,1);,
                                 cent_hat(i,j,k,2) = ccent(i,j,k,2););
                }
            }
        }

        if (itracker(i,j,k,0) > 0) {
            int const active_index = Gpu::Atomic::Add(active_count, 1);
            active_cells[active_index] = Dim3{i,j,k};
        }
    });

    int const num_active_cells = active_count_v.dataValue();

    // Need nrs in bxg5 in order to compute alpha in bxg4
    // nrs captures how many neighborhoods (r,s,t) is in
    // Need alpha in bxg4 in order to define cent_hat in bxg3
    // (which allows us to compute slopes in bxg1 because some slopes go out by 2)

    // Scatter neighborhood counts and construct local alpha only for active cells.
    amrex::ParallelFor(num_active_cells,
    [=] AMREX_GPU_DEVICE (int active_index) noexcept
    {
        Dim3 const cell = active_cells[active_index];
        int const i = cell.x;
        int const j = cell.y;
        int const k = cell.z;

        // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
        for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
        {
            int const r = i+imap[itracker(i,j,k,i_nbor)];
            int const s = j+jmap[itracker(i,j,k,i_nbor)];
            int const t = k+kmap[itracker(i,j,k,i_nbor)];
            if (domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))) &&
                bxg5.contains(IntVect(AMREX_D_DECL(r,s,t))))
            {
                amrex::Gpu::Atomic::Add(&nrs(r,s,t),1.0_rt);
            }
        }

        if (bxg4.contains(IntVect(AMREX_D_DECL(i,j,k))) &&
            !flag(i,j,k).isCovered())
        {
            Real vol_of_nbors = 0.;
            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int const r = i+imap[itracker(i,j,k,i_nbor)];
                int const s = j+jmap[itracker(i,j,k,i_nbor)];
                int const t = k+kmap[itracker(i,j,k,i_nbor)];
                vol_of_nbors += vfrac(r,s,t);
            }
            alpha(i,j,k,1) = (target_vol - vfrac(i,j,k)) / vol_of_nbors;
        }
    });
    //

    // ****************************************************************************************
    // DEFINING ALPHA
    // ****************************************************************************************
    //
    // Must loop over (i,j,k) in bxg5 in order to define alpha(r,s,t) in bxg4
    //
    amrex::ParallelFor(num_active_cells,
    [=] AMREX_GPU_DEVICE (int active_index) noexcept
    {
        Dim3 const cell = active_cells[active_index];
        int const i = cell.x;
        int const j = cell.y;
        int const k = cell.z;

        if (!flag(i,j,k).isCovered())
        {
            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];
                if ( bxg4.contains(IntVect(AMREX_D_DECL(r,s,t))) ) {
                    amrex::Gpu::Atomic::Add(&alpha(r,s,t,0),-(alpha(i,j,k,1)/nrs(r,s,t)));
                }
            }
        }
    });

    // Tracker-zero cells do not need nbhd_vol: their normalized neighborhood
    // state reduces to U_in. Compute volumes and centroid corrections only for
    // active cells.
    amrex::ParallelFor(num_active_cells,
    [=] AMREX_GPU_DEVICE (int active_index) noexcept
    {
        Dim3 const cell = active_cells[active_index];
        int const i = cell.x;
        int const j = cell.y;
        int const k = cell.z;

        if (bxg3.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
            bool const is_uncovered = !flag(i,j,k).isCovered();
            Real const cell_vfrac = vfrac(i,j,k);
            bool const weight_centroid = cell_vfrac > 0.0 &&
                domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k)));

            Real local_nbhd_vol = is_uncovered ? alpha(i,j,k,0) * cell_vfrac : 0.0;

            AMREX_D_TERM(Real local_cent_hat_x = 0.0;,
                         Real local_cent_hat_y = 0.0;,
                         Real local_cent_hat_z = 0.0;);

            if (weight_centroid) {
                Real const self_weight = alpha(i,j,k,0) * cell_vfrac;
                AMREX_D_TERM(local_cent_hat_x = ccent(i,j,k,0) * self_weight;,
                             local_cent_hat_y = ccent(i,j,k,1) * self_weight;,
                             local_cent_hat_z = ccent(i,j,k,2) * self_weight;);
            }

            if (is_uncovered || weight_centroid) {
                // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
                for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
                {
                    int const ii = imap[itracker(i,j,k,i_nbor)]; int const r = i+ii;
                    int const jj = jmap[itracker(i,j,k,i_nbor)]; int const s = j+jj;
                    int const kk = kmap[itracker(i,j,k,i_nbor)]; int const t = k+kk;
                    Real const weight = alpha(i,j,k,1) * vfrac(r,s,t) / nrs(r,s,t);

                    if (is_uncovered) {
                        local_nbhd_vol += weight;
                    }

                    if (weight_centroid) {
                        AMREX_D_TERM(local_cent_hat_x += (ccent(r,s,t,0) + Real(ii)) * weight;,
                                     local_cent_hat_y += (ccent(r,s,t,1) + Real(jj)) * weight;,
                                     local_cent_hat_z += (ccent(r,s,t,2) + Real(kk)) * weight;);
                    }
                }
            }

            if (is_uncovered) {
                nbhd_vol(i,j,k) = local_nbhd_vol;
            }

            if (weight_centroid) {
                AMREX_D_TERM(cent_hat(i,j,k,0) = local_cent_hat_x / local_nbhd_vol;,
                             cent_hat(i,j,k,1) = local_cent_hat_y / local_nbhd_vol;,
                             cent_hat(i,j,k,2) = local_cent_hat_z / local_nbhd_vol;);
            }
        }
    });
}

}
