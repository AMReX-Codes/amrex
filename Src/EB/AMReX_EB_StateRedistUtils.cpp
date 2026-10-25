/**
 * \file AMReX_EB_StateRedistUtils.cpp
 */

#include <AMReX_EB_Redistribution.H>
#include <AMReX_EB_StateRedistTracker_K.H>

namespace amrex {

using state_redist_detail::imap;
using state_redist_detail::jmap;
using state_redist_detail::kmap;

void
MakeStateRedistUtils ( Box const& bx,
                       Array4<EBCellFlag const> const& flag,
                       Array4<Real const> const& vfrac,
                       Array4<Real const> const& ccent,
                       Array4<int const  > const& itracker,
                       Array4<int       > const& nrs,
                       Array4<Real      > const& alpha,
                       Array4<Real      > const& nbhd_vol,
                       Array4<Real      > const& cent_hat,
                       Dim3 const* active_cells,
                       int num_active_cells,
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

    // Ghost extents propagate nrs through alpha and cent_hat to slopes in bxg1.

    // Scatter neighborhood counts and construct local alpha only for active cells.
    amrex::ParallelFor(num_active_cells,
    [=] AMREX_GPU_DEVICE (int active_index) noexcept
    {
        Dim3 const cell = active_cells[active_index];
        int const i = cell.x;
        int const j = cell.y;
        int const k = cell.z;

        for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
        {
            int const r = i+imap[itracker(i,j,k,i_nbor)];
            int const s = j+jmap[itracker(i,j,k,i_nbor)];
            int const t = k+kmap[itracker(i,j,k,i_nbor)];
            if (domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))) &&
                bxg5.contains(IntVect(AMREX_D_DECL(r,s,t))))
            {
                amrex::Gpu::Atomic::Add(&nrs(r,s,t),1);
            }
        }

        if (bxg4.contains(IntVect(AMREX_D_DECL(i,j,k))) &&
            !flag(i,j,k).isCovered())
        {
            Real vol_of_nbors = 0.;
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

    // Accumulate active-cell alpha corrections in neighboring cells.
    amrex::ParallelFor(num_active_cells,
    [=] AMREX_GPU_DEVICE (int active_index) noexcept
    {
        Dim3 const cell = active_cells[active_index];
        int const i = cell.x;
        int const j = cell.y;
        int const k = cell.z;

        if (!flag(i,j,k).isCovered())
        {
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];
                if ( bxg4.contains(IntVect(AMREX_D_DECL(r,s,t))) ) {
                    amrex::Gpu::Atomic::Add(&alpha(r,s,t,0),
                                           -(alpha(i,j,k,1) / static_cast<Real>(nrs(r,s,t))));
                }
            }
        }
    });

    // MakeITracker appends only positive-volume cells in domain_per_grown, so
    // both the volume and centroid accumulations are unconditional here.
    amrex::ParallelFor(num_active_cells,
    [=] AMREX_GPU_DEVICE (int active_index) noexcept
    {
        Dim3 const cell = active_cells[active_index];
        int const i = cell.x;
        int const j = cell.y;
        int const k = cell.z;

        if (bxg3.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
            Real const cell_vfrac = vfrac(i,j,k);
            Real local_nbhd_vol = alpha(i,j,k,0) * cell_vfrac;

            Real const self_weight = local_nbhd_vol;
            AMREX_D_TERM(Real local_cent_hat_x = ccent(i,j,k,0) * self_weight;,
                         Real local_cent_hat_y = ccent(i,j,k,1) * self_weight;,
                         Real local_cent_hat_z = ccent(i,j,k,2) * self_weight;);

            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int const ii = imap[itracker(i,j,k,i_nbor)]; int const r = i+ii;
                int const jj = jmap[itracker(i,j,k,i_nbor)]; int const s = j+jj;
                int const kk = kmap[itracker(i,j,k,i_nbor)]; int const t = k+kk;
                Real const weight = alpha(i,j,k,1) * vfrac(r,s,t) /
                    static_cast<Real>(nrs(r,s,t));

                local_nbhd_vol += weight;
                AMREX_D_TERM(local_cent_hat_x += (ccent(r,s,t,0) + Real(ii)) * weight;,
                             local_cent_hat_y += (ccent(r,s,t,1) + Real(jj)) * weight;,
                             local_cent_hat_z += (ccent(r,s,t,2) + Real(kk)) * weight;);
            }

            nbhd_vol(i,j,k) = local_nbhd_vol;
            AMREX_D_TERM(cent_hat(i,j,k,0) = local_cent_hat_x / local_nbhd_vol;,
                         cent_hat(i,j,k,1) = local_cent_hat_y / local_nbhd_vol;,
                         cent_hat(i,j,k,2) = local_cent_hat_z / local_nbhd_vol;);
        }
    });
}

}
