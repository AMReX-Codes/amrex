/**
 * \file AMReX_StateRedistUtils.cpp
 */

#include <AMReX_EB_Redistribution.H>

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
    //
#if (AMREX_SPACEDIM == 2)
    Array<int,9> imap{0,-1, 0, 1,-1, 1,-1, 0, 1};
    Array<int,9> jmap{0,-1,-1,-1, 0, 0, 1, 1, 1};
    Array<int,9> kmap{0, 0, 0, 0, 0, 0, 0, 0, 0};
#else
    Array<int,27>    imap{0,-1, 0, 1,-1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
    Array<int,27>    jmap{0,-1,-1,-1, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1,-1,-1,-1, 0, 0, 0, 1, 1, 1};
    Array<int,27>    kmap{0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#endif

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

    // ****************************************************************************************
    // DEFINING NRS
    // ****************************************************************************************
    //
    // Need nrs in bxg5 in order to compute alpha in bxg4
    // nrs captures how many neighborhoods (r,s,t) is in
    //
    amrex::ParallelFor(bxg5,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // Everyone is in their own neighborhood at least
        nrs(i,j,k) = 1.;
    });
    amrex::ParallelFor(bxg5,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
        for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
        {
            int r = i+imap[itracker(i,j,k,i_nbor)];
            int s = j+jmap[itracker(i,j,k,i_nbor)];
            int t = k+kmap[itracker(i,j,k,i_nbor)];
            if ( domain_per_grown.contains(IntVect(AMREX_D_DECL(r,s,t))) &&
                 bxg5.contains(IntVect(AMREX_D_DECL(r,s,t))) )
            {
                amrex::Gpu::Atomic::Add(&nrs(r,s,t),1.0_rt);
            }
        }
    });
    //

    // ****************************************************************************************
    // DEFINING ALPHA
    // ****************************************************************************************
    //
    // Need alpha in bxg4 in order to define cent_hat in bxg3
    // (which allows us to compute slopes in bxg1 because some slopes go out by 2)
    //
    amrex::ParallelFor(bxg4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        alpha(i,j,k,0) = 1.;
        alpha(i,j,k,1) = 1.;
    });
    amrex::ParallelFor(bxg4,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (!flag(i,j,k).isCovered())
        {
            Real vol_of_nbors = 0.;

            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];
                vol_of_nbors += vfrac(r,s,t);
            }

            if (itracker(i,j,k,0) > 0) {
                alpha(i,j,k,1) = (target_vol - vfrac(i,j,k)) / vol_of_nbors;
            }

        } else {
            alpha(i,j,k,0) = 0.;
            alpha(i,j,k,1) = 0.;
        }
    });

    //
    // Must loop over (i,j,k) in bxg5 in order to define alpha(r,s,t) in bxg4
    //
    amrex::ParallelFor(bxg5,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
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

    // ****************************************************************************************
    // DEFINING NBHD_VOL
    // ****************************************************************************************
    //
    // Redefine nbhd_vol in bxg3 in order to use it to define cent_hat in bxg3
    // (which allows us to compute slopes in bxg1)
    // To define nbhd_vol in bxg3 we also need alpha(i,j,k,0/1) in bxg3,
    //    and vfrac and nrs in bxg4
    //
    amrex::ParallelFor(bxg3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (!flag(i,j,k).isCovered())
        {
            nbhd_vol(i,j,k)  = alpha(i,j,k,0) * vfrac(i,j,k);

            // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
            for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
            {
                int r = i+imap[itracker(i,j,k,i_nbor)];
                int s = j+jmap[itracker(i,j,k,i_nbor)];
                int t = k+kmap[itracker(i,j,k,i_nbor)];
                amrex::Gpu::Atomic::Add(&nbhd_vol(i,j,k),alpha(i,j,k,1) * vfrac(r,s,t) / nrs(r,s,t));
            }
        } else {
            nbhd_vol(i,j,k) = 0.;
        }
    });

    // ****************************************************************************************
    // DEFINING CENT_HAT
    // ****************************************************************************************
    //
    // Need cent_hat(xhat,yhat,zhat) in bxg3 to compute slopes in bxg1
    //
    amrex::ParallelFor(bxg3,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        if (vfrac(i,j,k) > 0.0)
        {
            AMREX_D_TERM(cent_hat(i,j,k,0) = ccent(i,j,k,0);,
                         cent_hat(i,j,k,1) = ccent(i,j,k,1);,
                         cent_hat(i,j,k,2) = ccent(i,j,k,2););

            if ( itracker(i,j,k,0) > 0 &&
                 domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k))) )
            {
                AMREX_D_TERM(cent_hat(i,j,k,0) = ccent(i,j,k,0) * alpha(i,j,k,0) *vfrac(i,j,k);,
                             cent_hat(i,j,k,1) = ccent(i,j,k,1) * alpha(i,j,k,0) *vfrac(i,j,k);,
                             cent_hat(i,j,k,2) = ccent(i,j,k,2) * alpha(i,j,k,0) *vfrac(i,j,k););

                // This loops over the neighbors of (i,j,k), and doesn't include (i,j,k) itself
                for (int i_nbor = 1; i_nbor <= itracker(i,j,k,0); i_nbor++)
                {
                    int ii = imap[itracker(i,j,k,i_nbor)]; int r = i+ii;
                    int jj = jmap[itracker(i,j,k,i_nbor)]; int s = j+jj;
                    int kk = kmap[itracker(i,j,k,i_nbor)]; int t = k+kk;

                    AMREX_D_TERM(cent_hat(i,j,k,0) += (ccent(r,s,t,0) + ii) * alpha(i,j,k,1) * vfrac(r,s,t) / nrs(r,s,t);,
                                 cent_hat(i,j,k,1) += (ccent(r,s,t,1) + jj) * alpha(i,j,k,1) * vfrac(r,s,t) / nrs(r,s,t);,
                                 cent_hat(i,j,k,2) += (ccent(r,s,t,2) + kk) * alpha(i,j,k,1) * vfrac(r,s,t) / nrs(r,s,t););
                }

                AMREX_D_TERM(cent_hat(i,j,k,0) /= nbhd_vol(i,j,k);,
                             cent_hat(i,j,k,1) /= nbhd_vol(i,j,k);,
                             cent_hat(i,j,k,2) /= nbhd_vol(i,j,k););
            }
        } else {

                AMREX_D_TERM(cent_hat(i,j,k,0) = eb_covered_val;,
                             cent_hat(i,j,k,1) = eb_covered_val;,
                             cent_hat(i,j,k,2) = eb_covered_val;);
        }
    });
}

}
