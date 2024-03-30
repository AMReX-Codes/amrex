/**
 * \file AMReX_EB_RedistributionApply.cpp
 * @{
 *
 */

#include <AMReX_BCRec.H>
#include <AMReX_EB_Redistribution.H>
#include <AMReX_EB_utils.H>

namespace amrex {

void ApplyRedistribution ( Box const& bx, int ncomp,
                           Array4<Real      > const& dUdt_out,
                           Array4<Real      > const& dUdt_in,
                           Array4<Real const> const& U_in,
                           Array4<Real> const& scratch,
                           Array4<EBCellFlag const> const& flag,
                           AMREX_D_DECL(Array4<Real const> const& apx,
                                        Array4<Real const> const& apy,
                                        Array4<Real const> const& apz),
                           Array4<amrex::Real const> const& vfrac,
                           AMREX_D_DECL(Array4<Real const> const& fcx,
                                        Array4<Real const> const& fcy,
                                        Array4<Real const> const& fcz),
                           Array4<Real const> const& ccc,
                           amrex::BCRec  const* d_bcrec_ptr,
                           Geometry const& lev_geom, Real dt,
                           std::string const& redistribution_type,
                           bool use_wts_in_divnc,
                           int srd_max_order,
                           amrex::Real target_volfrac,
                           Array4<Real const> const& srd_update_scale)
{
    int as_crse = 0;
    int as_fine = 0;

    // These values are dummies because they won't be used if as_crse == as_fine == 0
    int level_mask_not_covered = -1;
    Real fac_for_deltaR        =  1.0;
    int icomp                  = 0;

    ApplyMLRedistribution (bx, ncomp, dUdt_out, dUdt_in, U_in, scratch, flag,
                           AMREX_D_DECL(apx,apy,apz), vfrac,
                           AMREX_D_DECL(fcx,fcy,fcz), ccc,
                           d_bcrec_ptr, lev_geom, dt,
                           redistribution_type,
                           as_crse, Array4<Real>(), Array4<int const>(),
                           as_fine, Array4<Real>(), Array4<int const>(),
                           level_mask_not_covered,
                           fac_for_deltaR, use_wts_in_divnc, icomp,
                           srd_max_order, target_volfrac, srd_update_scale);
}

void
ApplyMLRedistribution ( Box const& bx, int ncomp,
                        Array4<Real      > const& dUdt_out,
                        Array4<Real      > const& dUdt_in,
                        Array4<Real const> const& U_in,
                        Array4<Real      > const& scratch,
                        Array4<EBCellFlag const> const& flag,
                        AMREX_D_DECL(Array4<Real const> const& apx,
                                     Array4<Real const> const& apy,
                                     Array4<Real const> const& apz),
                        Array4<amrex::Real const> const& vfrac,
                        AMREX_D_DECL(Array4<Real const> const& fcx,
                                     Array4<Real const> const& fcy,
                                     Array4<Real const> const& fcz),
                        Array4<Real const> const& ccc,
                        amrex::BCRec  const* d_bcrec_ptr,
                        Geometry const& lev_geom, Real dt,
                        std::string const& redistribution_type,
                        int as_crse,
                        Array4<Real            > const& rr_drho_crse,
                        Array4<int        const> const& rr_flag_crse,
                        int as_fine,
                        Array4<Real            > const& dm_as_fine,
                        Array4<int        const> const& levmsk,
                        int level_mask_not_covered,
                        Real fac_for_deltaR,
                        bool use_wts_in_divnc,
                        int icomp,
                        int srd_max_order,
                        amrex::Real target_volfrac,
                        Array4<Real const> const& srd_update_scale)
{
    // redistribution_type = "NoRedist";       // no redistribution
    // redistribution_type = "FluxRedist"      // flux_redistribute
    // redistribution_type = "StateRedist";    // (weighted) state redistribute

    // amrex::Print() <<" Redistribution::ApplyML " << redistribution_type << '\n';

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            dUdt_out(i,j,k,n) = 0.;
        });

    if (redistribution_type == "FluxRedist")
    {
        // This is the re-redistribution version of FluxRedistribution
        amrex_flux_redistribute (bx, dUdt_out, dUdt_in, scratch, vfrac, flag,
                                 as_crse, rr_drho_crse, rr_flag_crse,
                                 as_fine, dm_as_fine, levmsk, lev_geom,
                                 use_wts_in_divnc, level_mask_not_covered, icomp, ncomp,
                                 fac_for_deltaR*dt);

    } else if (redistribution_type == "StateRedist") {

        Box const& bxg1 = grow(bx,1);
        Box const& bxg3 = grow(bx,3);
        Box const& bxg4 = grow(bx,4);
        Box const& bxg5 = grow(bx,5);

#if (AMREX_SPACEDIM == 2)
        // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors

        // itracker(i,j,n) holds the identifier for (r,s), the nth neighbor of (i,j)
        IArrayBox itracker(bxg5,4,The_Async_Arena());
#else
        // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors

        // itracker(i,j,k,n) holds the identifier for (r,s,t), the nth neighbor of (i,j,k)
        IArrayBox itracker(bxg5,8,The_Async_Arena());
#endif
        FArrayBox nrs_fab(bxg5,1,The_Async_Arena());
        FArrayBox alpha_fab(bxg4,2,The_Async_Arena());

        // Total volume of all cells in my nbhd
        FArrayBox nbhd_vol_fab(bxg3,1,The_Async_Arena());

        // Centroid of my nbhd
        FArrayBox cent_hat_fab(bxg3,AMREX_SPACEDIM,The_Async_Arena());

        Array4<int> itr = itracker.array();
        Array4<int const> itr_const = itracker.const_array();

        Array4<Real      > nrs       = nrs_fab.array();
        Array4<Real const> nrs_const = nrs_fab.const_array();

        Array4<Real      > alpha       = alpha_fab.array();
        Array4<Real const> alpha_const = alpha_fab.const_array();

        Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
        Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

        Array4<Real      > cent_hat       = cent_hat_fab.array();
        Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

        Box domain_per_grown = lev_geom.Domain();
        AMREX_D_TERM(if (lev_geom.isPeriodic(0)) { domain_per_grown.grow(0,1); },
                     if (lev_geom.isPeriodic(1)) { domain_per_grown.grow(1,1); },
                     if (lev_geom.isPeriodic(2)) { domain_per_grown.grow(2,1); })

        // At any external Dirichlet domain boundaries we need to set dUdt_in to 0
        //    in the cells just outside the domain because those values will be used
        //    in the slope computation in state redistribution.  We assume here that
        //    the ext_dir values of U_in itself have already been set.
        if (!domain_per_grown.contains(bxg1)) {
            amrex::ParallelFor(bxg1,ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!domain_per_grown.contains(IntVect(AMREX_D_DECL(i,j,k)))) {
                        dUdt_in(i,j,k,n) = 0.;
                    }
                });
        }

        amrex::ParallelFor(Box(scratch), ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                const Real scale = (srd_update_scale) ? srd_update_scale(i,j,k) : Real(1.0);
                scratch(i,j,k,n) = U_in(i,j,k,n) + dt * dUdt_in(i,j,k,n) / scale;
            }
        );

        MakeITracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, target_volfrac);

        MakeStateRedistUtils(bx, flag, vfrac, ccc, itr, nrs, alpha, nbhd_vol, cent_hat,
                             lev_geom, target_volfrac);

        MLStateRedistribute(bx, ncomp, dUdt_out, scratch, flag, vfrac,
                            AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                            itr_const, nrs_const, alpha_const, nbhd_vol_const,
                            cent_hat_const, lev_geom,
                            as_crse, rr_drho_crse, rr_flag_crse,
                            as_fine, dm_as_fine, levmsk,
                            level_mask_not_covered, fac_for_deltaR, srd_max_order);

        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                // Only update the values which actually changed -- this makes
                // the results insensitive to tiling -- otherwise cells that aren't
                // changed but are in a tile on which StateRedistribute gets called
                // will have precision-level changes due to adding/subtracting U_in
                // and multiplying/dividing by dt.   Here we test on whether (i,j,k)
                // has at least one neighbor and/or whether (i,j,k) is in the
                // neighborhood of another cell -- if either of those is true the
                // value may have changed

                if (itr(i,j,k,0) > 0 || nrs(i,j,k) > 1.)
                {
                   const Real scale = (srd_update_scale) ? srd_update_scale(i,j,k) : Real(1.0);

                   dUdt_out(i,j,k,n) = scale * (dUdt_out(i,j,k,n) - U_in(i,j,k,n)) / dt;

                }
                else
                {
                   dUdt_out(i,j,k,n) = dUdt_in(i,j,k,n);
                }
            }
        );

    } else if (redistribution_type == "NoRedist") {
        amrex::ParallelFor(bx, ncomp,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                dUdt_out(i,j,k,n) = dUdt_in(i,j,k,n);
            }
        );

    } else {
        amrex::Error("Not a legit redist_type in ApplyML");
    }
}

void
ApplyInitialRedistribution ( Box const& bx, int ncomp,
                             Array4<Real      > const& U_out,
                             Array4<Real      > const& U_in,
                             Array4<EBCellFlag const> const& flag,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                          amrex::Array4<amrex::Real const> const& apy,
                                          amrex::Array4<amrex::Real const> const& apz),
                             amrex::Array4<amrex::Real const> const& vfrac,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                          amrex::Array4<amrex::Real const> const& fcy,
                                          amrex::Array4<amrex::Real const> const& fcz),
                             amrex::Array4<amrex::Real const> const& ccc,
                             amrex::BCRec  const* d_bcrec_ptr,
                             Geometry const& lev_geom,
                             std::string const& redistribution_type,
                             int srd_max_order,
                             amrex::Real target_volfrac)
{
    if (redistribution_type != "StateRedist") {
    std::string msg = "ApplyInitialRedistribution: Shouldn't be here with redist type "+redistribution_type;
        amrex::Error(msg);
    }

    // amrex::Print() <<" Redistribution::ApplyInitial " << redistribution_type << '\n';

    Box const& bxg3 = grow(bx,3);
    Box const& bxg4 = grow(bx,4);
    Box const& bxg5 = grow(bx,5);

#if (AMREX_SPACEDIM == 2)
        // We assume that in 2D a cell will only need at most 3 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors

        // itracker(i,j,n) holds the identifier for (r,s), the nth neighbor of (i,j)
        IArrayBox itracker(bxg5,4,The_Async_Arena());
#else
        // We assume that in 3D a cell will only need at most 7 neighbors to merge with, and we
        //    use the first component of this for the number of neighbors

        // itracker(i,j,k,n) holds the identifier for (r,s,t), the nth neighbor of (i,j,k)
        IArrayBox itracker(bxg5,8,The_Async_Arena());
#endif
    FArrayBox nrs_fab(bxg5,1,The_Async_Arena());
    FArrayBox alpha_fab(bxg4,2,The_Async_Arena());

    // Total volume of all cells in my nbhd
    FArrayBox nbhd_vol_fab(bxg3,1,The_Async_Arena());

    // Centroid of my nbhd
    FArrayBox cent_hat_fab(bxg3,AMREX_SPACEDIM,The_Async_Arena());

    Array4<int> itr = itracker.array();
    Array4<int const> itr_const = itracker.const_array();

    Array4<Real      > nrs       = nrs_fab.array();
    Array4<Real const> nrs_const = nrs_fab.const_array();

    Array4<Real      > alpha       = alpha_fab.array();
    Array4<Real const> alpha_const = alpha_fab.const_array();

    Array4<Real      > nbhd_vol       = nbhd_vol_fab.array();
    Array4<Real const> nbhd_vol_const = nbhd_vol_fab.const_array();

    Array4<Real      > cent_hat       = cent_hat_fab.array();
    Array4<Real const> cent_hat_const = cent_hat_fab.const_array();

    amrex::ParallelFor(bx,ncomp,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {
        U_out(i,j,k,n) = 0.;
    });

    MakeITracker(bx, AMREX_D_DECL(apx, apy, apz), vfrac, itr, lev_geom, target_volfrac);

    MakeStateRedistUtils(bx, flag, vfrac, ccc, itr, nrs, alpha, nbhd_vol, cent_hat,
                            lev_geom, target_volfrac);

    StateRedistribute(bx, ncomp, U_out, U_in, flag, vfrac,
                         AMREX_D_DECL(fcx, fcy, fcz), ccc,  d_bcrec_ptr,
                      itr_const, nrs_const, alpha_const, nbhd_vol_const,
                      cent_hat_const, lev_geom, srd_max_order);
}

}
