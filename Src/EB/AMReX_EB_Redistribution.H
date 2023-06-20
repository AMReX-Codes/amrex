#ifndef AMREX_EB_REDIST_H_
#define AMREX_EB_REDIST_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB2.H>

#ifdef AMREX_USE_FLOAT
        constexpr amrex::Real eb_covered_val = amrex::Real(1.e20);
#else
        constexpr amrex::Real eb_covered_val = amrex::Real(1.e40);
#endif

namespace amrex {
#if (AMREX_SPACEDIM > 1)

    void single_level_redistribute (amrex::MultiFab& div_tmp_in, amrex::MultiFab& div_out,
                                    int div_comp, int ncomp, const amrex::Geometry& geom);

    void single_level_weighted_redistribute (amrex::MultiFab& div_tmp_in, amrex::MultiFab& div_out,
                                             const amrex::MultiFab& weights,
                                             int div_comp, int ncomp, const amrex::Geometry& geom);

#if (__cplusplus >= 201402L)
    [[deprecated]]
#endif
    inline void single_level_redistribute (int lev, amrex::MultiFab& div_tmp_in, amrex::MultiFab& div_out,
                                           int div_comp, int ncomp, const amrex::Vector<amrex::Geometry>& geom)
    { single_level_redistribute(div_tmp_in, div_out, div_comp, ncomp, geom[lev]); }

#if (__cplusplus >= 201402L)
    [[deprecated]]
#endif
    inline void single_level_weighted_redistribute (int lev, amrex::MultiFab& div_tmp_in, amrex::MultiFab& div_out,
                                                    const amrex::MultiFab& weights, int div_comp, int ncomp,
                                                    const amrex::Vector<amrex::Geometry>& geom)
    { single_level_weighted_redistribute(div_tmp_in, div_out, weights, div_comp, ncomp, geom[lev]); }

    // Flux redistribution scheme for case where C/F interface does *not* cross EB
    void apply_flux_redistribution ( const amrex::Box& bx,
                                     amrex::Array4<amrex::Real      > const& div,
                                     amrex::Array4<amrex::Real const> const& divc,
                                     amrex::Array4<amrex::Real const> const& wt,
                                     int icomp, int ncomp,
                                     amrex::Array4<amrex::EBCellFlag const> const& flags,
                                     amrex::Array4<amrex::Real       const> const& vfrac,
                                     const amrex::Geometry & geom);

    // Flux redistribution scheme for case where C/F interface crosses EB
    void amrex_flux_redistribute (const amrex::Box& bx,
                                  amrex::Array4<amrex::Real            > const& dqdt,
                                  amrex::Array4<amrex::Real       const> const& divc,
                                  amrex::Array4<amrex::Real       const> const& wt,
                                  amrex::Array4<amrex::Real       const> const& vfrac,
                                  amrex::Array4<amrex::EBCellFlag const> const& flag,
                                  int as_crse,
                                  amrex::Array4<amrex::Real     > const& rr_drho_crse,
                                  amrex::Array4<int        const> const& rr_flag_crse,
                                  int as_fine,
                                  amrex::Array4<amrex::Real     > const& dm_as_fine,
                                  amrex::Array4<int        const> const& levmsk,
                                  const amrex::Geometry & geom,
                                  bool use_wts_in_divnc,
                                  int level_mask_not_covered,
                                  int icomp, int ncomp, amrex::Real dt);
}

    // Interface to redistribution schemes that only calls single-level routines
    void ApplyRedistribution ( amrex::Box const& bx, int ncomp,
                 amrex::Array4<amrex::Real>       const& dUdt_out,
                 amrex::Array4<amrex::Real>       const& dUdt_in,
                 amrex::Array4<amrex::Real const> const& U_in,
                 amrex::Array4<amrex::Real> const& scratch,
                 amrex::Array4<amrex::EBCellFlag const> const& flag,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                              amrex::Array4<amrex::Real const> const& apy,
                              amrex::Array4<amrex::Real const> const& apz),
                 amrex::Array4<amrex::Real const> const& vfrac,
                 AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                              amrex::Array4<amrex::Real const> const& fcy,
                              amrex::Array4<amrex::Real const> const& fcz),
                 amrex::Array4<amrex::Real const> const& ccent,
                 amrex::BCRec  const* d_bcrec_ptr,
                 amrex::Geometry const& geom,
                 amrex::Real dt, std::string redistribution_type,
                 const int srd_max_order = 2,
                 amrex::Real target_volfrac = 0.5,
                 amrex::Array4<amrex::Real const> const& update_scale={});

    // Interface to redistribution schemes that calls multi-level routines
    void ApplyMLRedistribution (
        amrex::Box const& bx, int ncomp,
        amrex::Array4<amrex::Real      > const& dUdt_out,
        amrex::Array4<amrex::Real      > const& dUdt_in,
        amrex::Array4<amrex::Real const> const& U_in,
        amrex::Array4<amrex::Real      > const& wt,
        amrex::Array4<amrex::EBCellFlag const> const& flag,
        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                     amrex::Array4<amrex::Real const> const& apy,
                     amrex::Array4<amrex::Real const> const& apz),
        amrex::Array4<amrex::Real const> const& vfrac,
        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                     amrex::Array4<amrex::Real const> const& fcy,
                     amrex::Array4<amrex::Real const> const& fcz),
        amrex::Array4<amrex::Real const> const& ccc,
        amrex::BCRec  const* d_bcrec_ptr,
        amrex::Geometry const& lev_geom, amrex::Real dt,
        std::string redistribution_type,
        int as_crse,
        amrex::Array4<amrex::Real            > const& rr_drho_crse,
        amrex::Array4<int        const> const& rr_flag_crse,
        int as_fine,
        amrex::Array4<amrex::Real            > const& dm_as_fine,
        amrex::Array4<int        const> const& levmsk,
        int level_mask_not_covered,
        bool use_wts_in_divnc = false,
        int icomp = 0,
        const int srd_max_order = 2,
        amrex::Real target_volfrac = 0.5,
        amrex::Array4<amrex::Real const> const& update_scale={});

    void ApplyInitialRedistribution (
        amrex::Box const& bx, int ncomp,
        amrex::Array4<amrex::Real                  > const& U_out,
        amrex::Array4<amrex::Real                  > const& U_in,
        amrex::Array4<amrex::EBCellFlag const> const& flag,
        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                     amrex::Array4<amrex::Real const> const& apy,
                     amrex::Array4<amrex::Real const> const& apz),
        amrex::Array4<amrex::Real const> const& vfrac,
        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                     amrex::Array4<amrex::Real const> const& fcy,
                     amrex::Array4<amrex::Real const> const& fcz),
        amrex::Array4<amrex::Real const> const& ccent,
        amrex::BCRec  const* d_bcrec_ptr,
        amrex::Geometry& geom, std::string redistribution_type,
        const int srd_max_order = 2,
        amrex::Real target_volfrac = 0.5);

    void StateRedistribute ( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& dUdt_out,
                             amrex::Array4<amrex::Real> const& dUdt_in,
                             amrex::Array4<amrex::EBCellFlag const> const& flag,
                             amrex::Array4<amrex::Real const> const& vfrac,
                             AMREX_D_DECL(amrex::Array4<amrex::Real const> const& fcx,
                                          amrex::Array4<amrex::Real const> const& fcy,
                                          amrex::Array4<amrex::Real const> const& fcz),
                             amrex::Array4<amrex::Real const> const& ccent,
                             amrex::BCRec  const* d_bcrec_ptr,
                             amrex::Array4<int const> const& itracker,
                             amrex::Array4<amrex::Real const> const& nrs,
                             amrex::Array4<amrex::Real const> const& alpha,
                             amrex::Array4<amrex::Real const> const& nbhd_vol,
                             amrex::Array4<amrex::Real const> const& cent_hat,
                             amrex::Geometry const& geom,
                             const int max_order = 2);

    void MakeITracker ( amrex::Box const& bx,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                     amrex::Array4<amrex::Real const> const& apy,
                                     amrex::Array4<amrex::Real const> const& apz),
                        amrex::Array4<amrex::Real const> const& vfrac,
                        amrex::Array4<int> const& itracker,
                        amrex::Geometry const& geom,
                        amrex::Real target_volfrac);

    void MakeStateRedistUtils ( amrex::Box const& bx,
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
                                amrex::Array4<amrex::Real const> const& vfrac,
                                amrex::Array4<amrex::Real const> const& ccent,
                                amrex::Array4<        int const> const& itracker,
                                amrex::Array4<amrex::Real> const& nrs,
                                amrex::Array4<amrex::Real> const& nbhd_vol,
                                amrex::Array4<amrex::Real> const& alpha,
                                amrex::Array4<amrex::Real> const& cent_hat,
                                amrex::Geometry const& geom,
                                amrex::Real target_volfrac);

#endif
/** @} */
#endif
