#ifndef AMREX_EB_REDIST_H_
#define AMREX_EB_REDIST_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_BCRec.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MultiCutFab.H>
#include <AMReX_EB2.H>

/**
 * \file AMReX_EB_Redistribution.H
 *
 * Redistribution kernels used to stabilize EB finite-volume updates.
 */

namespace amrex {

#ifdef AMREX_USE_FLOAT
    //! Sentinel value used when marking EB-covered cells in float builds.
    static constexpr amrex::Real eb_covered_val = amrex::Real(1.e20);
#else
    //! Sentinel value used when marking EB-covered cells in double builds.
    static constexpr amrex::Real eb_covered_val = amrex::Real(1.e40);
#endif

    /**
     * \brief Redistribute the provisional divergence within a single level (no coarse/fine coupling).
     *
     * \param div_tmp_in  Input divergence (will be modified).
     * \param div_out     Destination MultiFab for redistributed divergence.
     * \param div_comp    Component to start from.
     * \param ncomp       Number of components to process.
     * \param geom        Geometry used for metrics.
     */
    void single_level_redistribute (amrex::MultiFab& div_tmp_in, amrex::MultiFab& div_out,
                                    int div_comp, int ncomp, const amrex::Geometry& geom);

    /**
     * \brief Weighted single-level redistribution (e.g., when weighting by density).
     *
     * \param div_tmp_in        Provisional divergence to be redistributed (modified in-place).
     * \param div_out           Destination MultiFab receiving the corrected divergence.
     * \param weights           Per-cell weights controlling distribution.
     * \param div_comp          Starting component index.
     * \param ncomp             Number of components to process.
     * \param geom              Geometry describing the current level.
     * \param use_wts_in_divnc  Whether to weight non-conservative divergence terms.
     */
    void single_level_weighted_redistribute (amrex::MultiFab& div_tmp_in, amrex::MultiFab& div_out,
                                             const amrex::MultiFab& weights,
                                             int div_comp, int ncomp, const amrex::Geometry& geom,
                                             bool use_wts_in_divnc);

    /**
     * \brief Flux redistribution scheme for coarse/fine interfaces that do not cross the EB.
     *
     * \param bx     Tile box being updated.
     * \param div    Output coarse divergence increments.
     * \param divc   Conservative divergence on coarse cells.
     * \param wt     Optional weights used for redistribution.
     * \param icomp  Component offset.
     * \param ncomp  Number of components.
     * \param flag_arr EB cut-cell flags.
     * \param vfrac  Volume fractions.
     * \param geom   Geometry describing the coarse level.
     * \param use_wts_in_divnc Whether to weight the non-conservative term.
     */
    void apply_flux_redistribution ( const amrex::Box& bx,
                                     amrex::Array4<amrex::Real      > const& div,
                                     amrex::Array4<amrex::Real const> const& divc,
                                     amrex::Array4<amrex::Real const> const& wt,
                                     int icomp, int ncomp,
                                     amrex::Array4<amrex::EBCellFlag const> const& flag_arr,
                                     amrex::Array4<amrex::Real       const> const& vfrac,
                                     const amrex::Geometry & geom,
                                     bool use_wts_in_divnc);

    /**
     * \brief Flux redistribution scheme for coarse/fine interfaces that cross the EB.
     *
     * Handles both coarse-as-fine and fine-as-coarse coupling along with neighbor masks.
     */
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

    /**
     * \brief Interface for applying single-level redistribution schemes.
     *
     * Chooses the requested algorithm (e.g., state redistribution vs. flux redistribution)
     * based on \p redistribution_type and operates solely on data from one AMR level.
     */
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
                 amrex::Array4<amrex::Real const> const& ccc,
                 amrex::BCRec  const* d_bcrec_ptr,
                 amrex::Geometry const& lev_geom,
                 amrex::Real dt, std::string const& redistribution_type,
                 bool use_wts_in_divnc = false,
                 int srd_max_order = 2,
                 amrex::Real target_volfrac = 0.5_rt,
                 amrex::Array4<amrex::Real const> const& update_scale={});

    /**
     * \brief Multi-level redistribution interface that couples coarse and fine levels.
     *
     * Parameters mirror ApplyRedistribution but also include buffers for coarse/fine overlap.
     */
    void ApplyMLRedistribution (
        amrex::Box const& bx, int ncomp,
        amrex::Array4<amrex::Real      > const& dUdt_out,
        amrex::Array4<amrex::Real      > const& dUdt_in,
        amrex::Array4<amrex::Real const> const& U_in,
        amrex::Array4<amrex::Real      > const& scratch,
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
        std::string const& redistribution_type,
        int as_crse,
        amrex::Array4<amrex::Real            > const& rr_drho_crse,
        amrex::Array4<int        const> const& rr_flag_crse,
        int as_fine,
        amrex::Array4<amrex::Real            > const& dm_as_fine,
        amrex::Array4<int        const> const& levmsk,
        int level_mask_not_covered,
        amrex::Real fac_for_deltaR = 1.0_rt,
        bool use_wts_in_divnc = false,
        int icomp = 0,
        int srd_max_order = 2,
        amrex::Real target_volfrac = 0.5_rt,
        amrex::Array4<amrex::Real const> const& update_scale={});

    /**
     * \brief Redistribute immediately after initialization (before time stepping).
     */
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
        amrex::Array4<amrex::Real const> const& ccc,
        amrex::BCRec  const* d_bcrec_ptr,
        amrex::Geometry const& geom, std::string const& redistribution_type,
        int srd_max_order = 2,
        amrex::Real target_volfrac = 0.5_rt);

    /**
     * \brief Perform state redistribution (SRD) within a single AMR level.
     */
    void StateRedistribute ( amrex::Box const& bx, int ncomp,
                             amrex::Array4<amrex::Real> const& U_out,
                             amrex::Array4<amrex::Real> const& U_in,
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
                             int max_order = 2);

    /**
     * \brief Multi-level state redistribution that couples coarse/fine tiles.
     */
    void MLStateRedistribute ( amrex::Box const& bx, int ncomp,
                               amrex::Array4<amrex::Real> const& U_out,
                               amrex::Array4<amrex::Real> const& U_in,
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
                               int as_crse,
                               Array4<Real            > const& drho_as_crse,
                               Array4<int        const> const& flag_as_crse,
                               int as_fine,
                               Array4<Real            > const& dm_as_fine,
                               Array4<int        const> const& levmsk,
                               int is_ghost_cell,
                               amrex::Real fac_for_deltaR,
                               int max_order = 2);

    /**
     * \brief Build the irregular-cell tracker used by SRD algorithms.
     */
    void MakeITracker ( amrex::Box const& bx,
                        AMREX_D_DECL(amrex::Array4<amrex::Real const> const& apx,
                                     amrex::Array4<amrex::Real const> const& apy,
                                     amrex::Array4<amrex::Real const> const& apz),
                        amrex::Array4<amrex::Real const> const& vfrac,
                        amrex::Array4<int> const& itracker,
                        amrex::Geometry const& geom,
                        amrex::Real target_volfrac);

    /**
     * \brief Construct helper arrays (normals, neighbor volumes, etc.) for SRD.
     */
    void MakeStateRedistUtils ( amrex::Box const& bx,
                                amrex::Array4<amrex::EBCellFlag const> const& flag,
                                amrex::Array4<amrex::Real const> const& vfrac,
                                amrex::Array4<amrex::Real const> const& ccent,
                                amrex::Array4<        int const> const& itracker,
                                amrex::Array4<amrex::Real> const& nrs,
                                amrex::Array4<amrex::Real> const& alpha,
                                amrex::Array4<amrex::Real> const& nbhd_vol,
                                amrex::Array4<amrex::Real> const& cent_hat,
                                amrex::Geometry const& geom,
                                amrex::Real target_volfrac);

}

#endif
