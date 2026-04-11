#ifndef AMREX_EBMULTIFABUTIL_H_
#define AMREX_EBMULTIFABUTIL_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BCRec.H>

/**
 * \file AMReX_EBMultiFabUtil.H
 *
 * Utility routines for manipulating MultiFabs that carry EB data.
 */

namespace amrex
{
    /**
     * \brief Fill all covered cells with a single value \p val.
     */
    void EB_set_covered (MultiFab& mf,                                               Real   val);
    /**
     * \brief Fill `ncomp` components (starting at `icomp`) of covered cells with \p val.
     *
     * \param mf    MultiFab to modify.
     * \param icomp First component to fill.
     * \param ncomp Number of components to fill.
     * \param ngrow Number of ghost cells to include.
     * \param val   Value assigned to covered cells.
     */
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp, int ngrow,              Real   val);
    /**
     * \brief Fill covered cells with per-component constants provided in \p vals.
     *
     * \param mf    MultiFab to modify.
     * \param icomp First component to fill.
     * \param ncomp Number of components to fill.
     * \param vals  Per-component values applied to covered cells.
     */
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp,            const Vector<Real>& vals);
    /**
     * \brief Fill covered cells (including \p ngrow ghost cells) with per-component values.
     *
     * \param mf    MultiFab to modify.
     * \param icomp First component to fill.
     * \param ncomp Number of components to fill.
     * \param ngrow Number of ghost cells included.
     * \param vals  Per-component values applied to covered cells.
     */
    void EB_set_covered (MultiFab& mf, int icomp, int ncomp, int ngrow, const Vector<Real>& vals);

    /**
     * \brief Fill covered faces in all directions with constant \p val.
     */
    void EB_set_covered_faces (const Array<MultiFab*,AMREX_SPACEDIM>& umac, Real val);
    /**
     * \brief Fill covered faces with per-component values.
     */
    void EB_set_covered_faces (const Array<MultiFab*,AMREX_SPACEDIM>& umac, int scomp, int ncomp, const Vector<Real>& vals);

    /**
     * \brief Volume-weighted average-down from fine to coarse using EB volume fractions.
     *
     * \param S_fine     Fine-level MultiFab.
     * \param S_crse     Coarse-level MultiFab to receive averaged data.
     * \param vol_fine   Fine cell volumes.
     * \param vfrac_fine Fine EB volume fractions.
     * \param scomp      Source component.
     * \param ncomp      Number of components.
     * \param ratio      Refinement ratio.
     */
    void EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, const MultiFab& vol_fine,
                          const MultiFab& vfrac_fine, int scomp, int ncomp, const IntVect& ratio);

    //! Convenience average-down that assumes unit volumes (with EB volume fraction weights).
    void EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp,
                          int ratio);
    //! Same as above but with anisotropic ratio.
    void EB_average_down (const MultiFab& S_fine, MultiFab& S_crse, int scomp, int ncomp,
                          const IntVect& ratio);

    /**
     * \brief Average face-centered data from fine to coarse grids weighted by area fractions.
     *
     * \param fine   Fine-level face-centered data (one MultiFab per direction).
     * \param crse   Coarse-level destinations (one per direction).
     * \param ratio  Refinement ratio between levels.
     * \param ngcrse Number of coarse ghost cells to fill.
     */
    void EB_average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                                const Array<      MultiFab*,AMREX_SPACEDIM>& crse,
                                int ratio, int ngcrse);
    //! Variant that accepts anisotropic refinement ratio.
    void EB_average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                                const Array<      MultiFab*,AMREX_SPACEDIM>& crse,
                                const IntVect& ratio, int ngcrse);

    /**
     * \brief Average face-centered data weighted by area fractions while respecting periodicity in \p crse_geom.
     */
    void EB_average_down_faces (const Array<const MultiFab*,AMREX_SPACEDIM>& fine,
                                const Array<      MultiFab*,AMREX_SPACEDIM>& crse,
                                const IntVect& ratio, const Geometry& crse_geom);

    //! Average EB boundary data from fine to coarse faces.
    void EB_average_down_boundaries (const MultiFab& fine, MultiFab& crse,
                                     int ratio, int ngcrse);
    //! Anisotropic refinement version of EB_average_down_boundaries().
    void EB_average_down_boundaries (const MultiFab& fine, MultiFab& crse,
                                     const IntVect& ratio, int ngcrse);

    /**
     * \brief Compute divergence of face-centered velocity \p umac.
     *
     * \param divu   Output MultiFab receiving the divergence.
     * \param umac   Face-centered velocities (array of pointers).
     * \param geom   Geometry describing the domain.
     * \param already_on_centroids True if fluxes already live on face centroids.
     */
    void EB_computeDivergence (MultiFab& divu, const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                               const Geometry& geom, bool already_on_centroids);

    /**
     * \brief Compute divergence including EB-normal velocities stored in \p vel_eb.
     *
     * \param divu    Output MultiFab receiving the divergence.
     * \param umac    Face-centered velocities (array of pointers).
     * \param geom    Geometry describing the domain.
     * \param already_on_centroids True if face fluxes already live on centroids.
     * \param vel_eb  EB-normal velocities defined on cut cells.
     */
    void EB_computeDivergence (MultiFab& divu, const Array<MultiFab const*,AMREX_SPACEDIM>& umac,
                               const Geometry& geom, bool already_on_centroids,
                               const MultiFab& vel_eb);

    /**
     * \brief Average face-centered values to cell centers.
     */
    void EB_average_face_to_cellcenter (MultiFab& ccmf, int dcomp,
                                        const Array<MultiFab const*,AMREX_SPACEDIM>& fmf);

    /**
     * \brief Interpolate cell-centered data to cell centroids.
     */
    void EB_interp_CC_to_Centroid (MultiFab& cent, const MultiFab& cc, int scomp, int dcomp, int ncomp, const Geometry& geom);

    /**
     * \brief Interpolate cell-centered data to face centroids for each direction.
     */
    void EB_interp_CC_to_FaceCentroid (const MultiFab& cc,
                                       AMREX_D_DECL( MultiFab& fc_x,
                                                     MultiFab& fc_y,
                                                     MultiFab& fc_z),
                                       int scomp, int dcomp, int nc,
                                       const Geometry& geom,
                                       const amrex::Vector<amrex::BCRec>& a_bcs);

    //! Interpolate cell centroids to face centroids (array overload).
    void EB_interp_CellCentroid_to_FaceCentroid (const MultiFab& phi_centroid,
                                                 const Array<MultiFab *,AMREX_SPACEDIM>& phi_faces,
                                                 int scomp, int dcomp, int nc,
                                                 const Geometry& geom,
                                                 const amrex::Vector<amrex::BCRec>& a_bcs);

    //! Same as above but using Vector container for faces.
    void EB_interp_CellCentroid_to_FaceCentroid (const MultiFab& phi_centroid,
                                                 const Vector<MultiFab*>& phi_faces,
                                                 int scomp, int dcomp, int nc,
                                                 const Geometry& geom,
                                                 const amrex::Vector<amrex::BCRec>& a_bcs);

    //! Component-wise overload taking explicit face references.
    void EB_interp_CellCentroid_to_FaceCentroid (const MultiFab& phi_centroid,
                                                 AMREX_D_DECL( MultiFab& phi_xface,
                                                               MultiFab& phi_yface,
                                                               MultiFab& phi_zface),
                                                 int scomp, int dcomp, int nc,
                                                 const Geometry& geom,
                                                 const amrex::Vector<amrex::BCRec>& a_bcs);
}

#endif
