/* Copyright 2019 Remi Lehe
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARALLELIZATION_COMM_H_
#define WARPX_PARALLELIZATION_COMM_H_

#include <AMReX_MultiFab.H>

/** \brief Fills the values of the current on the coarse patch by
 *  averaging the values of the current of the fine patch (on the same level).
 *  Also fills the guards of the coarse patch.
 *
 * \param[in]  fine             fine patches to interpolate from
 * \param[out] coarse           coarse patches to interpolate to
 * \param[in]  refinement_ratio integer ratio between the two
 */
void
interpolateCurrentFineToCoarse (
      std::array< amrex::MultiFab const *, 3 > const & fine,
      std::array< amrex::MultiFab       *, 3 > const & coarse,
      int const refinement_ratio);

/** \brief Fills the values of the charge density on the coarse patch by
 *  averaging the values of the charge density of the fine patch (on the same level).
 *
 * \param[in]  fine             fine patches to interpolate from
 * \param[out] coarse           coarse patches to interpolate to
 * \param[in]  refinement_ratio integer ratio between the two
 */
void
interpolateDensityFineToCoarse (
     const amrex::MultiFab& fine,
     amrex::MultiFab& coarse,
     int const refinement_ratio);

#endif // WARPX_PARALLELIZATION_COMM_H_
