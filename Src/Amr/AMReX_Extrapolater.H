#ifndef AMREX_EXTRAPOLATER_H_
#define AMREX_EXTRAPOLATER_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

/**
 * \file AMReX_Extrapolater.H
 * \brief Utilities for first-order cell-value extrapolation to domain boundaries.
 */

namespace amrex::Extrapolater
{
    // finebnd: boundary cells covered by fine cells (including periodically shifted fine cells)
    // crsebnd: boundary cells not covered by fine cells
    // physbnd: boundary cells outside the domain (excluding periodic boundaries)
    // interior: interior cells
    const int finebnd = 1;
    const int crsebnd = 0;
    const int physbnd = 0;
    const int interior = 1;

    /**
     * \brief Extrapolate cell-centered data outward one cell to fill boundary conditions.
     *
     * \param mf    MultiFab whose ghost cells will be filled.
     * \param geom  Geometry describing the level (used for domain bounds).
     * \param scomp Starting component index in \p mf.
     * \param ncomp Number of components to extrapolate.
     *
     * \note FillBoundary (with periodicity) should be called on \p mf before invoking this routine.
     */
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp);

    /**
     * \brief Overload that extrapolates into a caller-provided number of grow cells.
     *
     * \param mf    MultiFab whose ghost cells will be filled.
     * \param geom  Geometry describing the level.
     * \param scomp Starting component index in \p mf.
     * \param ncomp Number of components to extrapolate.
     * \param nGrow Number of ghost cells to update (instead of using mf.nGrow()).
     */
    void FirstOrderExtrap (MultiFab& mf, const Geometry& geom, int scomp, int ncomp, int nGrow);

}

#endif
