#ifndef AMREX_BCUTIL_H_
#define AMREX_BCUTIL_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_BCUtil.H
 * \brief Utility function for filling physical-domain ghost cells on cell-centered MultiFabs.
 */

#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BCRec.H>

namespace amrex
{
    /**
     * \brief Fill ghost cells outside the physical domain for cell-centered data
     *        (excluding periodic boundaries).
     *
     * It only fills BCType::foextrap, BCType::hoextrap, BCType::hoextrapcc,
     * BCType::reflect_even, and BCType::reflect_odd.  It does not fill
     * BCType::ext_dir and BCType::ext_dir_cc (i.e., external Dirichlet).
     * If you have BCType::ext_dir or BCType::ext_dir_cc, you can copy,
     * paste and rename this function, and then modify it.
     *
     * \param phi Cell-centered field whose ghost cells are updated in-place.
     * \param geom Geometry describing the physical domain extents.
     * \param bc Per-component boundary descriptors controlling the fill behavior.
     */
    void FillDomainBoundary (MultiFab& phi, const Geometry& geom, const Vector<BCRec>& bc);
}

#endif
