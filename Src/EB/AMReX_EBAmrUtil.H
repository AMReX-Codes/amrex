#ifndef AMREX_EBAMRUTIL_H_
#define AMREX_EBAMRUTIL_H_
#include <AMReX_Config.H>

#include <AMReX_TagBox.H>
#include <AMReX_MultiFab.H>

/**
 * \file AMReX_EBAmrUtil.H
 *
 * Tagging helpers used when adapting AMR meshes for EB simulations.
 */

namespace amrex {

    /**
     * \brief Tag every cell intersected by the EB interface.
     *
     * \param tags  Output tagging array (set to TagBox::SET).
     * \param state State MultiFab containing EB data (typically EB flags).
     */
    void TagCutCells (TagBoxArray& tags, const MultiFab& state);

    /**
     * \brief Tag cells whose volume fraction differs from unity by more than \p tol.
     *
     * \param tags    Output tagging array.
     * \param volfrac Volume-fraction MultiFab (usually from EB geometry).
     * \param tol     Threshold beyond which cells are tagged.
     */
    void TagVolfrac  (TagBoxArray& tags, const MultiFab& volfrac, Real tol = 0.000001_rt);

}

#endif
