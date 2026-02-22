#ifndef AMREX_EB_SUPPORT_H_
#define AMREX_EB_SUPPORT_H_
#include <AMReX_Config.H>

#include <AMReX_Config.H>

/**
 * \file AMReX_EBSupport.H
 *
 * Enumerates how much embedded-boundary metadata should be allocated.
 */

namespace amrex {

    //! Requested EB metadata level for factories and utilities.
    enum struct EBSupport : int {
        none     = 0,
        basic    = 1,     //!< EBCellFlag
        volume   = 2,     //!< + volume fraction
        full     = 3      //!< + area fraction, boundary centroids and face centroids
    };

}

#endif
