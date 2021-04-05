#ifndef AMREX_EB_SUPPORT_H_
#define AMREX_EB_SUPPORT_H_
#include <AMReX_Config.H>

namespace amrex {

    enum struct EBSupport : int {
        none     = 0,
        basic    = 1,     //!< EBCellFlag
        volume   = 2,     //!< + volume fraction
        full     = 3      //!< + area fraction, boundary centroids and face centroids
    };

}

#endif
