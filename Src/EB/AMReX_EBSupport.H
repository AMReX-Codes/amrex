#ifndef AMREX_EB_SUPPORT_H_
#define AMREX_EB_SUPPORT_H_

namespace amrex {

    enum class EBSupport {
        none,
        basic,     // EBCellFlag
        geometry   // volume fraction, area fraction, centroids
    };

}

#endif
