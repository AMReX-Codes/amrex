#ifndef AMREX_BASEFAB_UTILITY_H_
#define AMREX_BASEFAB_UTILITY_H_

#include <AMReX_BaseFab.H>

namespace amrex {

template <class Tto, class Tfrom>
AMREX_GPU_HOST_DEVICE
void
cast (Box const& bx, BaseFab<Tto>& tofab, BaseFab<Tfrom> const& fromfab,
      int tfcomp, int ffcomp, int ncomp)
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto tdata =   tofab.view(lo,tfcomp);
    const auto fdata = fromfab.view(lo,ffcomp);

    for (int n = 0; n < ncomp; ++n) {
        for         (int k = 0; k < len.z; ++k) {
            for     (int j = 0; j < len.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = 0; i < len.x; ++i) {
                    tdata(i,j,k,n) = static_cast<Tto>(fdata(i,j,k,n));
                }
            }
        }
    }
}

}

#endif
