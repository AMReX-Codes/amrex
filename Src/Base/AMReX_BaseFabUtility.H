#ifndef AMREX_BASEFAB_UTILITY_H_
#define AMREX_BASEFAB_UTILITY_H_

#include <AMReX_BaseFab.H>

namespace amrex {

template <class Tto, class Tfrom>
AMREX_GPU_HOST_DEVICE
void
cast (BaseFab<Tto>& tofab, BaseFab<Tfrom> const& fromfab,
      Box const& bx, SrcComp scomp, DestComp dcomp, NumComps ncomp) noexcept
{
    const auto len = length(bx);
    const auto lo  = lbound(bx);
    const auto tdata =   tofab.view(lo,dcomp.i);
    const auto fdata = fromfab.view(lo,scomp.i);

    for (int n = 0; n < ncomp.n; ++n) {
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
