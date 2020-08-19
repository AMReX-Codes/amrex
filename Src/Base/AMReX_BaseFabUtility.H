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
    auto const& tdata = tofab.array();
    auto const& fdata = fromfab.const_array();
    amrex::LoopConcurrent(bx, ncomp.n, [=] (int i, int j, int k, int n) noexcept
    {
        tdata(i,j,k,n+dcomp.i) = static_cast<Tto>(fdata(i,j,k,n+scomp.i));
    });
}

}

#endif
