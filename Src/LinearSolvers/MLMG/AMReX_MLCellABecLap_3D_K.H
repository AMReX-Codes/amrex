#ifndef AMREX_MLCELLABECLAP_3D_K_H_
#define AMREX_MLCELLABECLAP_3D_K_H_
#include <AMReX_Config.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int coarsen_overset_mask (Box const& bx, Array4<int> const& cmsk, Array4<int const> const& fmsk)
    noexcept
{
    int nerrors = 0;
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    for (int i = lo.x; i <= hi.x; ++i) {
        int ii = 2*i;
        int jj = 2*j;
        int kk = 2*k;
        cmsk(i,j,k) = fmsk(ii,jj  ,kk  ) + fmsk(ii+1,jj  ,kk  )
            +         fmsk(ii,jj+1,kk  ) + fmsk(ii+1,jj+1,kk  )
            +         fmsk(ii,jj  ,kk+1) + fmsk(ii+1,jj  ,kk+1)
            +         fmsk(ii,jj+1,kk+1) + fmsk(ii+1,jj+1,kk+1);
        if (cmsk(i,j,k) == 8) {
            cmsk(i,j,k) = 1;
        } else if (cmsk(i,j,k) != 0) {
            ++nerrors;
        }
    }}}
    return nerrors;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void coarsen_overset_mask (int i, int j, int k, Array4<int> const& cmsk,
                           Array4<int const> const& fmsk) noexcept
{
    int ii = 2*i;
    int jj = 2*j;
    int kk = 2*k;
    cmsk(i,j,k) = fmsk(ii,jj  ,kk  ) + fmsk(ii+1,jj  ,kk  )
        +         fmsk(ii,jj+1,kk  ) + fmsk(ii+1,jj+1,kk  )
        +         fmsk(ii,jj  ,kk+1) + fmsk(ii+1,jj  ,kk+1)
        +         fmsk(ii,jj+1,kk+1) + fmsk(ii+1,jj+1,kk+1);
    if (cmsk(i,j,k) == 8) {
        cmsk(i,j,k) = 1;
    }
}

}

#endif
