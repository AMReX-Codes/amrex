#ifndef AMREX_MG_2D_K_H_
#define AMREX_MG_2D_K_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mg_cc_interp (int i, int j, int /*k*/, int n,
                   Array4<Real> const& f, Array4<Real const> const& c) noexcept
{
    int i2 = 2*i;
    int j2 = 2*j;
    int i2p1 = i2+1;
    int j2p1 = j2+1;
    Real cv = c(i,j,0,n);
    f(i2  ,j2  ,0,n) += cv;
    f(i2p1,j2  ,0,n) += cv;
    f(i2  ,j2p1,0,n) += cv;
    f(i2p1,j2p1,0,n) += cv;
}

}

#endif
