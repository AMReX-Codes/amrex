#ifndef AMREX_MG_1D_K_H_
#define AMREX_MG_1D_K_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mg_cc_interp (int i, int /*j*/, int /*k*/, int n,
                   Array4<Real> const& f, Array4<Real const> const& c) noexcept
{
    int i2 = 2*i;
    int i2p1 = i2+1;
    Real cv = c(i,0,0,n);
    f(i2  ,0,0,n) += cv;
    f(i2p1,0,0,n) += cv;
}

}

#endif
