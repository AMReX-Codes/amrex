#ifndef AMREX_MLPOISSON_2D_K_H_
#define AMREX_MLPOISSON_2D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int /*j*/, int /*k*/, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx,
                      Real const* AMREX_RESTRICT re,
                      int rlo)
{
    y(i,0,0) = dhx * (re[i-rlo]*x(i-1,0,0)
                      - (re[i-rlo]+re[i+1-rlo])*x(i,0,0)
                      + re[i+1-rlo]*x(i+1,0,0));
}

}

#endif
