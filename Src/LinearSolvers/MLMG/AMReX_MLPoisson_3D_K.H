#ifndef AMREX_MLPOISSON_3D_K_H_
#define AMREX_MLPOISSON_3D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlpoisson_adotx (int i, int j, int k, Array4<Real> const& y,
                      Array4<Real const> const& x,
                      Real dhx, Real dhy, Real dhz)
{
    y(i,j,k) = dhx * (x(i-1,j,k) - 2.0*x(i,j,k) + x(i+1,j,k))
        +      dhy * (x(i,j-1,k) - 2.0*x(i,j,k) + x(i,j+1,k))
        +      dhz * (x(i,j,k-1) - 2.0*x(i,j,k) + x(i,j,k+1));
}

}

#endif
