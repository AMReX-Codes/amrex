#ifndef WARPX_K_H_
#define WARPX_K_H_

#include "cstdlib"
#include <AMReX_FArrayBox.H>

using namespace amrex;

AMREX_GPU_HOST_DEVICE AMREX_INLINE

void warpx_apply_filter_2d (int j, int k, int l, Array4<Real> const& dst, Array4<const Real> const& src,
                            Vector<Real> const& stencil_x, Vector<Real> const& stencil_z,
                            int nx, int nz)
{
    for (int ix=-nx+1; ix<nx; ix++){
    for (int iz=-nz+1; iz<nz; iz++){
        dst(j,k,0) = src(j,k,0) + stencil_x[abs(ix)]*stencil_z[abs(iz)]*src(j+ix,k+iz,0);
    }
    }
}
#endif
