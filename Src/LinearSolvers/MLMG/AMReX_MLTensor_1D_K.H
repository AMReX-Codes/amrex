#ifndef AMREX_MLTENSOR_1D_K_H_
#define AMREX_MLTENSOR_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mltensor_adotx_add_extra (Box const& box, Array4<Real> const& Ax,
                               Array4<Real const> const& vel,
                               Array4<Real const> const& gradeta,
                               GpuArray<Real,AMREX_SPACEDIM> const& dxinv)
{

}

}

#endif
