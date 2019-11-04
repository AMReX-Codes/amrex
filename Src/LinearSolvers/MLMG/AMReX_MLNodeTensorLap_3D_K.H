#ifndef AMREX_MLNODETENSORLAP_3D_K_H_
#define AMREX_MLNODETENSORLAP_3D_K_H_

namespace amrex {

namespace {

}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_interpadd (int i, int j, int k, Array4<Real> const& fine,
                          Array4<Real const> const& crse, Array4<int const> const& msk) noexcept
{
    if (!msk(i,j,k)) {
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_adotx (Box const& b, Array4<Real> const& y, Array4<Real const> const& x,
                      Array4<int const> const& msk, GpuArray<Real,6> const& s,
                      GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
}

}

#endif
