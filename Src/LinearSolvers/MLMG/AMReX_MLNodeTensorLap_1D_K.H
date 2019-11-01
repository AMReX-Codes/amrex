#ifndef AMREX_MLNODETENSORLAP_1D_K_H_
#define AMREX_MLNODETENSORLAP_1D_K_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_interpadd (int i, int, int, Array4<Real> const& fine,
                          Array4<Real const> const& crse, Array4<int const> const& msk) noexcept
{}



}

#endif
