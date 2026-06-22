#ifndef AMREX_MLNODETENSORLAP_1D_K_H_
#define AMREX_MLNODETENSORLAP_1D_K_H_
#include <AMReX_Config.H>

namespace amrex {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_interpadd (int, int, int, Array4<T> const&,
                          Array4<T const> const&, Array4<int const> const&) noexcept
{}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlndtslap_semi_interpadd (int, int, int, Array4<T> const&,
                               Array4<T const> const&, Array4<int const> const&,
                               int) noexcept
{
}

}

#endif
