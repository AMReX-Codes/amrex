#ifndef AMREX_COMPENSATEDSUM_H_
#define AMREX_COMPENSATEDSUM_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex
{

template <typename T, typename U>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compensatedSet (T& sum, T& compensation, U value) noexcept
{
    sum = static_cast<T>(value);
    compensation = T(0);
}

}

#endif
