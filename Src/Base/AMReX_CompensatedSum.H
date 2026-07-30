#ifndef AMREX_COMPENSATEDSUM_H_
#define AMREX_COMPENSATEDSUM_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

namespace amrex
{

/**
 * \brief Add a value to a running sum with Kahan compensation.
 *
 * Callers own both the accumulated sum and the compensation storage. Both
 * should be initialized before the first call, usually to zero.
 */
template <typename T, typename U>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compensatedAdd (T& sum, T& compensation, U value) noexcept
{
    T const y = static_cast<T>(value) - compensation;
    T const t = sum + y;
    compensation = (t - sum) - y;
    sum = t;
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compensatedReset (T& sum, T& compensation) noexcept
{
    sum = T(0);
    compensation = T(0);
}

template <typename T, typename U>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compensatedSet (T& sum, T& compensation, U value) noexcept
{
    sum = static_cast<T>(value);
    compensation = T(0);
}

}

#endif
