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
 * should be initialized before the first call, usually to zero. T should be
 * a floating-point type.
 *
 * Kahan compensation relies on (t - sum) - y evaluating to the rounding error
 * lost by sum + y, rather than to the exact-arithmetic value of 0.
 * Optimizations that assume floating-point addition is associative (e.g.
 * -ffast-math, -fassociative-math) are free to simplify that expression to 0,
 * silently turning this into a plain, uncompensated sum. The volatile
 * qualifiers below force the intermediate values to be materialized, creating
 * barriers the optimizer cannot reassociate across even when the caller is
 * compiled with fast-math flags.
 *
 * All three intermediates must be volatile. Marking only t and c still lets
 * some compilers substitute y into c - y and reassociate the expression into a
 * different rounding path under fast-math.
 *
 * \param[in,out] sum Running sum updated by \p value.
 * \param[in,out] compensation Kahan compensation term updated to track the
 *                             rounding error from the addition.
 * \param[in] value Value to add to \p sum.
 */
template <typename T, typename U>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compensatedAdd (T& sum, T& compensation, U value) noexcept
{
    T volatile y = static_cast<T>(value) - compensation;
    T volatile t = sum + y;
    T volatile c = t - sum;
    compensation = c - y;
    sum = t;
}

}

#endif
