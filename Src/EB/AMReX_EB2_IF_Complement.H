#ifndef AMREX_EB2_IF_COMPLEMENT_H_
#define AMREX_EB2_IF_COMPLEMENT_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>

/**
 * \file AMReX_EB2_IF_Complement.H
 *
 * Boolean complement helper that flips the sign of any implicit function.
 */

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F>
//! Implicit function wrapper that swaps interior/exterior.
class ComplementIF
{
public:

    //! Store the child implicit function to be complemented.
    ComplementIF (F a_f) : m_f(std::move(a_f)) {}

    //! Evaluate the complement in array form.
    [[nodiscard]] AMREX_FORCE_INLINE Real operator() (const RealArray& p) const noexcept
    {
        return -m_f(p);
    }


    //! Evaluate the complement.
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
        requires (IsGPUable<F>::value)
    {
        return -m_f(AMREX_D_DECL(x,y,z));
    }

protected:

    F m_f;
};

template <class F>
requires (IsGPUable<F>::value)
struct IsGPUable<ComplementIF<F>> : std::true_type {};

template <class F>
//! Helper that constructs a complement operator.
constexpr ComplementIF<std::decay_t<F>>
makeComplement (F&& f)
{
    return ComplementIF<std::decay_t<F>>(std::forward<F>(f));
}

}

#endif
