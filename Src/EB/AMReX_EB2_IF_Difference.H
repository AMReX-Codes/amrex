#ifndef AMREX_EB2_IF_DIFFERENCE_H_
#define AMREX_EB2_IF_DIFFERENCE_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>
#include <algorithm>
#include <utility>

/**
 * \file AMReX_EB2_IF_Difference.H
 *
 * Constructive solid-geometry helper that subtracts one implicit function from another.
 */

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F, class G>
//! Implicit function that keeps points inside \p F but outside \p G.
class DifferenceIF
{
public:

    //! Store the minuend `a_f` and subtrahend `a_g` implicit functions.
    DifferenceIF (F a_f, G a_g)
        : m_f(std::move(a_f)),
          m_g(std::move(a_g))
        {}

    //! Evaluate the set difference at \p p.
    [[nodiscard]] AMREX_FORCE_INLINE Real operator() (const RealArray& p) const noexcept
    {
        Real r1 = m_f(p);
        Real r2 = m_g(p);
        return amrex::min(r1, -r2);
    }

    //! Evaluate the set difference at given physical coordinates.
    template <class U=F, class V=G>
    requires (IsGPUable<U>::value && IsGPUable<V>::value)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        Real r1 = m_f(AMREX_D_DECL(x,y,z));
        Real r2 = m_g(AMREX_D_DECL(x,y,z));
        return amrex::min(r1, -r2);
    }

protected:

    F m_f;
    G m_g;
};

template <class F, class G>
requires (IsGPUable<F>::value && IsGPUable<G>::value)
struct IsGPUable<DifferenceIF<F,G>> : std::true_type {};

template <class F, class G>
//! Helper that constructs a difference operator.
constexpr DifferenceIF<std::decay_t<F>,
                       std::decay_t<G>>
makeDifference (F&& f, G&& g)
{
    return DifferenceIF<std::decay_t<F>,
                        std::decay_t<G>>
        (std::forward<F>(f), std::forward<G>(g));
}

}

#endif
