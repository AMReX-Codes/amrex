#ifndef AMREX_EB2_IF_TRANSLATION_H_
#define AMREX_EB2_IF_TRANSLATION_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>

#include <type_traits>

/**
 * \file AMReX_EB2_IF_Translation.H
 *
 * Translate an implicit function by a constant offset.
 */

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F>
//! Translate an implicit surface by subtracting an offset vector.
class TranslationIF
{
public:

    //! Construct a translated implicit function.
    TranslationIF (F a_f, const RealArray& a_offset)
        : m_f(std::move(a_f)),
          m_offset(makeXDim3(a_offset))
        {}

    //! Evaluate the translated function on an array of coordinates.
    [[nodiscard]] AMREX_FORCE_INLINE Real operator() (const RealArray& p) const noexcept
    {
        return m_f({AMREX_D_DECL(p[0]-m_offset.x,
                                 p[1]-m_offset.y,
                                 p[2]-m_offset.z)});
    }

    //! GPU-friendly overload for AMREX_D_DECL coordinates.
    template <class U=F>
    requires (IsGPUable<U>::value)
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return m_f(AMREX_D_DECL(x-m_offset.x,
                                y-m_offset.y,
                                z-m_offset.z));
    }

protected:

    F m_f;
    XDim3 m_offset;
};

template <class F>
requires(IsGPUable<F>::value)
struct IsGPUable<TranslationIF<F>> : std::true_type {};

template <class F>
constexpr TranslationIF<std::decay_t<F>>
translate (F&&f, const RealArray& offset)
{
    return TranslationIF<std::decay_t<F>>(std::forward<F>(f),offset);
}

}

#endif
