#ifndef AMREX_EB2_IF_TRANSLATION_H_
#define AMREX_EB2_IF_TRANSLATION_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>

#include <type_traits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F>
class TranslationIF
{
public:

    TranslationIF (F a_f, const RealArray& a_offset)
        : m_f(std::move(a_f)),
          m_offset(makeXDim3(a_offset))
        {}

    [[nodiscard]] inline Real operator() (const RealArray& p) const noexcept
    {
        return m_f({AMREX_D_DECL(p[0]-m_offset.x,
                                 p[1]-m_offset.y,
                                 p[2]-m_offset.z)});
    }

    template <class U=F, std::enable_if_t<IsGPUable<U>::value,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
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
struct IsGPUable<TranslationIF<F>, std::enable_if_t<IsGPUable<F>::value>>
    : std::true_type {};

template <class F>
constexpr TranslationIF<std::decay_t<F>>
translate (F&&f, const RealArray& offset)
{
    return TranslationIF<std::decay_t<F>>(std::forward<F>(f),offset);
}

}

#endif
