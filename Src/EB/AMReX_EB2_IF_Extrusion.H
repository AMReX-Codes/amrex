#ifndef AMREX_EB2_IF_EXTRUSION_H_
#define AMREX_EB2_IF_EXTRUSION_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F>
class ExtrusionIF
{
public:

    ExtrusionIF (F a_f, int direction)
        : m_f(std::move(a_f)),
          m_direction(direction)
        {}

    [[nodiscard]] inline Real operator() (const RealArray& p) const
    {
        RealArray x = p;
        x[m_direction] = 0.0;
        return m_f(x);
    }

    template <class U=F, std::enable_if_t<IsGPUable<U>::value,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        switch (m_direction)
        {
        case 0:
            return m_f(AMREX_D_DECL(0.0, y, z));
        case 1:
            return m_f(AMREX_D_DECL(x, 0.0, z));
        default:
            return m_f(AMREX_D_DECL(x, y, 0.0));
        }
    }

protected:

    F m_f;
    int m_direction;
};

template <class F>
struct IsGPUable<ExtrusionIF<F>, std::enable_if_t<IsGPUable<F>::value>>
    : std::true_type {};

template <class F>
constexpr ExtrusionIF<std::decay_t<F>>
extrude (F&&f, int direction)
{
    return ExtrusionIF<std::decay_t<F>>(std::forward<F>(f),direction);
}

}

#endif
