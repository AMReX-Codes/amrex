#ifndef AMREX_EB2_IF_LATHE_H_
#define AMREX_EB2_IF_LATHE_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>
#include <cmath>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F>
class LatheIF
{
public:

    LatheIF (F a_f) : m_f(std::move(a_f)) {}

    [[nodiscard]] inline Real operator() (const RealArray& p) const noexcept
    {
        Real r = std::hypot(p[0],p[1]);
#if (AMREX_SPACEDIM == 2)
        return m_f({r,0.0});
#else
        return m_f({r,p[2],0.0});
#endif
    }

    template <class U=F, std::enable_if_t<IsGPUable<U>::value,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        Real r = std::hypot(x,y);
#if (AMREX_SPACEDIM == 2)
        return m_f(r,0.0);
#else
        return m_f(r,z,0.0);
#endif
    }

protected:

    F m_f;
};

template <class F>
struct IsGPUable<LatheIF<F>, std::enable_if_t<IsGPUable<F>::value>>
    : std::true_type {};

template <class F>
constexpr LatheIF<std::decay_t<F>>
lathe (F&& f)
{
    return LatheIF<std::decay_t<F>>(std::forward<F>(f));
}

}

#endif
