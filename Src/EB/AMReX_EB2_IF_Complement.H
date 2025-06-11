#ifndef AMREX_EB2_IF_COMPLEMENT_H_
#define AMREX_EB2_IF_COMPLEMENT_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F>
class ComplementIF
{
public:

    ComplementIF (F a_f) : m_f(std::move(a_f)) {}

    [[nodiscard]] inline Real operator() (const RealArray& p) const noexcept
    {
        return -m_f(p);
    }

    template<class U=F, std::enable_if_t<IsGPUable<U>::value,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return -m_f(AMREX_D_DECL(x,y,z));
    }

protected:

    F m_f;
};

template <class F>
struct IsGPUable<ComplementIF<F>, std::enable_if_t<IsGPUable<F>::value>>
    : std::true_type {};

template <class F>
constexpr ComplementIF<std::decay_t<F>>
makeComplement (F&& f)
{
    return ComplementIF<std::decay_t<F>>(std::forward<F>(f));
}

}

#endif
