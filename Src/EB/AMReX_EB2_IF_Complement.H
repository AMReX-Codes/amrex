#ifndef AMREX_EB2_IF_COMPLEMENT_H_
#define AMREX_EB2_IF_COMPLEMENT_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F>
class ComplementIF
{
public:

    ComplementIF (F&& a_f) : m_f(std::move(a_f)) {}
    ComplementIF (F const& a_f) : m_f(a_f) {}
          
    ComplementIF (const ComplementIF& rhs) = default;
    ComplementIF (ComplementIF&& rhs) = default;
    ComplementIF& operator= (const ComplementIF& rhs) = delete;
    ComplementIF& operator= (ComplementIF&& rhs) = delete;

    inline Real operator() (const RealArray& p) const noexcept
    {
        return -m_f(p);
    }

    template<class U=F, class = typename std::enable_if<IsGPUable<U>::value>::type >
    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return -m_f(AMREX_D_DECL(x,y,z));
    }

protected:

    F m_f;
};

template <class F>
struct IsGPUable<ComplementIF<F>, typename std::enable_if<IsGPUable<F>::value>::type>
    : std::true_type {};

template <class F>
constexpr ComplementIF<typename std::decay<F>::type>
makeComplement (F&& f)
{
    return ComplementIF<typename std::decay<F>::type>(std::forward<F>(f));
}

}}

#endif

