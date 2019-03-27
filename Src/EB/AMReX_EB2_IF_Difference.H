#ifndef AMREX_EB2_IF_DIFFERENCE_H_
#define AMREX_EB2_IF_DIFFERENCE_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>
#include <algorithm>
#include <utility>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F, class G>
class DifferenceIF
{
public:

    DifferenceIF (F&& a_f, G&& a_g)
        : m_f(std::move(a_f)),
          m_g(std::move(a_g))
        {}
    DifferenceIF (F const& a_f, G const& a_g)
        : m_f(a_f),
          m_g(a_g)
    {}

    DifferenceIF (const DifferenceIF& rhs) = default;
    DifferenceIF (DifferenceIF&& rhs) = default;
    DifferenceIF& operator= (const DifferenceIF& rhs) = delete;
    DifferenceIF& operator= (DifferenceIF&& rhs) = delete;

    inline Real operator() (const RealArray& p) const noexcept
    {
        Real r1 = m_f(p);
        Real r2 = m_g(p);
        return amrex::min(r1, -r2);
    }

    template <class U=F, class V=G,
              class = typename std::enable_if<IsGPUable<U>::value &&
                                              IsGPUable<V>::value>::type>
    AMREX_GPU_HOST_DEVICE inline
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
struct IsGPUable<DifferenceIF<F,G>, typename std::enable_if<IsGPUable<F>::value &&
                                                            IsGPUable<G>::value>::type>
    : std::true_type {};

template <class F, class G>
constexpr DifferenceIF<typename std::decay<F>::type,
                       typename std::decay<G>::type>
makeDifference (F&& f, G&& g)
{
    return DifferenceIF<typename std::decay<F>::type,
                        typename std::decay<G>::type>
        (std::forward<F>(f), std::forward<G>(g));
}

}}

#endif
