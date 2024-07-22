#ifndef AMREX_EB2_IF_DIFFERENCE_H_
#define AMREX_EB2_IF_DIFFERENCE_H_
#include <AMReX_Config.H>

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <type_traits>
#include <algorithm>
#include <utility>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex::EB2 {

template <class F, class G>
class DifferenceIF
{
public:

    DifferenceIF (F a_f, G a_g)
        : m_f(std::move(a_f)),
          m_g(std::move(a_g))
        {}

    [[nodiscard]] inline Real operator() (const RealArray& p) const noexcept
    {
        Real r1 = m_f(p);
        Real r2 = m_g(p);
        return amrex::min(r1, -r2);
    }

    template <class U=F, class V=G,
              std::enable_if_t<IsGPUable<U>::value &&
                                      IsGPUable<V>::value, int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
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
struct IsGPUable<DifferenceIF<F,G>, std::enable_if_t<IsGPUable<F>::value &&
                                                            IsGPUable<G>::value>>
    : std::true_type {};

template <class F, class G>
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
