#ifndef AMREX_EB2_IF_DIFFERENCE_H_
#define AMREX_EB2_IF_DIFFERENCE_H_

#include <AMReX_Array.H>

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

    ~DifferenceIF () {}

    DifferenceIF (const DifferenceIF& rhs) = default;
    DifferenceIF (DifferenceIF&& rhs) = default;
    DifferenceIF& operator= (const DifferenceIF& rhs) = delete;
    DifferenceIF& operator= (DifferenceIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        Real r1 = m_f(p);
        Real r2 = m_g(p);
        return std::min(m_f(p), -m_g(p));
    }

protected:

    F m_f;
    G m_g;
};

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
