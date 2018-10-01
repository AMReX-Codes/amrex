#ifndef AMREX_EB2_IF_LATHE_H_
#define AMREX_EB2_IF_LATHE_H_

#include <AMReX_Array.H>

#include <type_traits>
#include <cmath>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F>
class LatheIF
{
public:

    LatheIF (F&& a_f) : m_f(std::move(a_f)) {}
    LatheIF (F const& a_f) : m_f(a_f) {}
          
    ~LatheIF () {}

    LatheIF (const LatheIF& rhs) = default;
    LatheIF (LatheIF&& rhs) = default;
    LatheIF& operator= (const LatheIF& rhs) = delete;
    LatheIF& operator= (LatheIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        Real r = std::hypot(p[0],p[1]);
#if (AMREX_SPACEDIM == 2)
        return m_f({r,0.0});
#else
        return m_f({r,p[2],0.0});
#endif
    }

protected:

    F m_f;
};

template <class F>
constexpr LatheIF<typename std::decay<F>::type>
lathe (F&& f)
{
    return LatheIF<typename std::decay<F>::type>(std::forward<F>(f));
}

}}

#endif
