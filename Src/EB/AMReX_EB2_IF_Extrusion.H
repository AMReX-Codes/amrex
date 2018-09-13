#ifndef AMREX_EB2_IF_EXTRUSION_H_
#define AMREX_EB2_IF_EXTRUSION_H_

#include <AMReX_Array.H>

#include <type_traits>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

template <class F>
class ExtrusionIF
{
public:

    ExtrusionIF (F&& a_f, int direction)
        : m_f(std::move(a_f)),
          m_direction(direction)
        {}
    ExtrusionIF (F const& a_f, int direction)
        : m_f(a_f),
          m_direction(direction)
    {}

    ~ExtrusionIF () {}

    ExtrusionIF (const ExtrusionIF& rhs) = default;
    ExtrusionIF (ExtrusionIF&& rhs) = default;
    ExtrusionIF& operator= (const ExtrusionIF& rhs) = delete;
    ExtrusionIF& operator= (ExtrusionIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        RealArray x = p;
        x[m_direction] = 0.0;
        return m_f(x);
    }

protected:

    F m_f;
    int m_direction;
};

template <class F>
constexpr ExtrusionIF<typename std::decay<F>::type>
extrude (F&&f, int direction)
{
    return ExtrusionIF<typename std::decay<F>::type>(std::forward<F>(f),direction);
}

}}

#endif
