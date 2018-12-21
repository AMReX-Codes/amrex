#ifndef AMREX_EB2_IF_TORUS_H_
#define AMREX_EB2_IF_TORUS_H_

#include <AMReX_Array.H>
#include <cmath>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class TorusIF
{
public:

    // inside: is the fluid inside the sphere?
    TorusIF (Real a_large_radius, Real a_small_radius, const RealArray& a_center, bool a_inside)
        : m_large_radius(a_large_radius),
          m_small_radius(a_small_radius),
          m_center(a_center), m_inside(a_inside),
          m_large_radius2(a_large_radius*a_large_radius),
          m_small_radius2(a_small_radius*a_small_radius),
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    ~TorusIF () {}

    TorusIF (const TorusIF& rhs) noexcept = default;
    TorusIF (TorusIF&& rhs) noexcept = default;
    TorusIF& operator= (const TorusIF& rhs) = delete;
    TorusIF& operator= (TorusIF&& rhs) = delete;

    Real operator() (const RealArray& p) const {
        Real d2 = std::sqrt((p[0]-m_center[0])*(p[0]-m_center[0]) +
                           +(p[1]-m_center[1])*(p[1]-m_center[1]));
#if (AMREX_SPACEDIM == 2)
        return m_sign*((m_large_radius-d2)*(m_large_radius-d2)
                      - m_small_radius2);
#else
        return m_sign*((m_large_radius-d2)*(m_large_radius-d2)
                      +(p[2]-m_center[2])*(p[2]-m_center[2])
                      - m_small_radius2);
#endif
    }

protected:

    Real      m_large_radius;
    Real      m_small_radius;
    RealArray m_center;
    bool      m_inside;
    //
    Real      m_large_radius2;  // precomputed radius squared
    Real      m_small_radius2;  // precomputed radius squared
    Real      m_sign;
};

}}

#endif
