#ifndef AMREX_EB2_SPHERE_IF_H_
#define AMREX_EB2_SPHERE_IF_H_

#include <AMReX_SPACE.H>
#include <AMReX_Array.H>

namespace amrex { namespace EB2 {

class SphereIF
{
public:

    // inside: is the fluid inside the sphere?
    SphereIF (Real a_radius, const RealArray& a_center, bool a_inside)
        : m_radius(a_radius), m_center(a_center), m_inside(a_inside),
          m_radius2(a_radius*a_radius),
          m_sign( a_inside ? -1.0 : 1.0 )
        {}

    ~SphereIF () {}

    SphereIF (const SphereIF& rhs) noexcept = default;
    SphereIF (SphereIF&& rhs) noexcept = default;
    SphereIF& operator= (const SphereIF& rhs) = delete;
    SphereIF& operator= (SphereIF&& rhs) = delete;

    // >0: body; =0: boundary; <0: fluid
    Real operator() (const RealArray& p) const {
        Real d2 = AMREX_D_TERM((  p[0]-m_center[0])*(p[0]-m_center[0]),
                               * (p[1]-m_center[1])*(p[1]-m_center[1]),
                               * (p[2]-m_center[2])*(p[2]-m_center[2]));
        return m_sign*(d2-m_radius2);
    }

protected:
  
    Real      m_radius;
    RealArray m_center;
    bool      m_inside;
    Real      m_radius2;  // precomputed radius squared
    Real      m_sign;
};

}}

#endif
