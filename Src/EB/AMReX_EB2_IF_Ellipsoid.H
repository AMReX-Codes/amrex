#ifndef AMREX_EB2_IF_ELLIPSOID_H_
#define AMREX_EB2_IF_ELLIPSOID_H_

#include <AMReX_Array.H>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class EllipsoidIF
{
public:

    // inside: is the fluid inside the ellipsoid?
    EllipsoidIF (const RealArray& a_radii, const RealArray& a_center, bool a_inside)
        : m_radii(a_radii), m_center(a_center), m_inside(a_inside),
          m_radii2_inv{AMREX_D_DECL(1.0/(a_radii[0]*a_radii[0]),
                                    1.0/(a_radii[1]*a_radii[1]),
                                    1.0/(a_radii[2]*a_radii[2]))},
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    ~EllipsoidIF () {}

    EllipsoidIF (const EllipsoidIF& rhs) noexcept = default;
    EllipsoidIF (EllipsoidIF&& rhs) noexcept = default;
    EllipsoidIF& operator= (const EllipsoidIF& rhs) = delete;
    EllipsoidIF& operator= (EllipsoidIF&& rhs) = delete;

    Real operator() (const RealArray& p) const {
        Real d2 = AMREX_D_TERM(  (p[0]-m_center[0])*(p[0]-m_center[0]) * m_radii2_inv[0],
                               + (p[1]-m_center[1])*(p[1]-m_center[1]) * m_radii2_inv[1],
                               + (p[2]-m_center[2])*(p[2]-m_center[2]) * m_radii2_inv[2]);
        return m_sign*(d2-1.0);
    }

protected:
  
    RealArray m_radii;
    RealArray m_center;
    bool      m_inside;
    //
    RealArray m_radii2_inv;  // precomputed radii squared
    Real      m_sign;
};

}}

#endif
