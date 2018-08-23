#ifndef AMREX_EB2_IF_PLANE_H_
#define AMREX_EB2_IF_PLANE_H_

#include <AMReX_Array.H>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class PlaneIF
{
public:

    PlaneIF (const RealArray& a_point, const RealArray& a_normal, bool a_inside = true)
        : m_point(a_point), m_normal(a_normal), m_inside(a_inside),
          m_sign( a_inside ? 1.0 : -1.0 )
    {}

    ~PlaneIF () {}

    PlaneIF (const PlaneIF& rhs) noexcept = default;
    PlaneIF (PlaneIF&& rhs) noexcept = default;
    PlaneIF& operator= (const PlaneIF& rhs) = delete;
    PlaneIF& operator= (PlaneIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        return AMREX_D_TERM( (p[0]-m_point[0])*m_normal[0]*m_sign,
                            +(p[1]-m_point[1])*m_normal[1]*m_sign,
                            +(p[2]-m_point[2])*m_normal[2]*m_sign );
    }

protected:

    RealArray m_point;
    RealArray m_normal;

    bool      m_inside;
    Real      m_sign;

};

}}

#endif
