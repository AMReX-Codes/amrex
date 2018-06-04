#ifndef AMREX_EB2_IF_PLANE_H_
#define AMREX_EB2_IF_PLANE_H_

#include <AMReX_Array.H>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class PlaneIF
{
public:

    PlaneIF (const RealArray& a_point, const RealArray& a_normal)
        : m_point(a_point), m_normal(a_normal)
    {}

    ~PlaneIF () {}

    PlaneIF (const PlaneIF& rhs) noexcept = default;
    PlaneIF (PlaneIF&& rhs) noexcept = default;
    PlaneIF& operator= (const PlaneIF& rhs) = delete;
    PlaneIF& operator= (PlaneIF&& rhs) = delete;

    Real operator() (const RealArray& p) const
    {
        return AMREX_D_TERM( (p[0]-m_point[0])*m_normal[0],
                            +(p[1]-m_point[1])*m_normal[1],
                            +(p[2]-m_point[2])*m_normal[2] );
    }

protected:

    RealArray m_point;
    RealArray m_normal;

};

}}

#endif
