#ifndef AMREX_EB2_IF_PLANE_H_
#define AMREX_EB2_IF_PLANE_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class PlaneIF
    : GPUable
{
public:

    PlaneIF (const RealArray& a_point, const RealArray& a_normal, bool a_inside = true)
        : m_point(makeXDim3(a_point)),
          m_normal(makeXDim3(a_normal)),
          m_sign( a_inside ? 1.0 : -1.0 )
    {}

    PlaneIF (const PlaneIF& rhs) noexcept = default;
    PlaneIF (PlaneIF&& rhs) noexcept = default;
    PlaneIF& operator= (const PlaneIF& rhs) = delete;
    PlaneIF& operator= (PlaneIF&& rhs) = delete;

    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return AMREX_D_TERM( (x-m_point.x)*m_normal.x*m_sign,
                            +(y-m_point.y)*m_normal.y*m_sign,
                            +(z-m_point.z)*m_normal.z*m_sign );
    }

    inline Real operator() (const RealArray& p) const noexcept
    {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

protected:

    XDim3 m_point;
    XDim3 m_normal;
    //
    Real  m_sign;

};

}}

#endif
