#ifndef AMREX_EB2_IF_ELLIPSOID_H_
#define AMREX_EB2_IF_ELLIPSOID_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class EllipsoidIF
    : public GPUable
{
public:

    // inside: is the fluid inside the ellipsoid?
    EllipsoidIF (const RealArray& a_radii, const RealArray& a_center, bool a_inside)
        : m_radii(makeXDim3(a_radii)),
          m_center(makeXDim3(a_center)),
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    EllipsoidIF (const EllipsoidIF& rhs) noexcept = default;
    EllipsoidIF (EllipsoidIF&& rhs) noexcept = default;
    EllipsoidIF& operator= (const EllipsoidIF& rhs) = delete;
    EllipsoidIF& operator= (EllipsoidIF&& rhs) = delete;


    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept {
        Real d2 = AMREX_D_TERM(  (x-m_center.x)*(x-m_center.x) / (m_radii.x*m_radii.x),
                               + (y-m_center.y)*(y-m_center.y) / (m_radii.y*m_radii.y),
                               + (z-m_center.z)*(z-m_center.z) / (m_radii.z*m_radii.z));
        return m_sign*(d2-1.0);
    }

    inline Real operator() (const RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

protected:

    XDim3 m_radii;
    XDim3 m_center;
    //
    Real  m_sign;
};

}}

#endif
