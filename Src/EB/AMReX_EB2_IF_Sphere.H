#ifndef AMREX_EB2_IF_SPHERE_H_
#define AMREX_EB2_IF_SPHERE_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class SphereIF
    : public GPUable
{
public:

    // inside: is the fluid inside the sphere?
    SphereIF (Real a_radius, const RealArray& a_center, bool a_inside)
        : m_radius(a_radius),
          m_center(makeXDim3(a_center)),
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    SphereIF (const SphereIF& rhs) noexcept = default;
    SphereIF (SphereIF&& rhs) noexcept = default;
    SphereIF& operator= (const SphereIF& rhs) = delete;
    SphereIF& operator= (SphereIF&& rhs) = delete;

    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept {
        Real d2 = AMREX_D_TERM(  (x-m_center.x)*(x-m_center.x),
                               + (y-m_center.y)*(y-m_center.y),
                               + (z-m_center.z)*(z-m_center.z));
        return m_sign*(d2-m_radius*m_radius);
    }

    inline Real operator() (const RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

protected:
  
    Real  m_radius;
    XDim3 m_center;
    //
    Real  m_sign;
};

}}

#endif
