#ifndef AMREX_EB2_IF_TORUS_H_
#define AMREX_EB2_IF_TORUS_H_

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <cmath>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

namespace amrex { namespace EB2 {

class TorusIF
    : GPUable
{
public:

    // inside: is the fluid inside the sphere?
    TorusIF (Real a_large_radius, Real a_small_radius, const RealArray& a_center, bool a_inside)
        : m_large_radius(a_large_radius),
          m_small_radius(a_small_radius),
          m_center(makeXDim3(a_center)),
          m_sign( a_inside ? 1.0 : -1.0 )
        {}

    TorusIF (const TorusIF& rhs) noexcept = default;
    TorusIF (TorusIF&& rhs) noexcept = default;
    TorusIF& operator= (const TorusIF& rhs) = delete;
    TorusIF& operator= (TorusIF&& rhs) = delete;

    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept {
        Real d = std::hypot(x-m_center.x, y-m_center.y);
#if (AMREX_SPACEDIM == 2)
        return m_sign*((m_large_radius-d)*(m_large_radius-d)
                       - m_small_radius*m_small_radius);
#else
        return m_sign*((m_large_radius-d)*(m_large_radius-d)
                       +(z-m_center.z)*(z-m_center.z)
                       - m_small_radius*m_small_radius);
#endif
    }

    inline Real operator() (const RealArray& p) const noexcept {
        return this->operator()(AMREX_D_DECL(p[0],p[1],p[2]));
    }

protected:

    Real      m_large_radius;
    Real      m_small_radius;
    XDim3     m_center;
    //
    Real      m_small_radius2;  // precomputed radius squared
    Real      m_sign;
};

}}

#endif
