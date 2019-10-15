#ifndef MY_EB_H_
#define MY_EB_H_

#include <AMReX_Array.H>
#include <AMReX_EB2_IF_Base.H>

#include <cmath>
#include <algorithm>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class FlowerIF
    : amrex::GPUable
{
public:

    FlowerIF (amrex::Real a_radius, amrex::Real a_delta, int a_npetals,
              const amrex::RealArray& a_center, bool a_inside)
        : m_r(a_radius),
          m_dr(a_delta),
          m_npetals(a_npetals),
          m_center(amrex::makeXDim3(a_center)),
          m_inside(a_inside),
          m_sign(a_inside ? 1.0 : -1.0)
        {}

    ~FlowerIF () {}
    
    FlowerIF (const FlowerIF& rhs) noexcept = default;
    FlowerIF (FlowerIF&& rhs) noexcept = default;
    FlowerIF& operator= (const FlowerIF& rhs) = delete;
    FlowerIF& operator= (FlowerIF&& rhs) = delete;

    AMREX_GPU_HOST_DEVICE inline
    amrex::Real operator() (AMREX_D_DECL(amrex::Real x, amrex::Real y, amrex::Real z))
        const noexcept
    {
        amrex::Real posx = x - m_center.x;
        amrex::Real posy = y - m_center.y;
        amrex::Real r = std::hypot(posx, posy);
        amrex::Real theta = std::atan2(posy, posx);
        return m_sign*(r - m_r - m_dr * std::cos(m_npetals*theta));
    }

    inline amrex::Real operator() (const amrex::RealArray& p) const noexcept
    {
        return this->operator() (AMREX_D_DECL(p[0], p[1], p[2]));
    }

protected:
    amrex::Real      m_r;
    amrex::Real      m_dr;
    amrex::Real      m_npetals;
    amrex::XDim3     m_center;
    bool             m_inside;
    //
    amrex::Real      m_sign;
};

#endif
