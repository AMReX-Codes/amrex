#ifndef MY_EB_H_
#define MY_EB_H_

#include <AMReX_Array.H>
#include <cmath>
#include <algorithm>

// For all implicit functions, >0: body; =0: boundary; <0: fluid

class FlowerIF
{
public:

    FlowerIF (amrex::Real a_radius, amrex::Real a_delta, int a_npetals,
              const amrex::RealArray& a_center, bool a_inside)
        : m_r(a_radius), m_dr(a_delta), m_npetals(a_npetals), m_center(a_center), m_inside(a_inside),
          m_sign(a_inside ? 1.0 : -1.0)
        {}

    ~FlowerIF () {}
    
    FlowerIF (const FlowerIF& rhs) noexcept = default;
    FlowerIF (FlowerIF&& rhs) noexcept = default;
    FlowerIF& operator= (const FlowerIF& rhs) = delete;
    FlowerIF& operator= (FlowerIF&& rhs) = delete;

    amrex::Real operator() (const amrex::RealArray& p) const
    {
        amrex::RealArray pos{AMREX_D_DECL(p[0]-m_center[0],
                                   p[1]-m_center[1],
                                   p[2]-m_center[2])};
        amrex::Real r = std::hypot(pos[0],pos[1]);
        amrex::Real theta = std::atan2(pos[1], pos[0]);
        return m_sign*(r - m_r - m_dr * std::cos(m_npetals*theta));
    }

protected:
    amrex::Real      m_r;
    amrex::Real      m_dr;
    amrex::Real      m_npetals;
    amrex::RealArray m_center;
    bool             m_inside;
    //
    amrex::Real      m_sign;
};

#endif
