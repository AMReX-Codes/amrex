#ifndef WARPX_CONST_H_
#define WARPX_CONST_H_

#include <AMReX_REAL.H>

// Math constants
namespace MathConst
{
    static constexpr amrex::Real pi = 3.14159265358979323846;
}

// Physical constants
namespace PhysConst
{
    static constexpr amrex::Real c     = 299792458.;
    static constexpr amrex::Real ep0   = 8.854187817e-12;
    static constexpr amrex::Real mu0   = 1.2566370614359173e-06;
    static constexpr amrex::Real q_e   = 1.6021764620000001e-19;
    static constexpr amrex::Real m_e   = 9.10938291e-31;
    static constexpr amrex::Real m_p   = 1.6726231000000001e-27;
    static constexpr amrex::Real hbar  = 1.05457180010e-34;
    static constexpr amrex::Real alpha = mu0/(4*MathConst::pi)*q_e*q_e*c/hbar;
    static constexpr amrex::Real r_e   = 1./(4*MathConst::pi*ep0) * q_e*q_e/(m_e*c*c);
}

#endif
