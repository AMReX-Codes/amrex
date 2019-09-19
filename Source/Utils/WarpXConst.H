#ifndef WARPX_CONST_H_
#define WARPX_CONST_H_

#include <AMReX_REAL.H>

// Math constants
namespace MathConst
{
    static constexpr amrex::Real pi = 3.14159265358979323846;
}

// Physical constants. Values are the 2018 CODATA recommended values
// https://physics.nist.gov/cuu/Constants/index.html
namespace PhysConst
{
    static constexpr amrex::Real c     = 299792458.;
    static constexpr amrex::Real ep0   = 8.8541878128e-12;
    static constexpr amrex::Real mu0   = 1.25663706212e-06;
    static constexpr amrex::Real q_e   = 1.602176634e-19;
    static constexpr amrex::Real m_e   = 9.1093837015e-31;
    static constexpr amrex::Real m_p   = 1.67262192369e-27;
    static constexpr amrex::Real hbar  = 1.054571817e-34;
    static constexpr amrex::Real alpha = mu0/(4*MathConst::pi)*q_e*q_e*c/hbar;
    static constexpr amrex::Real r_e   = 1./(4*MathConst::pi*ep0) * q_e*q_e/(m_e*c*c);
}

#endif
