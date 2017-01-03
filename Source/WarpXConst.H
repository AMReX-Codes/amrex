#ifndef WARPX_CONST_H_
#define WARPX_CONST_H_

#include <AMReX_REAL.H>

// Physical constants
namespace PhysConst
{
    constexpr Real c   = 299792458.;
    constexpr Real mu0 = 1.2566370614359173e-06;
    constexpr Real q_e = 1.6021764620000001e-19;
    constexpr Real m_e = 9.10938291e-31;
    constexpr Real m_p = 1.6726231000000001e-27;
};

#endif
