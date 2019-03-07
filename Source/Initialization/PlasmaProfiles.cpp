#include <PlasmaInjector.H>
#include <cmath>
#include <iostream>
#include <WarpXConst.H>

using namespace amrex;

Real PredefinedDensityProfile::getDensity(Real x, Real y, Real z) const {
    Real n;
    if ( which_profile == "parabolic_channel" ) {
        n = ParabolicChannel(x,y,z);
    }
    return n;
}

///
/// plateau between linear upramp and downramp, and parab transverse profile
///
Real PredefinedDensityProfile::ParabolicChannel(Real x, Real y, Real z) const {
    //  params = [ramp_up   plateau   ramp_down   rc       n0]
    Real ramp_up   = params[0];
    Real plateau   = params[1];
    Real ramp_down = params[2];
    Real rc        = params[3];
    Real n0        = params[4];
    Real n;
    Real kp = PhysConst::q_e/PhysConst::c*sqrt( n0/(PhysConst::m_e*PhysConst::ep0) );

    if        (z>=0               and z<ramp_up                  ) {
        n = z/ramp_up;
    } else if (z>=ramp_up         and z<ramp_up+plateau          ) {
        n = 1;
    } else if (z>=ramp_up+plateau and z<ramp_up+plateau+ramp_down) {
        n = 1-(z-ramp_up-plateau)/ramp_down;
    } else {
        n = 0;
    }
    n *= n0*(1+4*(x*x+y*y)/(kp*kp*std::pow(rc,4)));
    return n;
}
