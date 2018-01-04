#include <PlasmaInjector.H>

using namespace amrex;

///
/// This "custom" density profile just does constant
///
Real CustomDensityProfile::getDensity(Real x, Real y, Real z) const {
  const Real on_axis_density = params[0];
  const Real plasma_zmin = params[1];
  const Real plasma_zmax = params[2];
  const Real plasma_lramp_start = params[3];
  const Real plasma_lramp_end = params[4];
  const Real plasma_rcap = params[5];
  const Real plasma_rdownramp = params[6];
  const Real plasma_rchannel = params[7];
  static const Real re = 2.8178403227e-15; // Electron classical radius
  static const Real pi = 3.14159265359;

  Real r2 = x*x + y*y;
  Real r = std::sqrt( r2 );

  // Transverse part of the profile
  Real nr;
  if (r<plasma_rcap) {
      nr = 1. + 1./(pi*on_axis_density*re) * r2/pow(plasma_rchannel, 4);
  } else {
      nr = 1. + 1./(pi*on_axis_density*re) *
            pow(plasma_rcap, 2)/pow(plasma_rchannel, 4) *
            (plasma_rcap+plasma_rdownramp-r)/plasma_rdownramp;
  }
  // Longitudinal part of the profile
  Real nz;
  if (z<plasma_zmin) {
      nz = 0;
  } else if (z<plasma_zmin+plasma_lramp_start) {
      nz = (z-plasma_zmin)/plasma_lramp_start;
  } else if (z<plasma_zmax-plasma_lramp_end) {
      nz = 1.;
  } else if (z<plasma_zmax){
      nz = -(z-plasma_zmax)/plasma_lramp_end;
  } else {
      nz = 0;
  }
  // Combine and saturate profile
  Real n = nr*nz;
  if (n > 4.) {
      n = 4.;
  }

  return on_axis_density*n;
}
