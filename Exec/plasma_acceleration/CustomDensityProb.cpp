#include <PlasmaInjector.H>

#include <iostream>

using namespace amrex;

///
/// This "custom" density profile just does constant
///
Real CustomDensityProfile::getDensity(Real x, Real y, Real z) const {
  return params[0];
}
