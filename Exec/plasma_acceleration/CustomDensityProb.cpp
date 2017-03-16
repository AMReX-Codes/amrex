#include <PlasmaInjector.H>

#include <iostream>

using namespace amrex;

///
/// This "custom" density profile just does constant
///
Real CustomDensityProfile::getDensity(Real x, Real y, Real z) const {
  std::cout << "Using custom getDensity" << std::endl;
  return params[0];
}
