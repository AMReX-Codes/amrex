#include <PlasmaInjector.H>

#include <iostream>

Real CustomPlasmaInjector::getDensity(Real x, Real y, Real z) {
  static_assert(false,
		"If running with a custom density profile, you must supply a CustomDensityProb.cpp file");
  return 0.0;
}
