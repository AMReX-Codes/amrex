#include <PlasmaInjector.H>

#include <iostream>

using namespace amrex;

///
/// This "custom" momentum distribution just does 0 momentum
///
void CustomMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
  u[0] = 0;
  u[1] = 0;
  u[2] = 0;
}
