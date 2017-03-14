#include "PlasmaInjector.H"

#include <iostream>

using namespace amrex;

PlasmaInjector::PlasmaInjector() {
  ParmParse pp("plasma");

  pp.get("xmin", xmin);
  pp.get("ymin", ymin);
  pp.get("zmin", zmin);
  pp.get("xmax", xmax);
  pp.get("ymax", ymax);
  pp.get("zmax", zmax);

  pp.get("density", density);
}

bool PlasmaInjector::insideBounds(Real x, Real y, Real z) {
  if (x >= xmax || x < xmin ||
      y >= ymax || y < ymin ||
      z >= zmax || z < zmin ) return true;
  return false;
}

Real ConstantPlasmaInjector::getDensity(Real x, Real y, Real z) {
  return density;
}

DoubleRampPlasmaInjector::DoubleRampPlasmaInjector() 
  : PlasmaInjector()
{
  ParmParse pp("plasma");
  pp.get("ramp_length", ramp_length);
  pp.get("plateau_length", plateau_length);
}

Real DoubleRampPlasmaInjector::getDensity(Real x, Real y, Real z) {
  return density;
}

