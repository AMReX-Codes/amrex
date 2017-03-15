#include "PlasmaInjector.H"

#include <sstream>

using namespace amrex;

PlasmaInjector::PlasmaInjector(int ispecies, const std::string& name) 
    : species_id(ispecies), species_name(name)
{
    ParmParse pp(species_name);

    pp.get("xmin", xmin);
    pp.get("ymin", ymin);
    pp.get("zmin", zmin);
    pp.get("xmax", xmax);
    pp.get("ymax", ymax);
    pp.get("zmax", zmax);
    
    pp.get("num_particles_per_cell", num_particles_per_cell);
    pp.get("density", density);
    pp.get("gamma", gamma);
}

bool PlasmaInjector::insideBounds(Real x, Real y, Real z) {
  if (x >= xmax || x < xmin ||
      y >= ymax || y < ymin ||
      z >= zmax || z < zmin ) return false;
  return true;
}

ConstantPlasmaInjector::ConstantPlasmaInjector(int ispecies, const std::string& name) 
    : PlasmaInjector(ispecies, name)
{
}

Real ConstantPlasmaInjector::getDensity(Real x, Real y, Real z) {
  return density;
}

CustomPlasmaInjector::CustomPlasmaInjector(int ispecies, const std::string& name) 
    : PlasmaInjector(ispecies, name)
{
}

DoubleRampPlasmaInjector::DoubleRampPlasmaInjector(int ispecies, const std::string& name) 
    : PlasmaInjector(ispecies, name)
{
  ParmParse pp(species_name);
  pp.get("ramp_length", ramp_length);
  pp.get("plateau_length", plateau_length);
}

Real DoubleRampPlasmaInjector::getDensity(Real x, Real y, Real z) {
  return density;
}


