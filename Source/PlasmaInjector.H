#ifndef PLASMA_INJECTOR_H_
#define PLASMA_INJECTOR_H_

#include "AMReX_Real.H"
#include "AMReX_ParmParse.H"

class PlasmaInjector
{
public:
  PlasmaInjector();

  virtual amrex::Real getDensity(amrex::Real x, amrex::Real y, amrex::Real z) = 0;

  bool insideBounds(amrex::Real x, amrex::Real y, amrex::Real z);

protected:

  amrex::Real xmin, xmax;
  amrex::Real ymin, ymax;
  amrex::Real zmin, zmax;

  amrex::Real density;

};

class ConstantPlasmaInjector : public PlasmaInjector {
public:
  amrex::Real getDensity( amrex::Real x, amrex::Real y, amrex::Real z);
};

class CustomPlasmaInjector : public PlasmaInjector {
public:
  amrex::Real getDensity( amrex::Real x, amrex::Real y, amrex::Real z);
};

class DoubleRampPlasmaInjector : public PlasmaInjector {
public:
  DoubleRampPlasmaInjector();
  amrex::Real getDensity( amrex::Real x, amrex::Real y, amrex::Real z);

protected:
  amrex::Real plateau_length;
  amrex::Real ramp_length;
};

#endif
