#ifndef PLASMA_INJECTOR_H_
#define PLASMA_INJECTOR_H_

#include "Real.H"
#include "ParmParse.H"

class PlasmaInjector
{
public:
  PlasmaInjector();

  virtual Real getDensity(Real x, Real y, Real z) = 0;

  bool insideBounds(Real x, Real y, Real z);

protected:

  Real xmin, xmax;
  Real ymin, ymax;
  Real zmin, zmax;

  Real density;

};

class ConstantPlasmaInjector : public PlasmaInjector {
public:
  Real getDensity( Real x, Real y, Real z);
};

class CustomPlasmaInjector : public PlasmaInjector {
public:
  Real getDensity( Real x, Real y, Real z);
};

class DoubleRampPlasmaInjector : public PlasmaInjector {
public:
  DoubleRampPlasmaInjector();
  Real getDensity( Real x, Real y, Real z);

protected:
  Real plateau_length;
  Real ramp_length;
};

#endif
