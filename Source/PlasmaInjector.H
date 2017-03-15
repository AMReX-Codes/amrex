#ifndef PLASMA_INJECTOR_H_
#define PLASMA_INJECTOR_H_

#include "AMReX_Real.H"
#include "AMReX_ParmParse.H"

class PlasmaInjector
{
public:
    PlasmaInjector(int ispecies, const std::string& name);

    virtual amrex::Real getDensity(amrex::Real x, amrex::Real y, amrex::Real z) = 0;

    bool insideBounds(amrex::Real x, amrex::Real y, amrex::Real z);

    int numParticlesPerCell() {return num_particles_per_cell;};

    amrex::Real getGamma() {return gamma;}

protected:

    amrex::Real xmin, xmax;
    amrex::Real ymin, ymax;
    amrex::Real zmin, zmax;
    
    amrex::Real density;
    int num_particles_per_cell;

    amrex::Real gamma;

    int species_id;
    std::string species_name;
};

class ConstantPlasmaInjector : public PlasmaInjector {
public:
    ConstantPlasmaInjector(int ispecies, const std::string& name);
    amrex::Real getDensity( amrex::Real x, amrex::Real y, amrex::Real z);
};

class CustomPlasmaInjector : public PlasmaInjector {
public:
    CustomPlasmaInjector(int ispecies, const std::string& name);
    amrex::Real getDensity( amrex::Real x, amrex::Real y, amrex::Real z);
};

class DoubleRampPlasmaInjector : public PlasmaInjector {
public:
    DoubleRampPlasmaInjector(int ispecies, const std::string& name);
    amrex::Real getDensity( amrex::Real x, amrex::Real y, amrex::Real z);
    
protected:
    amrex::Real plateau_length;
    amrex::Real ramp_length;
};

#endif
