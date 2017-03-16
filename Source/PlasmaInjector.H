#ifndef PLASMA_INJECTOR_H_
#define PLASMA_INJECTOR_H_

#include <random>

#include "AMReX_Real.H"
#include "AMReX_ParmParse.H"

class PlasmaMomentumDistribution
{
public:
    virtual void getMomentum(amrex::Real* u) = 0;
};

class ConstantMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    ConstantMomentumDistribution(amrex::Real ux,
                                 amrex::Real uy,
                                 amrex::Real uz);
    void getMomentum(amrex::Real* u);
    
private:
    amrex::Real _ux;
    amrex::Real _uy;
    amrex::Real _uz;
};

class GaussianRandomMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    GaussianRandomMomentumDistribution(amrex::Real ux_m,
                                       amrex::Real uy_m, 
                                       amrex::Real uz_m,
                                       amrex::Real u_th);
    void getMomentum(amrex::Real* u);
private:
    amrex::Real _ux_m;
    amrex::Real _uy_m;
    amrex::Real _uz_m;    
    amrex::Real _u_th;    

    std::normal_distribution<amrex::Real> momentum_distribution;
    std::default_random_engine generator;
};

class PlasmaInjector
{

public:
    PlasmaInjector(int ispecies, const std::string& name);

    virtual amrex::Real getDensity(amrex::Real x, amrex::Real y, amrex::Real z) = 0;

    bool insideBounds(amrex::Real x, amrex::Real y, amrex::Real z);

    int numParticlesPerCell() {return num_particles_per_cell;};

    void getMomentum(amrex::Real* u);

    amrex::Real getCharge() {return charge;}
    amrex::Real getMass() {return mass;}

protected:

    amrex::Real mass, charge;

    amrex::Real xmin, xmax;
    amrex::Real ymin, ymax;
    amrex::Real zmin, zmax;
    
    amrex::Real density;
    int num_particles_per_cell;

    int species_id;
    std::string species_name;

    std::unique_ptr<PlasmaMomentumDistribution> mom_dist;
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
