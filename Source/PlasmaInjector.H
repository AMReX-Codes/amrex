#ifndef PLASMA_INJECTOR_H_
#define PLASMA_INJECTOR_H_

#include <random>

#include "AMReX_Real.H"
#include "AMReX_ParmParse.H"

///
/// PlasmaDensityDistribution describes how the charge density
/// is set in particle initialization. Subclasses must define a
/// getDensity function that describes the charge density as a
/// function of x, y, and z.
///
class PlasmaDensityDistribution
{
public:
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const = 0;
};

///
/// This describes a constant density distribution.
///
class ConstantDensityDistribution : public PlasmaDensityDistribution
{
public:
    ConstantDensityDistribution(amrex::Real density);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;

private:
    amrex::Real _density;
};

///
/// This describes a custom density distribution. Users can supply
/// a getDensity function.
///
///
class CustomDensityDistribution : public PlasmaDensityDistribution
{
public:
    CustomDensityDistribution(amrex::Real density);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;
};

///
/// PlasmaMomentumDistribution describes how the particle momenta
/// are set. Subclasses must define a getMomentum method that fills
/// a u with the 3 components of the particle momentum
///
class PlasmaMomentumDistribution
{
public:
    virtual void getMomentum(amrex::Real* u) = 0;
};

///
/// This is a constant momentum distribution - all particles will
/// have the same ux, uy, and uz
///
class ConstantMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    ConstantMomentumDistribution(amrex::Real ux,
                                 amrex::Real uy,
                                 amrex::Real uz);
    virtual void getMomentum(amrex::Real* u) override;
    
private:
    amrex::Real _ux;
    amrex::Real _uy;
    amrex::Real _uz;
};

///
/// This is a Gaussian Random momentum distribution. 
/// Particles will get random momenta, drawn from a normal.
/// ux_m, ux_y, and ux_z describe the mean components in the x, y, and z
/// directions, while u_th is the standard deviation of the random 
/// component.
///
class GaussianRandomMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    GaussianRandomMomentumDistribution(amrex::Real ux_m,
                                       amrex::Real uy_m, 
                                       amrex::Real uz_m,
                                       amrex::Real u_th);
    virtual void getMomentum(amrex::Real* u) override;
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

    amrex::Real getDensity(amrex::Real x, amrex::Real y, amrex::Real z);

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
    std::string injection_style;
    
    std::unique_ptr<PlasmaDensityDistribution>  rho_dist;
    std::unique_ptr<PlasmaMomentumDistribution> mom_dist;
};

#endif
