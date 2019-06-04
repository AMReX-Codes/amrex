#ifndef PLASMA_INJECTOR_H_
#define PLASMA_INJECTOR_H_

#include <array>

#include "AMReX_REAL.H"
#include <AMReX_Vector.H>
#include <WarpXConst.H>
#include <WarpXParser.H>
#include "AMReX_ParmParse.H"
#include "AMReX_Utility.H"

enum class predefined_profile_flag { Null, parabolic_channel };

///
/// PlasmaDensityProfile describes how the charge density
/// is set in particle initialization. Subclasses must define a
/// getDensity function that describes the charge density as a
/// function of x, y, and z.
///
class PlasmaDensityProfile
{
public:
    virtual ~PlasmaDensityProfile() {};
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const = 0;
protected:
    std::string _species_name;
};

///
/// This describes a constant density distribution.
///
class ConstantDensityProfile : public PlasmaDensityProfile
{
public:
    ConstantDensityProfile(amrex::Real _density);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;

private:
    amrex::Real _density;
};

///
/// This describes a custom density distribution. Users can supply
/// in their problem directory.
///
///
class CustomDensityProfile : public PlasmaDensityProfile
{
public:
    CustomDensityProfile(const std::string& species_name);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;
private:
    amrex::Vector<amrex::Real> params;
};

///
/// This describes predefined density distributions.
///
class PredefinedDensityProfile : public PlasmaDensityProfile
{
public:
    PredefinedDensityProfile(const std::string& species_name);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;
    amrex::Real ParabolicChannel(amrex::Real x,
                                 amrex::Real y,
                                 amrex::Real z) const;
private:
    predefined_profile_flag which_profile = predefined_profile_flag::Null;
    amrex::Vector<amrex::Real> params;
};

///
/// This describes a density function parsed in the input file. 
///
class ParseDensityProfile : public PlasmaDensityProfile
{
public:
    ParseDensityProfile(const std::string _parse_density_function);
    virtual amrex::Real getDensity(amrex::Real x,
                                   amrex::Real y,
                                   amrex::Real z) const override;
private:
    std::string _parse_density_function;
    WarpXParser parser_density;
};

///
/// PlasmaMomentumDistribution describes how the particle momenta
/// are set. Subclasses must define a getMomentum method that fills
/// a u with the 3 components of the particle momentum
///
class PlasmaMomentumDistribution
{
public:
    using vec3 = std::array<amrex::Real, 3>;
    virtual ~PlasmaMomentumDistribution() {};
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) = 0;
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
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;

private:
    amrex::Real _ux;
    amrex::Real _uy;
    amrex::Real _uz;
};

///
/// This describes a custom momentum distribution. Users can supply
/// in their problem directory.
///
///
class CustomMomentumDistribution : public PlasmaMomentumDistribution
{
public:
    CustomMomentumDistribution(const std::string& species_name);
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;

private:
    amrex::Vector<amrex::Real> params;
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
                                       amrex::Real ux_th,
                                       amrex::Real uy_th,
                                       amrex::Real uz_th);
    virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;
private:
    amrex::Real _ux_m;
    amrex::Real _uy_m;
    amrex::Real _uz_m;
    amrex::Real _ux_th;
    amrex::Real _uy_th;
    amrex::Real _uz_th;
};

///
/// This is a radially expanding momentum distribution
/// Particles will have a radial momentum proportional to their 
/// radius, with proportionality constant u_over_r
class RadialExpansionMomentumDistribution : public PlasmaMomentumDistribution
{
public:
  RadialExpansionMomentumDistribution( amrex::Real u_over_r );
  virtual void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z) override;
private:
    amrex::Real _u_over_r;
};

///
/// This describes a momentum distribution function parsed in the input file. 
///
class ParseMomentumFunction : public PlasmaMomentumDistribution
{
public:
    ParseMomentumFunction(const std::string _parse_momentum_function_ux,
                          const std::string _parse_momentum_function_uy,
                          const std::string _parse_momentum_function_uz);
    virtual void getMomentum(vec3& u, 
                             amrex::Real x,
                             amrex::Real y,
                             amrex::Real z) override;
private:
    std::string _parse_momentum_function_ux;
    std::string _parse_momentum_function_uy;
    std::string _parse_momentum_function_uz;
    WarpXParser parser_ux;
    WarpXParser parser_uy;
    WarpXParser parser_uz;
};


///
/// PlasmaParticlePosition describes how particles are initialized
/// into each cell box. Subclasses must define a
/// getPositionUnitBox function that returns the position of
/// particle number i_part in a unitary box.
///
class PlasmaParticlePosition{
public:
  using vec3 = std::array<amrex::Real, 3>;
  virtual ~PlasmaParticlePosition() {};
    virtual void getPositionUnitBox(vec3& r, int i_part, int ref_fac=1) = 0;
};

///
/// Particles are initialized with a random uniform
/// distribution inside each cell
///
class RandomPosition : public PlasmaParticlePosition{
public:
    RandomPosition(int num_particles_per_cell);
    virtual void getPositionUnitBox(vec3& r, int i_part, int ref_fac=1) override;
private:
    amrex::Real _x;
    amrex::Real _y;
    amrex::Real _z;
    int _num_particles_per_cell;
};

///
/// Particles are regularly distributed inside each cell. The user provides
/// a 3d (resp. 2d) vector num_particles_per_cell_each_dim that contains
/// the number of particles per cell along each dimension.
///
class RegularPosition : public PlasmaParticlePosition{
public:
  RegularPosition(const amrex::Vector<int>& num_particles_per_cell_each_dim);
    virtual void getPositionUnitBox(vec3& r, int i_part, int ref_fac=1) override;
private:
  amrex::Real _x;
  amrex::Real _y;
  amrex::Real _z;
  amrex::Vector<int> _num_particles_per_cell_each_dim;
};

///
/// The PlasmaInjector class parses and stores information about the plasma
/// type used in the particle container. This information is used to create the
/// particles on initialization and whenever the window moves.
///
class PlasmaInjector
{

public:

    using vec3 = std::array<amrex::Real, 3>;

    PlasmaInjector();

    PlasmaInjector(int ispecies, const std::string& name);

    amrex::Real getDensity(amrex::Real x, amrex::Real y, amrex::Real z);

    bool insideBounds(amrex::Real x, amrex::Real y, amrex::Real z);

    int num_particles_per_cell;

    amrex::Vector<int> num_particles_per_cell_each_dim;

    void getMomentum(vec3& u, amrex::Real x, amrex::Real y, amrex::Real z);

    void getPositionUnitBox(vec3& r, int i_part, int ref_fac=1);

    amrex::Real getCharge() {return charge;}
    amrex::Real getMass() {return mass;}

    bool doInjection() { return part_pos != NULL;}

    bool add_single_particle = false;
    amrex::Vector<amrex::Real> single_particle_pos;
    amrex::Vector<amrex::Real> single_particle_vel;
    amrex::Real single_particle_weight;

    bool gaussian_beam = false;
    amrex::Real x_m;
    amrex::Real y_m;
    amrex::Real z_m;
    amrex::Real x_rms;
    amrex::Real y_rms;
    amrex::Real z_rms;
    amrex::Real q_tot;
    long npart;
    int do_symmetrize = 0;

    bool radially_weighted = true;

    std::string str_density_function;
    std::string str_momentum_function_ux;
    std::string str_momentum_function_uy;
    std::string str_momentum_function_uz;

    amrex::Real xmin, xmax;
    amrex::Real ymin, ymax;
    amrex::Real zmin, zmax;

protected:

    amrex::Real mass, charge;

    amrex::Real density;

    int species_id;
    std::string species_name;

    std::unique_ptr<PlasmaDensityProfile> rho_prof;
    std::unique_ptr<PlasmaMomentumDistribution> mom_dist;
    std::unique_ptr<PlasmaParticlePosition> part_pos;
    
    void parseDensity(amrex::ParmParse pp);
    void parseMomentum(amrex::ParmParse pp);
    
};

#endif
