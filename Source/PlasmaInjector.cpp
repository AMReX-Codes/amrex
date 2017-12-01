#include "PlasmaInjector.H"

#include <sstream>

#include <WarpXConst.H>
#include <AMReX.H>

using namespace amrex;

namespace {
    void StringParseAbortMessage(const std::string& var,
                                 const std::string& name) {
        std::stringstream stringstream;
        std::string string;
        stringstream << var << " string '" << name << "' not recognized.";
        string = stringstream.str();
        amrex::Abort(string.c_str());
    }

    Real parseChargeName(const ParmParse pp, const std::string& name) {
        Real result;
        if (name == "q_e") {
            return PhysConst::q_e;
        } else if (pp.query("charge", result)) {
            return result;
        } else {
            StringParseAbortMessage("Charge", name);
            return 0.0;
        }
    }

    Real parseChargeString(const ParmParse pp, const std::string& name) {
        if(name.substr(0, 1) == "-")
            return -1.0 * parseChargeName(pp, name.substr(1, name.size() - 1));
        return parseChargeName(pp, name);
    }

    Real parseMassString(const ParmParse pp, const std::string& name) {
        Real result;
        if (name == "m_e") {
            return PhysConst::m_e;
        } else if (name == "m_p"){
            return PhysConst::m_p;
        } else if (pp.query("mass", result)) {
            return result;
        } else {
            StringParseAbortMessage("Mass", name);
            return 0.0;
        }
    }
}

ConstantDensityProfile::ConstantDensityProfile(Real density)
    : _density(density)
{}

Real ConstantDensityProfile::getDensity(Real x, Real y, Real z) const
{
    return _density;
}

CustomDensityProfile::CustomDensityProfile(const std::string& species_name)
{
    ParmParse pp(species_name);
    pp.getarr("custom_profile_params", params);
}

ConstantMomentumDistribution::ConstantMomentumDistribution(Real ux,
                                                           Real uy,
                                                           Real uz)
    : _ux(ux), _uy(uy), _uz(uz)
{}

void ConstantMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
    u[0] = _ux;
    u[1] = _uy;
    u[2] = _uz;
}

GaussianRandomMomentumDistribution::GaussianRandomMomentumDistribution(Real ux_m,
                                                                       Real uy_m,
                                                                       Real uz_m,
                                                                       Real u_th)
    : _ux_m(ux_m), _uy_m(uy_m), _uz_m(uz_m), _u_th(u_th)
{
}

void GaussianRandomMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
    Real ux_th = amrex::RandomNormal(0.0, _u_th);
    Real uy_th = amrex::RandomNormal(0.0, _u_th);
    Real uz_th = amrex::RandomNormal(0.0, _u_th);

    u[0] = _ux_m + ux_th;
    u[1] = _uy_m + uy_th;
    u[2] = _uz_m + uz_th;
}
RadialExpansionMomentumDistribution::RadialExpansionMomentumDistribution(Real u_over_r) : _u_over_r( u_over_r )
{
}

void RadialExpansionMomentumDistribution::getMomentum(vec3& u, Real x, Real y, Real z) {
  u[0] = _u_over_r * x;
  u[1] = _u_over_r * y;
  u[2] = _u_over_r * z;
}

RandomPosition::RandomPosition(int num_particles_per_cell):
  _num_particles_per_cell(num_particles_per_cell)
{}

void RandomPosition::getPositionUnitBox(vec3& r, int i_part){
    r[0] = amrex::Random();
    r[1] = amrex::Random();
    r[2] = amrex::Random();
}

RegularPosition::RegularPosition(const amrex::Vector<int>& num_particles_per_cell_each_dim)
    : _num_particles_per_cell_each_dim(num_particles_per_cell_each_dim)
{}

void RegularPosition::getPositionUnitBox(vec3& r, int i_part)
{
  int nx = _num_particles_per_cell_each_dim[0];
  int ny = _num_particles_per_cell_each_dim[1];
  int nz = _num_particles_per_cell_each_dim[2];

  int ix_part = i_part/(ny * nz);
  int iy_part = (i_part % (ny * nz)) % ny;
  int iz_part = (i_part % (ny * nz)) / ny;

  r[0] = (0.5+ix_part)/nx;
  r[1] = (0.5+iy_part)/ny;
  r[2] = (0.5+iz_part)/nz;
}

PlasmaInjector::PlasmaInjector(int ispecies, const std::string& name)
    : species_id(ispecies), species_name(name)
{
    ParmParse pp(species_name);

    // parse charge and mass
    std::string charge_s;
    pp.get("charge", charge_s);
    std::transform(charge_s.begin(),
                   charge_s.end(),
                   charge_s.begin(),
                   ::tolower);
    charge = parseChargeString(pp, charge_s);

    std::string mass_s;
    pp.get("mass", mass_s);
    std::transform(mass_s.begin(),
                   mass_s.end(),
                   mass_s.begin(),
                   ::tolower);
    mass = parseMassString(pp, mass_s);

    // parse injection style
    std::string part_pos_s;
    pp.get("injection_style", part_pos_s);
    std::transform(part_pos_s.begin(),
                   part_pos_s.end(),
                   part_pos_s.begin(),
                   ::tolower);
    if (part_pos_s == "python") {
        return;
    } else if (part_pos_s == "singleparticle") {
        pp.getarr("single_particle_pos", single_particle_pos, 0, 3);
        pp.getarr("single_particle_vel", single_particle_vel, 0, 3);
        for (auto& x : single_particle_vel) {
            x *= PhysConst::c;
        }
        pp.get("single_particle_weight", single_particle_weight);
        add_single_particle = true;
        return;
    } else if (part_pos_s == "gaussian_beam") {
        pp.get("x_m", x_m);
        pp.get("y_m", y_m);
        pp.get("z_m", z_m);
        pp.get("x_rms", x_rms);
        pp.get("y_rms", y_rms);
        pp.get("z_rms", z_rms);
        pp.get("q_tot", q_tot);
        pp.get("npart", npart);
        gaussian_beam = true;
    }
    else if (part_pos_s == "nrandompercell") {
        pp.query("num_particles_per_cell", num_particles_per_cell);
        part_pos.reset(new RandomPosition(num_particles_per_cell));
    } else if (part_pos_s == "nuniformpercell") {
        num_particles_per_cell_each_dim.resize(3);
        pp.getarr("num_particles_per_cell_each_dim", num_particles_per_cell_each_dim);
#if ( BL_SPACEDIM == 2 )
        num_particles_per_cell_each_dim[2] = 1;
#endif
        part_pos.reset(new RegularPosition(num_particles_per_cell_each_dim));
        num_particles_per_cell = num_particles_per_cell_each_dim[0] *
                                 num_particles_per_cell_each_dim[1] *
                                 num_particles_per_cell_each_dim[2];
    } else {
        StringParseAbortMessage("Injection style", part_pos_s);
    }

    // parse plasma boundaries
    xmin = std::numeric_limits<amrex::Real>::lowest();
    ymin = std::numeric_limits<amrex::Real>::lowest();
    zmin = std::numeric_limits<amrex::Real>::lowest();

    xmax = std::numeric_limits<amrex::Real>::max();
    ymax = std::numeric_limits<amrex::Real>::max();
    zmax = std::numeric_limits<amrex::Real>::max();

    pp.query("xmin", xmin);
    pp.query("ymin", ymin);
    pp.query("zmin", zmin);
    pp.query("xmax", xmax);
    pp.query("ymax", ymax);
    pp.query("zmax", zmax);

    // parse density information
    std::string rho_prof_s;
    pp.get("profile", rho_prof_s);
    std::transform(rho_prof_s.begin(),
                   rho_prof_s.end(),
                   rho_prof_s.begin(),
                   ::tolower);
    if (rho_prof_s == "constant") {
        pp.get("density", density);
        rho_prof.reset(new ConstantDensityProfile(density));
    } else if (rho_prof_s == "custom") {
        rho_prof.reset(new CustomDensityProfile(species_name));
    } else {
        StringParseAbortMessage("Density profile type", rho_prof_s);
    }

    // parse momentum information
    std::string mom_dist_s;
    pp.get("momentum_distribution_type", mom_dist_s);
    std::transform(mom_dist_s.begin(),
                   mom_dist_s.end(),
                   mom_dist_s.begin(),
                   ::tolower);
    if (mom_dist_s == "constant") {
        Real ux = 0.;
        Real uy = 0.;
        Real uz = 0.;
        pp.query("ux", ux);
        pp.query("uy", uy);
        pp.query("uz", uz);
        mom_dist.reset(new ConstantMomentumDistribution(ux, uy, uz));
    } else if (mom_dist_s == "gaussian") {
        Real ux_m = 0.;
        Real uy_m = 0.;
        Real uz_m = 0.;
        Real u_th = 0.;
        pp.query("ux_m", ux_m);
        pp.query("uy_m", uy_m);
        pp.query("uz_m", uz_m);
        pp.query("u_th", u_th);
        mom_dist.reset(new GaussianRandomMomentumDistribution(ux_m, uy_m, uz_m, u_th));
    } else if (mom_dist_s == "radial_expansion") {
        Real u_over_r = 0.;
	pp.query("u_over_r", u_over_r);
        mom_dist.reset(new RadialExpansionMomentumDistribution(u_over_r));
    } else {
        StringParseAbortMessage("Momentum distribution type", mom_dist_s);
    }
}

void PlasmaInjector::getPositionUnitBox(vec3& r, int i_part) {
  return part_pos->getPositionUnitBox(r, i_part);
}

void PlasmaInjector::getMomentum(vec3& u, Real x, Real y, Real z) {
    mom_dist->getMomentum(u, x, y, z);
    u[0] *= PhysConst::c;
    u[1] *= PhysConst::c;
    u[2] *= PhysConst::c;
}

bool PlasmaInjector::insideBounds(Real x, Real y, Real z) {
  if (x >= xmax || x < xmin ||
      y >= ymax || y < ymin ||
      z >= zmax || z < zmin ) return false;
  return true;
}

Real PlasmaInjector::getDensity(Real x, Real y, Real z) {
    return rho_prof->getDensity(x, y, z);
}
