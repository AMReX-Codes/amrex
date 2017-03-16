#include "PlasmaInjector.H"

#include <sstream>

#include <WarpXConst.H>
#include <AMReX.H>

using namespace amrex;

void StringParseAbortMessage(const std::string& var,
                             const std::string& name) {
    std::stringstream stringstream;
    std::string string;
    stringstream << var << " string '" << name << "' not recognized."; 
    string = stringstream.str();
    amrex::Abort(string.c_str());    
}

Real parseChargeName(const std::string& name) {
    if (name == "q_e") {
        return PhysConst::q_e;
    } else {
        StringParseAbortMessage("Charge", name);
        return 0.0;
    }
}

Real parseChargeString(const std::string& name) {
    if(name.substr(0, 1) == "-")
        return -1.0 * parseChargeName(name.substr(1, name.size() - 1)); 
    return parseChargeName(name);
}

Real parseMassString(const std::string& name) {
    if (name == "m_e") {
        return PhysConst::m_e;
    } else if (name == "m_p"){
        return PhysConst::m_p;
    } else {
        StringParseAbortMessage("Mass", name);
        return 0.0;
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

void ConstantMomentumDistribution::getMomentum(Real* u) {
    u[0] = _ux;
    u[1] = _uy;
    u[2] = _uz;
}

GaussianRandomMomentumDistribution::GaussianRandomMomentumDistribution(Real ux_m,
                                                                       Real uy_m,
                                                                       Real uz_m,
                                                                       Real u_th) 
    : _ux_m(ux_m), _uy_m(uy_m), _uz_m(uz_m), _u_th(u_th), 
      momentum_distribution(0.0, u_th)
{
}

void GaussianRandomMomentumDistribution::getMomentum(Real* u) {
    Real ux_th = momentum_distribution(generator);
    Real uy_th = momentum_distribution(generator);
    Real uz_th = momentum_distribution(generator);
    u[0] = _ux_m + ux_th;
    u[1] = _uy_m + uy_th;
    u[2] = _uz_m + uz_th;
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
    charge = parseChargeString(charge_s);

    std::string mass_s;
    pp.get("mass", mass_s);
    std::transform(mass_s.begin(), 
                   mass_s.end(), 
                   mass_s.begin(), 
                   ::tolower);
    mass = parseMassString(mass_s);

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
        Real density;
        pp.get("density", density);
        rho_prof.reset(new ConstantDensityProfile(density));
    } else if (rho_prof_s == "custom") {
        rho_prof.reset(new CustomDensityProfile(species_name));
    } else {
        StringParseAbortMessage("Density profile type", rho_prof_s);
    }
    pp.get("num_particles_per_cell", num_particles_per_cell);
    
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
    } else {
        StringParseAbortMessage("Momentum distribution type", mom_dist_s);
    }

    // get injection style
    pp.get("injection_style", injection_style);
    std::transform(injection_style.begin(), 
                   injection_style.end(), 
                   injection_style.begin(), 
                   ::tolower);
    if (injection_style == "nrandomnormal" or
        injection_style == "nrandomuniformpercell" or
        injection_style == "ndiagpercell") {
        return;
    } else {
        StringParseAbortMessage("Injection style", injection_style);
    }
}

void PlasmaInjector::getMomentum(Real* u) {
    mom_dist->getMomentum(u);
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
