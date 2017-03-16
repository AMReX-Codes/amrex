#include "PlasmaInjector.H"

#include <sstream>

#include <WarpXConst.H>
#include <AMReX.H>

using namespace amrex;

Real parseChargeName(const std::string& name) {
    if (name == "q_e") {
        return PhysConst::q_e;
    } else {
        std::stringstream stringstream;
        std::string string;
        stringstream << "Charge string '" << name << "' not recognized."; 
        string = stringstream.str();
        amrex::Abort(string.c_str());
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
        std::stringstream stringstream;
        std::string string;
        stringstream << "Mass string '" << name << "' not recognized."; 
        string = stringstream.str();
        amrex::Abort(string.c_str());
        return 0.0;
    }
}

ConstantMomentumDistribution::ConstantMomentumDistribution(Real ux,
                                                           Real uy,
                                                           Real uz) 
    : _ux(ux), _uy(uy), _uz(uz)
{
}

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
    
    pp.get("num_particles_per_cell", num_particles_per_cell);
    pp.get("density", density);
    
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
        std::stringstream stringstream;
        std::string string;
        stringstream << "Momentum distribution type '" << mom_dist_s << "' not recognized."; 
        string = stringstream.str();
        amrex::Abort(string.c_str());        
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
