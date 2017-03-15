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
        
    pp.get("xmin", xmin);
    pp.get("ymin", ymin);
    pp.get("zmin", zmin);
    pp.get("xmax", xmax);
    pp.get("ymax", ymax);
    pp.get("zmax", zmax);
    
    pp.get("num_particles_per_cell", num_particles_per_cell);
    pp.get("density", density);
    
    ux = 0.;
    uy = 0.;
    uz = 0.;
    pp.query("ux", ux);
    pp.query("uy", uy);
    pp.query("uz", uz);
}

void PlasmaInjector::getMomentum(Real* u) {
    u[0] = ux * PhysConst::c;
    u[1] = uy * PhysConst::c;
    u[2] = uz * PhysConst::c;
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


