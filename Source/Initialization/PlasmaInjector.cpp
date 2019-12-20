#include "PlasmaInjector.H"

#include <WarpXConst.H>
#include <WarpX_f.H>
#include <WarpX.H>

#include <AMReX.H>

#include <sstream>
#include <functional>

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

    Real parseChargeName(const ParmParse& pp, const std::string& name) {
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

    Real parseChargeString(const ParmParse& pp, const std::string& name) {
        if(name.substr(0, 1) == "-")
            return -1.0 * parseChargeName(pp, name.substr(1, name.size() - 1));
        return parseChargeName(pp, name);
    }

    Real parseMassString(const ParmParse& pp, const std::string& name) {
        Real result;
        if (name == "m_e") {
            return PhysConst::m_e;
        } else if (name == "m_p"){
            return PhysConst::m_p;
        } else if (name == "inf"){
            return std::numeric_limits<double>::infinity();
        } else if (pp.query("mass", result)) {
            return result;
        } else {
            StringParseAbortMessage("Mass", name);
            return 0.0;
        }
    }
}

PlasmaInjector::PlasmaInjector () {}

PlasmaInjector::PlasmaInjector (int ispecies, const std::string& name)
    : species_id(ispecies), species_name(name)
{
    ParmParse pp(species_name);

    pp.query("radially_weighted", radially_weighted);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(radially_weighted, "ERROR: Only radially_weighted=true is supported");

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

    pp.query("density_min", density_min);
    pp.query("density_max", density_max);

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
    num_particles_per_cell_each_dim.assign(3, 0);
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
        pp.query("do_symmetrize", do_symmetrize);
        gaussian_beam = true;
        parseMomentum(pp);
    }
    // Depending on injection type at runtime, initialize inj_pos
    // so that inj_pos->getPositionUnitBox calls
    // InjectorPosition[Random or Regular].getPositionUnitBox.
    else if (part_pos_s == "nrandompercell") {
        pp.query("num_particles_per_cell", num_particles_per_cell);
        // Construct InjectorPosition with InjectorPositionRandom.
        inj_pos.reset(new InjectorPosition((InjectorPositionRandom*)nullptr,
                                           xmin, xmax, ymin, ymax, zmin, zmax));
        parseDensity(pp);
        parseMomentum(pp);
    } else if (part_pos_s == "nuniformpercell") {
        // Note that for RZ, three numbers are expected, r, theta, and z.
        // For 2D, only two are expected. The third is overwritten with 1.
        num_particles_per_cell_each_dim.assign(3, 1);
        pp.getarr("num_particles_per_cell_each_dim", num_particles_per_cell_each_dim);
#if WARPX_DIM_XZ
        num_particles_per_cell_each_dim[2] = 1;
#endif
        // Construct InjectorPosition from InjectorPositionRegular.
        inj_pos.reset(new InjectorPosition((InjectorPositionRegular*)nullptr,
                                           xmin, xmax, ymin, ymax, zmin, zmax,
                                           Dim3{num_particles_per_cell_each_dim[0],
                                                num_particles_per_cell_each_dim[1],
                                                num_particles_per_cell_each_dim[2]}));
        num_particles_per_cell = num_particles_per_cell_each_dim[0] *
                                 num_particles_per_cell_each_dim[1] *
                                 num_particles_per_cell_each_dim[2];
        parseDensity(pp);
        parseMomentum(pp);
    } else {
        StringParseAbortMessage("Injection style", part_pos_s);
    }
}

namespace {
WarpXParser makeParser (std::string const& parse_function)
{
    WarpXParser parser(parse_function);
    parser.registerVariables({"x","y","z"});

    ParmParse pp("my_constants");
    std::set<std::string> symbols = parser.symbols();
    symbols.erase("x");
    symbols.erase("y");
    symbols.erase("z"); // after removing variables, we are left with constants
    for (auto it = symbols.begin(); it != symbols.end(); ) {
        Real v;
        if (pp.query(it->c_str(), v)) {
            parser.setConstant(*it, v);
            it = symbols.erase(it);
        } else {
            ++it;
        }
    }
    for (auto const& s : symbols) { // make sure there no unknown symbols
        amrex::Abort("PlasmaInjector::makeParser: Unknown symbol "+s);
    }

    return parser;
}
}

// Depending on injection type at runtime, initialize inj_rho
// so that inj_rho->getDensity calls
// InjectorPosition[Constant or Custom or etc.].getDensity.
void PlasmaInjector::parseDensity (ParmParse& pp)
{
    // parse density information
    std::string rho_prof_s;
    pp.get("profile", rho_prof_s);
    std::transform(rho_prof_s.begin(), rho_prof_s.end(),
                   rho_prof_s.begin(), ::tolower);
    if (rho_prof_s == "constant") {
        pp.get("density", density);
        // Construct InjectorDensity with InjectorDensityConstant.
        inj_rho.reset(new InjectorDensity((InjectorDensityConstant*)nullptr, density));
    } else if (rho_prof_s == "custom") {
        // Construct InjectorDensity with InjectorDensityCustom.
        inj_rho.reset(new InjectorDensity((InjectorDensityCustom*)nullptr, species_name));
    } else if (rho_prof_s == "predefined") {
        // Construct InjectorDensity with InjectorDensityPredefined.
        inj_rho.reset(new InjectorDensity((InjectorDensityPredefined*)nullptr,species_name));
    } else if (rho_prof_s == "parse_density_function") {
        std::vector<std::string> f;
        pp.getarr("density_function(x,y,z)", f);
        for (auto const& s : f) {
            str_density_function += s;
        }
        // Construct InjectorDensity with InjectorDensityParser.
        inj_rho.reset(new InjectorDensity((InjectorDensityParser*)nullptr,
                                          makeParser(str_density_function)));
    } else {
        StringParseAbortMessage("Density profile type", rho_prof_s);
    }
}

// Depending on injection type at runtime, initialize inj_mom
// so that inj_mom->getMomentum calls
// InjectorMomentum[Constant or Custom or etc.].getMomentum.
void PlasmaInjector::parseMomentum (ParmParse& pp)
{
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
        // Construct InjectorMomentum with InjectorMomentumConstant.
        inj_mom.reset(new InjectorMomentum((InjectorMomentumConstant*)nullptr, ux,uy, uz));
    } else if (mom_dist_s == "custom") {
        // Construct InjectorMomentum with InjectorMomentumCustom.
        inj_mom.reset(new InjectorMomentum((InjectorMomentumCustom*)nullptr, species_name));
    } else if (mom_dist_s == "gaussian") {
        Real ux_m = 0.;
        Real uy_m = 0.;
        Real uz_m = 0.;
        Real ux_th = 0.;
        Real uy_th = 0.;
        Real uz_th = 0.;
        pp.query("ux_m", ux_m);
        pp.query("uy_m", uy_m);
        pp.query("uz_m", uz_m);
        pp.query("ux_th", ux_th);
        pp.query("uy_th", uy_th);
        pp.query("uz_th", uz_th);
        // Construct InjectorMomentum with InjectorMomentumGaussian.
        inj_mom.reset(new InjectorMomentum((InjectorMomentumGaussian*)nullptr,
                                           ux_m, uy_m, uz_m, ux_th, uy_th, uz_th));
    } else if (mom_dist_s == "maxwell_boltzmann"){
        Real beta = 0.;
        Real theta = 10.;
        int dir = 0;
        std::string direction = "x";
        pp.query("beta", beta);
        if(beta < 0){
            amrex::Abort("Please enter a positive beta value. Drift direction is set with <s_name>.bulk_vel_dir = 'x' or '+x', '-x', 'y' or '+y', etc.");
        }
        pp.query("theta", theta);
        pp.query("bulk_vel_dir", direction);
        if(direction[0] == '-'){
            beta = -beta;
        }
        if((direction == "x" || direction[1] == 'x') ||
           (direction == "X" || direction[1] == 'X')){
            dir = 0;
        } else if ((direction == "y" || direction[1] == 'y') ||
                   (direction == "Y" || direction[1] == 'Y')){
            dir = 1;
        } else if ((direction == "z" || direction[1] == 'z') ||
                   (direction == "Z" || direction[1] == 'Z')){
            dir = 2;
        } else{
            std::stringstream stringstream;
            stringstream << "Cannot interpret <s_name>.bulk_vel_dir input '" << direction << "'. Please enter +/- x, y, or z with no whitespace between the sign and other character.";
            direction = stringstream.str();
            amrex::Abort(direction.c_str());
        }
        // Construct InjectorMomentum with InjectorMomentumBoltzmann.
        inj_mom.reset(new InjectorMomentum((InjectorMomentumBoltzmann*)nullptr, theta, beta, dir));
    } else if (mom_dist_s == "maxwell_juttner"){
        Real beta = 0.;
        Real theta = 10.;
        int dir = 0;
        std::string direction = "x";
        pp.query("beta", beta);
        if(beta < 0){
            amrex::Abort("Please enter a positive beta value. Drift direction is set with <s_name>.bulk_vel_dir = 'x' or '+x', '-x', 'y' or '+y', etc.");
        }
        pp.query("theta", theta);
        pp.query("bulk_vel_dir", direction);
        if(direction[0] == '-'){
            beta = -beta;
        }
        if((direction == "x" || direction[1] == 'x') ||
           (direction == "X" || direction[1] == 'X')){
            dir = 0;
        } else if ((direction == "y" || direction[1] == 'y') ||
                   (direction == "Y" || direction[1] == 'Y')){
            dir = 1;
        } else if ((direction == "z" || direction[1] == 'z') ||
                   (direction == "Z" || direction[1] == 'Z')){
            dir = 2;
        } else{
            std::stringstream stringstream;
            stringstream << "Cannot interpret <s_name>.bulk_vel_dir input '" << direction << "'. Please enter +/- x, y, or z with no whitespace between the sign and other character.";
            direction = stringstream.str();
            amrex::Abort(direction.c_str());
        }
        // Construct InjectorMomentum with InjectorMomentumJuttner.
        inj_mom.reset(new InjectorMomentum((InjectorMomentumJuttner*)nullptr, theta, beta, dir));
    } else if (mom_dist_s == "radial_expansion") {
        Real u_over_r = 0.;
        pp.query("u_over_r", u_over_r);
        // Construct InjectorMomentum with InjectorMomentumRadialExpansion.
        inj_mom.reset(new InjectorMomentum
                      ((InjectorMomentumRadialExpansion*)nullptr, u_over_r));
    } else if (mom_dist_s == "parse_momentum_function") {
        std::vector<std::string> f;
        pp.getarr("momentum_function_ux(x,y,z)", f);
        for (auto const& s : f) {
            str_momentum_function_ux += s;
        }
        f.clear();
        pp.getarr("momentum_function_uy(x,y,z)", f);
        for (auto const& s : f) {
            str_momentum_function_uy += s;
        }
        f.clear();
        pp.getarr("momentum_function_uz(x,y,z)", f);
        for (auto const& s : f) {
            str_momentum_function_uz += s;
        }
        // Construct InjectorMomentum with InjectorMomentumParser.
        inj_mom.reset(new InjectorMomentum((InjectorMomentumParser*)nullptr,
                                           makeParser(str_momentum_function_ux),
                                           makeParser(str_momentum_function_uy),
                                           makeParser(str_momentum_function_uz)));
    } else {
        StringParseAbortMessage("Momentum distribution type", mom_dist_s);
    }
}

XDim3 PlasmaInjector::getMomentum (Real x, Real y, Real z) const noexcept
{
    return inj_mom->getMomentum(x, y, z); // gamma*beta
}

bool PlasmaInjector::insideBounds (Real x, Real y, Real z) const noexcept
{
    return (x < xmax and x >= xmin and
            y < ymax and y >= ymin and
            z < zmax and z >= zmin);
}

InjectorPosition*
PlasmaInjector::getInjectorPosition ()
{
    return inj_pos.get();
}

InjectorDensity*
PlasmaInjector::getInjectorDensity ()
{
    return inj_rho.get();
}

InjectorMomentum*
PlasmaInjector::getInjectorMomentum ()
{
    return inj_mom.get();
}

