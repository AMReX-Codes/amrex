#include <LaserProfiles.H>

#include <WarpX_Complex.H>

using namespace amrex;
using namespace WarpXLaserProfiles;

void
FieldFunctionLaserProfile::init (
    const amrex::ParmParse& ppl,
    const amrex::ParmParse& ppc,
    CommonLaserParameters params)
{
    // Parse the properties of the parse_field_function profile
    ppl.get("field_function(X,Y,t)", m_params.field_function);
    m_parser.define(m_params.field_function);
    m_parser.registerVariables({"X","Y","t"});

    std::set<std::string> symbols = m_parser.symbols();
    symbols.erase("X");
    symbols.erase("Y");
    symbols.erase("t"); // after removing variables, we are left with constants
    for (auto it = symbols.begin(); it != symbols.end(); ) {
        Real v;
        if (ppc.query(it->c_str(), v)) {
            m_parser.setConstant(*it, v);
            it = symbols.erase(it);
        } else {
            ++it;
        }
    }
    for (auto const& s : symbols) { // make sure there no unknown symbols
        amrex::Abort("Laser Profile: Unknown symbol "+s);
    }
}

void
FieldFunctionLaserProfile::fill_amplitude (
    const int np, Real const * AMREX_RESTRICT const Xp, Real const * AMREX_RESTRICT const Yp,
    Real t, Real * AMREX_RESTRICT const amplitude)
{
    for (int i = 0; i < np; ++i) {
        amplitude[i] = m_parser.eval(Xp[i], Yp[i], t);
    }
}
