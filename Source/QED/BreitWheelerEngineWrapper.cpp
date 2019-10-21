#include "BreitWheelerEngineWrapper.H"

#include "QedTableParserHelperFunctions.H"

#include <AMReX_Print.H>

#include <utility>
#include <vector>

using namespace std;
using namespace QedUtils;
using namespace amrex;

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the PICSAR library

// Factory class =============================

BreitWheelerEngine::BreitWheelerEngine (){}

//Builds the functor to initialize the optical depth
BreitWheelerGetOpticalDepth
BreitWheelerEngine::build_optical_depth_functor ()
{
    return BreitWheelerGetOpticalDepth();
}

//Builds the functor to evolve the optical depth
BreitWheelerEvolveOpticalDepth
BreitWheelerEngine::build_evolve_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return BreitWheelerEvolveOpticalDepth(&innards);
}

//Builds the functor to generate the pairs
BreitWheelerGeneratePairs
BreitWheelerEngine::build_pair_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return BreitWheelerGeneratePairs(&innards);
}

bool BreitWheelerEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}

/* \brief Reads lookup tables from 'file' on disk */
bool
BreitWheelerEngine::init_lookup_tables_from_raw_data (
    const vector<char>& raw_data)
{
    const char* p_data = raw_data.data();
    const char* const p_last = &raw_data.back();
    bool is_ok;

    //Header (control parameters)
    tie(is_ok, innards.ctrl.chi_phot_min, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_phot_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_phot_tdndt_min, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_phot_tdndt_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_phot_tdndt_max, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_phot_tdndt_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_phot_tdndt_how_many, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_phot_tdndt_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_phot_tpair_min, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_phot_tpair_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_phot_tpair_max, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_phot_tpair_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_phot_tpair_how_many, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_phot_tpair_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_frac_tpair_how_many, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_frac_tpair_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    //___________________________

    //Data
    vector<Real> tndt_coords(innards.ctrl.chi_phot_tdndt_how_many);
    vector<Real> tndt_data(innards.ctrl.chi_phot_tdndt_how_many);
    vector<Real> cum_tab_coords1(innards.ctrl.chi_phot_tpair_how_many);
    vector<Real> cum_tab_coords2(innards.ctrl.chi_frac_tpair_how_many);
    vector<Real> cum_tab_data(innards.ctrl.chi_phot_tpair_how_many*
        innards.ctrl.chi_frac_tpair_how_many);

    tie(is_ok, tndt_coords, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_coords.size(), p_last);
    if(!is_ok) return false;
    innards.TTfunc_coords.assign(tndt_coords.begin(), tndt_coords.end());

    tie(is_ok, tndt_data, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_data.size(), p_last);
    if(!is_ok) return false;
    innards.TTfunc_data.assign(tndt_data.begin(), tndt_data.end());

    tie(is_ok, cum_tab_coords1, p_data) =
        parse_raw_data_vec<Real>(
            p_data, cum_tab_coords1.size(), p_last);
    if(!is_ok) return false;
    innards.cum_distrib_coords_1.assign(
        cum_tab_coords1.begin(), cum_tab_coords1.end());

    tie(is_ok, cum_tab_coords2, p_data) =
        parse_raw_data_vec<Real>(
            p_data, cum_tab_coords2.size(), p_last);
    if(!is_ok) return false;
    innards.cum_distrib_coords_2.assign(
        cum_tab_coords2.begin(), cum_tab_coords2.end());

    tie(is_ok, cum_tab_data, p_data) =
        parse_raw_data_vec<Real>(
            p_data, cum_tab_data.size(), p_last);
    if(!is_ok) return false;
    innards.cum_distrib_data.assign(
        cum_tab_data.begin(), cum_tab_data.end());

    //___________________________
    lookup_tables_initialized = true;

    return true;
}

//Initializes the Lookup tables using the default settings
//provided by the library
void BreitWheelerEngine::compute_lookup_tables_default ()
{
    //A control parameters structure
    //with the default values provided by the library
    PicsarBreitWheelerCtrl ctrl_default;

    compute_lookup_tables(ctrl_default);

    lookup_tables_initialized = true;
}

// Computes the Lookup tables using user-defined settings
void BreitWheelerEngine::compute_custom_lookup_tables (
    PicsarBreitWheelerCtrl ctrl)
{
    compute_lookup_tables(ctrl);

    lookup_tables_initialized = true;
}


/* \brief Writes lookup tables on disk in 'file'
 *  return false if it fails. */
std::vector<char> BreitWheelerEngine::export_lookup_tables_data () const
{
    if(!lookup_tables_initialized)
        return std::vector<char>{};

    //TODO
    return std::vector<char>{};
}



//Private function which actually computes the lookup tables
void BreitWheelerEngine::compute_lookup_tables (
    PicsarBreitWheelerCtrl ctrl)
{
#ifdef WARPX_QED_TABLE_GEN
    table_builder.compute_table(ctrl, innards);
#else
    amrex::Print() <<
        "Error: use QED_TABLE_GEN=TRUE to enable table generation!\n";
#endif
}

//============================================
