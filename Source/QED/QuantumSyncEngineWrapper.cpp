#include "QuantumSyncEngineWrapper.H"

#include "QedTableParserHelperFunctions.H"

#include <AMReX_Print.H>

#include <utility>
#include <vector>

using namespace std;
using namespace QedUtils;
using namespace amrex;

//This file provides a wrapper aroud the quantum_sync engine
//provided by the PICSAR library

// Factory class =============================

QuantumSynchrotronEngine::QuantumSynchrotronEngine (){}

//Builds the functor to evolve the optical depth
QuantumSynchrotronGetOpticalDepth
QuantumSynchrotronEngine::build_optical_depth_functor ()
{
    return QuantumSynchrotronGetOpticalDepth();
}

//Builds the functor to evolve the optical depth
QuantumSynchrotronEvolveOpticalDepth QuantumSynchrotronEngine::build_evolve_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return QuantumSynchrotronEvolveOpticalDepth(&innards);
}

bool QuantumSynchrotronEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}

/* \brief Reads lookup tables from 'file' on disk */
bool
QuantumSynchrotronEngine::init_lookup_tables_from_raw_data (
    const vector<char>& raw_data)
{
    const char* p_data = raw_data.data();
    const char* const p_last = &raw_data.back();
    bool is_ok;

    //Header (control parameters)
    tie(is_ok, innards.ctrl.chi_part_min, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_part_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_part_tdndt_min, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_part_tdndt_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_part_tdndt_max, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_part_tdndt_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_part_tdndt_how_many, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_part_tdndt_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_part_tem_min, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_part_tem_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_part_tem_max, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_part_tem_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.chi_part_tem_how_many, p_data) =
        parse_raw_data<decltype(innards.ctrl.chi_part_tem_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, innards.ctrl.prob_tem_how_many, p_data) =
        parse_raw_data<decltype(innards.ctrl.prob_tem_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    //___________________________

    //Data
    vector<Real> tndt_coords(innards.ctrl.chi_part_tdndt_min);
    vector<Real> tndt_data(innards.ctrl.chi_part_tdndt_min);
    vector<Real> cum_tab_coords1(innards.ctrl.chi_part_tem_how_many);
    vector<Real> cum_tab_coords2(innards.ctrl.prob_tem_how_many);
    vector<Real> cum_tab_data(innards.ctrl.chi_part_tem_how_many*
        innards.ctrl.prob_tem_how_many);

    tie(is_ok, tndt_coords, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_coords.size(), p_last);
    if(!is_ok) return false;
    innards.KKfunc_coords.assign(tndt_coords.begin(), tndt_coords.end());

    tie(is_ok, tndt_data, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_data.size(), p_last);
    if(!is_ok) return false;
    innards.KKfunc_data.assign(tndt_data.begin(), tndt_data.end());

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
void QuantumSynchrotronEngine::compute_lookup_tables_default ()
{
    //A control parameters structure
    //with the default values provided by the library
    PicsarQuantumSynchrotronCtrl ctrl_default;

    compute_lookup_tables(ctrl_default);

    lookup_tables_initialized = true;
}

// Computes the Lookup tables using user-defined settings
void QuantumSynchrotronEngine::compute_custom_lookup_tables (
    PicsarQuantumSynchrotronCtrl ctrl)
{
    compute_lookup_tables(ctrl);

    lookup_tables_initialized = true;
}

/* \brief Writes lookup tables on disk in 'file'
 *  return false if it fails. */
std::vector<char> QuantumSynchrotronEngine::export_lookup_tables_data () const
{
    if(!lookup_tables_initialized)
        return std::vector<char>{};

    //TODO
    return std::vector<char>{};
}

//Private function which actually computes the lookup tables
void QuantumSynchrotronEngine::compute_lookup_tables (
    PicsarQuantumSynchrotronCtrl ctrl)
{
#ifdef WARPX_QED_TABLE_GEN
    table_builder.compute_table(ctrl, innards);
#else
    amrex::Print() <<
        "Error: use QED_TABLE_GEN=TRUE to enable table generation!\n";
#endif
}


//============================================
