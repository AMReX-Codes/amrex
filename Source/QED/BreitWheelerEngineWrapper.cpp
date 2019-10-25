#include "BreitWheelerEngineWrapper.H"

#include "QedTableParserHelperFunctions.H"

#include <utility>

using namespace std;
using namespace QedUtils;
using namespace amrex;

//This file provides a wrapper aroud the breit_wheeler engine
//provided by the PICSAR library

// Factory class =============================

BreitWheelerEngine::BreitWheelerEngine (){}

BreitWheelerGetOpticalDepth
BreitWheelerEngine::build_optical_depth_functor ()
{
    return BreitWheelerGetOpticalDepth();
}

BreitWheelerEvolveOpticalDepth
BreitWheelerEngine::build_evolve_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return BreitWheelerEvolveOpticalDepth(innards);
}

BreitWheelerGeneratePairs
BreitWheelerEngine::build_pair_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return BreitWheelerGeneratePairs(innards);
}

bool BreitWheelerEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}

bool
BreitWheelerEngine::init_lookup_tables_from_raw_data (
    const Vector<char>& raw_data)
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
    Vector<Real> tndt_coords(innards.ctrl.chi_phot_tdndt_how_many);
    Vector<Real> tndt_data(innards.ctrl.chi_phot_tdndt_how_many);
    Vector<Real> cum_tab_coords1(innards.ctrl.chi_phot_tpair_how_many);
    Vector<Real> cum_tab_coords2(innards.ctrl.chi_frac_tpair_how_many);
    Vector<Real> cum_tab_data(innards.ctrl.chi_phot_tpair_how_many*
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

Vector<char> BreitWheelerEngine::export_lookup_tables_data () const
{
   Vector<char> res{};

    if(!lookup_tables_initialized)
        return res;

    auto add_data_to_buf = [&res](auto data){
        res.insert(res.end(),
        reinterpret_cast<const char*>(&data),
        reinterpret_cast<const char*>(&data) +
        sizeof(data));
    };

    add_data_to_buf(innards.ctrl.chi_phot_min);
    add_data_to_buf(innards.ctrl.chi_phot_tdndt_min);
    add_data_to_buf(innards.ctrl.chi_phot_tdndt_max);
    add_data_to_buf(innards.ctrl.chi_phot_tdndt_how_many);
    add_data_to_buf(innards.ctrl.chi_phot_tpair_min);
    add_data_to_buf(innards.ctrl.chi_phot_tpair_max);
    add_data_to_buf(innards.ctrl.chi_phot_tpair_how_many);
    add_data_to_buf(innards.ctrl.chi_frac_tpair_how_many);

    auto add_data_to_buf_vv = [&res](const auto* data, size_t how_many){
        res.insert(res.end(),
        reinterpret_cast<const char*>(data),
        reinterpret_cast<const char*>(data) +
        sizeof(*data)*how_many);
    };

    add_data_to_buf_vv(innards.TTfunc_coords.data(),
        innards.TTfunc_coords.size());
    add_data_to_buf_vv(innards.TTfunc_data.data(),
        innards.TTfunc_data.size());
    add_data_to_buf_vv(innards.cum_distrib_coords_1.data(),
        innards.cum_distrib_coords_1.size());
    add_data_to_buf_vv(innards.cum_distrib_coords_2.data(),
        innards.cum_distrib_coords_2.size());
    add_data_to_buf_vv(innards.cum_distrib_data.data(),
        innards.cum_distrib_data.size());

    return res;
}

PicsarBreitWheelerCtrl
BreitWheelerEngine::get_default_ctrl() const
{
    return PicsarBreitWheelerCtrl();
}

void BreitWheelerEngine::compute_lookup_tables (
    PicsarBreitWheelerCtrl ctrl)
{
#ifdef WARPX_QED_TABLE_GEN
    table_builder.compute_table(ctrl, innards);
    lookup_tables_initialized = true;
#endif
}

//============================================
