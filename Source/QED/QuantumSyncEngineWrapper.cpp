#include "QuantumSyncEngineWrapper.H"

#include "QedTableParserHelperFunctions.H"

#include <utility>

using namespace std;
using namespace QedUtils;
using namespace amrex;

//This file provides a wrapper aroud the quantum_sync engine
//provided by the PICSAR library

// Factory class =============================

QuantumSynchrotronEngine::QuantumSynchrotronEngine (){}

QuantumSynchrotronGetOpticalDepth
QuantumSynchrotronEngine::build_optical_depth_functor ()
{
    return QuantumSynchrotronGetOpticalDepth();
}

QuantumSynchrotronEvolveOpticalDepth QuantumSynchrotronEngine::build_evolve_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return QuantumSynchrotronEvolveOpticalDepth(innards);
}

QuantumSynchrotronGeneratePhotonAndUpdateMomentum QuantumSynchrotronEngine::build_phot_em_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return QuantumSynchrotronGeneratePhotonAndUpdateMomentum(innards);

}

bool QuantumSynchrotronEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}

bool
QuantumSynchrotronEngine::init_lookup_tables_from_raw_data (
    const Vector<char>& raw_data)
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
    Vector<Real> tndt_coords(innards.ctrl.chi_part_tdndt_min);
    Vector<Real> tndt_data(innards.ctrl.chi_part_tdndt_min);
    Vector<Real> cum_tab_coords1(innards.ctrl.chi_part_tem_how_many);
    Vector<Real> cum_tab_coords2(innards.ctrl.prob_tem_how_many);
    Vector<Real> cum_tab_data(innards.ctrl.chi_part_tem_how_many*
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

Vector<char> QuantumSynchrotronEngine::export_lookup_tables_data () const
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

    add_data_to_buf(innards.ctrl.chi_part_min);
    add_data_to_buf(innards.ctrl.chi_part_tdndt_min);
    add_data_to_buf(innards.ctrl.chi_part_tdndt_max);
    add_data_to_buf(innards.ctrl.chi_part_tdndt_how_many);
    add_data_to_buf(innards.ctrl.chi_part_tem_min);
    add_data_to_buf(innards.ctrl.chi_part_tem_max);
    add_data_to_buf(innards.ctrl.chi_part_tem_how_many);
    add_data_to_buf(innards.ctrl.prob_tem_how_many);

    auto add_data_to_buf_vv = [&res](const auto* data, size_t how_many){
        res.insert(res.end(),
        reinterpret_cast<const char*>(data),
        reinterpret_cast<const char*>(data) +
        sizeof(*data)*how_many);
    };

    add_data_to_buf_vv(innards.KKfunc_coords.data(),
        innards.KKfunc_coords.size());
    add_data_to_buf_vv(innards.KKfunc_data.data(),
        innards.KKfunc_data.size());
    add_data_to_buf_vv(innards.cum_distrib_coords_1.data(),
        innards.cum_distrib_coords_1.size());
    add_data_to_buf_vv(innards.cum_distrib_coords_2.data(),
        innards.cum_distrib_coords_2.size());
    add_data_to_buf_vv(innards.cum_distrib_data.data(),
        innards.cum_distrib_data.size());

    return res;
}

PicsarQuantumSynchrotronCtrl
QuantumSynchrotronEngine::get_default_ctrl() const
{
    return PicsarQuantumSynchrotronCtrl();
}

void QuantumSynchrotronEngine::compute_lookup_tables (
    PicsarQuantumSynchrotronCtrl ctrl)
{
#ifdef WARPX_QED_TABLE_GEN
    table_builder.compute_table(ctrl, innards);
    lookup_tables_initialized = true;
#endif
}

//============================================
