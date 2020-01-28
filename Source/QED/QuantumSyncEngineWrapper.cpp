/* Copyright 2019 Luca Fedeli, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "QuantumSyncEngineWrapper.H"

#include "QedTableParserHelperFunctions.H"
#include "QuantumSyncDummyTable.H"

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
    AMREX_ALWAYS_ASSERT(m_lookup_tables_initialized);

    return QuantumSynchrotronEvolveOpticalDepth(m_innards);
}

QuantumSynchrotronGeneratePhotonAndUpdateMomentum QuantumSynchrotronEngine::build_phot_em_functor ()
{
    AMREX_ALWAYS_ASSERT(m_lookup_tables_initialized);

    return QuantumSynchrotronGeneratePhotonAndUpdateMomentum(m_innards);

}

bool QuantumSynchrotronEngine::are_lookup_tables_initialized () const
{
    return m_lookup_tables_initialized;
}

bool
QuantumSynchrotronEngine::init_lookup_tables_from_raw_data (
    const Vector<char>& raw_data)
{
    const char* p_data = raw_data.data();
    const char* const p_last = &raw_data.back();
    bool is_ok;

    //Header (control parameters)
    tie(is_ok, m_innards.ctrl.chi_part_min, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_part_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_part_tdndt_min, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_part_tdndt_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_part_tdndt_max, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_part_tdndt_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_part_tdndt_how_many, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_part_tdndt_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_part_tem_min, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_part_tem_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_part_tem_max, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_part_tem_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_part_tem_how_many, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_part_tem_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.prob_tem_how_many, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.prob_tem_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    //___________________________

    //Data
    Vector<Real> tndt_coords(m_innards.ctrl.chi_part_tdndt_how_many);
    Vector<Real> tndt_data(m_innards.ctrl.chi_part_tdndt_how_many);
    Vector<Real> cum_tab_coords1(m_innards.ctrl.chi_part_tem_how_many);
    Vector<Real> cum_tab_coords2(m_innards.ctrl.prob_tem_how_many);
    Vector<Real> cum_tab_data(m_innards.ctrl.chi_part_tem_how_many*
        m_innards.ctrl.prob_tem_how_many);

    tie(is_ok, tndt_coords, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_coords.size(), p_last);
    if(!is_ok) return false;
    m_innards.KKfunc_coords.assign(tndt_coords.begin(), tndt_coords.end());

    tie(is_ok, tndt_data, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_data.size(), p_last);
    if(!is_ok) return false;
    m_innards.KKfunc_data.assign(tndt_data.begin(), tndt_data.end());

    tie(is_ok, cum_tab_coords1, p_data) =
        parse_raw_data_vec<Real>(
            p_data, cum_tab_coords1.size(), p_last);
    if(!is_ok) return false;
    m_innards.cum_distrib_coords_1.assign(
        cum_tab_coords1.begin(), cum_tab_coords1.end());

    tie(is_ok, cum_tab_coords2, p_data) =
        parse_raw_data_vec<Real>(
            p_data, cum_tab_coords2.size(), p_last);
    if(!is_ok) return false;
    m_innards.cum_distrib_coords_2.assign(
        cum_tab_coords2.begin(), cum_tab_coords2.end());

    tie(is_ok, cum_tab_data, p_data) =
        parse_raw_data_vec<Real>(
            p_data, cum_tab_data.size(), p_last);
    if(!is_ok) return false;
    m_innards.cum_distrib_data.assign(
        cum_tab_data.begin(), cum_tab_data.end());

    //___________________________
    m_lookup_tables_initialized = true;

    return true;
}

void QuantumSynchrotronEngine::init_dummy_tables()
{
    m_innards.ctrl = QedUtils::QuantumSyncEngineInnardsDummy.ctrl;
    m_innards.KKfunc_coords.assign(
        QedUtils::QuantumSyncEngineInnardsDummy.KKfunc_coords.begin(),
        QedUtils::QuantumSyncEngineInnardsDummy.KKfunc_coords.end());
    m_innards.KKfunc_data.assign(
        QedUtils::QuantumSyncEngineInnardsDummy.KKfunc_data.begin(),
        QedUtils::QuantumSyncEngineInnardsDummy.KKfunc_data.end());
    m_innards.cum_distrib_coords_1.assign(
        QedUtils::QuantumSyncEngineInnardsDummy.cum_distrib_coords_1.begin(),
        QedUtils::QuantumSyncEngineInnardsDummy.cum_distrib_coords_1.end());
    m_innards.cum_distrib_coords_2.assign(
        QedUtils::QuantumSyncEngineInnardsDummy.cum_distrib_coords_2.begin(),
        QedUtils::QuantumSyncEngineInnardsDummy.cum_distrib_coords_2.end());
    m_innards.cum_distrib_data.assign(
        QedUtils::QuantumSyncEngineInnardsDummy.cum_distrib_data.begin(),
        QedUtils::QuantumSyncEngineInnardsDummy.cum_distrib_data.end());

    m_lookup_tables_initialized = true;
}

Vector<char> QuantumSynchrotronEngine::export_lookup_tables_data () const
{
    Vector<char> res{};

    if(!m_lookup_tables_initialized)
        return res;

    add_data_to_vector_char(&m_innards.ctrl.chi_part_min, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_part_tdndt_min, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_part_tdndt_max, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_part_tdndt_how_many, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_part_tem_min, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_part_tem_max, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_part_tem_how_many, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.prob_tem_how_many, 1, res);

    add_data_to_vector_char(m_innards.KKfunc_coords.data(),
        m_innards.KKfunc_coords.size(), res);
    add_data_to_vector_char(m_innards.KKfunc_data.data(),
        m_innards.KKfunc_data.size(), res);
    add_data_to_vector_char(m_innards.cum_distrib_coords_1.data(),
        m_innards.cum_distrib_coords_1.size(), res);
    add_data_to_vector_char(m_innards.cum_distrib_coords_2.data(),
        m_innards.cum_distrib_coords_2.size(), res);
    add_data_to_vector_char(m_innards.cum_distrib_data.data(),
        m_innards.cum_distrib_data.size(), res);

    return res;
}

PicsarQuantumSynchrotronCtrl
QuantumSynchrotronEngine::get_default_ctrl() const
{
    return PicsarQuantumSynchrotronCtrl();
}

const PicsarQuantumSynchrotronCtrl&
QuantumSynchrotronEngine::get_ref_ctrl() const
{
    return m_innards.ctrl;
}

void QuantumSynchrotronEngine::compute_lookup_tables (
    PicsarQuantumSynchrotronCtrl ctrl)
{
#ifdef WARPX_QED_TABLE_GEN
    m_table_builder.compute_table(ctrl, m_innards);
    m_lookup_tables_initialized = true;
#endif
}

//============================================
