/* Copyright 2019 Luca Fedeli, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "BreitWheelerEngineWrapper.H"

#include "QedTableParserHelperFunctions.H"
#include "BreitWheelerDummyTable.H"

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
    AMREX_ALWAYS_ASSERT(m_lookup_tables_initialized);

    return BreitWheelerEvolveOpticalDepth(m_innards);
}

BreitWheelerGeneratePairs
BreitWheelerEngine::build_pair_functor ()
{
    AMREX_ALWAYS_ASSERT(m_lookup_tables_initialized);

    return BreitWheelerGeneratePairs(m_innards);
}

bool BreitWheelerEngine::are_lookup_tables_initialized () const
{
    return m_lookup_tables_initialized;
}

bool
BreitWheelerEngine::init_lookup_tables_from_raw_data (
    const Vector<char>& raw_data)
{
    const char* p_data = raw_data.data();
    const char* const p_last = &raw_data.back();
    bool is_ok;

    //Header (control parameters)
    tie(is_ok, m_innards.ctrl.chi_phot_min, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_phot_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_phot_tdndt_min, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_phot_tdndt_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_phot_tdndt_max, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_phot_tdndt_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_phot_tdndt_how_many, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_phot_tdndt_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_phot_tpair_min, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_phot_tpair_min)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_phot_tpair_max, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_phot_tpair_max)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_phot_tpair_how_many, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_phot_tpair_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    tie(is_ok, m_innards.ctrl.chi_frac_tpair_how_many, p_data) =
        parse_raw_data<decltype(m_innards.ctrl.chi_frac_tpair_how_many)>(
            p_data, p_last);
    if(!is_ok) return false;

    //___________________________

    //Data
    Vector<Real> tndt_coords(m_innards.ctrl.chi_phot_tdndt_how_many);
    Vector<Real> tndt_data(m_innards.ctrl.chi_phot_tdndt_how_many);
    Vector<Real> cum_tab_coords1(m_innards.ctrl.chi_phot_tpair_how_many);
    Vector<Real> cum_tab_coords2(m_innards.ctrl.chi_frac_tpair_how_many);
    Vector<Real> cum_tab_data(m_innards.ctrl.chi_phot_tpair_how_many*
        m_innards.ctrl.chi_frac_tpair_how_many);

    tie(is_ok, tndt_coords, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_coords.size(), p_last);
    if(!is_ok) return false;
    m_innards.TTfunc_coords.assign(tndt_coords.begin(), tndt_coords.end());

    tie(is_ok, tndt_data, p_data) =
        parse_raw_data_vec<Real>(
            p_data, tndt_data.size(), p_last);
    if(!is_ok) return false;
    m_innards.TTfunc_data.assign(tndt_data.begin(), tndt_data.end());

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

void BreitWheelerEngine::init_dummy_tables()
{
    m_innards.ctrl = QedUtils::BreitWheelerEngineInnardsDummy.ctrl;
    m_innards.TTfunc_coords.assign(
        QedUtils::BreitWheelerEngineInnardsDummy.TTfunc_coords.begin(),
        QedUtils::BreitWheelerEngineInnardsDummy.TTfunc_coords.end());
    m_innards.TTfunc_data.assign(
        QedUtils::BreitWheelerEngineInnardsDummy.TTfunc_data.begin(),
        QedUtils::BreitWheelerEngineInnardsDummy.TTfunc_data.end());
    m_innards.cum_distrib_coords_1.assign(
        QedUtils::BreitWheelerEngineInnardsDummy.cum_distrib_coords_1.begin(),
        QedUtils::BreitWheelerEngineInnardsDummy.cum_distrib_coords_1.end());
    m_innards.cum_distrib_coords_2.assign(
        QedUtils::BreitWheelerEngineInnardsDummy.cum_distrib_coords_2.begin(),
        QedUtils::BreitWheelerEngineInnardsDummy.cum_distrib_coords_2.end());
    m_innards.cum_distrib_data.assign(
        QedUtils::BreitWheelerEngineInnardsDummy.cum_distrib_data.begin(),
        QedUtils::BreitWheelerEngineInnardsDummy.cum_distrib_data.end());

    m_lookup_tables_initialized = true;
}

Vector<char> BreitWheelerEngine::export_lookup_tables_data () const
{
   Vector<char> res{};

    if(!m_lookup_tables_initialized)
        return res;

    add_data_to_vector_char(&m_innards.ctrl.chi_phot_min, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_phot_tdndt_min, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_phot_tdndt_max, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_phot_tdndt_how_many, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_phot_tpair_min, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_phot_tpair_max, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_phot_tpair_how_many, 1, res);
    add_data_to_vector_char(&m_innards.ctrl.chi_frac_tpair_how_many, 1, res);

    add_data_to_vector_char(m_innards.TTfunc_coords.data(),
        m_innards.TTfunc_coords.size(), res);
    add_data_to_vector_char(m_innards.TTfunc_data.data(),
        m_innards.TTfunc_data.size(), res);
    add_data_to_vector_char(m_innards.cum_distrib_coords_1.data(),
        m_innards.cum_distrib_coords_1.size(), res);
    add_data_to_vector_char(m_innards.cum_distrib_coords_2.data(),
        m_innards.cum_distrib_coords_2.size(), res);
    add_data_to_vector_char(m_innards.cum_distrib_data.data(),
        m_innards.cum_distrib_data.size(), res);

    return res;
}

PicsarBreitWheelerCtrl
BreitWheelerEngine::get_default_ctrl() const
{
    return PicsarBreitWheelerCtrl();
}

const PicsarBreitWheelerCtrl&
BreitWheelerEngine::get_ref_ctrl() const
{
    return m_innards.ctrl;
}

void BreitWheelerEngine::compute_lookup_tables (
    PicsarBreitWheelerCtrl ctrl)
{
#ifdef WARPX_QED_TABLE_GEN
    m_table_builder.compute_table(ctrl, m_innards);
    m_lookup_tables_initialized = true;
#endif
}

//============================================
