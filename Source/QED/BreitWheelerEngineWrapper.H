#ifndef WARPX_breit_wheeler_engine_wrapper_h_
#define WARPX_breit_wheeler_engine_wrapper_h_

#include "QedWrapperCommons.H"
#include "BreitWheelerEngineInnards.H"

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

//BW ENGINE from PICSAR
//Include only essential parts of the library in this file
#define PXRMP_CORE_ONLY
#include <breit_wheeler_engine.hpp>

#ifdef WARPX_QED_TABLE_GEN
    #include "BreitWheelerEngineTableBuilder.H"
#endif

#include <string>

using PicsarBreitWheelerEngine = picsar::multi_physics::
    breit_wheeler_engine<amrex::Real, QedUtils::DummyStruct>;

using PicsarBreitWheelerCtrl =
    picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real>;

// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

/**
 * \brief Functor to initialize the optical depth of photons for the
 * Breit-Wheeler process
 */
class BreitWheelerGetOpticalDepth
{
public:
    BreitWheelerGetOpticalDepth ()
    {};

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real operator() () const
    {
        return PicsarBreitWheelerEngine::
            internal_get_optical_depth(amrex::Random());
    }
};
//____________________________________________

// Evolution of the optical depth (returns true if
// an event occurs)
class BreitWheelerEvolveOpticalDepth
{
public:
    BreitWheelerEvolveOpticalDepth(
        BreitWheelerEngineInnards* t_innards):
        p_innards{t_innards}{};

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    bool operator()(
    amrex::Real px, amrex::Real py, amrex::Real pz,
    amrex::Real ex, amrex::Real ey, amrex::Real ez,
    amrex::Real bx, amrex::Real by, amrex::Real bz,
    amrex::Real dt, amrex::Real& opt_depth) const
    {
        bool has_event_happend{false};
        amrex::Real dummy_lambda{1.0};
        amrex::Real unused_event_time{0.0};

        const auto table = picsar::multi_physics::lookup_1d<amrex::Real>
            (p_innards->TTfunc_data.size(),
            p_innards->TTfunc_coords.data(),
            p_innards->TTfunc_data.data());

        PicsarBreitWheelerEngine::
        internal_evolve_opt_depth_and_determine_event(
            px, py, pz,
            ex, ey, ez,
            bx, by, bz,
            dt, opt_depth,
            has_event_happend, unused_event_time,
            dummy_lambda,
            table,
            p_innards->ctrl);

        return has_event_happend;
    }

private:
    BreitWheelerEngineInnards* p_innards;
};

// Generates an electron-positron pair via the Breit-Wheeler process
// (returns false if errors occur)
class BreitWheelerGeneratePairs
{
public:
    BreitWheelerGeneratePairs(
        BreitWheelerEngineInnards* t_innards):
        p_innards{t_innards}{};

    template <size_t sampling>
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    bool operator()(
    amrex::Real px, amrex::Real py, amrex::Real pz,
    amrex::Real ex, amrex::Real ey, amrex::Real ez,
    amrex::Real bx, amrex::Real by, amrex::Real bz,
    amrex::Real weight,
    amrex::Real* e_px, amrex::Real* e_py, amrex::Real* e_pz,
    amrex::Real* p_px, amrex::Real* p_py, amrex::Real* p_pz,
    amrex::Real* e_weight, amrex::Real* p_weight) const
    {
        amrex::Real dummy_lambda{1.0};
        picsar::multi_physics::picsar_array<amrex::Real, sampling>
            rand_zero_one_minus_epsi;
        for(auto& el : rand_zero_one_minus_epsi)
            el = amrex::Random();

        const auto table = picsar::multi_physics::lookup_2d<amrex::Real>
            (p_innards->cum_distrib_coords_1.size(),
            p_innards->cum_distrib_coords_1.data(),
            p_innards->cum_distrib_coords_2.size(),
            p_innards->cum_distrib_coords_2.data(),
            p_innards->cum_distrib_data.data());

        bool stat = PicsarBreitWheelerEngine::
        internal_generate_breit_wheeler_pairs(
            px, py, pz,
            ex, ey, ez,
            bx, by, bz,
            weight, sampling,
            e_px, e_py, e_pz,
            p_px, p_py, p_pz,
            e_weight, p_weight,
            dummy_lambda,
            table,
            p_innards->ctrl,
            rand_zero_one_minus_epsi.data());

        return stat;
    }

private:
    BreitWheelerEngineInnards* p_innards;
};

// Factory class =============================

/**
 * \brief Wrapper for the Breit Wheeler engine of the PICSAR library
 */
class BreitWheelerEngine
{
public:
    BreitWheelerEngine ();

    /**
     * \brief Builds the functor to initialize the optical depth
     */
    BreitWheelerGetOpticalDepth build_optical_depth_functor ();

    /* \brief Builds the functor to evolve the optical depth */
    BreitWheelerEvolveOpticalDepth build_evolve_functor ();

    /* \brief Builds the functor to generate the pairs */
    BreitWheelerGeneratePairs build_pair_functor ();

    /* \brief Checks if lookup tables are properly initialized */
    bool are_lookup_tables_initialized () const;

    /* \brief Reads lookup tables from 'file' on disk */
    bool init_lookup_tables_from_raw_data (const amrex::Vector<char>& raw_data);

    /* \brief Writes lookup tables on disk in 'file'
     *  return false if it fails. */
    amrex::Vector<char> export_lookup_tables_data () const;

    /* \brief Computes the Lookup tables using user-defined settings */
    void compute_lookup_tables (PicsarBreitWheelerCtrl ctrl);

    PicsarBreitWheelerCtrl get_default_ctrl() const;

private:
    bool lookup_tables_initialized = false;

    BreitWheelerEngineInnards innards;

#ifdef WARPX_QED_TABLE_GEN
    BreitWheelerEngineTableBuilder table_builder;
#endif

};

//============================================

#endif //WARPX_breit_wheeler_engine_wrapper_H_
