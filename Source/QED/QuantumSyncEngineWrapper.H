#ifndef WARPX_quantum_sync_engine_wrapper_h_
#define WARPX_quantum_sync_engine_wrapper_h_

#include "QedWrapperCommons.H"
#include "QuantumSyncEngineInnards.H"

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Gpu.H>

//QS ENGINE from PICSAR
//Include only essential parts of the library in this file
#define PXRMP_CORE_ONLY
#include <quantum_sync_engine.hpp>

#ifdef WARPX_QED_TABLE_GEN
    #include "QuantumSyncEngineTableBuilder.H"
#endif

#include <string>

using PicsarQuantumSynchrotronEngine = picsar::multi_physics::
    quantum_synchrotron_engine<amrex::Real, QedUtils::DummyStruct>;

using PicsarQuantumSynchrotronCtrl =
    picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real>;


// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

/**
*  Functor to initialize the optical depth of leptons for the
*  Quantum Synchrotron process
*/
class QuantumSynchrotronGetOpticalDepth
{
public:
    QuantumSynchrotronGetOpticalDepth ()
    {};

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    amrex::Real operator() () const
    {
        return PicsarQuantumSynchrotronEngine::
            internal_get_optical_depth(amrex::Random());
    }
};
//____________________________________________

// Evolution of the optical depth (returns true if
// an event occurs)
class QuantumSynchrotronEvolveOpticalDepth
{
public:
    QuantumSynchrotronEvolveOpticalDepth(
        QuantumSynchrotronEngineInnards* t_innards):
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

        const auto table = picsar::multi_physics
        ::lookup_1d<amrex::Real>(
            p_innards->KKfunc_data.size(),
            p_innards->KKfunc_coords.data(),
            p_innards->KKfunc_data.data());

        PicsarQuantumSynchrotronEngine::
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
    QuantumSynchrotronEngineInnards* p_innards;
};

// Generates a photon via the Quantum Synchrotron process
// and updates momentum accordingly (returns false if errors occur)
class QuantumSynchrotronGeneratePhotonAndUpdateMomentum
{
public:
    QuantumSynchrotronGeneratePhotonAndUpdateMomentum(
        QuantumSynchrotronEngineInnards* t_innards):
        p_innards{t_innards}{};

    template <size_t sampling>
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    bool operator()(
    amrex::Real* px, amrex::Real* py, amrex::Real* pz,
    amrex::Real ex, amrex::Real ey, amrex::Real ez,
    amrex::Real bx, amrex::Real by, amrex::Real bz,
    amrex::Real weight,
    amrex::Real* g_px, amrex::Real* g_py, amrex::Real* g_pz,
    amrex::Real* g_weight) const
    {
        amrex::Real dummy_lambda{1.0};
        amrex::GpuArray<amrex::Real, sampling>
            rand_zero_one_minus_epsi;
        for(auto& el : rand_zero_one_minus_epsi)
            el = amrex::Random();

        const auto table = picsar::multi_physics::lookup_2d<amrex::Real>
            (p_innards->cum_distrib_coords_1.size(),
            p_innards->cum_distrib_coords_1.data(),
            p_innards->cum_distrib_coords_2.size(),
            p_innards->cum_distrib_coords_2.data(),
            p_innards->cum_distrib_data.data());

        bool stat = PicsarQuantumSynchrotronEngine::
        internal_generate_photons_and_update_momentum(
            *px, *py, *pz,
            ex, ey, ez,
            bx, by, bz,
            weight, sampling,
            *g_px, *g_py, *g_pz,
            *g_weight,
            dummy_lambda,
            table,
            p_innards->ctrl,
            rand_zero_one_minus_epsi.data());

        return stat;
    }

private:
    QuantumSynchrotronEngineInnards* p_innards;
};

// Factory class =============================

/**
 * \brief Wrapper for the Quantum Synchrotron engine of the PICSAR library
 */
class QuantumSynchrotronEngine
{
public:
    QuantumSynchrotronEngine ();

    /**
     * \brief Builds the functor to initialize the optical depth
     */
    QuantumSynchrotronGetOpticalDepth build_optical_depth_functor ();

    /* \brief Builds the functor to evolve the optical depth */
    QuantumSynchrotronEvolveOpticalDepth build_evolve_functor ();

    /* \brief Checks if lookup tables are properly initialized */
    bool are_lookup_tables_initialized () const;

      /* \brief Reads lookup tables from 'file' on disk */
    bool init_lookup_tables_from_raw_data (const amrex::Vector<char>& raw_data);

    /* \brief Writes lookup tables on disk in 'file'
     *  return false if it fails. */
    amrex::Vector<char> export_lookup_tables_data () const;

    /* \brief Computes the Lookup tables using user-defined settings */
    void compute_lookup_tables (PicsarQuantumSynchrotronCtrl ctrl);

    PicsarQuantumSynchrotronCtrl get_default_ctrl() const;

private:
    bool lookup_tables_initialized = false;

    QuantumSynchrotronEngineInnards innards;

#ifdef WARPX_QED_TABLE_GEN
    QuantumSynchrotronEngineTableBuilder table_builder;
#endif

};

//============================================

#endif //WARPX_quantum_sync_engine_wrapper_h_
