#ifndef WARPX_quantum_sync_engine_wrapper_h_
#define WARPX_quantum_sync_engine_wrapper_h_

#include "QedWrapperCommons.H"

#include <AMReX_Array.H>

#include<string>
#include<fstream>

//QS ENGINE from PICSAR
#define PXRMP_CORE_ONLY
#include "quantum_sync_engine.hpp"

using WarpXQuantumSynchrotronWrapper =
    picsar::multi_physics::quantum_synchrotron_engine<amrex::Real, DummyStruct>;

using WarpXQuantumSynchrotronWrapperCtrl =
    picsar::multi_physics::quantum_synchrotron_engine_ctrl<amrex::Real>;

// Struct to hold engine data ================

struct QuantumSynchrotronEngineInnards
{
    // Control parameters
    WarpXQuantumSynchrotronWrapperCtrl ctrl;

    //Lookup table data
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_coords;
    amrex::Gpu::ManagedDeviceVector<amrex::Real> KKfunc_data;

    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_1;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_coords_2;
    amrex::Gpu::ManagedVector<amrex::Real> cum_distrib_data;
    //______
};

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
        return WarpXQuantumSynchrotronWrapper::
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
        QuantumSynchrotronEngineInnards* _innards):
        innards{_innards}{};

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
            innards->KKfunc_data.size(),
            innards->KKfunc_coords.data(),
            innards->KKfunc_data.data());

        WarpXQuantumSynchrotronWrapper::
        internal_evolve_opt_depth_and_determine_event(
            px, py, pz,
            ex, ey, ez,
            bx, by, bz,
            dt, opt_depth,
            has_event_happend, unused_event_time,
            dummy_lambda,
            table,
            innards->ctrl);

        return has_event_happend;
    }

private:
    QuantumSynchrotronEngineInnards* innards;
};

// Generates a photon via the Quantum Synchrotron process
// and updates momentum accordingly (returns false if errors occur)
class QuantumSynchrotronGeneratePhotonAndUpdateMomentum
{
public:
    QuantumSynchrotronGeneratePhotonAndUpdateMomentum(
        QuantumSynchrotronEngineInnards* _innards):
        innards{_innards}{};

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
        picsar::multi_physics::picsar_array<amrex::Real, sampling>
            rand_zero_one_minus_epsi;
        for(auto& el : rand_zero_one_minus_epsi)
            el = amrex::Random();

        const auto table = picsar::multi_physics::lookup_2d<amrex::Real>
            (innards->cum_distrib_coords_1.size(),
            innards->cum_distrib_coords_1.data(),
            innards->cum_distrib_coords_2.size(),
            innards->cum_distrib_coords_2.data(),
            innards->cum_distrib_data.data());

        bool stat = WarpXQuantumSynchrotronWrapper::
        internal_generate_photons_and_update_momentum(
            *px, *py, *pz,
            ex, ey, ez,
            bx, by, bz,
            weight, sampling,
            *g_px, *g_py, *g_pz,
            *g_weight,
            dummy_lambda,
            table,
            innards->ctrl,
            rand_zero_one_minus_epsi.data());

        return stat;
    }

private:
    QuantumSynchrotronEngineInnards* innards;
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
    void read_lookup_tables (std::string file);

private:
    bool lookup_tables_initialized = false;

    QuantumSynchrotronEngineInnards innards;

#ifdef WARPX_QED_TABLE_GEN
public:
    /* \brief Computes the Lookup tables using the default settings
     *  provided by the PICSAR library */
    void compute_lookup_tables_default ();

    /* \brief Computes the Lookup tables using user-defined settings */
    void compute_custom_lookup_tables (WarpXQuantumSynchrotronWrapperCtrl ctrl);
        /* \brief Writes lookup tables on disk in 'file'
     *  return false if it fails. */
    bool write_lookup_tables (std::string file) const;

private:
    //Private function which actually computes the lookup tables
    void computes_lookup_tables (
        WarpXQuantumSynchrotronWrapperCtrl ctrl);

#endif
};

//============================================

#endif //WARPX_quantum_sync_engine_wrapper_h_
