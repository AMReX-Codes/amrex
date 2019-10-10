#include "QuantumSyncEngineWrapper.h"
//This file provides a wrapper aroud the quantum_sync engine
//provided by the PICSAR library

// Functors ==================================

// Initialization of the optical depth
AMREX_GPU_DEVICE
amrex::Real
QuantumSynchrotronGetOpticalDepth::operator() () const
{
    return WarpXQuantumSynchrotronWrapper::
        internal_get_optical_depth(amrex::Random());
}
//____________________________________________

// Evolution of the optical depth (returns true if
// an event occurs)
AMREX_GPU_DEVICE
bool QuantumSynchrotronEvolveOpticalDepth::operator()(
    amrex::Real px, amrex::Real py, amrex::Real pz, 
    amrex::Real ex, amrex::Real ey, amrex::Real ez,
    amrex::Real bx, amrex::Real by, amrex::Real bz,
    amrex::Real dt, amrex::Real& opt_depth) const
{
    bool has_event_happend = false;
    amrex::Real dummy_lambda = 1.0;
    amrex::Real unused_event_time = 0.0;

    const auto table = picsar::multi_physics::lookup_1d<amrex::Real>
        (innards->KKfunc_data.size(),
         innards->KKfunc_coords.data(),
         innards->KKfunc_data.data());

    WarpXQuantumSynchrotronWrapper::internal_evolve_opt_depth_and_determine_event(
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

//Initializes the Lookup tables using the default settings 
//provided by the library
void QuantumSynchrotronEngine::computes_lookup_tables_default ()
{
    //A control parameters structure
    //with the default values provided by the library
    WarpXQuantumSynchrotronWrapperCtrl ctrl_default;

    computes_lookup_tables(ctrl_default);

    lookup_tables_initialized = true;
}

bool QuantumSynchrotronEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}

//Private function which actually computes the lookup tables
void QuantumSynchrotronEngine::computes_lookup_tables (
    WarpXQuantumSynchrotronWrapperCtrl ctrl)
{
    //Lambda is not actually used if S.I. units are enabled
    WarpXQuantumSynchrotronWrapper qs_engine(std::move(DummyStruct()), 1.0, ctrl);
    
    qs_engine.compute_dN_dt_lookup_table();
    //qs_engine.compute_cumulative_pair_table();

    auto qs_innards_picsar = qs_engine.export_innards();  

    //Copy data in a GPU-friendly data-structure
    innards.ctrl = qs_innards_picsar.qs_ctrl;
    innards.KKfunc_coords.assign(qs_innards_picsar.KKfunc_table_coords_ptr, 
        qs_innards_picsar.KKfunc_table_coords_ptr + 
        qs_innards_picsar.KKfunc_table_coords_how_many);
    innards.KKfunc_data.assign(qs_innards_picsar.KKfunc_table_data_ptr,
        qs_innards_picsar.KKfunc_table_data_ptr + 
        qs_innards_picsar.KKfunc_table_data_how_many);
    //____
}

//============================================
