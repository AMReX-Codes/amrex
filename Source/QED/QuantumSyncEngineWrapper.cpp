#include "QuantumSyncEngineWrapper.H"
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
