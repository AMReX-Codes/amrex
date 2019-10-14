#include "BreitWheelerEngineWrapper.H"
//This file provides a wrapper aroud the breit_wheeler engine
//provided by the PICSAR library

// Factory class =============================

BreitWheelerEngine::BreitWheelerEngine (){}

//Builds the functor to initialize the optical depth
BreitWheelerGetOpticalDepth BreitWheelerEngine::build_optical_depth_functor ()
{
    return BreitWheelerGetOpticalDepth();
}

//Builds the functor to evolve the optical depth
BreitWheelerEvolveOpticalDepth BreitWheelerEngine::build_evolve_functor ()
{
    AMREX_ALWAYS_ASSERT(lookup_tables_initialized);

    return BreitWheelerEvolveOpticalDepth(&innards);
}


//Initializes the Lookup tables using the default settings
//provided by the library
void BreitWheelerEngine::computes_lookup_tables_default ()
{
    //A control parameters structure
    //with the default values provided by the library
    WarpXBreitWheelerWrapperCtrl ctrl_default;

    computes_lookup_tables(ctrl_default);

    lookup_tables_initialized = true;
}

bool BreitWheelerEngine::are_lookup_tables_initialized () const
{
    return lookup_tables_initialized;
}



/* \brief Writes lookup tables on disk in 'file'
 *  return false if it fails. */
bool BreitWheelerEngine::write_lookup_tables (
        std::string file) const
{
    if(!lookup_tables_initialized)
        return false;    

    return true;
}

//Private function which actually computes the lookup tables
void BreitWheelerEngine::computes_lookup_tables (
    WarpXBreitWheelerWrapperCtrl ctrl)
{
    //Lambda is not actually used if S.I. units are enabled
    WarpXBreitWheelerWrapper bw_engine(std::move(DummyStruct()), 1.0, ctrl);

    bw_engine.compute_dN_dt_lookup_table();
    //bw_engine.compute_cumulative_pair_table();

    auto bw_innards_picsar = bw_engine.export_innards();

    //Copy data in a GPU-friendly data-structure
    innards.ctrl = bw_innards_picsar.bw_ctrl;
    innards.TTfunc_coords.assign(bw_innards_picsar.TTfunc_table_coords_ptr,
        bw_innards_picsar.TTfunc_table_coords_ptr +
        bw_innards_picsar.TTfunc_table_coords_how_many);
    innards.TTfunc_data.assign(bw_innards_picsar.TTfunc_table_data_ptr,
        bw_innards_picsar.TTfunc_table_data_ptr +
        bw_innards_picsar.TTfunc_table_data_how_many);
    //____
}

//============================================
