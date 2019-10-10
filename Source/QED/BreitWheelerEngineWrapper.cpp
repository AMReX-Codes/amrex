#include "BreitWheelerEngineWrapper.h"
//This file provides a wrapper aroud the breit_wheeler engine
//provided by the PICSAR library

// Functors ==================================

// Initialization of the optical depth
AMREX_GPU_DEVICE
amrex::Real
BreitWheelerGetOpticalDepth::operator() () const
{
    return WarpXBreitWheelerWrapper::
        internal_get_optical_depth(amrex::Random());
}
//____________________________________________

// Evolution of the optical depth (returns true if
// an event occurs)
AMREX_GPU_DEVICE
bool BreitWheelerEvolveOpticalDepth::operator()(
    amrex::Real px, amrex::Real py, amrex::Real pz, 
    amrex::Real ex, amrex::Real ey, amrex::Real ez,
    amrex::Real bx, amrex::Real by, amrex::Real bz,
    amrex::Real dt, amrex::Real& opt_depth) const
{
    bool has_event_happend = false;
    amrex::Real dummy_lambda = 1.0;
    amrex::Real unused_event_time = 0.0;

    const auto table = picsar::multi_physics::lookup_1d<amrex::Real>
        (innards->TTfunc_data.size(),
         innards->TTfunc_coords.data(),
         innards->TTfunc_data.data());

    WarpXBreitWheelerWrapper::internal_evolve_opt_depth_and_determine_event(
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


// Writes lookup tables on disk in 'folder'
// return false if it fails. */
bool BreitWheelerEngine::write_lookup_tables (
        std::string folder) const
{
    if(!lookup_tables_initialized)
        return false;

    auto all_data = make_tuple(
        std::ref(innards.ctrl.chi_phot_min),
        std::ref(innards.ctrl.chi_phot_tdndt_min),
        std::ref(innards.ctrl.chi_phot_tdndt_max),
        std::ref(innards.ctrl.chi_phot_tdndt_how_many),
        std::ref(innards.ctrl.chi_phot_tpair_min),
        std::ref(innards.ctrl.chi_phot_tpair_max),
        std::ref(innards.ctrl.chi_phot_tpair_how_many),
        std::ref(innards.ctrl.chi_frac_tpair_how_many),
        std::ref(innards.TTfunc_coords),
        std::ref(innards.TTfunc_data)); 

    
       

    char* data_dump =  new char(buf_size);

    size_t count = 0;
    auto copy_and_advance = [&count] (char* source, char*dest, size_t size) {
        count += size;
    };    
    
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
