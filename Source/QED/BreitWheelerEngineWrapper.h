#ifndef WARPX_breit_wheeler_engine_wrapper_h_
#define WARPX_breit_wheeler_engine_wrapper_h_

#include "QedWrapperCommons.h"

#include<string>
#include<vector>
#include<tuple>

//BW ENGINE from PICSAR
#include "breit_wheeler_engine.hpp"

using WarpXBreitWheelerWrapper =
    picsar::multi_physics::breit_wheeler_engine<amrex::Real, DummyStruct>;

using WarpXBreitWheelerWrapperCtrl =
    picsar::multi_physics::breit_wheeler_engine_ctrl<amrex::Real>;

// Struct to hold engine data ================

struct BreitWheelerEngineInnards
{
    // Control parameters
    WarpXBreitWheelerWrapperCtrl ctrl;

    //Lookup table data
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_coords;
    amrex::Gpu::ManagedVector<amrex::Real> TTfunc_data; 
    //______
};

// Functors ==================================

// These functors provide the core elementary functions of the library
// Can be included in GPU kernels

/* \brief Functor to initialize the optical depth of photons for the
*   Breit-Wheeler process */
class BreitWheelerGetOpticalDepth
{
public:
    BreitWheelerGetOpticalDepth ()
    {};

    AMREX_GPU_DEVICE
    amrex::Real operator() () const;
};
//____________________________________________

// Evolution of the optical depth (returns true if
// an event occurs)
class BreitWheelerEvolveOpticalDepth
{
public:
    BreitWheelerEvolveOpticalDepth(
        BreitWheelerEngineInnards* _innards):
        innards{_innards}{};

    AMREX_GPU_DEVICE
    bool operator()(
    amrex::Real px, amrex::Real py, amrex::Real pz, 
    amrex::Real ex, amrex::Real ey, amrex::Real ez,
    amrex::Real bx, amrex::Real by, amrex::Real bz,
    amrex::Real dt, amrex::Real& opt_depth) const;

private:
    BreitWheelerEngineInnards* innards;
};

// Factory class =============================

/* \brief Wrapper for the Breit Wheeler engine of the PICSAR library */
class BreitWheelerEngine
{
public:
    BreitWheelerEngine ();

    /* \brief Builds the functor to initialize the optical depth */
    BreitWheelerGetOpticalDepth build_optical_depth_functor ();
    
    /* \brief Builds the functor to evolve the optical depth */
    BreitWheelerEvolveOpticalDepth build_evolve_functor ();

    /* \brief Computes the Lookup tables using the default settings 
     *  provided by the PICSAR library */
    void computes_lookup_tables_default ();

    /* \brief Checks if lookup tables are properly initialized */
    bool are_lookup_tables_initialized () const;

    /* \brief Writes lookup tables on disk in 'file'
     *  return false if it fails. */
    bool write_lookup_tables (std::string file) const;   

private:
    bool lookup_tables_initialized = false;

    BreitWheelerEngineInnards innards;

    //Private function which actually computes the lookup tables
    void computes_lookup_tables (
        WarpXBreitWheelerWrapperCtrl ctrl);
};

//============================================

#endif //WARPX_breit_wheeler_engine_wrapper_H_