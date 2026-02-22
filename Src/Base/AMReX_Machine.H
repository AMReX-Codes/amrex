#ifndef AMREX_MACHINE_H
#define AMREX_MACHINE_H
#include <AMReX_Config.H>

#include <string>

/**
 * \file AMReX_Machine.H
 * \brief Helpers for querying machine/host metadata.
 */

namespace amrex::Machine {

//! Gather host metadata (hostname etc.). Called automatically by amrex::Initialize().
void Initialize ();

//! Release cached machine information.
void Finalize ();

//! Return the cached machine/host name for display/logging.
std::string const& name ();

}

#endif
