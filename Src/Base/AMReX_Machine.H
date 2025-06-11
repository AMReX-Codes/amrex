#ifndef AMREX_MACHINE_H
#define AMREX_MACHINE_H
#include <AMReX_Config.H>

#include <string>

namespace amrex::Machine {

void Initialize (); //!< called in amrex::Initialize()

void Finalize ();

std::string const& name ();

}

#endif
