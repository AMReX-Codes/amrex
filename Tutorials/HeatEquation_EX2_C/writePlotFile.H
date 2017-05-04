#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>

void
writePlotFile (const std::string& dir, 
               const amrex::MultiFab&    mf, 
               const amrex::Geometry&    geom,
	       const amrex::Real&        time);
