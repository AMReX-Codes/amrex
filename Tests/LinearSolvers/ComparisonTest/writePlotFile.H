#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

void writePlotFile (const std::string& dir, 
		    const amrex::Array<amrex::MultiFab*>& soln,
		    const amrex::Array<amrex::MultiFab*>& exac, 
		    const amrex::Array<amrex::MultiFab*>& alph,
		    const amrex::Array<amrex::MultiFab*>& beta, 
		    const amrex::Array<amrex::MultiFab*>& rhs, 
		    const std::vector<amrex::Geometry>& geom, 
		    const std::vector<amrex::BoxArray>& grids,
		    int nsoln, int iCpp, int iF90, int iHyp);
