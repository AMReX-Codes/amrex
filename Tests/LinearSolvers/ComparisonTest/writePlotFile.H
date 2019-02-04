#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

void writePlotFile (const std::string& dir, 
		    const amrex::Vector<amrex::MultiFab*>& soln,
		    const amrex::Vector<amrex::MultiFab*>& exac, 
		    const amrex::Vector<amrex::MultiFab*>& alph,
		    const amrex::Vector<amrex::MultiFab*>& beta, 
		    const amrex::Vector<amrex::MultiFab*>& rhs, 
		    const std::vector<amrex::Geometry>& geom, 
		    const std::vector<amrex::BoxArray>& grids,
		    int nsoln, int iCpp, int iHyp);
