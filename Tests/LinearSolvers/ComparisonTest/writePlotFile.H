#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

void writePlotFile (const std::string& dir, 
		    const Array<MultiFab*>& soln,
		    const Array<MultiFab*>& exac, 
		    const Array<MultiFab*>& alph,
		    const Array<MultiFab*>& beta, 
		    const Array<MultiFab*>& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int nsoln, int iCpp, int iF90, int iHyp);
