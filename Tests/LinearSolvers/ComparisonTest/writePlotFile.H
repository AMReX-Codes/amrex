#include <Utility.H>
#include <Geometry.H>
#include <PArray.H>
#include <MultiFab.H>

void writePlotFile (const std::string& dir, 
		    const PArray<MultiFab>& soln, const PArray<MultiFab>& exac, 
		    const PArray<MultiFab>& alph, const PArray<MultiFab>& beta, 
		    const PArray<MultiFab>& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int nsoln, int iCpp, int iF90, int iHyp);
