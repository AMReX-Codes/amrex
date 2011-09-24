#include <Utility.H>

void writePlotFile (const std::string& dir, 
		    const PArray<MultiFab>& soln, const PArray<MultiFab>& exac, 
		    const PArray<MultiFab>& alph, const PArray<MultiFab>& beta, 
		    const PArray<MultiFab>& rhs, Geometry geom[], BoxArray bas[],
		    int nsoln, int iF90, int iCpp);
