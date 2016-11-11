#include <Utility.H>
#include <Geometry.H>
#include <MultiFab.H>

void writePlotFile (const std::string& dir, 
		    const Array<std::unique_ptr<MultiFab> >& soln,
		    const Array<std::unique_ptr<MultiFab> >& exac, 
		    const Array<std::unique_ptr<MultiFab> >& alph,
		    const Array<std::unique_ptr<MultiFab> >& beta, 
		    const Array<std::unique_ptr<MultiFab> >& rhs, 
		    const std::vector<Geometry>& geom, 
		    const std::vector<BoxArray>& grids,
		    int nsoln, int iCpp, int iF90, int iHyp);
