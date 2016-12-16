#include <fstream>
#include <iomanip>

#include <AMReX_Geometry.H>
#include <AMReX_VisMF.H>

void
writePlotFile (const std::string& dir, 
               const MultiFab&    mf, 
               const Geometry&    geom,
	       const Real&        time);
