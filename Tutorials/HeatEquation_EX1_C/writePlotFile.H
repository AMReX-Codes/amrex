#include <fstream>
#include <iomanip>

#include <Geometry.H>
#include <VisMF.H>

void
writePlotFile (const std::string& dir, 
               const MultiFab&    mf, 
               const Geometry&    geom,
	       const Real&        time);
