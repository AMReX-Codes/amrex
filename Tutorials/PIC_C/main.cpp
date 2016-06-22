#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MGT_Solver.H>
#include <mg_cpp_f.h>
#include <stencil_types.H>
#include <MultiFabUtil.H>

#include "Particles.H"

// declare routines below
void single_level(int nlevs, int nx, int ny, int nz, int max_grid_size, int nppc, bool verbose);
void    two_level(int nlevs, int nx, int ny, int nz, int max_grid_size, int nppc, bool verbose);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    const Real strt_total = ParallelDescriptor::second();

    int    nx, ny, nz;
    double hx, hy, hz;

    ParmParse pp;

    pp.get("nx", nx);
    pp.get("ny", ny);
    pp.get("nz", nz);

    int max_level;
    int max_grid_size;

    pp.get("max_grid_size", max_grid_size);
    pp.get("max_level", max_level);

    int nlevs = max_level + 1;

    // Number of particles per cell
    int nppc = -1;
    pp.get("nppc", nppc);

    if (nppc < 1 && ParallelDescriptor::IOProcessor())
       BoxLib::Abort("Must specify at least one particle per cell");

    bool verbose = false;
    pp.query("verbose", verbose);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
       std::cout << "                              " << std::endl;
       std::cout << "Number of levels             : " << nlevs << std::endl;
       std::cout << "Number of particles per cell : " << nppc  << std::endl;
       std::cout << "Size of domain               : " << nx << " " << ny << " " << nz << std::endl;
    }

    Real strt_single, end_single;
    if (nlevs == 1) {
       single_level(nlevs,nx,ny,nz,max_grid_size,nppc,verbose);
    } else if (nlevs == 2) {
          two_level(nlevs,nx,ny,nz,max_grid_size,nppc,verbose);
    } else {
       BoxLib::Abort("Right now we only take max_level = 0 or 1");
    }

    Real end_total = ParallelDescriptor::second() - strt_total;

    ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
    if (verbose && ParallelDescriptor::IOProcessor()) 
       std::cout << "Total Time           : " << end_total << '\n' << '\n';
;

    BoxLib::Finalize();
}
