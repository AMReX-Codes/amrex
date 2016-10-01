#include <iostream>

#include <BoxLib.H>
#include <MultiFab.H>
#include <MultiFabUtil.H>
#include <BLFort.H>
#include <MacBndry.H>
#include <MultiFabUtil.H>

#include "Particles.H"

// declare routines below
void single_level(int nlevs, int nx, int ny, int nz, int max_grid_size, int order, bool verbose);

int main(int argc, char* argv[])
{
    BoxLib::Initialize(argc,argv);

    const Real strt_total = ParallelDescriptor::second();

    int    nx, ny, nz;

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
    int nppc = 1;

    if (nppc < 1 && ParallelDescriptor::IOProcessor())
       BoxLib::Abort("Must specify at least one particle per cell");

    // Order of charge deposition routine
    int order = 1;
    pp.query("order", order);

    if ( (order < 1 || order > 3) && ParallelDescriptor::IOProcessor())
       BoxLib::Abort("The order must be 1 <= order <= 3");

    bool verbose = false;
    pp.query("verbose", verbose);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
       std::cout << "                              " << std::endl;
       std::cout << "Number of levels             : " << nlevs << std::endl;
       std::cout << "Number of particles per cell : " << nppc  << std::endl;
       std::cout << "Size of domain               : " << nx << " " << ny << " " << nz << std::endl;
       std::cout << "Order of charge deposition   : " << order << std::endl;
    }

    single_level(nlevs,nx,ny,nz,max_grid_size,order,verbose);

    Real end_total = ParallelDescriptor::second() - strt_total;

    ParallelDescriptor::ReduceRealMax(end_total ,ParallelDescriptor::IOProcessorNumber());
    if (verbose && ParallelDescriptor::IOProcessor()) 
       std::cout << "Total Time                     : " << end_total << '\n';

;

    BoxLib::Finalize();
}
