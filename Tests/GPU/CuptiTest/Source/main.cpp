#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#ifdef AMREX_USE_CUPTI
#include <AMReX_ActivityTraceAsync.H>
#endif // AMREX_USE_CUPTI

#include "myfunc.H"

using namespace amrex;

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // Initialize the CUPTI trace
#ifdef AMREX_USE_CUPTI
    initTrace();
#endif // AMREX_USE_CUPTI
    
    main_main();
    
    amrex::Finalize();
    return 0;
}

void main_main ()
{

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    // Set parameters
    // Number of cells on each side of a square (or cubic) domain
    n_cell = 12;

    // The domain is broken into boxes of size max_grid_size
    max_grid_size = 7;

    // Run for 2 steps
    nsteps = 2;

    // Make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
	
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size);

	// This defines the physical box, [-1,1] in each direction
        RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array 
    int Nghost = 1;
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;
  
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // Allocate a MultiFab
    MultiFab mf(ba, dm, Ncomp, Nghost);

    GpuArray<Real,AMREX_SPACEDIM> dx = geom.CellSizeArray();

    init_mf(mf, geom);

    for (int n = 1; n <= nsteps; ++n)
    {
        doRandomSleep(mf); 
        
        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";
    }
}
