
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

#include "myfunc.H"
#include "myfunc_F.H"

// additional header files needed by MUI
#include <mpi.h>
#include <lib_mpi_split.h>
#include <mui.h>

using namespace amrex;
using namespace mui;

int main (int argc, char* argv[])
{
    MPI_Comm comm = mui::mpi_split_by_app( argc, argv );
    amrex::Initialize(argc,argv,true,comm);
    
    main_main();
    
    amrex::Finalize();
    return 0;
}

void main_main ()
{
    static_assert(AMREX_SPACEDIM == 3, "3D only");

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = amrex::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size_3d, plot_int, verbosity;
    Vector<int> is_periodic(AMREX_SPACEDIM,1);  // periodic in all direction by default

    // inputs parameters
    {
        // ParmParse is way of reading inputs from the inputs file
        ParmParse pp;

        // We need to get n_cell from the inputs file - this is the number of cells on each side of 
        //   a square (or cubic) domain.
        pp.get("n_cell",n_cell);

        // The domain is broken into boxes of size max_grid_size
        pp.get("max_grid_size_3d",max_grid_size_3d);

        // Default plot_int to -1, allow us to set it to something else in the inputs file
        //  If plot_int < 0 then no plot files will be writtenq
        plot_int = -1;
        pp.query("plot_int",plot_int);

        // Default verbosity to 10, allow us to set it to something else in the inputs file
        verbosity = 0;
        pp.query("verbosity",verbosity);

        pp.queryarr("is_periodic", is_periodic);
    }

    // make BoxArray and Geometry
    BoxArray ba;
    Geometry geom;
    {
        IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
        Box domain(dom_lo, dom_hi);

        // Initialize the boxarray "ba" from the single box "bx"
        ba.define(domain);
        // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
        ba.maxSize(max_grid_size_3d);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }
  
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);
    
    if (verbosity > 0) {
      Print() << "3D:" << std::endl << ba << std::endl;
      Print() << "3D:" << std::endl << dm << std::endl;
    }

    // Note: must have no ghost cells for direct pointer access to work
    MultiFab phi(ba, dm, 1, 0);

    int myrank = ParallelDescriptor::MyProc();  // Can also use: MPI_Comm_rank(comm, &myrank);
    int root_proc = 0;

    // Initialize phi by calling a Fortran routine.
    // MFIter = MultiFab Iterator
    for ( MFIter mfi(phi); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();

        init_phi(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(phi[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi());
    }

    // time
    int n = 0;
    Real time = 0.0;

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
    {
        const std::string& pltfile_orig = amrex::Concatenate("plt3D_orig_",n,5);
        WriteSingleLevelPlotfile(pltfile_orig, phi, {"phi_original"}, geom, time, 0);
    }

    // define an interface.  Note on the MUI side, the string "FHD-KMC-coupling" must match
    // note you could use an uniface2d, however, uniface3d allows targeted receive of data on 3D side, 
    // because multiple 3D grids contain the same (x,y) coordinates, where (x,y,z) is unique
    mui::uniface3d uniface( "mpi://FHD-side/FHD-KMC-coupling" );
    mui::sampler_exact3d<double> r;
    mui::sampler_exact3d<int> s;
    mui::chrono_sampler_exact3d t;

    //////////////////////////////////////////////////////////////
    // Send 2D slice (in x & y, z index = 0) to 2D executable

    for ( MFIter mfi(phi); mfi.isValid(); ++mfi )
      {

	// Pull box from MFIter loop of MultiFab
	const Box& bx = mfi.validbox();
	
	// Determine box dimensions & low indices
	int nx = bx.size()[0];
	int ny = bx.size()[1];
	int nz = bx.size()[2];
	int low_x = bx.smallEnd(0);
	int low_y = bx.smallEnd(1);
	int low_z = bx.smallEnd(2);
	
	if (verbosity > 0)
	  Print() << "3D: nx, ny = " << nx << ", " << ny << std::endl;
	
	// Declare single index to iterate through each box
	int local_indx = 0;

	// Integer position vector used as unique identifier in MUI send & receive 
	Vector< int > pos;

	// Temporary variable to store data to send
	double x_send;

	// z-index of slice can be specified, in this case set to "lower most" face
        int k = 0;

	if (k >= low_z && k <= low_z+nz) {  // Only send if on correct face
	  
	  // annouce send span
	  geometry::box3d send_region( {(double)low_x, (double)low_y, (double)k}, 
				       {(double)(low_x + nx-1), (double)(low_y + ny-1), (double)k} );
	  printf( "send region for rank %d: %d %d - %d %d\n", myrank, low_x, low_y, 
		  low_x+nx, low_y+ny );
	  uniface.announce_send_span( 0, 1, send_region );

	  for(int j=0; j<ny; j++) {
	    for(int i=0; i<nx; i++) {
	      // Send 

	      // Determine global position in domain
	      pos = {(int)(low_x+i), (int)(low_y+j)};
	      point3d loc( (double)pos[0], (double)pos[1], (double)k);
	      
	      // Access MultiFab values directly (0th component), belonging to MFI grid
	      // Note: make sure no ghost cells
	      x_send = phi[mfi].dataPtr(0)[k*nx*ny + local_indx];

	      // Send data point, specifying unique global location as identifier
	      uniface.push( "channel1", loc, x_send);
	      
	      if (verbosity > 0) {
		printf("SENDER %d, step %d: channel1 at (%d,%d) is %f\n", 
		       ParallelDescriptor::MyProc(), n, pos[0], pos[1], x_send);
	      }

	      local_indx++;
	    }
	  }
	}
      }
    
    // Make sure to "commit" pushed data
    uniface.commit( n );

    //////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////
    // Receive modified data from 2D executable

    // Wait for other executable
    // "uniface.commit( n )" command in other code signals current process to continue
    uniface.barrier( n );

    for ( MFIter mfi(phi); mfi.isValid(); ++mfi )
      {
	
	// Pull box from MFIter loop of MultiFab
	const Box& bx = mfi.validbox();
	
	// Determine box dimensions & low indices
	int nx = bx.size()[0];
	int ny = bx.size()[1];
	int nz = bx.size()[2];
	int low_x = bx.smallEnd(0);
	int low_y = bx.smallEnd(1);
	int low_z = bx.smallEnd(2);
	
	if (verbosity > 0)
	  Print() << "3D: nx, ny = " << nx << ", " << ny << std::endl;
	
	// Declare single index to iterate through each box
	int local_indx = 0;

	// Integer position vector used as unique identifier in MUI send & receive 
	Vector< int > pos;

	// Temporary variable to store incoming data
	double x_recv;

	// z-index of slice can be specified, in this case set to "lower most" face
        int k = 0;

	if (k >= low_z && k <= low_z+nz) {  // Only recv if on correct face
	  
	  // announce 'span' for smart sending
	  geometry::box3d recv_region( {(double)low_x, (double)low_y, (double)k}, 
				       {(double)(low_x + nx-1), (double)(low_y + ny-1), (double)k} );
	  uniface.announce_recv_span( 0, 1, recv_region);
	  
	  for(int j=0; j<ny; j++) {
	    for(int i=0; i<nx; i++) {
	      // Receive 

	      // Determine global position in domain
	      pos = {(int)(low_x+i), (int)(low_y+j)};
	      point3d loc( (double)pos[0], (double)pos[1], (double)k);

	      // Receive data point, specifying unique global location as identifier
	      x_recv = uniface.fetch( "channel1", loc, n, r, t );
	      
	      if (verbosity > 0) {
		printf("RETURNED %d, step %d: channel2 at (%d,%d) is %f\n", 
		       myrank, n, pos[0], pos[1], x_recv);
	      }

	      // Access MultiFab values directly (0th component), belonging to MFI grid
	      // Note: make sure no ghost cells
	      phi[mfi].dataPtr(0)[k*nx*ny + local_indx] = x_recv;

	      local_indx++;
	    }
	  }
	}
      }
    
    // Clear data after it is received to free-up storage
    uniface.forget( n );

    //////////////////////////////////////////////////////////////

    // Write a plotfile of the current data (plot_int was defined in the inputs file)
    if (plot_int > 0 && n%plot_int == 0)
      {
	const std::string& pltfile = amrex::Concatenate("plt3D_mod_",n,5);
	WriteSingleLevelPlotfile(pltfile, phi, {"phi_modified"}, geom, time, n);
      }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = amrex::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
