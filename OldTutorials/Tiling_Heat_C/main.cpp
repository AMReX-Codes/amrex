
#include <AMReX_BLFort.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>

#include <writePlotFile.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

extern "C"
{
  void advance_phi(const int* lo, const int* hi,
		   const BL_FORT_FAB_ARG(phiold),
		   const BL_FORT_FAB_ARG(phinew),
		   const int& ncomp, const Real* dx, const Real& dt);

  void advance_phi2(const int* lo, const int* hi,
		    const BL_FORT_FAB_ARG(phiold),
		    const BL_FORT_FAB_ARG(phinew),
		    const int& ncomp, const Real* dx, const Real& dt);

  void init_phi(const int* lo, const int* hi,
		BL_FORT_FAB_ARG(phi),
		const int& ncomp, const Real* dx, const Real* prob_lo, const Real* prob_hi);
}

static Real kernel_time = 0;
static Real FB_time = 0;
static int do_tiling = 1;

void advance (MultiFab* old_phi, MultiFab* new_phi, Real* dx, Real dt, Geometry geom)
{
    int Ncomp = old_phi->nComp();

    // Fill the ghost cells of each grid from the other grids
    Real t1 = ParallelDescriptor::second();
    old_phi->FillBoundary(geom.periodicity());
    FB_time += ParallelDescriptor::second() - t1;

    Real t0 = ParallelDescriptor::second();
    
    if (do_tiling) {
#ifdef _OPENMP
#pragma omp parallel
#endif
	for ( MFIter mfi(*old_phi,true); mfi.isValid(); ++mfi )
	{
	    const Box& bx = mfi.tilebox();
	
	    advance_phi(bx.loVect(), bx.hiVect(),
			BL_TO_FORTRAN((*old_phi)[mfi]),
			BL_TO_FORTRAN((*new_phi)[mfi]),
			Ncomp,dx, dt);
	}
    } else {
	for ( MFIter mfi(*old_phi); mfi.isValid(); ++mfi )
	{
	    const Box& bx = mfi.validbox();
	
	    advance_phi2(bx.loVect(), bx.hiVect(),
			 BL_TO_FORTRAN((*old_phi)[mfi]),
			 BL_TO_FORTRAN((*new_phi)[mfi]),
			 Ncomp,dx, dt);
	}
    }
    
    kernel_time += ParallelDescriptor::second() - t0;
}



Real compute_dt (Real dx)
{
    return 0.9*dx*dx / (2.0*BL_SPACEDIM);
}

int
main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    std::cout << std::setprecision(15);

    // ParmParse is way of reading inputs from the inputs file
    ParmParse pp;
    
    int verbose = 0;
    pp.query("verbose", verbose);
    
    // We need to get n_cell from the inputs file - this is the number of cells on each side of 
    //   a square (or cubic) domain.
    int n_cell;
    pp.get("n_cell",n_cell);

    int max_grid_size;
    pp.get("max_grid_size",max_grid_size);

    // Default plot_int to 1, allow us to set it to something else in the inputs file
    //  If plot_int < 0 then no plot files will be written
    int plot_int = 1;
    pp.query("plot_int",plot_int);

    // Default nsteps to 0, allow us to set it to something else in the inputs file
    int nsteps   = 0;
    pp.query("nsteps",nsteps);

    pp.query("do_tiling", do_tiling);

    // Define a single box covering the domain
    IntVect dom_lo(0,0,0);
    IntVect dom_hi(n_cell-1,n_cell-1,n_cell-1);
    Box domain(dom_lo,dom_hi);

    // Initialize the boxarray "bs" from the single box "bx"
    BoxArray bs(domain);

    // Break up boxarray "bs" into chunks no larger than "max_grid_size" along a direction
    bs.maxSize(max_grid_size);

    // This defines the physical size of the box.  Right now the box is [-1,1] in each direction.
    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
	real_box.setLo(n,-1.0);
	real_box.setHi(n, 1.0);
    }

    // This says we are using Cartesian coordinates
    int coord = 0;
    
    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1; 
    
    // This defines a Geometry object which is useful for writing the plotfiles  
    Geometry geom(domain,&real_box,coord,is_per);

    // This defines the mesh spacing
    Real dx[BL_SPACEDIM];
    for ( int n=0; n<BL_SPACEDIM; n++ )
	dx[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/domain.length(n);

    // Nghost = number of ghost cells for each array 
    int Nghost = 1;

    // Ncomp = number of components for each array
    int Ncomp  = 1;
    pp.query("ncomp", Ncomp);

    DistributionMapping dm{bs};

    // Allocate space for the old_phi and new_phi -- we define old_phi and new_phi as
    Vector < std::unique_ptr<MultiFab> > phis(2);
    phis[0].reset(new MultiFab(bs, dm, Ncomp, Nghost));
    phis[1].reset(new MultiFab(bs, dm, Ncomp, Nghost));
    MultiFab* old_phi = phis[0].get();
    MultiFab* new_phi = phis[1].get();

    // Initialize both to zero (just because)
    old_phi->setVal(0.0);
    new_phi->setVal(0.0);

    // Initialize phi by calling a Fortran routine.
    // MFIter = MultiFab Iterator
#ifdef _OPENMP
#pragma omp parallel
#endif
    for ( MFIter mfi(*new_phi,true); mfi.isValid(); ++mfi )
    {
	const Box& bx = mfi.tilebox();

	init_phi(bx.loVect(),bx.hiVect(), 
		 BL_TO_FORTRAN((*new_phi)[mfi]),Ncomp,
		 dx,geom.ProbLo(),geom.ProbHi());
    }

    // Call the compute_dt routine to return a time step which we will pass to advance
    Real dt = compute_dt(dx[0]);

    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0) {
	int n = 0;
	const std::string& pltfile = amrex::Concatenate("plt",n,5);
	writePlotFile(pltfile, *new_phi, geom);
    }

    Real adv_start_time = ParallelDescriptor::second();

    for (int n = 1; n <= nsteps; n++)
    {
	// Swap the pointers so we don't have to allocate and de-allocate data
	std::swap(old_phi, new_phi);

	// new_phi = old_phi + dt * (something)
	advance(old_phi, new_phi, dx, dt, geom); 

	// Tell the I/O Processor to write out which step we're doing
	if (verbose && ParallelDescriptor::IOProcessor())
	    std::cout << "Advanced step " << n << std::endl;

	// Write a plotfile of the current data (plot_int was defined in the inputs file)
	if (plot_int > 0 && n%plot_int == 0) {
	    const std::string& pltfile = amrex::Concatenate("plt",n,5);
	    writePlotFile(pltfile, *new_phi, geom);
	}
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real advance_time = ParallelDescriptor::second() - adv_start_time;
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);
    ParallelDescriptor::ReduceRealMax(advance_time,IOProc);
    ParallelDescriptor::ReduceRealMax(kernel_time,IOProc);
    ParallelDescriptor::ReduceRealMax(FB_time,IOProc);
    
    // Tell the I/O Processor to write out the "run time"
    if (ParallelDescriptor::IOProcessor()) {
	std::cout << "Kernel    time = " << kernel_time << std::endl;
	std::cout << "FB        time = " << FB_time << std::endl;
	std::cout << "Advance   time = " << advance_time << std::endl;
	std::cout << "Total run time = " << stop_time << std::endl;
    }
    
    // Say goodbye to MPI, etc...
    amrex::Finalize();
}
