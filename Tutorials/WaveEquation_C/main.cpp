#include <fstream>
#include <iomanip>

#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>
#include <writePlotFile.H>

#include <AMReX_ArrayLim.H>

using namespace amrex;

extern "C"
{
    void init_data(Real* data, const int* lo, const int* hi,
		   const int* Ncomp,  const int* ng, 
		   const Real* dx, const Real* prob_lo, const Real* prob_hi);

    void advance(Real* old_data , Real* new_data, const int* lo, const int* hi, 
		 const int* Ncomp, const int* ng, const Real* dx, const Real* dt);
}

static
void advance(MultiFab* old_data, MultiFab* new_data, Real* dx, Real dt)
{
  // Fill the ghost cells of each grid from the other grids
  old_data->FillBoundary();

  int Ncomp = old_data->nComp();

  int ng = old_data->nGrow();

  // Advance the solution one grid at a time
  for ( MFIter mfi(*new_data); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.validbox();

    advance((*old_data)[mfi].dataPtr(),
	    (*new_data)[mfi].dataPtr(),
	    bx.loVect(),bx.hiVect(),&Ncomp, &ng,
	    &(dx[0]),&dt);
  }
}

static
Real compute_dt(Real dx)
{
   return 0.1 * dx;
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

  // We need to get n_cell from the inputs file - this is the number of cells on each side of 
  //   a square (or cubic) domain.
  int n_cell;
  pp.get("n_cell",n_cell);

  // Default Nsteps to 0, allow us to set it to something else in the inputs file
  int max_grid_size;
  pp.get("max_grid_size",max_grid_size);

  // Default plot_int to 1, allow us to set it to something else in the inputs file
  //  If plot_int < 0 then no plot files will be written
  int plot_int = 1;
  pp.query("plot_int",plot_int);

  // Default Nsteps to 0, allow us to set it to something else in the inputs file
  int Nsteps   = 0;
  pp.query("Nsteps",Nsteps);

  // Default verbose to 0, allow us to set it to 1 i the inputs file
  int verbose   = 0;
  pp.query("verbose",verbose);

  // Define a single box covering the domain
#if (BL_SPACEDIM == 2)
  IntVect dom_lo(0,0);
  IntVect dom_hi(n_cell-1,n_cell-1);
#else
  IntVect dom_lo(0,0,0);
  IntVect dom_hi(n_cell-1,n_cell-1,n_cell-1);
#endif
  Box domain(dom_lo,dom_hi);

  // Initialize the boxarray "bs" from the single box "bx"
  BoxArray bs(domain);

  // Break up boxarray "bs" into chunks no larger than "max_grid_size" along a direction
  bs.maxSize(max_grid_size);

  // This defines the physical size of the box.  Right now the box is [-1,1] in each direction.
  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++)
  {
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
  int Nghost = 6;

  // Ncomp = number of components for each array
  int Ncomp  = 2;

  // Make sure we can fill the ghost cells from the adjacent grid
  if (Nghost > max_grid_size)
    std::cout <<  "NGHOST < MAX_GRID_SIZE --  grids are too small! " << std::endl;

  DistributionMapping dm{bs};

  // Allocate space for the old_data and new_data -- we define old_data and new_data as
  //   pointers to the MultiFabs
  MultiFab* old_data = new MultiFab(bs, dm, Ncomp, Nghost);
  MultiFab* new_data = new MultiFab(bs, dm, Ncomp, Nghost);

  // Initialize both to zero (just because)
  old_data->setVal(0.0);
  new_data->setVal(0.0);

  // Initialize the old_data by calling a Fortran routine.
  // MFIter = MultiFab Iterator
  for ( MFIter mfi(*new_data); mfi.isValid(); ++mfi )
  {
    const Box& bx = mfi.validbox();

    init_data((*new_data)[mfi].dataPtr(),
	      bx.loVect(),bx.hiVect(), &Ncomp, &Nghost,
	      dx,geom.ProbLo(),geom.ProbHi());
  }

  // Call the compute_dt routine to return a time step which we will pass to advance
  Real dt = compute_dt(dx[0]);

  // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
  if (plot_int > 0)
  {
     int n = 0;
     const std::string& pltfile = amrex::Concatenate("plt",n,5);
     writePlotFile(pltfile, *new_data, geom);
  }

  for (int n = 1; n <= Nsteps; n++)
  {
     // Swap the pointers so we don't have to allocate and de-allocate data
     std::swap(old_data, new_data);

     // new_data = old_data + dt * (something)
     advance(old_data, new_data, dx, dt); 

     // Tell the I/O Processor to write out which step we're doing (verbose is defined in the inputs file)
     if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Advanced step " << n << std::endl;

     // Write a plotfile of the current data (plot_int was defined in the inputs file)
     if (plot_int > 0 && n%plot_int == 0)
     {
        const std::string& pltfile = amrex::Concatenate("plt",n,5);
        writePlotFile(pltfile, *new_data, geom);
     }
  }

  // Call the timer again and compute the maximum difference between the start time and stop time
  //   over all processors
  Real stop_time = ParallelDescriptor::second() - strt_time;
  const int IOProc = ParallelDescriptor::IOProcessorNumber();
  ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

  // Tell the I/O Processor to write out the "run time"
  if (ParallelDescriptor::IOProcessor())
     std::cout << "Run time = " << stop_time << std::endl;
  
  // Say goodbye to MPI, etc...
  amrex::Finalize();

}
