#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "myfunc.H"
#include "myfunc_F.H"  // includes advance.cpp
#include "AMReX_SDCstruct.H"
int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    
    main_main();
    
    amrex::Finalize();
    return 0;
}

void main_main ()
{
    Real a;  // advection coef.
    Real d;  // diffusion coef.
    Real r;  // reaction coef. 

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, Nsteps, plot_int;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    // ParmParse is way of reading inputs from the inputs file
    ParmParse pp;
    
    // We need to get n_cell from the inputs file - this is the number of cells on each side of 
    //   a square (or cubic) domain.
    pp.get("n_cell",n_cell);
    
    // The domain is broken into boxes of size max_grid_size
    pp.get("max_grid_size",max_grid_size);
    
    // Default plot_int to -1, allow us to set it to something else in the inputs file
    //  If plot_int < 0 then no plot files will be writtenq
    plot_int = -1;
    pp.query("plot_int",plot_int);

    // Set  plot_err = 1 to  output the error to plot files instead of the solution
    int plot_err = 1;
    
    // Read in number of steps and final time
    pp.query("Nsteps",Nsteps);
    Real Tfin=0.0;
    pp.query("Tfin",Tfin);    
    Real dt = Tfin/Nsteps;  // Set the time step

    // read in BC; see Src/Base/AMReX_BC_TYPES.H for supported types
    pp.queryarr("bc_lo", bc_lo);
    pp.queryarr("bc_hi", bc_hi);
    
    //  Read in the coefficients for A-D-R
    pp.query("a",a);
    pp.query("d",d);
    pp.query("r",r);
    
    // determine whether boundary conditions are periodic
    Vector<int> is_periodic(AMREX_SPACEDIM,0);
    for (int idim=0; idim < AMREX_SPACEDIM; ++idim) {
        if (bc_lo[idim] == INT_DIR && bc_hi[idim] == INT_DIR) {
            is_periodic[idim] = 1;
        }
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
        ba.maxSize(max_grid_size);

       // This defines the physical box, [-1,1] in each direction.
        RealBox real_box({AMREX_D_DECL(-1.0,-1.0,-1.0)},
                         {AMREX_D_DECL( 1.0, 1.0, 1.0)});

        // This defines a Geometry object
        geom.define(domain,&real_box,CoordSys::cartesian,is_periodic.data());
    }

    // Nghost = number of ghost cells for each array 
    int Nghost = 2;
    
    // Ncomp = number of components for each array
    int Ncomp  = 1;

    // time = starting time in the simulation
    Real time = 0.0;
  
    // How Boxes are distrubuted among MPI processes
    DistributionMapping dm(ba);

    // We allocate two phi multifabs; one will store the old state, the other the new.
    MultiFab phi_old(ba, dm, Ncomp, Nghost);
    MultiFab phi_new(ba, dm, Ncomp, Nghost);

    // Initialize phi_new by calling a Fortran routine (init_phi_2d.f90).
    // MFIter = MultiFab Iterator
    for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
    {
        const Box& bx = mfi.validbox();
        init_phi(BL_TO_FORTRAN_BOX(bx),
                 BL_TO_FORTRAN_ANYD(phi_new[mfi]),
                 geom.CellSize(), geom.ProbLo(), geom.ProbHi());
    }
    
    // Set up BCRec; see Src/Base/AMReX_BC_TYPES.H for supported types
    Vector<BCRec> bc(phi_old.nComp());
    for (int n = 0; n < phi_old.nComp(); ++n)
      {
	for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	  {
	    // lo-side BCs
	    if (bc_lo[idim] == INT_DIR) {
	      bc[n].setLo(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
	    }
	    else {
	      amrex::Abort("Invalid bc_lo");
	    }
	    
	    // hi-side BCs
	    if (bc_hi[idim] == INT_DIR) {
	      bc[n].setHi(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
	    }
	    else {
	      amrex::Abort("Invalid bc_hi");
	    }
	  } 
      }

    // Build the flux multifabs
    std::array<MultiFab, AMREX_SPACEDIM> flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        flux[dir].define(edge_ba, dm, 1, 0);
    }

    // Make an SDC structure
    int Nnodes=5;  // Default to 8th order
    int Npieces=3; // Default is full MISDC
    int Nsweeps=2*Nnodes-2;  //  This will give highest formal accuracy for Lobatto nodes
    pp.get("Nnodes",Nnodes);
    pp.get("Npieces",Npieces);
    //    pp.get("Nsweeps",Nsweeps);  //  Uncomment to adjust Nsweeps          

    //  Build the structure
    SDCstruct SDCmats(Nnodes,Npieces,phi_old);
    SDCmats.Nsweeps =Nsweeps;  // Number of SDC sweeps per time step
    
    const Real* dx = geom.CellSize();
    
    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    if (plot_int > 0)
      {
	if (plot_err == 1)  // Turn the solution into the error
	  {
	    MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);
	    for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
	      {
		const Box& bx = mfi.validbox();
		err_phi(BL_TO_FORTRAN_BOX(bx),
			BL_TO_FORTRAN_ANYD(phi_new[mfi]),
			geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&a,&d,&r,&time);
	      }
	  }
	int n = 0;
	const std::string& pltfile = amrex::Concatenate("plt",n,5);
	WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);
	if (plot_err == 1)  // Put the solution back
	  MultiFab::Copy(phi_new, phi_old, 0, 0, 1, 0);	
      }

  // Set an assorment of solver and parallization options and parameters
  // see AMReX_MLLinOp.H for the defaults, accessors, and mutators
  LPInfo info;
  
  // Implicit solve using MLABecLaplacian class
  MLABecLaplacian mlabec({geom}, {ba}, {dm}, info);
  
  // order of stencil
  int linop_maxorder = 2;
  mlabec.setMaxOrder(linop_maxorder);
  
  // build array of boundary conditions needed by MLABecLaplacian
  // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
  std::array<LinOpBCType,AMREX_SPACEDIM> mgbc_lo;
  std::array<LinOpBCType,AMREX_SPACEDIM> mgbc_hi;
  
  for (int n = 0; n < phi_old.nComp(); ++n) 
    {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	{
	  // lo-side BCs
	  if (bc[n].lo(idim) == BCType::int_dir) {
	    mgbc_lo[idim] = LinOpBCType::Periodic;
	  }
	  else {
	    amrex::Abort("Invalid bc_lo");
	  }
	  
	  // hi-side BCs
	  if (bc[n].hi(idim) == BCType::int_dir) {
	    mgbc_hi[idim] = LinOpBCType::Periodic;
	  }
	  else {
	    amrex::Abort("Invalid bc_hi");
	  }
	}
    }
  
  // tell the solver what the domain boundary conditions are
  mlabec.setDomainBC(mgbc_lo, mgbc_hi);
  
  // scaling factors
  Real ascalar = 1.0;
  Real bscalar = 1.0;
  mlabec.setScalars(ascalar, bscalar);
  
  // Set up coefficient matrices
  MultiFab acoef(ba, dm, 1, 0);
  
  // fill in the acoef MultiFab and load this into the solver
  acoef.setVal(1.0);
  mlabec.setACoeffs(0, acoef);
  
  // bcoef lives on faces so we make an array of face-centered MultiFabs
  //   then we will in face_bcoef MultiFabs and load them into the solver.
    std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
      {
        const BoxArray& bamg = amrex::convert(acoef.boxArray(),
  					  IntVect::TheDimensionVector(idim));
        face_bcoef[idim].define(bamg, acoef.DistributionMap(), 1, 0);
	face_bcoef[idim].setVal(1.0);	      	
      }
    mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));
  
  // build an MLMG solver
  MLMG mlmg(mlabec);
  
  // set solver parameters
  int max_iter = 100;
  mlmg.setMaxIter(max_iter);
  int max_fmg_iter = 0;
  mlmg.setMaxFmgIter(max_fmg_iter);
  int verbose = 1;
  mlmg.setVerbose(verbose);
  int cg_verbose = 0;
  mlmg.setCGVerbose(cg_verbose);
  

  //  Do the time stepp[ing
  MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);
  for (int n = 1; n <= Nsteps; ++n)
    {
      
      // Do an SDC step
      SDC_advance(phi_old, phi_new,flux, dt, geom, bc, mlmg,mlabec,SDCmats,a,d,r); 

      time = time + dt;
      MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);    
      
      
      if (plot_err == 1)  // Turn the solution into the error
	for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
	  {
	    const Box& bx = mfi.validbox();
	    err_phi(BL_TO_FORTRAN_BOX(bx),
		    BL_TO_FORTRAN_ANYD(phi_new[mfi]),
		    geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&a,&d,&r,&time);
	  }

        // Tell the I/O Processor to write out which step we're doing
        amrex::Print() << "Advanced step " << n << "\n";

        // Write a plotfile of the current data (plot_int was defined in the inputs file)
        if (plot_int > 0 && n%plot_int == 0)
        {
            const std::string& pltfile = amrex::Concatenate("plt",n,5);
            WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, n);
        }
    }

    // Call the timer again and compute the maximum difference between the start time and stop time
    //   over all processors
    Real stop_time = ParallelDescriptor::second() - strt_time;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

    // Tell the I/O Processor to write out the "run time"
    amrex::Print() << "Run time = " << stop_time << std::endl;
}
