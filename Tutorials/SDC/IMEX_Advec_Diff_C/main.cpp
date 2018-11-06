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
    // What time is it now?  We'll use this to compute total run time.
    Real strt_time = ParallelDescriptor::second();

    // AMREX_SPACEDIM: number of dimensions
    int n_cell, max_grid_size, nsteps, plot_int;
    Vector<int> bc_lo(AMREX_SPACEDIM,0);
    Vector<int> bc_hi(AMREX_SPACEDIM,0);

    // inputs parameters
    {
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

        // Default nsteps to 0, allow us to set it to something else in the inputs file
        nsteps = 10;
        pp.query("nsteps",nsteps);

        // read in BC; see Src/Base/AMReX_BC_TYPES.H for supported types
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);
    }

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

    // we allocate two phi multifabs; one will store the old state, the other the new.
    MultiFab phi_old(ba, dm, Ncomp, Nghost);
    MultiFab phi_new(ba, dm, Ncomp, Nghost);

    // Initialize phi_new by calling a Fortran routine.
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
            else if (bc_lo[idim] == FOEXTRAP) {
                bc[n].setLo(idim, BCType::foextrap); // first-order extrapolation 
            }
            else if (bc_lo[idim] == EXT_DIR) {
                bc[n].setLo(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_lo");
            }

            // hi-side BCs
            if (bc_hi[idim] == INT_DIR) {
                bc[n].setHi(idim, BCType::int_dir);  // periodic uses "internal Dirichlet"
            }
            else if (bc_hi[idim] == FOEXTRAP) {
                bc[n].setHi(idim, BCType::foextrap); // first-order extrapolation (homogeneous Neumann)
            }
            else if (bc_hi[idim] == EXT_DIR) {
                bc[n].setHi(idim, BCType::ext_dir);  // external Dirichlet
            }
            else {
                amrex::Abort("Invalid bc_hi");
            }

        }
    }
    // build the flux multifabs
    std::array<MultiFab, AMREX_SPACEDIM> flux;
    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
    {
        // flux(dir) has one component, zero ghost cells, and is nodal in direction dir
        BoxArray edge_ba = ba;
        edge_ba.surroundingNodes(dir);
        flux[dir].define(edge_ba, dm, 1, 0);
    }

    // Make an SDC structure
    int Nnodes=5;
    int Npieces=2;
    SDCstruct SDCmats(Nnodes,Npieces,phi_old);

    // Compute the time step
    // Implicit time step is imFactor*(explicit time step)
    const Real* dx = geom.CellSize();
    const int imFactor = pow(10, AMREX_SPACEDIM-1);
    Real Tfin = 0.5;
    Real dt = imFactor*0.9*dx[0]*dx[0] / (2.0*AMREX_SPACEDIM);
    dt = Tfin/nsteps;
    // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
    //    if (plot_int > 0)
    //    {
    //        int n = 0;
    //        const std::string& pltfile = amrex::Concatenate("plt",n,5);
    //        WriteSingleLevelPlotfile(pltfile, phi_new, {"phi"}, geom, time, 0);
    //    }

    // assorment of solver and parallization options and parameters
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
	  else if (bc[n].lo(idim) == BCType::foextrap) {
	    mgbc_lo[idim] = LinOpBCType::Neumann;
	  }
	  else if (bc[n].lo(idim) == BCType::ext_dir) {
	    mgbc_lo[idim] = LinOpBCType::Dirichlet;
	  }
	    else {
	      amrex::Abort("Invalid bc_lo");
	    }
	  
	  // hi-side BCs
	  if (bc[n].hi(idim) == BCType::int_dir) {
	    mgbc_hi[idim] = LinOpBCType::Periodic;
	  }
	  else if (bc[n].hi(idim) == BCType::foextrap) {
	    mgbc_hi[idim] = LinOpBCType::Neumann;
	  }
	  else if (bc[n].hi(idim) == BCType::ext_dir) {
	    mgbc_hi[idim] = LinOpBCType::Dirichlet;
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
  // then we will in face_bcoef MultiFabs and load them into the solver.
  std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      const BoxArray& bamg = amrex::convert(acoef.boxArray(),
					  IntVect::TheDimensionVector(idim));
      face_bcoef[idim].define(bamg, acoef.DistributionMap(), 1, 0);
      face_bcoef[idim].setVal(dt);
    }
  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));
  
  // build an MLMG solver
  MLMG mlmg(mlabec);
  
  // set solver parameters
  int max_iter = 100;
  mlmg.setMaxIter(max_iter);
  int max_fmg_iter = 0;
  mlmg.setMaxFmgIter(max_fmg_iter);
  int verbose = 2;
  mlmg.setVerbose(verbose);
  int cg_verbose = 0;
  mlmg.setCGVerbose(cg_verbose);
  
  // relative and absolute tolerances for linear solve
  const Real tol_rel = 1.e-10;
  const Real tol_abs = 0.0;

    MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);
    SDCmats.Nsweeps =8;
    for (int n = 1; n <= nsteps; ++n)
    {

        // Do and SDC step
      SDC_advance(phi_old, phi_new,flux, dt, geom, bc, mlmg,mlabec,SDCmats); 
        time = time + dt;
	MultiFab::Copy(phi_old, phi_new, 0, 0, 1, 0);    

	// Turn the solution into the error
	int plot_err = 1;
	if (plot_err == 1) 
	  for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
	    {
	      const Box& bx = mfi.validbox();
	      err_phi(BL_TO_FORTRAN_BOX(bx),
		      BL_TO_FORTRAN_ANYD(phi_new[mfi]),
		      geom.CellSize(), geom.ProbLo(), geom.ProbHi(),&time);
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
