
#include "myfunc.H"
#include "myfunc_F.H"

#include <AMReX_BCUtil.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MultiFabUtil.H>



void sweep(MultiFab& phi_old,
	   MultiFab& phi_new,
	   std::array<MultiFab, AMREX_SPACEDIM>& flux,
	   Vector<MultiFab>& phi_sdc,
	   Vector<MultiFab>& res_sdc,
	   Vector<Vector<MultiFab> >& f_sdc,
	   Real dt,
	   const Geometry& geom,
	   const BoxArray& grids, 
	   const DistributionMapping& dmap, 
	   const Vector<BCRec>& bc,
           SDCstuff sdcmats)
{
  /*  We use an MLABecLaplacian operator:

      (ascalar*acoef - bscalar div bcoef grad) phi = RHS
      
      for an implicit discretization of the heat equation
      
      (I - div dt grad) phi^{n+1} = phi^n
  */
  Real dt_m,nudt_m;
  Real qij;

  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  
  MultiFab::Copy(phi_sdc[0],phi_old, 0, 0, 1, 0);
  // Fill the ghost cells of each grid from the other grids
  // includes periodic domain boundaries
  phi_sdc[0].FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries
  FillDomainBoundary(phi_sdc[0], geom, bc);
  
  //  Compute the first function value
  int sdc_m=0;
  for (int nf = 0; nf < SDC_NPIECES; nf++)
    {
      for ( MFIter mfi(phi_sdc[sdc_m]); mfi.isValid(); ++mfi )
	{
	  const Box& bx = mfi.validbox();
	  
	  compute_f(BL_TO_FORTRAN_BOX(bx),
		    BL_TO_FORTRAN_BOX(domain_bx),
		    BL_TO_FORTRAN_ANYD(phi_sdc[sdc_m][mfi]),
		    BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		    BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		    BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
		    BL_TO_FORTRAN_ANYD(f_sdc[nf][sdc_m][mfi]),
		    dx,&nf);
	}
    }


  // Copy first function value to all nodes
  for (int nf = 0; nf < SDC_NPIECES; nf++)
    {
      
      for (int sdc_n = 1; sdc_n < SDC_NNODES; sdc_n++)
	{
	  MultiFab::Copy(f_sdc[nf][sdc_n],f_sdc[nf][0], 0, 0, 1, 0);
	}
    }
  // assorment of solver and parallization options and parameters
  // see AMReX_MLLinOp.H for the defaults, accessors, and mutators
  LPInfo info;
  
  // Implicit solve using MLABecLaplacian class
  MLABecLaplacian mlabec({geom}, {grids}, {dmap}, info);
  
  // order of stencil
  int linop_maxorder = 2;
  mlabec.setMaxOrder(linop_maxorder);
  
  // build array of boundary conditions needed by MLABecLaplacian
  // see Src/Boundary/AMReX_LO_BCTYPES.H for supported types
  std::array<LinOpBCType,AMREX_SPACEDIM> bc_lo;
  std::array<LinOpBCType,AMREX_SPACEDIM> bc_hi;
  
  for (int n = 0; n < phi_sdc[0].nComp(); ++n) 
    {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	{
	  // lo-side BCs
	  if (bc[n].lo(idim) == BCType::int_dir) {
	    bc_lo[idim] = LinOpBCType::Periodic;
	  }
	  else if (bc[n].lo(idim) == BCType::foextrap) {
	    bc_lo[idim] = LinOpBCType::Neumann;
	  }
	  else if (bc[n].lo(idim) == BCType::ext_dir) {
	    bc_lo[idim] = LinOpBCType::Dirichlet;
	  }
	    else {
	      amrex::Abort("Invalid bc_lo");
	    }
	  
	  // hi-side BCs
	  if (bc[n].hi(idim) == BCType::int_dir) {
	    bc_hi[idim] = LinOpBCType::Periodic;
	  }
	  else if (bc[n].hi(idim) == BCType::foextrap) {
	    bc_hi[idim] = LinOpBCType::Neumann;
	  }
	  else if (bc[n].hi(idim) == BCType::ext_dir) {
	    bc_hi[idim] = LinOpBCType::Dirichlet;
	  }
	  else {
	    amrex::Abort("Invalid bc_hi");
	  }
	}
    }
  
  // tell the solver what the domain boundary conditions are
  mlabec.setDomainBC(bc_lo, bc_hi);
  
  // set the boundary conditions
  mlabec.setLevelBC(0, &phi_sdc[0]);
  
  // scaling factors
  Real ascalar = 1.0;
  Real bscalar = 1.0;
  mlabec.setScalars(ascalar, bscalar);
  
  // Set up coefficient matrices
  MultiFab acoef(grids, dmap, 1, 0);
  
  // fill in the acoef MultiFab and load this into the solver
  acoef.setVal(1.0);
  mlabec.setACoeffs(0, acoef);
  
  // bcoef lives on faces so we make an array of face-centered MultiFabs
  // then we will in face_bcoef MultiFabs and load them into the solver.
  std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      const BoxArray& ba = amrex::convert(acoef.boxArray(),
					  IntVect::TheDimensionVector(idim));
      face_bcoef[idim].define(ba, acoef.DistributionMap(), 1, 0);
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
  
  for (int k=1; k <= sdcmats.nsweeps; ++k)
    {
      amrex::Print() << "sweep " << k << "\n";
      // Compute the quadrature term
      for (int sdc_m = 0; sdc_m < SDC_NNODES-1; sdc_m++)
	{
	  for ( MFIter mfi(res_sdc[sdc_m]); mfi.isValid(); ++mfi )
	    {
	      const Box& bx = mfi.validbox();
	      res_sdc[sdc_m].setVal(0.0);
	      for (int sdc_n = 0; sdc_n < SDC_NNODES; sdc_n++)
		{
		  for (int nf = 0; nf < SDC_NPIECES; nf++)
		    {
		      qij = dt*(sdcmats.qmats[0][sdc_m][sdc_n]-sdcmats.qmats[nf+1][sdc_m][sdc_n]);
		      res_sdc[sdc_m][mfi].saxpy(qij,f_sdc[nf][sdc_n][mfi],bx,bx,0,0,1);
		    }
		}
	    }
	}
      
      //  Substep over SDC nodes
      for (int sdc_m = 0; sdc_m < SDC_NNODES-1; sdc_m++)
	{
	  dt_m = dt*(sdcmats.qnodes[sdc_m+1]-sdcmats.qnodes[sdc_m]);
	  
	  // use phi_new as rhs
	  MultiFab::Copy(phi_new,phi_sdc[0], 0, 0, 1, 0);
	  for ( MFIter mfi(phi_new); mfi.isValid(); ++mfi )
	    {
	      const Box& bx = mfi.validbox();
	      phi_new[mfi].saxpy(1.0,res_sdc[sdc_m][mfi],bx,bx,0,0,1);
	      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
		{
		  for (int nf = 0; nf < SDC_NPIECES; nf++)
		    {
		      qij = dt*sdcmats.qmats[nf+1][sdc_m][sdc_n];
		      phi_new[mfi].saxpy(qij,f_sdc[nf][sdc_n][mfi],bx,bx,0,0,1);
		    }
		}
	    }
	  
	  // Fill the ghost cells of each grid from the other grids
	  // includes periodic domain boundaries
	  phi_new.FillBoundary(geom.periodicity());
	  
	  // Fill non-periodic physical boundaries
	  FillDomainBoundary(phi_new, geom, bc);
	  
	  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	    {
	      //		const BoxArray& ba = amrex::convert(acoef.boxArray(),
	      //					    IntVect::TheDimensionVector(idim));
	      //		face_bcoef[idim].define(ba, acoef.DistributionMap(), 1, 0);
	      face_bcoef[idim].setVal(0.1*dt_m);
	    }
	  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));
	  
	  // set the boundary conditions
	  mlabec.setLevelBC(0, &phi_new);
	  
	  // get the best initial guess
	  MultiFab::Copy(phi_sdc[sdc_m+1],phi_new, 0, 0, 1, 0);
	  for ( MFIter mfi(phi_sdc[sdc_m+1]); mfi.isValid(); ++mfi )
	    {
	      const Box& bx = mfi.validbox();
	      qij = dt*sdcmats.qmats[2][sdc_m][sdc_m+1];
	      phi_sdc[sdc_m+1][mfi].saxpy(qij,f_sdc[1][sdc_m+1][mfi],bx,bx,0,0,1);
	    }
	  phi_sdc[sdc_m+1].FillBoundary(geom.periodicity());
	  phi_new.FillBoundary(geom.periodicity());
	  // Solve linear system
	  mlmg.solve({&phi_sdc[sdc_m+1]}, {&phi_new}, tol_rel, tol_abs);

	  // Compute the fluxes at node sdc_m+1
	  phi_sdc[sdc_m+1].FillBoundary(geom.periodicity());
	  for (int nf = 0; nf < SDC_NPIECES; nf++)
	    {
	      for ( MFIter mfi(phi_sdc[sdc_m+1]); mfi.isValid(); ++mfi )
		{
		  const Box& bx = mfi.validbox();
		  
		  compute_f(BL_TO_FORTRAN_BOX(bx),
			    BL_TO_FORTRAN_BOX(domain_bx),
			    BL_TO_FORTRAN_ANYD(phi_sdc[sdc_m+1][mfi]),
			    BL_TO_FORTRAN_ANYD(flux[0][mfi]),
			    BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
			    BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
			    BL_TO_FORTRAN_ANYD(f_sdc[nf][sdc_m+1][mfi]),
			    dx,&nf);
		}
	    }
	  

	  
	} // end SDC substep loop
    }  // end sweeps loop
    
  // Return the last node in phi_new
  MultiFab::Copy(phi_new, phi_sdc[SDC_NNODES-1], 0, 0, 1, 0);

}

