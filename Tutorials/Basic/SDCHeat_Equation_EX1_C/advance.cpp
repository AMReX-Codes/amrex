
#include "myfunc.H"
#include "myfunc_F.H"

void sweep (MultiFab& phi_old,
	    MultiFab& phi_new,
	    std::array<MultiFab, AMREX_SPACEDIM>& flux,
   	    Vector<MultiFab>& phi_sdc,
   	    Vector<MultiFab>& res_sdc,
	    Vector<Vector<MultiFab> >& f_sdc,
	    Vector<Vector<MultiFab> >& f_sdc_old,
	    Real dt,
	    const Geometry& geom,
	    Real qnodes [],
	    Real qmat [][SDC_NNODES],
  	    Real qmatE [][SDC_NNODES],
	    Real qmatI [][SDC_NNODES],
	    int nsweeps)


{

  int Ncomp = phi_old.nComp();
  int ng_p = phi_old.nGrow();
  int ng_f = flux[0].nGrow();
  Real dt_m;
  Real qij;
  const Real* dx = geom.CellSize();


  // Note that this simple example is not optimized.
  // The following two MFIter loops could be merged
  // and we do not have to use flux MultiFab.
  // 
  //  pf_quadrature(qtype_in, SDC_NNODES,SDC_NNODES, );
  const Box& domain_bx = geom.Domain();

  MultiFab::Copy(phi_sdc[0],phi_old, 0, 0, 1, 0);

  //  Compute the first function value
  int sdc_m=0;
  phi_sdc[sdc_m].FillBoundary(geom.periodicity());
  for ( MFIter mfi(phi_sdc[sdc_m]); mfi.isValid(); ++mfi )
	{
	  const Box& bx = mfi.validbox();
	  
	  compute_flux(BL_TO_FORTRAN_BOX(bx),
		       BL_TO_FORTRAN_BOX(domain_bx),
		       BL_TO_FORTRAN_ANYD(phi_sdc[sdc_m][mfi]),
		       BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		       BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		       BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
		       BL_TO_FORTRAN_ANYD(f_sdc[1][sdc_m][mfi]),
		       dx);
	} 
  // Copy first function value to all nodes
  for (int sdc_n = 1; sdc_n < SDC_NNODES; sdc_n++)
    MultiFab::Copy(f_sdc[1][sdc_n],f_sdc[1][0], 0, 0, 1, 0);
  

  for (int k; k < nsweeps; ++k)
    {
      // Compute the quadrature term
      for (int sdc_m = 0; sdc_m < SDC_NNODES-1; sdc_m++)
	{
	  for ( MFIter mfi(res_sdc[sdc_m]); mfi.isValid(); ++mfi )
	    {
	      res_sdc[sdc_m].setVal(0.0);
	      for (int sdc_n = 0; sdc_n < SDC_NNODES; sdc_n++)
		{
		  qij = dt*(qmat[sdc_m][sdc_n]-qmatE[sdc_m][sdc_n]);
		  res_sdc[sdc_m][mfi].saxpy(qij,f_sdc[1][sdc_n][mfi]);
		}
	    }
	}

      //  Substep over SDC nodes
      for (int sdc_m = 0; sdc_m < SDC_NNODES-1; sdc_m++)
	{
	  dt_m = dt*(qnodes[sdc_m+1]-qnodes[sdc_m]);
	  
	  // Fill the ghost cells of each grid from the other grids
	  // includes periodic domain boundaries
	  
	  
	  // Advance the solution one grid at a time
	  //      MultiFab::Copy(phi_sdc[sdc_m+1],phi_sdc[sdc_m], 0, 0, 1, 0);
	  //      for ( MFIter mfi(phi_sdc[sdc_m+1]); mfi.isValid(); ++mfi )
	  //      	{
	  //      	  const Box& bx = mfi.validbox();
	  //      	  phi_sdc[sdc_m+1][mfi].saxpy(dt_m,f_sdc[1][sdc_m][mfi]);
	  //      	}

	  // Add the linear combination of function values
	  MultiFab::Copy(phi_sdc[sdc_m+1],phi_sdc[0], 0, 0, 1, 0);
	  for ( MFIter mfi(phi_sdc[sdc_m+1]); mfi.isValid(); ++mfi )
	    {
	      const Box& bx = mfi.validbox();
	      phi_sdc[sdc_m+1][mfi].saxpy(1.0,res_sdc[sdc_m][mfi]);
	      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
		{
		  qij = dt*qmatE[sdc_m][sdc_n];
		  //      	      amrex::Print() << "sdc_n " << sdc_n << "\n";
		  //amrex::Print() << "qij " << qij << "\n";
		  //amrex::Print() << "dt_m " << dt_m << "\n";	      
		  
		  phi_sdc[sdc_m+1][mfi].saxpy(qij,f_sdc[1][sdc_n][mfi]);
		}
	    }
	  
	  // Compute the fluxes at node sdc_m+1
	  phi_sdc[sdc_m+1].FillBoundary(geom.periodicity());
	  for ( MFIter mfi(phi_sdc[sdc_m+1]); mfi.isValid(); ++mfi )
	    {
	      const Box& bx = mfi.validbox();
	      
	      compute_flux(BL_TO_FORTRAN_BOX(bx),
			   BL_TO_FORTRAN_BOX(domain_bx),
			   BL_TO_FORTRAN_ANYD(phi_sdc[sdc_m+1][mfi]),
			   BL_TO_FORTRAN_ANYD(flux[0][mfi]),
			   BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
			   BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
			   BL_TO_FORTRAN_ANYD(f_sdc[1][sdc_m+1][mfi]),
			   dx);
	    }
	  
	}  // end SDC substep loop
    }  // end sweeps loop
  // Return the last node in phi_new
  MultiFab::Copy(phi_new, phi_sdc[SDC_NNODES-1], 0, 0, 1, 0);
  //  MultiFab::Copy(phi_new, phi_sdc[0], 0, 0, 1, 0);  
}

