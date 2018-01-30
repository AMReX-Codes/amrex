
#include "myfunc.H"
#include "myfunc_F.H"

void sweep (MultiFab& phi_old,
	    MultiFab& phi_new,
	    std::array<MultiFab, AMREX_SPACEDIM>& flux,
	    std::array<MultiFab, SDC_NNODES>& phi_sdc,
	    std::array<MultiFab, SDC_NNODES>& f_sdc,
	    Real dt,
	    const Geometry& geom,
	    int nnodes,
	    Real qnodes [],
	    Real qmat [][SDC_NNODES],
  	    Real qmatE [][SDC_NNODES],
	    Real qmatI [][SDC_NNODES])

{

  int Ncomp = phi_old.nComp();
  int ng_p = phi_old.nGrow();
  int ng_f = flux[0].nGrow();
  Real dt_m = dt/(SDC_NNODES-1);
  const Real* dx = geom.CellSize();


  // Note that this simple example is not optimized.
  // The following two MFIter loops could be merged
  // and we do not have to use flux MultiFab.
  // 
  //  pf_quadrature(qtype_in, SDC_NNODES,SDC_NNODES, );
  const Box& domain_bx = geom.Domain();

  // Compute the residual
  for (int sdc_m = 0; sdc_m < SDC_NNODES-2; sdc_m++)
    {
      //      compute_resid()
    }

  //  Substep
  for (int sdc_m = 0; sdc_m < SDC_NNODES-1; sdc_m++)
    {
      // Fill the ghost cells of each grid from the other grids
      // includes periodic domain boundaries
      phi_sdc[sdc_m].FillBoundary(geom.periodicity());

      
      // Loop over SDC nodes
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
		       BL_TO_FORTRAN_ANYD(f_sdc[sdc_m][mfi]),
		       dx);
	}
      
      // Advance the solution one grid at a time
      for ( MFIter mfi(phi_sdc[sdc_m]); mfi.isValid(); ++mfi )
	{
	  const Box& bx = mfi.validbox();
	  
	  update_phi(BL_TO_FORTRAN_BOX(bx),
		     BL_TO_FORTRAN_ANYD(phi_sdc[sdc_m][mfi]),
		     BL_TO_FORTRAN_ANYD(phi_sdc[sdc_m+1][mfi]),
		     BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		     BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		     BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif
		     dx, dt_m);
	  
	  
	}
    }  // end SDC substep loop

  //  Compute the final function value
  phi_sdc[SDC_NNODES-1].FillBoundary(geom.periodicity());
  for ( MFIter mfi(phi_sdc[SDC_NNODES-1]); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();
      compute_flux(BL_TO_FORTRAN_BOX(bx),
		   BL_TO_FORTRAN_BOX(domain_bx),
		   BL_TO_FORTRAN_ANYD(phi_sdc[SDC_NNODES-1][mfi]),
		   BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		   BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		   BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
		   BL_TO_FORTRAN_ANYD(f_sdc[SDC_NNODES-1][mfi]),
		   dx);
    }
      
}
