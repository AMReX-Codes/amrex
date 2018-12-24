
#include "myfunc.H"
#include "myfunc_F.H"
#include <AMReX_BCUtil.H>
#include <AMReX_MultiFabUtil.H>

void SDC_advance(MultiFab& phi_old,
		 MultiFab& phi_new,
		 std::array<MultiFab, AMREX_SPACEDIM>& flux,
		 Real dt,
		 const Geometry& geom,
		 const Vector<BCRec>& bc,
		 MLMG&  mlmg,
		 MLABecLaplacian& mlabec,
		 SDCstruct &SDC, Real a, Real d, Real r)
{

  /*  This is a multi-implicit SDC example time step for an 
      advection-diffusion-reaction equation of the form

      phi_t = A(phi)+D(phi)+R(phi)
      
      The advection is treated explicilty and the diffusion and reaction implicitly
      and uncoupled

      the constants a,d, and r control the strength of each term
  */
  Real qij;

  const BoxArray &ba=phi_old.boxArray();
  const DistributionMapping &dm=phi_old.DistributionMap();

  // Copy old phi into first SDC node
  MultiFab::Copy(SDC.sol[0],phi_old, 0, 0, 1, 0);

  // Fill the ghost cells of each grid from the other grids
  // includes periodic domain boundaries
  SDC.sol[0].FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries
  FillDomainBoundary(SDC.sol[0], geom, bc);
  
  //  Compute the first function value
  int sdc_m=0;
  SDC_feval(flux,geom,bc,SDC,a,d,r,sdc_m,-1);

  // Copy first function value to all nodes
  for (int sdc_n = 1; sdc_n < SDC.Nnodes; sdc_n++)
    {
      MultiFab::Copy(SDC.f[0][sdc_n],SDC.f[0][0], 0, 0, 1, 0);
      MultiFab::Copy(SDC.f[1][sdc_n],SDC.f[1][0], 0, 0, 1, 0);
      if (SDC.Npieces==3)
	MultiFab::Copy(SDC.f[2][sdc_n],SDC.f[2][0], 0, 0, 1, 0);      
    }


  //  Now do the actual sweeps
  for (int k=1; k <= SDC.Nsweeps; ++k)
    {
      amrex::Print() << "sweep " << k << "\n";

      //  Compute RHS integrals
      SDC.SDC_rhs_integrals(dt);
      
      //  Substep over SDC nodes
      for (int sdc_m = 0; sdc_m < SDC.Nnodes-1; sdc_m++)
	{
	  // use phi_new as rhs and fill it with terms at this iteration
	  SDC.SDC_rhs_k_plus_one(phi_new,dt,sdc_m);
	  
	  // get the best initial guess for implicit solve
	  MultiFab::Copy(SDC.sol[sdc_m+1],phi_new, 0, 0, 1, 0);
	  for ( MFIter mfi(SDC.sol[sdc_m+1]); mfi.isValid(); ++mfi )
	    {
	      //	      const Box& bx = mfi.validbox();
	      qij = dt*SDC.Qimp[sdc_m][sdc_m+1];
	      SDC.sol[sdc_m+1][mfi].saxpy(qij,SDC.f[1][sdc_m+1][mfi]);
	    }

	  // Solve for the first implicit part
	  SDC_fcomp(phi_new, geom, bc, SDC, mlmg, mlabec,dt,d,r,sdc_m+1,1);

	  if (SDC.Npieces==3)
	    {
	      // Build rhs for 2nd solve
	      MultiFab::Copy(phi_new, SDC.sol[sdc_m+1],0, 0, 1, 0);
	      
	      // Add in the part for the 2nd implicit term to rhs
	      SDC.SDC_rhs_misdc(phi_new,dt,sdc_m);
	      
	      // Solve for the second implicit part
	      SDC_fcomp(phi_new, geom, bc, SDC, mlmg, mlabec,dt,d,r,sdc_m+1,2);
	    }
	  // Compute the function values at node sdc_m+1
	  SDC_feval(flux,geom,bc,SDC,a,d,r,sdc_m+1,-1);	  

	} // end SDC substep loop
    }  // end sweeps loop
    
  // Return the last node in phi_new
  MultiFab::Copy(phi_new, SDC.sol[SDC.Nnodes-1], 0, 0, 1, 0);

}

void SDC_feval(std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
	       Real a,Real d,Real r,
	       int sdc_m,int npiece)
{
  /*  Evaluate explicitly the rhs terms of the equation at the SDC node "sdc_m".
      The input parameter "npiece" describes which term to do.  
      If npiece = -1, do all the pieces */
  const BoxArray &ba=SDC.sol[0].boxArray();
  const DistributionMapping &dm=SDC.sol[0].DistributionMap();

  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  int nlo,nhi;
  if (npiece < 0)
    {
      nlo=0;
      nhi=SDC.Npieces;
    }
  else
    {
      nlo=npiece;
      nhi=npiece+1;
    }

  SDC.sol[sdc_m].FillBoundary(geom.periodicity());
  for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();
      for (int n = nlo; n < nhi; n++)    
	{
	  SDC_feval_F(BL_TO_FORTRAN_BOX(bx),
		      BL_TO_FORTRAN_BOX(domain_bx),
		      BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
		      BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		      BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		      BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
		      BL_TO_FORTRAN_ANYD(SDC.f[n][sdc_m][mfi]),
		      dx,&a,&d,&r,&n);
	}
      
    }
}
void SDC_fcomp(MultiFab& rhs,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
	       MLMG &mlmg,
	       MLABecLaplacian& mlabec,	      
	       Real dt,Real d,Real r,
	       int sdc_m,int npiece)
{
  /*  Solve implicitly for the implicit terms of the equation at the SDC node "sdc_m".
      The input parameter "npiece" describes which term to do.  */

  const BoxArray &ba=SDC.sol[0].boxArray();
  const DistributionMapping &dm=SDC.sol[0].DistributionMap();

  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  Real qij;

    // relative and absolute tolerances for linear solve
  const Real tol_rel = 1.e-10;
  const Real tol_abs = 0.0;

  if (npiece == 1)  
    {
      // Do diffusion solve
      
      // Fill the ghost cells of each grid from the other grids
      // includes periodic domain boundaries
      rhs.FillBoundary(geom.periodicity());
      
      // Fill non-periodic physical boundaries
      FillDomainBoundary(rhs, geom, bc);

      //  Set diffusion scalar in solve
      qij = d*dt*SDC.Qimp[sdc_m-1][sdc_m];	      
      Real ascalar = 1.0;
      mlabec.setScalars(ascalar, qij);

      // set the boundary conditions
      mlabec.setLevelBC(0, &rhs);
      mlabec.setLevelBC(0, &SDC.sol[sdc_m]);  
      //  Do the multigrid solve
      mlmg.solve({&SDC.sol[sdc_m]}, {&rhs}, tol_rel, tol_abs);
      
    }
  else
    {  // Do reaction solve  y - qij*y*(1-y)*(y-1/2) = rhs

      //  make a flag to change how the reaction is done
      int nflag=1;  // Lazy approximation

      qij = r*dt*SDC.Qimp[sdc_m-1][sdc_m];	            
      for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
	{
	  const Box& bx = mfi.validbox();
	  SDC_fcomp_reaction_F(BL_TO_FORTRAN_BOX(bx),
			       BL_TO_FORTRAN_BOX(domain_bx),
			       BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
			       BL_TO_FORTRAN_ANYD(rhs[mfi]),		      
			       BL_TO_FORTRAN_ANYD(SDC.f[2][sdc_m][mfi]),
			       &qij,&nflag); 
      }
      
    }

}


