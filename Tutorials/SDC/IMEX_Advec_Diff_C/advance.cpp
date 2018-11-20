
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
		 SDCstruct &SDC)
{

  /*  We use an MLABecLaplacian operator:

      (ascalar*acoef - bscalar div bcoef grad) phi = RHS
      
      for an implicit discretization of the heat equation
      
      (I - div dt grad) phi^{n+1} = phi^n
  */
  Real dt_m;
  Real qij;
  int nf=1;

  const BoxArray &ba=phi_old.boxArray();
  const DistributionMapping &dm=phi_old.DistributionMap();

  const Box& domain_bx = geom.Domain();
  const Real* dx = geom.CellSize();
  // set the boundary conditions

  
  // relative and absolute tolerances for linear solve
  const Real tol_rel = 1.e-10;
  const Real tol_abs = 0.0;
  // Set up coefficient matrices
  MultiFab acoef(ba, dm, 1, 0);
  
  // fill in the acoef MultiFab and load this into the solver
  acoef.setVal(1.0);
  mlabec.setACoeffs(0, acoef);

  std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      const BoxArray& bamg = amrex::convert(acoef.boxArray(),
					  IntVect::TheDimensionVector(idim));
      face_bcoef[idim].define(bamg, acoef.DistributionMap(), 1, 0);
    }
  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));

  // Copy old phi into first SDC node
  MultiFab::Copy(SDC.sol[0],phi_old, 0, 0, 1, 0);

  // Fill the ghost cells of each grid from the other grids
  // includes periodic domain boundaries
  SDC.sol[0].FillBoundary(geom.periodicity());
  
  // Fill non-periodic physical boundaries
  FillDomainBoundary(SDC.sol[0], geom, bc);
  
  //  Compute the first function value
  int sdc_m=0;
  SDC_feval(flux,geom,bc,SDC,sdc_m,-1);

  // Copy first function value to all nodes
  for (int sdc_n = 1; sdc_n < SDC_NNODES; sdc_n++)
    {
      MultiFab::Copy(SDC.f[0][sdc_n],SDC.f[0][0], 0, 0, 1, 0);
      MultiFab::Copy(SDC.f[1][sdc_n],SDC.f[1][0], 0, 0, 1, 0);
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
	  dt_m = dt*(SDC.qnodes[sdc_m+1]-SDC.qnodes[sdc_m]);
	  
	  // use phi_new as rhs and fill it with terms at this iteration
	  SDC.SDC_rhs_k_plus_one(phi_new,dt,sdc_m);
	  
	  // get the best initial guess for implicit solve
	  MultiFab::Copy(SDC.sol[sdc_m+1],phi_new, 0, 0, 1, 0);
	  for ( MFIter mfi(SDC.sol[sdc_m+1]); mfi.isValid(); ++mfi )
	    {
	      const Box& bx = mfi.validbox();
	      qij = dt*SDC.qmats[2][sdc_m][sdc_m+1];
	      SDC.sol[sdc_m+1][mfi].saxpy(qij,SDC.f[1][sdc_m+1][mfi]);
	    }

	  // Solve for the implicit part
	  SDC_comp(phi_new, geom, bc, SDC, mlmg, mlabec,dt,sdc_m+1,1);

	  // Compute the function values at node sdc_m+1
	  SDC.sol[sdc_m+1].FillBoundary(geom.periodicity());
	  SDC_feval(flux,geom,bc,SDC,sdc_m+1,-1);	  

	} // end SDC substep loop
    }  // end sweeps loop
    
  // Return the last node in phi_new
  MultiFab::Copy(phi_new, SDC.sol[SDC_NNODES-1], 0, 0, 1, 0);

}

void SDC_feval(std::array<MultiFab, AMREX_SPACEDIM>& flux,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
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
  for (int n = nlo; n < nhi; n++)    
    for ( MFIter mfi(SDC.sol[sdc_m]); mfi.isValid(); ++mfi )
      {
      const Box& bx = mfi.validbox();
      compute_flux(BL_TO_FORTRAN_BOX(bx),
		   BL_TO_FORTRAN_BOX(domain_bx),
		   BL_TO_FORTRAN_ANYD(SDC.sol[sdc_m][mfi]),
		   BL_TO_FORTRAN_ANYD(flux[0][mfi]),
		   BL_TO_FORTRAN_ANYD(flux[1][mfi]),
#if (AMREX_SPACEDIM == 3)   
		   BL_TO_FORTRAN_ANYD(flux[2][mfi]),
#endif		       
		   BL_TO_FORTRAN_ANYD(SDC.f[n][sdc_m][mfi]),
		   dx,&n);
      }
  

}
void SDC_comp(MultiFab& rhs,
	       const Geometry& geom,
	       const Vector<BCRec>& bc,
	       SDCstruct &SDC,
	       MLMG &mlmg,
		 MLABecLaplacian& mlabec,	      
	      Real dt,
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
  // Set up coefficient matrices
  MultiFab acoef(ba, dm, 1, 0);
  
  // fill in the acoef MultiFab and load this into the solver
  acoef.setVal(1.0);
  mlabec.setACoeffs(0, acoef);

  std::array<MultiFab,AMREX_SPACEDIM> face_bcoef;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      const BoxArray& bamg = amrex::convert(acoef.boxArray(),
					  IntVect::TheDimensionVector(idim));
      face_bcoef[idim].define(bamg, acoef.DistributionMap(), 1, 0);
    }
  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));


  // Fill the ghost cells of each grid from the other grids
  // includes periodic domain boundaries
  rhs.FillBoundary(geom.periodicity());
	  
  // Fill non-periodic physical boundaries
  FillDomainBoundary(rhs, geom, bc);
	  
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
      qij = dt*SDC.qmats[2][sdc_m][sdc_m+1];	      
      face_bcoef[idim].setVal(0.1*qij);	      
    }
  mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(face_bcoef));
	  
  // set the boundary conditions
  mlabec.setLevelBC(0, &rhs);

  mlmg.solve({&SDC.sol[sdc_m]}, {&rhs}, tol_rel, tol_abs);  


}

