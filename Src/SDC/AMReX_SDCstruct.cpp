#include "AMReX_SDCstruct.H"


SDCstruct::SDCstruct(int Nnodes_in,int Npieces_in, MultiFab& sol_in)
{
  
  Nnodes=Nnodes_in;
  Npieces=Npieces_in;       
  
  qnodes= new Real[Nnodes];
  Nflags= new int[Nnodes];
  //  qmats = new Real[4][Nnodes-1][Nnodes];    


  //  Make the quadrature tables
  SDC_quadrature(&qtype, &Nnodes, &Nnodes,qnodes,Nflags, &qmats[0][0][0]);
  
  sol.resize(Nnodes);
  res.resize(Nnodes);
  f.resize(Npieces);
  if (Npieces == 3) 
    Ithree.resize(Nnodes);

  //  Assign  geomety and multifab info
  const BoxArray &ba=sol_in.boxArray();
  const DistributionMapping &dm=sol_in.DistributionMap();
  const int Nghost=sol_in.nGrow();
  const int Ncomp=sol_in.nComp();

  for (auto& v : f) v.resize(Nnodes);  
  for (int sdc_m = 0; sdc_m < Nnodes; sdc_m++)
    {
      sol[sdc_m].define(ba, dm, Ncomp, Nghost);
      res[sdc_m].define(ba, dm, Ncomp, Nghost);
      for (int i = 0; i < Npieces; i++)
	{
	  f[i][sdc_m].define(ba, dm, Ncomp, Nghost);
	}
      if (Npieces == 3)
	Ithree[sdc_m].define(ba, dm, Ncomp, Nghost);      
    }
  
}

void SDCstruct::SDC_rhs_integrals(Real dt)
{

  Real qij;

  // Compute the quadrature terms from last iteration
  for (int sdc_m = 0; sdc_m < Nnodes-1; sdc_m++)
    {
      for ( MFIter mfi(res[sdc_m]); mfi.isValid(); ++mfi )
	{
	  res[sdc_m].setVal(0.0);
	  if (Npieces == 3)
	    Ithree[sdc_m].setVal(0.0);
	  for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
	    {
	      qij = dt*(qmats[0][sdc_m][sdc_n]-qmats[1][sdc_m][sdc_n]);
	      res[sdc_m][mfi].saxpy(qij,f[0][sdc_n][mfi]);
	      qij = dt*(qmats[0][sdc_m][sdc_n]-qmats[2][sdc_m][sdc_n]);
	      res[sdc_m][mfi].saxpy(qij,f[1][sdc_n][mfi]);
	    }
	  if (Npieces == 3)
	    {
	      for (int sdc_n = 0; sdc_n < Nnodes; sdc_n++)
		{ //  // MISDC pieces
		  qij = dt*(qmats[0][sdc_m][sdc_n]);  // leave off -dt*Qtil and add it later
		  res[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		  // Compute seperate integral for f_3 piece		  
		  qij = -dt*(qmats[2][sdc_m][sdc_n]);  
		  Ithree[sdc_m][mfi].saxpy(qij,f[2][sdc_n][mfi]);
		}
	    }
	}
    }
}

void SDCstruct::SDC_rhs_k_plus_one(MultiFab& sol_new, Real dt,int sdc_m)
{
  //  Compute the rhs terms for the implicit solve
  Real qij;
  
  //  Copy first the initial value
  MultiFab::Copy(sol_new,sol[0], 0, 0, 1, 0);
  for ( MFIter mfi(sol_new); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();
      sol_new[mfi].saxpy(1.0,res[sdc_m][mfi]);
      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
	{
	  qij = dt*qmats[1][sdc_m][sdc_n];
	  sol_new[mfi].saxpy(qij,f[0][sdc_n][mfi]);
	}
      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
	{
	  qij = dt*qmats[2][sdc_m][sdc_n];
	  sol_new[mfi].saxpy(qij,f[1][sdc_n][mfi]);
	}
    }
  
}
void SDCstruct::SDC_rhs_misdc(MultiFab& sol_new, Real dt,int sdc_m)
{
  //  Add the terms to the rhs before the second implicit solve
  Real qij;
  
  for ( MFIter mfi(sol_new); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();
      sol_new[mfi].saxpy(1.0,Ithree[sdc_m][mfi]);
      for (int sdc_n = 0; sdc_n < sdc_m+1; sdc_n++)
	{
	  qij = dt*qmats[2][sdc_m][sdc_n];
	  sol_new[mfi].saxpy(qij,f[2][sdc_n][mfi]);
	}
    }
  
}
