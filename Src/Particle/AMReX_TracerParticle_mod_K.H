#ifndef AMREX_TRACERPARTICLE_MOD_K_H
#define AMREX_TRACERPARTICLE_MOD_K_H
#include <AMReX_TracerParticles.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <cmath>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include "AMReX_TracerParticles.H"


namespace amrex{
  
template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void cic_interpolate (const P& p,
		      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
		      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
		      const amrex::Array4<const double> &  uccarr,
		      amrex::Real * val){
  
  BL_ASSERT(val != 0);

#if (AMREX_SPACEDIM == 1)
  
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] - 0.5; //len

    int i = std::floor(lx); //cell

    amrex::Real xint = lx - i; //frac

    amrex::Real sx[] = {1.-xint, xint}; 

    for int d=0; d < AMREX_SPACEDIM; ++d)
      {
      for (int ii = 0; ii<=1; ++ii)
	{
        val[d] += sx[ii]*uccarr(i+ii,0,0,d)
	}
      }

  
#elif (AMREX_SPACEDIM == 2)

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] - 0.5; 
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] - 0.5;

    int i = std::floor(lx); 
    int j = std::floor(ly);

    amrex::Real xint = lx - i; 
    amrex::Real yint = ly - j;

    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};

    for (int d=0; d < AMREX_SPACEDIM; ++d)
      {
        val[d] = 0.0;
        for (int jj = 0; jj <= 1; ++jj)
	  {
            for (int ii = 0; ii <= 1; ++ii)
	      {
                val[d] += sx[ii]*sy[jj]*uccarr(i+ii,j+jj,0,d);
              }
          }
      }


#elif (AMREX_SPACEDIM == 3)

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] - 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] - 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] - 0.5;

    int i = std::floor(lx);
    int j = std::floor(ly);
    int k = std::floor(lz);

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;

    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};
    amrex::Real sz[] = {1.-zint, zint};

    for (int d=0; d < AMREX_SPACEDIM; ++d)
      {
      val[d] = 0.0;
      for (int kk = 0; kk<=1; ++kk)
	{
	for (int jj = 0; jj <= 1; ++jj)
	  {
          for (int ii = 0; ii <= 1; ++ii)
	    {
              val[d] += sx[ii]*sy[jj]*sz[kk]*uccarr(i+ii,j+jj,k+kk,d);
            }
	  }
        }
      }


#endif
}



template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void mac_interpolate (const P& p,
		      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
		      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
		      amrex::GpuArray<amrex::Array4<const double>,AMREX_SPACEDIM> const& p_uccarr,
		      amrex::Real * val)
{

#if (AMREX_SPACEDIM == 1)
  
  for (int d=0; d < AMREX_SPACEDIM; ++d)
    {
      amrex::Real lx = (p.m_rdata.pos[0]-plo[0])*dxi[0] - 0.5; //len
  
      int i = std::floor(lx); //cell
  
      amrex::Real xint = lx - i; //frac
  
      int ecellx = i;
  
      Real efracx = xint;
  

      if (d==0)
        {
          //cells
          ecellx = i +1;
          //fracs
          efracx = p.m_rdata.pos[0]-plo[0]*dxi[0] - floor((p.m_rdata.pos[0]-plo[0])*dxi[0])
        }

      if (efracx > 1) efracx =1;
      if (efracx < 0) efracx =0;
      
      amrex::Real sx[] = {1.-efracx, efracx};
      
      val[d] = 0.0;
      for (int ii = 0; ii <= 1; ++ii)
	{
	  val[d] += (p_uccarr[d])(ecellx+ii, 0 0,0)*sx[ii];
	}		
    }



    
#elif (AMREX_SPACEDIM == 2)
  
  for (int d=0; d < AMREX_SPACEDIM; ++d)
    {
      amrex::Real lx = (p.m_rdata.pos[0]-plo[0])*dxi[0] - 0.5; //len 
      amrex::Real ly = (p.m_rdata.pos[1]-plo[1])*dxi[1] - 0.5;

      int i = std::floor(lx); //cell
      int j = std::floor(ly);
  
      amrex::Real xint = lx - i; //frac
      amrex::Real yint = ly - j;
      int ecellx = i;
      int ecelly = j;
      Real efracx = xint;
      Real efracy = yint; 

      // This is a MAC grid, so in the x direction, y/z will be CC, 
      // but x will be on faces, adding 1 extra index
      if (d==0)
	{
	  //cells 
	  ecellx = i +1;
	  ecelly = j;
	  
	  //fracs
	  efracx = p.m_rdata.pos[0]-plo[0]*dxi[0] - std::floor((p.m_rdata.pos[0]-plo[0])*dxi[0]);
	  efracy = yint;
	}
      
      else if (d==1)
	{
	  //cells
	  ecellx = i;
	  ecelly = j +1;
	
	  //fracs
	  efracx = xint;
	  efracy = p.m_rdata.pos[1]-plo[1]*dxi[1] - std::floor((p.m_rdata.pos[1]-plo[1])*dxi[1]);
	}

      //The fraction must be between zero and one. 
      if (efracx > 1) efracx =1;
      if (efracy > 1) efracy =1;
      if (efracx < 0) efracx =0;
      if (efracy < 0) efracy =0;
      
      
      amrex::Real sx[] = {1.-efracx, efracx};
      amrex::Real sy[] = {1.-efracy, efracy};
      
    val[d] = 0.0;
    for (int jj = 0; jj <= 1; ++jj)
      {
	for (int ii = 0; ii <= 1; ++ii)
	  {
	    val[d] += (p_uccarr[d])(ecellx+ii, ecelly+jj, 0,0)*sx[ii]*sy[jj];
	  }
      }
    
    }


#elif (AMREX_SPACEDIM == 3)

  for (int d=0; d < AMREX_SPACEDIM; ++d)
    {
      amrex::Real lx = (p.m_rdata.pos[0]-plo[0])*dxi[0] - 0.5; //len                                                               
      amrex::Real ly = (p.m_rdata.pos[1]-plo[1])*dxi[1] - 0.5;
      amrex::Real lz = (p.m_rdata.pos[2]-plo[2])*dxi[2] - 0.5;

      
      int i = std::floor(lx); //cell                                                                                               
      int j = std::floor(ly);
      int k = std::floor(lz);
      
      amrex::Real xint = lx - i; //frac                                                                                            
      amrex::Real yint = ly - j;
      amrex::Real zint = lz - k;
      
      int  ecellx, ecelly, ecellz;
      Real efracx, efracy, efracz;

      // This is a MAC grid, so in the x direction, y/z will be CC, 
      // but x will be on faces, adding 1 extra index 
      if (d==0)
        {
          //cells
          ecellx = i +1;
          ecelly = j;
	  ecellz = k;

          //fracs
          efracx = p.m_rdata.pos[0]-plo[0]*dxi[0] - std::floor((p.m_rdata.pos[0]-plo[0])*dxi[0]);
          efracy = yint;
	  efracz = zint;
        }

      else if (d==1)
        {
          //cells
          ecellx = i;
          ecelly = j +1;
	  ecellz = k;

          //fracs
          efracx = xint;
          efracy = p.m_rdata.pos[1]-plo[1]*dxi[1] - std::floor((p.m_rdata.pos[1]-plo[1])*dxi[1]);
	  ecellz = zint;
	}

      else if (d==2)
        {
          //cells
          ecellx = i;
          ecelly = j;
          ecellz = k+1;

          //fracs
          efracx = xint;
	  efracy = yint;
          efracz = p.m_rdata.pos[2]-plo[2]*dxi[2] - std::floor((p.m_rdata.pos[2]-plo[2])*dxi[2]);
        }
      
      // The fraction must be between 0 and 1.
      if (efracx > 1) efracx =1;
      if (efracy > 1) efracy =1;
      if (efracz > 1) efracz =1;
      if (efracx < 0) efracx =0;
      if (efracy < 0) efracy =0;
      if (efracz < 0) efracz =0;

      amrex::Real sx[] = {1.-efracx, efracx};
      amrex::Real sy[] = {1.-efracy, efracy};
      amrex::Real sz[] = {1.-efracz, efracz};
      
      val[d] = 0.0;
      for (int kk = 0; kk <=1; ++kk)
	{
	  for (int jj = 0; jj <= 1; ++jj)
	    {
	      for (int ii = 0; ii <= 1; ++ii)
		{
		  val[d] += (p_uccarr[d])(ecellx+ii, ecelly+jj, ecellz+kk ,0)*sx[ii]*sy[jj]*sz[kk];
		}
	    }
	}
      
    } 
#endif
}
}
#endif
