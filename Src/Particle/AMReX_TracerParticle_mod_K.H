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
void amrex_interpolate_CIC_2 (const P& p,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
			      const amrex::Array4<const double> &  uccarr,
			      amrex::Real * val){
  
  BL_ASSERT(val != 0);



 #if (AMREX_SPACEDIM == 1)
  
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5; //len

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

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5; 
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;

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

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

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
void amrex_interpolate_MAC_2 (const P& p,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
			      amrex::GpuArray<amrex::Array4<const double>,2> * p_uccarr[2],
			      amrex::Real * val,
			      amrex::IntVect& index,
			      int & lev)
{

  //  if (p.m_idata.id <=0) continue;

  amrex::Real lx = (p.m_rdata.pos[0]-plo[0])*dxi[0] + 0.5; //len 
  amrex::Real ly = (p.m_rdata.pos[1]-plo[1])*dxi[1] + 0.5;
  
  int i = std::floor(lx); //cell
  int j = std::floor(ly);
  
  amrex::Real xint = lx - i;
  amrex::Real yint = ly - j;

  i = i+index[0];
  j = j+index[1];

  auto xint2 = p.m_rdata.pos[0]-plo[0]*dxi[0] - index[0];
  auto yint2 = p.m_rdata.pos[1]-plo[1]*dxi[1] - index[1];
   
  if (xint2 > 1) xint2 =1;
  if (yint2 > 1) yint2 =1;
  if (xint2 < 0) xint2 =0;
  if (yint2 < 0) yint2 =0;

  amrex::Real sx[] = {1.-xint2, xint2};
  amrex::Real sy[] = {1.-yint2, yint2};

  int d;
    for (d=0 ; d < AMREX_SPACEDIM; ++d)
      {
      val = 0.0;
      for (int jj = 0; jj <= 1; ++jj)
        {
	for (int ii = 0; ii <= 1; ++ii)
	  {
	  val += (*p_uccarr)[d](i+ii, j+jj, 0,d)*sx[ii]*sy[jj];
	  }
	}
      }
}
}
#endif

