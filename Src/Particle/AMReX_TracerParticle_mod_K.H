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

//template <typename P> 
//AMREX_GPU_HOST_DEVICE AMREX_INLINE
//void amrex_interpolate_CIC_2(const int p, 
//			     const GeometryData& geomdata, 
//			     array4<Real> const& uccar, 
//			     int idx, 
//			     int val, 
//			     int comp){

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void amrex_interpolate_CIC_2 (const P& p, 
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
			      const amrex::Array4<const double> &  uccarr,
			      amrex::Real * val,
			      const int ncomp)  
{
  BL_ASSERT(val != 0);
  
  amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
  amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
  
  int i = std::floor(lx);
  int j = std::floor(ly);
  
  amrex::Real xint = lx - i;
  amrex::Real yint = ly - j;
  
  amrex::Real sx[] = {1.-xint, xint};
  amrex::Real sy[] = {1.-yint, yint};

  for (int comp=0; comp < ncomp; ++comp) {
      val[comp] = 0.0;
      for (int jj = 0; jj <= 1; ++jj) { 
          for (int ii = 0; ii <= 1; ++ii) { 
              val[comp] += sx[ii]*sy[jj]*uccarr(i+ii,j+jj,0,comp);
          }
      }
  }
}
#endif
