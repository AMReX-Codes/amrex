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
			                        amrex::Real * val)
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

  for (int d=0; d < AMREX_SPACEDIM; ++d) {
      val[d] = 0.0;
      for (int jj = 0; jj <= 1; ++jj) {
          for (int ii = 0; ii <= 1; ++ii) {
              val[d] += sx[ii]*sy[jj]*uccarr(i+ii,j+jj,0,d);
          }
      }
  }
}



template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void amrex_interpolate_MAC_2 (const P& p,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                              amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
			                        const amrex::Array4<const double> &  uccarr,
			                        amrex::Real * val)
{

  if (p.m_idata.id <=0) continue;

  amrex::Real lx = (p.m_rdata.pos[0]-plo[0])*dxi[0] + Real(0.5);  //change to  multiply inverse dx
  amrex::Real ly = (p.m_rdata.pos[1]-plo[1])*dxi[1] + Real(0.5);
  amrex::Real lz = (p.m_rdata.pos[2]-plo[2])*dxi[2] + Real(0.5);

  int i = std::floor(lx);
  int j = std::floor(ly);
  int k = std::floor(lz);

  amrex::Real xint = lx - i;
  amrex::Real yint = ly - j;
  amrex::Real zint = lz - k;


  const IntVect& cc_cell = Index(p,lev);

  ecellx = cc_cell[0] + 1.;
  ecelly = cc_cell[1] + 1.;
  ecellz = cc_cell[2] + 1.;

  efracx = (p.m_rdata.pos[0] - plo[0])*dxi[0] - cc_cell[0];
  efracy = (p.m_rdata.pos[1] - plo[1])*dxi[1] - cc_cell[1];
  efracz = (p.m_rdata.pos[2] - plo[2])*dxi[2] - cc_cell[2];

  // If fraction is greater than 1 or less than 0, set it to 1 or 0.
  if (efracx > 1) efracx = 1.;
  if (efracy > 1) efracy = 1.;
  if (efracz > 1) efracz = 1.;

  if (efracx < 0) efracx = 0.;
  if (efracy < 0) efracy = 0.;
  if (efracz < 0) efracz = 0.;

  for (int d=0; d < AMREX_SPACEDIM; ++d) {
      val[d] = 0.0;
      for (int jj = 0; jj <= 1; ++jj) {
          for (int ii = 0; ii <= 1; ++ii) {
              val[d] += sx[ii]*sy[jj]*uccarr(i+ii,j+jj,0,d);






  }












}

#endif
