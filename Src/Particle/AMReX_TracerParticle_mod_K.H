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
		      const amrex::Array4<amrex::Real const> &  uccarr,
		      amrex::ParticleReal * val){

  AMREX_ASSERT(val != 0);

#if (AMREX_SPACEDIM == 1)

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] - 0.5; //len

    int const i = static_cast<int>(amrex::Math::floor(lx)); //cell

    amrex::Real xint = lx - i; //frac

    amrex::Real sx[] = {Real(1.0) - xint, xint};

    for (int d=0; d < AMREX_SPACEDIM; ++d)
    {
        val[d] = ParticleReal(0.0);
        for (int ii = 0; ii<=1; ++ii)
	{
            val[d] += static_cast<ParticleReal>(sx[ii]*uccarr(i+ii,0,0,d));
        }
    }


#elif (AMREX_SPACEDIM == 2)

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] - 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] - 0.5;

    int const i = static_cast<int>(amrex::Math::floor(lx));
    int const j = static_cast<int>(amrex::Math::floor(ly));

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;

    amrex::Real sx[] = {Real(1.0) - xint, xint};
    amrex::Real sy[] = {Real(1.0) - yint, yint};

    for (int d=0; d < AMREX_SPACEDIM; ++d)
      {
        val[d] = ParticleReal(0.0);
        for (int jj = 0; jj <= 1; ++jj)
	  {
            for (int ii = 0; ii <= 1; ++ii)
	      {
                  val[d] += static_cast<ParticleReal>(sx[ii]*sy[jj]*uccarr(i+ii,j+jj,0,d));
              }
          }
      }


#elif (AMREX_SPACEDIM == 3)

    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] - 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] - 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] - 0.5;

    int const i = static_cast<int>(amrex::Math::floor(lx));
    int const j = static_cast<int>(amrex::Math::floor(ly));
    int const k = static_cast<int>(amrex::Math::floor(lz));

    amrex::Real const xint = lx - i;
    amrex::Real const yint = ly - j;
    amrex::Real const zint = lz - k;

    amrex::Real sx[] = {Real(1.0) - xint, xint};
    amrex::Real sy[] = {Real(1.0) - yint, yint};
    amrex::Real sz[] = {Real(1.0) - zint, zint};

    for (int d=0; d < AMREX_SPACEDIM; ++d)
    {
        val[d] = ParticleReal(0.0);
        for (int kk = 0; kk<=1; ++kk)
	{
            for (int jj = 0; jj <= 1; ++jj)
            {
                for (int ii = 0; ii <= 1; ++ii)
                {
                    val[d] += static_cast<ParticleReal>(sx[ii]*sy[jj]*sz[kk]*uccarr(i+ii,j+jj,k+kk,d));
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
		      amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& p_uccarr,
		      amrex::ParticleReal * val)
{

#if (AMREX_SPACEDIM == 1)
    for (int d=0; d < AMREX_SPACEDIM; ++d)
    {
        amrex::Real lx = (p.pos(0)-plo[0])*dxi[0] - (d != 0)*0.5;

        int const i = static_cast<int>(amrex::Math::floor(lx));

        amrex::Real const xint = lx - i;

        amrex::Real sx[] = {Real(1.0) - xint, xint};

        val[d] = ParticleReal(0.0);
        for (int ii = 0; ii <= 1; ++ii)
        {
            val[d] += static_cast<ParticleReal>((p_uccarr[d])(i+ii, 0, 0, 0)*sx[ii]);
        }
    }

#elif (AMREX_SPACEDIM == 2)

  for (int d=0; d < AMREX_SPACEDIM; ++d)
  {
      amrex::Real lx = (p.pos(0)-plo[0])*dxi[0] - (d != 0)*0.5;
      amrex::Real ly = (p.pos(1)-plo[1])*dxi[1] - (d != 1)*0.5;

      int const i = static_cast<int>(amrex::Math::floor(lx));
      int const j = static_cast<int>(amrex::Math::floor(ly));

      amrex::Real const xint = lx - i;
      amrex::Real const yint = ly - j;

      amrex::Real sx[] = {Real(1.0) - xint, xint};
      amrex::Real sy[] = {Real(1.0) - yint, yint};

      val[d] = ParticleReal(0.0);
      for (int jj = 0; jj <= 1; ++jj)
      {
          for (int ii = 0; ii <= 1; ++ii)
          {
              val[d] += static_cast<ParticleReal>((p_uccarr[d])(i+ii, j+jj, 0, 0)*sx[ii]*sy[jj]);
          }
      }
  }


#elif (AMREX_SPACEDIM == 3)

  for (int d=0; d < AMREX_SPACEDIM; ++d)
  {
      amrex::Real lx = (p.pos(0)-plo[0])*dxi[0] - (d != 0)*0.5;
      amrex::Real ly = (p.pos(1)-plo[1])*dxi[1] - (d != 1)*0.5;
      amrex::Real lz = (p.pos(2)-plo[2])*dxi[2] - (d != 2)*0.5;

      int const i = static_cast<int>(amrex::Math::floor(lx));
      int const j = static_cast<int>(amrex::Math::floor(ly));
      int const k = static_cast<int>(amrex::Math::floor(lz));

      amrex::Real const xint = lx - i;
      amrex::Real const yint = ly - j;
      amrex::Real const zint = lz - k;

      amrex::Real sx[] = {Real(1.0) - xint, xint};
      amrex::Real sy[] = {Real(1.0) - yint, yint};
      amrex::Real sz[] = {Real(1.0) - zint, zint};

      val[d] = ParticleReal(0.0);
      for (int kk = 0; kk <=1; ++kk)
      {
	  for (int jj = 0; jj <= 1; ++jj)
          {
	      for (int ii = 0; ii <= 1; ++ii)
              {
		  val[d] += static_cast<ParticleReal>((p_uccarr[d])(i+ii, j+jj, k+kk ,0)*sx[ii]*sy[jj]*sz[kk]);
              }
          }
      }
  }
#endif
}
}
#endif
