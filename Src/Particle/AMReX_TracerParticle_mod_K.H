#ifndef AMREX_TRACERPARTICLE_MOD_K_H
#define AMREX_TRACERPARTICLE_MOD_K_H
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_IntVect.H>
#include <AMReX_TracerParticles.H>

#include <cmath>

namespace amrex{

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate (const P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      const amrex::Array4<amrex::Real const> &  uccarr,
                      amrex::ParticleReal * val, int M = AMREX_SPACEDIM)
{
    AMREX_ASSERT(val != nullptr);

    AMREX_D_TERM(amrex::Real lx = (Real(p.pos(0)) - plo[0]) * dxi[0] - Real(0.5);,
                 amrex::Real ly = (Real(p.pos(1)) - plo[1]) * dxi[1] - Real(0.5);,
                 amrex::Real lz = (Real(p.pos(2)) - plo[2]) * dxi[2] - Real(0.5));

    AMREX_D_TERM(int const i = static_cast<int>(amrex::Math::floor(lx));,
                 int const j = static_cast<int>(amrex::Math::floor(ly));,
                 int const k = static_cast<int>(amrex::Math::floor(lz)));

    AMREX_D_TERM(amrex::Real const xint = lx - static_cast<Real>(i);,
                 amrex::Real const yint = ly - static_cast<Real>(j);,
                 amrex::Real const zint = lz - static_cast<Real>(k));

#if AMREX_SPACEDIM > 2
    amrex::Real sz[] = {Real(1.0) - zint, zint};
#endif
#if AMREX_SPACEDIM > 1
    amrex::Real sy[] = {Real(1.0) - yint, yint};
#endif
    amrex::Real sx[] = {Real(1.0) - xint, xint};

    for (int d=0; d < M; ++d) {
        val[d] = ParticleReal(0.0);
#if AMREX_SPACEDIM > 2
        for (int kk = 0; kk<=1; ++kk) {
#endif
#if AMREX_SPACEDIM > 1
            for (int jj = 0; jj <= 1; ++jj) {
#endif
                for (int ii = 0; ii <= 1; ++ii) {
                    val[d] += static_cast<ParticleReal>(AMREX_D_TERM(sx[ii],*sy[jj],*sz[kk])*uccarr(IntVect(AMREX_D_DECL(i+ii,j+jj,k+kk)),d));
    AMREX_D_TERM(},},});
    }
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void cic_interpolate_mapped (const P& p,
                             const amrex::Array4<amrex::Real const>&  ucc_arr,
                             const amrex::Array4<amrex::Real const>&  loc_arr,
                             amrex::ParticleReal* val, int M = AMREX_SPACEDIM)
{
    AMREX_ASSERT(val != nullptr);

    AMREX_D_TERM(int i = p.idata(0);,
                 int j = p.idata(1);,
                 int k = p.idata(2));

    AMREX_D_TERM(amrex::Real lx = amrex::Real(p.pos(0)) - 0.5*(loc_arr(IntVect(AMREX_D_DECL(i-1,j,k)),0) + loc_arr(IntVect(AMREX_D_DECL(i,j,k)),0));,
                 amrex::Real ly = amrex::Real(p.pos(1)) - 0.5*(loc_arr(IntVect(AMREX_D_DECL(i,j-1,k)),1) + loc_arr(IntVect(AMREX_D_DECL(i,j,k)),1));,
                 amrex::Real lz = amrex::Real(p.pos(2)) - 0.5*(loc_arr(IntVect(AMREX_D_DECL(i,j,k-1)),2) + loc_arr(IntVect(AMREX_D_DECL(i,j,k)),2)));

    AMREX_D_TERM(int i0 = i + static_cast<int>(amrex::Math::floor(lx/(loc_arr(IntVect(AMREX_D_DECL(i+1,j,k)),0) - loc_arr(IntVect(AMREX_D_DECL(i,j,k)),0))));,
                 int j0 = j + static_cast<int>(amrex::Math::floor(ly/(loc_arr(IntVect(AMREX_D_DECL(i,j+1,k)),1) - loc_arr(IntVect(AMREX_D_DECL(i,j,k)),1))));,
                 int k0 = k + static_cast<int>(amrex::Math::floor(lz/(loc_arr(IntVect(AMREX_D_DECL(i,j,k+1)),2) - loc_arr(IntVect(AMREX_D_DECL(i,j,k)),2)))));

    AMREX_D_TERM(amrex::Real const xint = 2.0*lx/(loc_arr(IntVect(AMREX_D_DECL(i0+2,j,k)),0) - loc_arr(IntVect(AMREX_D_DECL(i0,j,k)),0));,
                 amrex::Real const yint = 2.0*ly/(loc_arr(IntVect(AMREX_D_DECL(i,j0+2,k)),1) - loc_arr(IntVect(AMREX_D_DECL(i,j0,k)),1));,
                 amrex::Real const zint = 2.0*lz/(loc_arr(IntVect(AMREX_D_DECL(i,j,k0+2)),2) - loc_arr(IntVect(AMREX_D_DECL(i,j,k0)),2)));

#if AMREX_SPACEDIM > 2
    amrex::Real sz[] = {Real(1.0) - zint, zint};
#endif
#if AMREX_SPACEDIM > 1
    amrex::Real sy[] = {Real(1.0) - yint, yint};
#endif
    amrex::Real sx[] = {Real(1.0) - xint, xint};

    for (int d=0; d < M; ++d) {
        val[d] = ParticleReal(0.0);
#if AMREX_SPACEDIM > 2
        for (int kk = 0; kk<=1; ++kk) {
#endif
#if AMREX_SPACEDIM > 1
            for (int jj = 0; jj <= 1; ++jj) {
#endif
                for (int ii = 0; ii <= 1; ++ii) {
                    val[d] += static_cast<ParticleReal>(AMREX_D_TERM(sx[ii],*sy[jj],*sz[kk])*ucc_arr(IntVect(AMREX_D_DECL(i0+ii,j0+jj,k0+kk)),d));
    AMREX_D_TERM(},},});
    }
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mac_interpolate (const P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& umac_arr,
                      amrex::ParticleReal * val)
{
    AMREX_ASSERT(val != nullptr);

    for (int d=0; d < AMREX_SPACEDIM; ++d) {

        AMREX_D_TERM(amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(d != 0)*Real(0.5);,
                     amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - static_cast<Real>(d != 1)*Real(0.5);,
                     amrex::Real lz = (Real(p.pos(2))-plo[2])*dxi[2] - static_cast<Real>(d != 2)*Real(0.5));

        AMREX_D_TERM(int const i = static_cast<int>(amrex::Math::floor(lx));,
                     int const j = static_cast<int>(amrex::Math::floor(ly));,
                     int const k = static_cast<int>(amrex::Math::floor(lz)));

        AMREX_D_TERM(amrex::Real const xint = lx - static_cast<Real>(i);,
                     amrex::Real const yint = ly - static_cast<Real>(j);,
                     amrex::Real const zint = lz - static_cast<Real>(k));

#if AMREX_SPACEDIM > 2
        amrex::Real sz[] = {Real(1.0) - zint, zint};
#endif
#if AMREX_SPACEDIM > 1
      amrex::Real sy[] = {Real(1.0) - yint, yint};
#endif
      amrex::Real sx[] = {Real(1.0) - xint, xint};

      val[d] = ParticleReal(0.0);
#if AMREX_SPACEDIM > 2
      for (int kk = 0; kk <=1; ++kk) {
#endif
#if AMREX_SPACEDIM > 1
          for (int jj = 0; jj <= 1; ++jj) {
#endif
              for (int ii = 0; ii <= 1; ++ii) {
                  val[d] += static_cast<ParticleReal>((umac_arr[d])(IntVect(AMREX_D_DECL(i+ii, j+jj, k+kk)), 0)*AMREX_D_TERM(sx[ii],*sy[jj],*sz[kk]));
    AMREX_D_TERM(},},});
    }
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mac_interpolate_mapped_z (const P& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& umac_arr,
#if AMREX_SPACEDIM == 3
                               const amrex::Array4<amrex::Real const>&  zloc_arr,
                               bool use_terrain,
#endif
                               amrex::ParticleReal * val)
{
    AMREX_ASSERT(val != nullptr);

    for (int d=0; d < AMREX_SPACEDIM; ++d) {

#if AMREX_SPACEDIM == 3
        if (use_terrain && d == 2) {continue;}
#endif

        AMREX_D_TERM(amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - static_cast<Real>(d != 0)*Real(0.5);,
                     amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - static_cast<Real>(d != 1)*Real(0.5);,
                     amrex::Real lz = (Real(p.pos(2))-plo[2])*dxi[2] - static_cast<Real>(d != 2)*Real(0.5));

        AMREX_D_TERM(int const i = static_cast<int>(amrex::Math::floor(lx));,
                     int const j = static_cast<int>(amrex::Math::floor(ly));,
                     int const k = static_cast<int>(amrex::Math::floor(lz)));

        AMREX_D_TERM(amrex::Real const xint = lx - static_cast<Real>(i);,
                     amrex::Real const yint = ly - static_cast<Real>(j);,
                     amrex::Real const zint = lz - static_cast<Real>(k));

#if AMREX_SPACEDIM > 2
        amrex::Real sz[] = {Real(1.0) - zint, zint};
#endif
#if AMREX_SPACEDIM > 1
      amrex::Real sy[] = {Real(1.0) - yint, yint};
#endif
      amrex::Real sx[] = {Real(1.0) - xint, xint};

      val[d] = ParticleReal(0.0);
#if AMREX_SPACEDIM > 2
      for (int kk = 0; kk <=1; ++kk) {
#endif
#if AMREX_SPACEDIM > 1
          for (int jj = 0; jj <= 1; ++jj) {
#endif
              for (int ii = 0; ii <= 1; ++ii) {
                  val[d] += static_cast<ParticleReal>((umac_arr[d])(IntVect(AMREX_D_DECL(i+ii, j+jj, k+kk)), 0)*AMREX_D_TERM(sx[ii],*sy[jj],*sz[kk]));
    AMREX_D_TERM(},},});
    }

#if AMREX_SPACEDIM == 3
    if (use_terrain) {
        amrex::Real lx = (Real(p.pos(0))-plo[0])*dxi[0] - Real(0.5);
        amrex::Real ly = (Real(p.pos(1))-plo[1])*dxi[1] - Real(0.5);

        int const i = static_cast<int>(amrex::Math::floor(lx));
        int const j = static_cast<int>(amrex::Math::floor(ly));
        int const k = p.idata(0);

        amrex::Real const xint = lx - static_cast<Real>(i);
        amrex::Real const yint = ly - static_cast<Real>(j);
        amrex::Real const zint = (Real(p.pos(2)) - zloc_arr(IntVect(i,j,k))) / (zloc_arr(i,j,k+1) - zloc_arr(i,j,k));

        amrex::Real sz[] = {Real(1.0) - zint, zint};
        amrex::Real sy[] = {Real(1.0) - yint, yint};
        amrex::Real sx[] = {Real(1.0) - xint, xint};

      val[2] = ParticleReal(0.0);
      for (int kk = 0; kk <=1; ++kk) {
          for (int jj = 0; jj <= 1; ++jj) {
              for (int ii = 0; ii <= 1; ++ii) {
                  val[2] += static_cast<ParticleReal>((umac_arr[2])(IntVect(i+ii, j+jj, k+kk), 0)*sx[ii]*sy[jj]*sz[kk]);
              }
          }
      }
    }
#endif
}

}  // namespace amrex
#endif  // include guard
