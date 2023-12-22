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
void cic_interpolate_mapped_z (const P& p,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                               const amrex::Array4<amrex::Real const>&  uccarr,
                               const amrex::Array4<amrex::Real const>&  zheight,
                               bool use_terrain, amrex::ParticleReal* val, int M = AMREX_SPACEDIM)
{
    AMREX_ASSERT(val != nullptr);

    // handle x and y normally
    amrex::Real lx = (Real(p.pos(0)) - plo[0]) * dxi[0] - Real(0.5);
    amrex::Real ly = (Real(p.pos(1)) - plo[1]) * dxi[1] - Real(0.5);

    int const i = static_cast<int>(amrex::Math::floor(lx));
    int const j = static_cast<int>(amrex::Math::floor(ly));

    amrex::Real const xint = lx - static_cast<Real>(i);
    amrex::Real const yint = ly - static_cast<Real>(j);

    amrex::Real sy[] = {Real(1.0) - yint, yint};
    amrex::Real sx[] = {Real(1.0) - xint, xint};

    // handle mapped z if in use
    int const k = p.idata(0);
    ParticleReal zlo, zhi;
    if (use_terrain) {
        zlo = zheight(i, j, k  );
        zhi = zheight(i, j, k+1);
    } else {
        zlo =  k    / dxi[2];
        zhi = (k+1) / dxi[2];
    }

    amrex::Real lz = Real((p.pos(2) - zlo) / (zhi - zlo)) - Real(0.5);
    amrex::Real const zint = lz;
    amrex::Real sz[] = {Real(1.0) - zint, zint};

    for (int d=0; d < M; ++d) {
        val[d] = ParticleReal(0.0);
        for (int kk = 0; kk<=1; ++kk) {
            for (int jj = 0; jj <= 1; ++jj) {
                for (int ii = 0; ii <= 1; ++ii) {
                    val[d] += static_cast<ParticleReal>(sx[ii]*sy[jj]*sz[kk]*uccarr(i+ii,j+jj,k+kk),d);
                }
            }
        }
    }
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mac_interpolate (const P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                      amrex::GpuArray<amrex::Array4<amrex::Real const>,AMREX_SPACEDIM> const& p_uccarr,
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
                  val[d] += static_cast<ParticleReal>((p_uccarr[d])(IntVect(AMREX_D_DECL(i+ii, j+jj, k+kk)), 0)*AMREX_D_TERM(sx[ii],*sy[jj],*sz[kk]));
    AMREX_D_TERM(},},});
    }
}

}  // namespace amrex
#endif  // include guard
