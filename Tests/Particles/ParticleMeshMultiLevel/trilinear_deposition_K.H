
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_AmrParticles.H>

struct TrilinearDeposition
{
  int start_part_comp;
  int start_mesh_comp;
  int   num_comp;

  AMREX_GPU_DEVICE
  void operator()
                  (const MyParticleContainer::ParticleType& p,
                   amrex::Array4<amrex::Real> const& rho,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) const noexcept
  {
      amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
      amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
      amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

      int i = static_cast<int>(amrex::Math::floor(lx));
      int j = static_cast<int>(amrex::Math::floor(ly));
      int k = static_cast<int>(amrex::Math::floor(lz));

      amrex::Real xint = lx - i;
      amrex::Real yint = ly - j;
      amrex::Real zint = lz - k;

      amrex::Real sx[] = {1.-xint, xint};
      amrex::Real sy[] = {1.-yint, yint};
      amrex::Real sz[] = {1.-zint, zint};

      for (int kk = 0; kk <= 1; ++kk) {
          for (int jj = 0; jj <= 1; ++jj) {
              for (int ii = 0; ii <= 1; ++ii) {
                  for (int icomp = 0; icomp < num_comp; icomp++)
                  {
                      amrex::Gpu::Atomic::AddNoRet(&rho(i+ii-1, j+jj-1, k+kk-1, start_mesh_comp+icomp),
                                                   sx[ii]*sy[jj]*sz[kk]*p.rdata(start_part_comp+icomp));
                  }
              }
          }
      }
  }
};
