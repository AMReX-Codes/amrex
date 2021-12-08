
#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Particles.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_ParticleInterpolators.H>

struct TrilinearDeposition
{
    int start_part_comp;
    int start_mesh_comp;
    int   num_comp;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void operator() (const MyParticleContainer::ParticleType& p,
                     amrex::Array4<amrex::Real> const& rho,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi) const noexcept
    {
        amrex::ParticleInterpolator::Linear interp(p, plo, dxi);

        interp.ParticleToMesh(p, rho, start_part_comp, start_mesh_comp, num_comp,
            [=] AMREX_GPU_DEVICE (const MyParticleContainer::ParticleType& part, int comp)
            {
                return part.rdata(comp);  // no weighting
            });
    }
};
