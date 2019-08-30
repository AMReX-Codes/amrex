#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_MFIter.H>

#include <limits>

namespace amrex
{

AMREX_GPU_HOST_DEVICE
int getTileIndex (const IntVect& iv, const Box& box, const bool a_do_tiling, 
                  const IntVect& a_tile_size, Box& tbx);

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getParticleCell (P const& p,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const Box& domain) noexcept
{
    IntVect iv = IntVect(
        AMREX_D_DECL(floor((p.pos(0)-plo[0])*dxi[0]),
                     floor((p.pos(1)-plo[1])*dxi[1]),
                     floor((p.pos(2)-plo[2])*dxi[2])));
    iv += domain.smallEnd();
    return iv;
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int getParticleGrid (P const& p, amrex::Array4<int> const& mask,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                     const Box& domain) noexcept
{
    if (p.id() < 0) return -1;    
    IntVect iv = getParticleCell(p, plo, dxi, domain);
    return mask(iv);
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void enforcePeriodic (P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& phi,
                      amrex::GpuArray<int,AMREX_SPACEDIM> const& is_per) noexcept
{
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (not is_per[idim]) continue;
        if (p.pos(idim) > phi[idim]) { 
            p.pos(idim) -= (phi[idim] - plo[idim]); 
        } else if (p.pos(idim) < plo[idim]) {
            p.pos(idim) += (phi[idim] - plo[idim]);
        }
    }
}

template <class PC>
bool
numParticlesOutOfRange (PC const& pc, int nGrow)
{
    return numParticlesOutOfRange(pc, 0, pc.finestLevel(), nGrow);
}

template <class PC>
bool
numParticlesOutOfRange (PC const& pc, int lev_min, int lev_max, int nGrow)
{
    BL_PROFILE("numParticlesOutOfRange()");

    using ParIter = typename PC::ParConstIterType;
    using ParticleType = typename PC::ParticleType;
    int num_wrong = 0;
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        const auto plo = pc.Geom(lev).ProbLoArray();
        const auto dxi = pc.Geom(lev).InvCellSizeArray();
        const Box domain = pc.Geom(lev).Domain();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion() and !system::regtest_reduction) reduction(+:num_wrong)
#endif
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            const auto& tile = pti.GetParticleTile();
            const auto np = tile.numParticles();
            const auto& aos = tile.GetArrayOfStructs();
            const auto pstruct = aos().dataPtr();
            
            Box box = pti.tilebox();
            box.grow(nGrow);

            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(np, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
            {
                const ParticleType& p = pstruct[i];
                IntVect iv = IntVect(
                        AMREX_D_DECL(floor((p.pos(0)-plo[0])*dxi[0]),
                                     floor((p.pos(1)-plo[1])*dxi[1]),
                                     floor((p.pos(2)-plo[2])*dxi[2])));
                iv += domain.smallEnd();
                return !box.contains(iv);
            });
            int hv = amrex::get<0>(reduce_data.value());
            num_wrong += hv;
        }
    }
    ParallelDescriptor::ReduceIntSum(num_wrong);

    return num_wrong;
}

}

#endif // include guard
