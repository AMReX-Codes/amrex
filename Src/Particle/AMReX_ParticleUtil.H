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
int getParticleGrid (P const& p, amrex::Array4<int> const& mask,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                     const Box& domain)
{
    if (p.id() < 0) return -1;
    IntVect iv = IntVect(
            AMREX_D_DECL(floor((p.pos(0)-plo[0])*dxi[0]),
                         floor((p.pos(1)-plo[1])*dxi[1]),
                         floor((p.pos(2)-plo[2])*dxi[2])));
    iv += domain.smallEnd();
    return mask(iv);
}
    
}

#endif // include guard
