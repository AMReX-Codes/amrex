#ifndef WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
#define WARPX_PARTICLES_SORTING_SORTINGUTILS_H_

#include <AMReX_Gpu.H>
#include <AMReX_CudaContainers.H>

// TODO: Add documentation
void fillWithConsecutiveIntegers( amrex::Gpu::ManagedDeviceVector<long>& v ) {
#ifdef AMREX_USE_GPU
    // On GPU: Use thrust
    thrust::sequence( v.begin(), v.end() );
#else
    // On CPU: Use std library
    std::iota( v.begin(), v.end(), 0L );
#endif
}

// TODO: Add documentation
template< typename PType >
AMREX_GPU_HOST_DEVICE AMREX_INLINE
amrex::IntVect findParticleIndex ( PType p,
                    const amrex::Real prob_lo_x,
                    const amrex::Real prob_lo_y,
                    const amrex::Real prob_lo_z,
                    const amrex::Real inv_dx,
                    const amrex::Real inv_dy,
                    const amrex::Real inv_dz,
                    const amrex::IntVect domain_small_end )
{
    amrex::IntVect iv;
    AMREX_D_TERM(iv[0]=static_cast<int>(floor((p.m_rdata.pos[0]-prob_lo_x)*inv_dx));,
                 iv[1]=static_cast<int>(floor((p.m_rdata.pos[1]-prob_lo_y)*inv_dy));,
                 iv[2]=static_cast<int>(floor((p.m_rdata.pos[2]-prob_lo_z)*inv_dz)););

    iv += domain_small_end;

    return iv;
};

#endif // WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
