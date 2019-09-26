#ifndef WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
#define WARPX_PARTICLES_SORTING_SORTINGUTILS_H_

#include <WarpXParticleContainer.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_Gpu.H>

// TODO: Add documentation
void fillWithConsecutiveIntegers( amrex::Gpu::ManagedDeviceVector<long>& v )
{
#ifdef AMREX_USE_GPU
    // On GPU: Use thrust
    thrust::sequence( v.begin(), v.end() );
#else
    // On CPU: Use std library
    std::iota( v.begin(), v.end(), 0L );
#endif
}

// TODO: Add documentation
template< typename ForwardIterator >
ForwardIterator stablePartition(ForwardIterator index_begin,
                               ForwardIterator index_end,
                               amrex::Gpu::ManagedDeviceVector<int>& predicate)
{
#ifdef AMREX_USE_GPU
    // On GPU: Use thrust
    int* AMREX_RESTRICT predicate_ptr = predicate.dataPtr();
    ForwardIterator sep = thrust::stable_partition(
                    thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                    index_begin, index_end,
                    AMREX_GPU_HOST_DEVICE
                    [predicate_ptr](long i) { return predicate_ptr[i]; }
                );
#else
    // On CPU: Use std library
    ForwardIterator sep = std::stable_partition(
                    index_begin, index_end,
                    [&predicate](long i) { return predicate[i]; }
                );
#endif
    return sep;
}

// TODO: Add documentation
template< typename ForwardIterator >
int iteratorDistance(ForwardIterator first,
                     ForwardIterator last)
{
#ifdef AMREX_USE_GPU
    // On GPU: Use thrust
    return thrust::distance( first, last );
#else
    return std::distance( first, last );
#endif
}

// TODO: Add documentation
class fillBufferFlag
{
    public:
        fillBufferFlag( WarpXParIter& pti, const amrex::iMultiFab* bmasks,
                        amrex::Gpu::ManagedDeviceVector<int>& inexflag,
                        const amrex::Geometry& geom ) {

            // Extract simple structure that can be used directly on the GPU
            m_particles = &(pti.GetArrayOfStructs()[0]);
            m_buffer_mask = (*bmasks)[pti].array();
            m_inexflag_ptr = inexflag.dataPtr();
            m_domain_small_end = geom.Domain().smallEnd();
            for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
                m_prob_lo[idim] = geom.ProbLo(idim);
                m_inv_cell_size[idim] = geom.InvCellSize(idim);
            }
        };


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator()( const long i ) const {
            // Select a particle
            auto& p = m_particles[i];
            // Find the index of the cell where this particle is located
            amrex::IntVect iv;
            for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
                iv[idim] = static_cast<int>(floor(
                    (p.m_rdata.pos[idim]-m_prob_lo[idim])*m_inv_cell_size[idim]
                ));
            }
            // Find the value of the buffer flag in this cell and
            // store it at the corresponding particle position in the array `inexflag`
            m_inexflag_ptr[i] = m_buffer_mask(iv);
        };

    private:
        amrex::Real m_prob_lo[AMREX_SPACEDIM];
        amrex::Real m_inv_cell_size[AMREX_SPACEDIM];
        amrex::IntVect m_domain_small_end;
        int* m_inexflag_ptr;
        WarpXParticleContainer::ParticleType* m_particles;
        amrex::Array4<const int> m_buffer_mask;
};

#endif // WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
