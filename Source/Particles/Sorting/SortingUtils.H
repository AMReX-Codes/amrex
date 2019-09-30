#ifndef WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
#define WARPX_PARTICLES_SORTING_SORTINGUTILS_H_

#include <WarpXParticleContainer.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_Gpu.H>

/* \brief Fill the elements of the input vector with consecutive integer,
 *        starting from 0
 *
 * \param[inout] v Vector of integers, to be filled by this routine
 */
void fillWithConsecutiveIntegers( amrex::Gpu::DeviceVector<long>& v )
{
#ifdef AMREX_USE_GPU
    // On GPU: Use thrust
    thrust::sequence( v.begin(), v.end() );
#else
    // On CPU: Use std library
    std::iota( v.begin(), v.end(), 0L );
#endif
}

/* \brief Find the indices that would reorder the elements of `predicate`
 * so that the elements with non-zero value precede the other elements
 *
 * \param[in] index_begin Point to the beginning of the vector which is
 *            to be filled with these indices
 * \param[in] index_begin Point to the end of the vector which is
 *            to be filled with these indices
 * \param[in] Vector that indicates the elements that need to be reordered first
 */
template< typename ForwardIterator >
ForwardIterator stablePartition(ForwardIterator index_begin,
                               ForwardIterator index_end,
                               amrex::Gpu::DeviceVector<int>& predicate)
{
#ifdef AMREX_USE_GPU
    // On GPU: Use thrust
    int* AMREX_RESTRICT predicate_ptr = predicate.dataPtr();
    ForwardIterator sep = thrust::stable_partition(
        thrust::cuda::par(amrex::Cuda::The_ThrustCachedAllocator()),
        index_begin, index_end,
        [predicate_ptr] AMREX_GPU_DEVICE (long i) { return predicate_ptr[i]; }
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

/* \brief Return the number of elements between `first` and `last`
 *
 * \param[in] fist Points to a position in a vector
 * \param[in] last Points to another position in a vector
 * \return The number of elements between `first` and `last`
 */
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

/* \brief Functor that fills the elements of the particle array `inexflag`
 *  with the value of the spatial array `bmasks`, at the corresponding particle position.
 *
 * This is done only for the elements from `start_index` to the end of `inexflag`
 *
 * \param[in] pti Contains information on the particle positions
 * \param[in] bmasks Spatial array, that contains a flag indicating
 *         whether each cell is part of the gathering/deposition buffers
 * \param[out] inexflag Vector to be filled with the value of `bmasks`
 * \param[in] geom Geometry object, necessary to locate particles within the array `bmasks`
 * \param[in] start_index Index that which elements start to be modified
 */
class fillBufferFlag
{
    public:
        fillBufferFlag( WarpXParIter& pti, const amrex::iMultiFab* bmasks,
                        amrex::Gpu::DeviceVector<int>& inexflag,
                        const amrex::Geometry& geom, long start_index=0 ) {

            // Extract simple structure that can be used directly on the GPU
            m_particles = &(pti.GetArrayOfStructs()[0]);
            m_buffer_mask = (*bmasks)[pti].array();
            m_inexflag_ptr = inexflag.dataPtr();
            m_domain_small_end = geom.Domain().smallEnd();
            for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
                m_prob_lo[idim] = geom.ProbLo(idim);
                m_inv_cell_size[idim] = geom.InvCellSize(idim);
            }
            m_start_index = start_index;
        };


        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void operator()( const long i ) const {
            // Select a particle
            auto& p = m_particles[i+m_start_index];
            // Find the index of the cell where this particle is located
            amrex::IntVect iv;
            for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
                iv[idim] = static_cast<int>(floor(
                    (p.m_rdata.pos[idim]-m_prob_lo[idim])*m_inv_cell_size[idim]
                ));
            }
            // Find the value of the buffer flag in this cell and
            // store it at the corresponding particle position in the array `inexflag`
            m_inexflag_ptr[i+m_start_index] = m_buffer_mask(iv);
        };

    private:
        amrex::Real m_prob_lo[AMREX_SPACEDIM];
        amrex::Real m_inv_cell_size[AMREX_SPACEDIM];
        amrex::IntVect m_domain_small_end;
        int* m_inexflag_ptr;
        WarpXParticleContainer::ParticleType* m_particles;
        amrex::Array4<const int> m_buffer_mask;
        long m_start_index;
};

/* \brief Functor that copies the elements of `src` into `dst`,
 *       while reordering them according to `indices`
 *
 * \param src Source vector
 * \param dst Destination vector
 * \param indices Array of indices that indicate how to reorder elements
 */
template <typename T>
class copyAndReorder
{
    public:
        copyAndReorder(
            amrex::Gpu::ManagedDeviceVector<T>& src,
            amrex::Gpu::ManagedDeviceVector<T>& dst,
            amrex::Gpu::DeviceVector<long>& indices ) {
            // Extract simple structure that can be used directly on the GPU
            m_src_ptr = src.dataPtr();
            m_dst_ptr = dst.dataPtr();
            m_indices_ptr = indices.dataPtr();
        };

        AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        void operator()( const long ip ) const {
            m_dst_ptr[ip] = m_src_ptr[ m_indices_ptr[ip] ];
        };

    private:
        T* m_src_ptr;
        T* m_dst_ptr;
        long* m_indices_ptr;
};

#endif // WARPX_PARTICLES_SORTING_SORTINGUTILS_H_
