#ifndef AMREX_PARTICLETRANSFORMATION_H_
#define AMREX_PARTICLETRANSFORMATION_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_ParticleTile.H>

namespace amrex
{

/**
 * \brief A general single particle copying routine that can run on the GPU. 
 *
 * \tparam NSR number of extra reals in the particle struct
 * \tparam NSI number of extra ints in the particle struct
 * \tparam NAR number of reals in the struct-of-arrays
 * \tparam NAI number of ints in the struct-of-arrays
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param src_i the index in the source to read from
 * \param dst_i the index in the destination to write to
 *
 */
template <int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void copyParticle(const      ParticleTileData<NSR, NSI, NAR, NAI>& dst, 
                  const ConstParticleTileData<NSR, NSI, NAR, NAI>& src, 
                  int src_i, int dst_i) noexcept
{
    AMREX_ASSERT(dst.m_num_runtime_real == src.m_num_runtime_real);
    AMREX_ASSERT(dst.m_num_runtime_int  == src.m_num_runtime_int );

    dst.m_aos[dst_i] = src.m_aos[src_i];
    for (int j = 0; j < NAR; ++j)
        dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
    for (int j = 0; j < dst.m_num_runtime_real; ++j)
        dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
    for (int j = 0; j < NAI; ++j)
        dst.m_idata[j][dst_i] = src.m_idata[j][src_i];
    for (int j = 0; j < dst.m_num_runtime_int; ++j)
        dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
}

/**
 * \brief Copy particles from src to dst. This version copies all the 
 * particles, writing them to the beginning of dst.
 * 
 * \tparam PTile the particle tile type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param f the function that will be applied to each particle
 *
 */    
template <typename PTile>
void copyParticles(PTile& dst, const PTile& src) noexcept
{
    auto np = src.numParticles();
    copyParticles(dst, src, 0, 0, np);
}

/**
 * \brief Copy particles from src to dst. This version copies n particles
 * starting at index src_start, writing the result starting at dst_start.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam N the size type, e.g. long
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 *
 */
template <typename PTile, typename IndexType, typename N,
          amrex::EnableIf_t<std::is_integral<IndexType>::value, int> foo = 0>
void copyParticles(PTile& dst, const PTile& src,
                   IndexType src_start, IndexType dst_start, N n) noexcept
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( n, i,
    {
        copyParticle(dst_data, src_data, src_start+i, dst_start+i);
    });
}
    
/**
 * \brief Apply the function f to all the particles in src, writing the
 * result to dst. This version does all the particles in src.
 * 
 * \tparam PTile the particle tile type
 * \tparam F a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param f the function that will be applied to each particle
 *
 */    
template <typename PTile, typename F>
void transformParticles(PTile& dst, const PTile& src, F&& f) noexcept
{
    auto np = src.numParticles();
    transformParticles(dst, src, 0, 0, np, std::forward<F>(f));
}

/**
 * \brief Apply the function f to particles in src, writing the
 * result to dst. This version applies the function to n particles
 * starting at index src_start, writing the result starting at dst_start.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam N the size type, e.g. long
 * \tparam F a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 * \param f the function that will be applied to each particle
 *
 */
template <typename PTile, typename IndexType, typename N, typename F,
          amrex::EnableIf_t<std::is_integral<IndexType>::value, int> foo = 0>
void transformParticles(PTile& dst, const PTile& src,
                        IndexType src_start, IndexType dst_start, N n, F&& f) noexcept
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( n, i,
    {
        f(dst_data, src_data, src_start+i, dst_start+i);
    });
}

/**
 * \brief Conditionally copy particles from src to dst based on the value of mask.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 *
 */
template <typename PTile, typename IndexType,
          amrex::EnableIf_t<std::is_integral<IndexType>::value, int> foo = 0>
IndexType filterParticles(PTile& dst, const PTile& src, const IndexType* mask) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    IndexType last_mask, last_offset;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);    
                               
    auto p_offsets = offsets.dataPtr();
    
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i]) copyParticle(dst_data, src_data, i, p_offsets[i]);
    });

    Gpu::streamSynchronize();
    return last_mask + last_offset;
}

/**
 * \brief Conditionally copy particles from src to dst based on a predicate.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam Pred a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 *
 */    
template <typename PTile, typename Pred>
auto filterParticles(PTile& dst, const PTile& src, Pred&& p) noexcept
    -> decltype(p(typename PTile::ConstParticleTileDataType(), 0))
{
    using IndexType = decltype(p(typename PTile::ConstParticleTileDataType(), 0));
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getConstParticleTileData();
    
    AMREX_HOST_DEVICE_FOR_1D(np, i,
    {
        p_mask[i] = p(src_data, i);
    });

    return filterParticles(dst, src, mask.dataPtr());
}

/**
 * \brief Conditionally copy particles from src to dst based on the value of mask.
 * A transformation will also be applied to the particles on copy.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename PTile, typename IndexType, typename F,
          amrex::EnableIf_t<std::is_integral<IndexType>::value, int> foo = 0>
IndexType filterAndTransformParticles(PTile& dst, const PTile& src, IndexType* mask, F&& f) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    IndexType last_mask, last_offset;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);    

    auto p_offsets = offsets.dataPtr();
    
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i]) f(dst_data, src_data, i, p_offsets[i]);
    });

    Gpu::streamSynchronize();
    return last_mask + last_offset;
}

/**
 * \brief Conditionally copy particles from src to dst based on a predicate.
 * A transformation will also be applied to the particles on copy.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam Pred a function object
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */    
template <typename PTile, typename Pred, typename F>
auto filterAndTransformParticles(PTile& dst, const PTile& src, Pred&& p, F&& f) noexcept
    -> decltype(p(typename PTile::ConstParticleTileDataType(), 0))
{
    using IndexType = decltype(p(typename PTile::ConstParticleTileDataType(), 0));
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getConstParticleTileData();
    
    AMREX_HOST_DEVICE_FOR_1D(np, i,
    {
        p_mask[i] = p(src_data, i);
    });

    return filterAndTransformParticles(dst, src, mask.dataPtr(), std::forward<F>(f));
}
    
/**
 * \brief Gather particles copies particles into contiguous order from an 
 * arbitrary order. Specifically, the particle at the index inds[i] in src
 * will be copied to the index i in dst.
 * 
 * \tparam PTile the particle tile type
 * \tparam N the size type, e.g. long
 * \tparam IndexType the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param np the number of particles
 * \param inds pointer to the permutation array
 *
 */
template <typename PTile, typename N, typename IndexType,
          amrex::EnableIf_t<std::is_integral<IndexType>::value, int> foo = 0>
void gatherParticles (PTile& dst, const PTile& src, N np, const IndexType* inds)
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();
          
    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        copyParticle(dst_data, src_data, inds[i], i);
    });
}

/**
 * \brief Scatter particles copies particles from contiguous order into an
 * arbitrary order. Specifically, the particle at the index i in src
 * will be copied to the index inds[i] in dst.
 * 
 * \tparam PTile the particle tile type
 * \tparam N the size type, e.g. long
 * \tparam IndexType the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param np the number of particles
 * \param inds pointer to the permutation array
 *
 */
template <typename PTile, typename N, typename IndexType,
          amrex::EnableIf_t<std::is_integral<IndexType>::value, int> foo = 0>
void scatterParticles (PTile& dst, const PTile& src, N np, const IndexType* inds)
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();
          
    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        copyParticle(dst_data, src_data, i, inds[i]);
    });
}

}

#endif // include guard
