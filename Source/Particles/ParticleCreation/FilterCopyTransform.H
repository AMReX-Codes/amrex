#ifndef FILTER_COPY_TRANSFORM_H_
#define FILTER_COPY_TRANSFORM_H_

#include <AMReX_GpuContainers.H>
#include <AMReX_TypeTraits.H>

/**
 * \brief Apply a filter, copy, and transform operation to the particles
 * in src, in that order, writing the result to dst, starting at dst_index.
 * The dst tile will be extended so all the particles will fit, if needed.
 *
 * Note that the transform function operates on both the src and the dst,
 * so both can be modified.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam TransFunc the transform function type
 * \tparam CopyFunc the copy function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param dst_index the location at which to starting writing the result to dst
 * \param copy defines the function that will be applied for the "copy" step
 * \param transform defines that transformation.
 *
 */
template <typename DstTile, typename SrcTile, typename Index,
          typename TransFunc, typename CopyFunc,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index filterCopyTransformParticles (DstTile& dst, SrcTile& src, Index* mask, Index dst_index,
                                    CopyFunc&& copy, TransFunc&& transform) noexcept
{
    using namespace amrex;

    auto np = src.numParticles();
    if (np == 0) return 0;

    Gpu::DeviceVector<Index> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    Index last_mask, last_offset;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);

    Index num_added = last_mask + last_offset;
    dst.resize(std::max(dst_index + num_added, dst.numParticles()));
    
    auto p_offsets = offsets.dataPtr();
    
    const auto src_data = src.getParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i])
        {
            copy(dst_data, src_data, i, p_offsets[i] + dst_index);
            transform(dst_data, src_data, i, p_offsets[i] + dst_index);
        }
    });

    Gpu::streamSynchronize();
    return num_added;
}

/**
 * \brief Apply a filter, copy, and transform operation to the particles
 * in src, in that order, writing the result to dst, starting at dst_index.
 * The dst tile will be extended so all the particles will fit, if needed.
 *
 * Note that the transform function operates on both the src and the dst,
 * so both can be modified.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam Filter the filter function type
 * \tparam TransFunc the transform function type
 * \tparam CopyFunc the copy function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param dst_index the location at which to starting writing the result to dst
 * \param filter a function returning true if that particle is to be copied and transformed
 * \param copy defines the function that will be applied for the "copy" step
 * \param transform defines that transformation.
 *
 */
template <typename DstTile, typename SrcTile, typename Index,
          typename PredFunc, typename TransFunc, typename CopyFunc>
Index filterCopyTransformParticles (DstTile& dst, SrcTile& src, Index dst_index,
                                    PredFunc&& filter, CopyFunc&& copy, TransFunc&& transform) noexcept
{
    using namespace amrex;

    auto np = src.numParticles();
    if (np == 0) return 0;

    Gpu::DeviceVector<Index> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D(np, i,
    {
        p_mask[i] = filter(src_data, i);
    });

    return filterCopyTransformParticles(dst, src, mask.dataPtr(), dst_index,
                                        std::forward<CopyFunc>(copy),
                                        std::forward<TransFunc>(transform));
}

#endif
