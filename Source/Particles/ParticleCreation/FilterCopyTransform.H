/* Copyright 2019-2020 Andrew Myers, Axel Huebl,
 * Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
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
 * \tparam N number of particles created in the dst(s) for each filtered src particle
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
 * \param copy callable that defines what will be done for the "copy" step.
 * \param transform callable that defines the transformation to apply on dst and src.
 *
 *        The form of the callable should model:
 *            template <typename DstData, typename SrcData>
 *            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 *            void operator() (DstData& dst, SrcData& src, int i_src, int i_dst);
 *
 *        where dst and src refer to the destination and source tiles and
 *        i_src and i_dst and the particle indices in each tile.
 *
 * \return num_added the number of particles that were written to dst.
 */
template <int N, typename DstTile, typename SrcTile, typename Index,
          typename TransFunc, typename CopyFunc,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index filterCopyTransformParticles (DstTile& dst, SrcTile& src, Index* mask, Index dst_index,
                                    CopyFunc&& copy, TransFunc&& transform) noexcept
{
    using namespace amrex;

    const auto np = src.numParticles();
    if (np == 0) return 0;

    Gpu::DeviceVector<Index> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    Index last_mask, last_offset;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);

    Gpu::streamSynchronize();
    const Index num_added = N * (last_mask + last_offset);
    dst.resize(std::max(dst_index + num_added, dst.numParticles()));

    const auto p_offsets = offsets.dataPtr();

    const auto src_data = src.getParticleTileData();
    const auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i])
        {
            for (int j = 0; j < N; ++j)
                copy(dst_data, src_data, i, N*p_offsets[i] + dst_index + j);
            transform(dst_data, src_data, i, N*p_offsets[i] + dst_index);
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
 * \tparam N number of particles created in the dst(s) for each filtered src particle
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
 * \param filter a callable returning true if that particle is to be copied and transformed
 * \param copy callable that defines what will be done for the "copy" step.
 * \param transform callable that defines the transformation to apply on dst and src.
 *
 *        The form of the callable should model:
 *            template <typename DstData, typename SrcData>
 *            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 *            void operator() (DstData& dst, SrcData& src, int i_src, int i_dst);
 *
 *        where dst and src refer to the destination and source tiles and
 *        i_src and i_dst and the particle indices in each tile.
 *
 * \return num_added the number of particles that were written to dst.
 */
template <int N, typename DstTile, typename SrcTile, typename Index,
          typename PredFunc, typename TransFunc, typename CopyFunc>
Index filterCopyTransformParticles (DstTile& dst, SrcTile& src, Index dst_index,
                                    PredFunc&& filter, CopyFunc&& copy, TransFunc&& transform) noexcept
{
    using namespace amrex;

    const auto np = src.numParticles();
    if (np == 0) return 0;

    Gpu::DeviceVector<Index> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D(np, i,
    {
        p_mask[i] = filter(src_data, i);
    });

    return filterCopyTransformParticles<N>(dst, src, mask.dataPtr(), dst_index,
                                                      std::forward<CopyFunc>(copy),
                                                      std::forward<TransFunc>(transform));
}

/**
 * \brief Apply a filter, copy, and transform operation to the particles
 * in src, in that order, writing the results to dst1 and dst2, starting
 * at dst1_index and dst2_index. The dst tiles will be extended so all the
 * particles will fit, if needed.
 *
 * Note that the transform function operates on all of src, dst1, and dst2,
 * so all of them can be modified.
 *
 * \tparam N number of particles created in the dst(s) for each filtered src particle
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam TransFunc the transform function type
 * \tparam CopyFunc1 the copy function type for src-->dst1
 * \tparam CopyFunc2 the copy function type for src-->dst2
 *
 * \param dst1 the destination tile
 * \param dst2 the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param dst1_index the location at which to starting writing the result to dst1
 * \param dst2_index the location at which to starting writing the result to dst2
 * \param copy1 callable that defines what will be done for the "copy" step for src-->dst1.
 * \param copy2 callable that defines what will be done for the "copy" step for src-->dst2.
 * \param transform callable that defines the transformation to apply on dst and src.
 *
 *        The form of the callable should model:
 *            template <typename DstData, typename SrcData>
 *            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 *            void operator() (DstData& dst, SrcData& src, int i_src, int i_dst);
 *
 *        where dst and src refer to the destination and source tiles and
 *        i_src and i_dst and the particle indices in each tile.
 *
 * \return num_added the number of particles that were written to dst.
 */
template <int N, typename DstTile, typename SrcTile, typename Index,
          typename TransFunc, typename CopyFunc1, typename CopyFunc2,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index filterCopyTransformParticles (DstTile& dst1, DstTile& dst2, SrcTile& src, Index* mask,
                                    Index dst1_index, Index dst2_index,
                                    CopyFunc1&& copy1, CopyFunc2&& copy2,
                                    TransFunc&& transform) noexcept
{
    using namespace amrex;

    auto np = src.numParticles();
    if (np == 0) return 0;

    Gpu::DeviceVector<Index> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    Index last_mask, last_offset;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);

    Gpu::streamSynchronize();
    const Index num_added = N*(last_mask + last_offset);
    dst1.resize(std::max(dst1_index + num_added, dst1.numParticles()));
    dst2.resize(std::max(dst2_index + num_added, dst2.numParticles()));

    auto p_offsets = offsets.dataPtr();

    const auto src_data  =  src.getParticleTileData();
    const auto dst1_data = dst1.getParticleTileData();
    const auto dst2_data = dst2.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i])
        {
            for (int j = 0; j < N; ++j)
            {
                copy1(dst1_data, src_data, i, N*p_offsets[i] + dst1_index + j);
                copy2(dst2_data, src_data, i, N*p_offsets[i] + dst2_index + j);
            }
            transform(dst1_data, dst2_data, src_data, i,
                      N*p_offsets[i] + dst1_index,
                      N*p_offsets[i] + dst2_index);
        }
    });

    Gpu::streamSynchronize();
    return num_added;
}

/**
 * \brief Apply a filter, copy, and transform operation to the particles
 * in src, in that order, writing the results to dst1 and dst2, starting
 * at dst1_index and dst2_index. The dst tiles will be extended so all the
 * particles will fit, if needed.
 *
 * Note that the transform function operates on all of src, dst1, and dst2,
 * so all of them can be modified.
 *
 * \tparam N number of particles created in the dst(s) for each filtered src particle
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam Filter the filter function type
 * \tparam TransFunc the transform function type
 * \tparam CopyFunc1 the copy function type for src-->dst1
 * \tparam CopyFunc2 the copy function type for src-->dst2
 *
 * \param dst1 the destination tile
 * \param dst2 the destination tile
 * \param src the source tile
 * \param dst1_index the location at which to starting writing the result to dst1
 * \param dst2_index the location at which to starting writing the result to dst2
 * \param filter a callable returning true if that particle is to be copied and transformed
 * \param copy1 callable that defines what will be done for the "copy" step for src-->dst1.
 * \param copy2 callable that defines what will be done for the "copy" step for src-->dst2.
 * \param transform callable that defines the transformation to apply on dst and src.
 *
 *        The form of the callable should model:
 *            template <typename DstData, typename SrcData>
 *            AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
 *            void operator() (DstData& dst, SrcData& src, int i_src, int i_dst);
 *
 *        where dst and src refer to the destination and source tiles and
 *        i_src and i_dst and the particle indices in each tile.
 *
 * \return num_added the number of particles that were written to dst.
 */
template <int N, typename DstTile, typename SrcTile, typename Index,
          typename PredFunc, typename TransFunc, typename CopyFunc1, typename CopyFunc2>
Index filterCopyTransformParticles (DstTile& dst1, DstTile& dst2, SrcTile& src,
                                    Index dst1_index, Index dst2_index,
                                    PredFunc&& filter, CopyFunc1&& copy1, CopyFunc2&& copy2,
                                    TransFunc&& transform) noexcept
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

    return filterCopyTransformParticles<N>(dst1, dst2, src, mask.dataPtr(),
                                        dst1_index, dst2_index,
                                        std::forward<CopyFunc1>(copy1),
                                        std::forward<CopyFunc2>(copy2),
                                        std::forward<TransFunc>(transform));
}

#endif
