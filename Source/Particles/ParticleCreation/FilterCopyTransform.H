#ifndef FILTER_COPY_TRANSFORM_H_
#define FILTER_COPY_TRANSFORM_H_

#include <AMReX_GpuContainers.H>
#include <AMReX_TypeTraits.H>

template <typename DstTile, typename SrcTile, typename Index, typename TransFunc, typename CopyFunc,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index filterCopyTransformParticles (DstTile& dst, SrcTile& src, Index* mask, int dst_index,
                                    CopyFunc&& f, TransFunc&& t) noexcept
{
    using namespace amrex;
    auto np = src.numParticles();
    Gpu::DeviceVector<Index> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    Index last_mask, last_offset;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);

    auto p_offsets = offsets.dataPtr();

    const auto src_data = src.getParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i])
        {
            f(dst_data, src_data, i, p_offsets[i] + dst_index);
            t(dst_data, src_data, i, p_offsets[i] + dst_index);
        }
    });

    Gpu::streamSynchronize();
    return last_mask + last_offset;
}

template <typename DstTile, typename SrcTile, typename PredFunc, typename TransFunc, typename CopyFunc>
int filterCopyTransformParticles (DstTile& dst, SrcTile& src, int dst_index,
                                  PredFunc&& p, CopyFunc&& f, TransFunc&& t) noexcept
{
    using namespace amrex;
    auto np = src.numParticles();
    Gpu::DeviceVector<int> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D(np, i,
    {
        p_mask[i] = p(src_data, i);
    });

    return filterCopyTransformParticles(dst, src, mask.dataPtr(), dst_index,
                                        std::forward<TransFunc>(t),
                                        std::forward<CopyFunc>(f));
}

#endif
