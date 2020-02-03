#ifndef SMART_COPY_H_
#define SMART_COPY_H_

#include <DefaultInitialization.H>

#include <AMReX_GpuContainers.H>
#include <AMReX_TypeTraits.H>

#include <map>
#include <string>

using NameMap = std::map<std::string, int>;
using PolicyVec = amrex::Gpu::DeviceVector<InitializationPolicy>;


struct SmartCopyTag
{
    std::vector<std::string> common_names;
    amrex::Gpu::DeviceVector<int> src_comps;
    amrex::Gpu::DeviceVector<int> dst_comps;

    int size () const noexcept { return common_names.size(); }
};

PolicyVec getPolicies (const NameMap& names) noexcept;
SmartCopyTag getSmartCopyTag (const NameMap& src, const NameMap& dst) noexcept;

struct SmartCopy
{
    int m_num_copy_real;
    const int* m_src_comps_r;
    const int* m_dst_comps_r;

    int m_num_copy_int;
    const int* m_src_comps_i;
    const int* m_dst_comps_i;

    const InitializationPolicy* m_policy_real;
    const InitializationPolicy* m_policy_int;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (DstData& dst, const SrcData& src, int i_src, int i_dst) const noexcept
    {
        // the particle struct is always copied over
        dst.m_aos[i_dst] = src.m_aos[i_src];

        // initialize the real components
        for (int j = 0; j < DstData::NAR; ++j)
            dst.m_rdata[j][i_dst] = initializeRealValue(m_policy_real[j]);
        for (int j = 0; j < dst.m_num_runtime_real; ++j)
            dst.m_rdata[j][i_dst] = initializeRealValue(m_policy_real[j+DstData::NAR]);

        // initialize the int components
        for (int j = 0; j < DstData::NAI; ++j)
            dst.m_idata[j][i_dst] = initializeIntValue(m_policy_int[j]);
        for (int j = 0; j < dst.m_num_runtime_int; ++j)
            dst.m_idata[j][i_dst] = initializeIntValue(m_policy_int[j+DstData::NAI]);

        // copy the shared real components
        for (int j = 0; j < m_num_copy_real; ++j)
        {
            int src_comp = (m_src_comps_r[j] < SrcData::NAR) ? m_src_comps_r[j] : m_src_comps_r[j] - SrcData::NAR;
            int dst_comp = (m_dst_comps_r[j] < DstData::NAR) ? m_dst_comps_r[j] : m_dst_comps_r[j] - DstData::NAR;

            amrex::ParticleReal* AMREX_RESTRICT dst_data = (m_dst_comps_r[j] < DstData::NAR) ? dst.m_rdata[dst_comp] : dst.m_runtime_rdata[dst_comp];
            const amrex::ParticleReal* AMREX_RESTRICT src_data = (m_src_comps_r[j] < SrcData::NAR) ? src.m_rdata[src_comp] : src.m_runtime_rdata[src_comp];

            dst_data[i_dst] = src_data[i_src];
        }

        // copy the shared int components
        for (int j = 0; j < m_num_copy_int; ++j)
        {
            int src_comp = (m_src_comps_i[j] < SrcData::NAI) ? m_src_comps_i[j] : m_src_comps_i[j] - SrcData::NAI;
            int dst_comp = (m_dst_comps_i[j] < DstData::NAI) ? m_dst_comps_i[j] : m_dst_comps_i[j] - DstData::NAI;

            int* AMREX_RESTRICT dst_data = (m_dst_comps_i[j] < DstData::NAI) ? dst.m_idata[dst_comp] : dst.m_runtime_idata[dst_comp];
            const int* AMREX_RESTRICT src_data = (m_src_comps_i[j] < SrcData::NAI) ? src.m_idata[src_comp] : src.m_runtime_idata[src_comp];

            dst_data[i_dst] = src_data[i_src];
        }
    }
};

class SmartCopyFactory
{
    SmartCopyTag m_tag_real;
    SmartCopyTag m_tag_int;
    PolicyVec m_policy_real;
    PolicyVec m_policy_int;
    bool m_defined;

public:
    template <class SrcPC, class DstPC>
    SmartCopyFactory (const SrcPC& src, const DstPC& dst) noexcept
        : m_defined(false)
    {
        m_tag_real = getSmartCopyTag(src.getParticleComps(),  dst.getParticleComps());
        m_tag_int  = getSmartCopyTag(src.getParticleiComps(), dst.getParticleiComps());

        m_policy_real = getPolicies(dst.getParticleComps());
        m_policy_int  = getPolicies(dst.getParticleiComps());

        m_defined = true;
    }

    SmartCopy getSmartCopy () const noexcept
    {
        AMREX_ASSERT(m_defined);
        return SmartCopy{m_tag_real.size(),
                         m_tag_real.src_comps.dataPtr(),
                         m_tag_real.dst_comps.dataPtr(),
                         m_tag_int.size(),
                         m_tag_int. src_comps.dataPtr(),
                         m_tag_int. dst_comps.dataPtr(),
                         m_policy_real.dataPtr(),
                         m_policy_int.dataPtr()};
    }

    bool isDefined () const noexcept { return m_defined; }
};

/**
 * \brief Conditionally copy particles from src to dst based on the value of mask.
 * A transformation will also be applied to the particles on copy.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename DstTile, typename SrcTile, typename Index, typename TransFunc, typename CopyFunc,
          amrex::EnableIf_t<std::is_integral<Index>::value, int> foo = 0>
Index filterAndTransformParticles (DstTile& dst, SrcTile& src, Index* mask,
                                   TransFunc&& t, CopyFunc&& f) noexcept
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
            f(dst_data, src_data, i, p_offsets[i]);
            t(dst_data, src_data, i, p_offsets[i]);
        }
    });

    Gpu::streamSynchronize();
    return last_mask + last_offset;
}

/**
 * \brief Conditionally copy particles from src to dst based on a predicate.
 * A transformation will also be applied to the particles on copy.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Pred a function object
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename DstTile, typename SrcTile, typename PredFunc, typename TransFunc, typename CopyFunc>
int filterAndTransformParticles (DstTile& dst, SrcTile& src,
                                 PredFunc&& p, TransFunc&& t, CopyFunc&& f) noexcept
{
    auto np = src.numParticles();
    amrex::Gpu::DeviceVector<int> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D(np, i,
    {
        p_mask[i] = p(src_data, i);
    });

    return filterAndTransformParticles(dst, src, mask.dataPtr(),
                                       std::forward<TransFunc>(t),
                                       std::forward<CopyFunc>(f));
}


#endif
