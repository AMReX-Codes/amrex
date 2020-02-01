#ifndef SMART_COPY_H_
#define SMART_COPY_H_

#include <AMReX_GpuContainers.H>

#include <map>
#include <string>

using NameMap = std::map<std::string, int>;

struct SmartCopyTag
{
    std::vector<std::string> common_names;
    amrex::Gpu::DeviceVector<int> src_comps;
    amrex::Gpu::DeviceVector<int> dst_comps;

    int size () const noexcept { return common_names.size(); }
};

SmartCopyTag getSmartCopyTag (const NameMap& src, const NameMap& dst);

struct SmartCopy
{
    int  m_num_real;
    const int* m_src_comps_r;
    const int* m_dst_comps_r;

    int  m_num_int;
    const int* m_src_comps_i;
    const int* m_dst_comps_i;

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void operator() (const DstData& dst, const SrcData& src, int i_src, int i_dst) const noexcept
    {
        dst.m_aos[i_dst] = src.m_aos[i_src];

        for (int j = 0; j < m_num_real; ++j)
        {
            int src_comp = (m_src_comps_r[j] < SrcData::NAR) ? m_src_comps_r[j] : m_src_comps_r[j] - SrcData::NAR;
            int dst_comp = (m_dst_comps_r[j] < DstData::NAR) ? m_dst_comps_r[j] : m_dst_comps_r[j] - DstData::NAR;

            amrex::ParticleReal* dst_data = (m_dst_comps_r[j] < DstData::NAR) ? dst.m_rdata : dst.m_runtime_rdata;
            amrex::ParticleReal* src_data = (m_src_comps_r[j] < SrcData::NAR) ? src.m_rdata : src.m_runtime_rdata;

            dst_data[dst_comp] = src_data[src_comp];            
        }

        for (int j = 0; j < m_num_int; ++j)
        {
            int src_comp = (m_src_comps_i[j] < SrcData::NAI) ? m_src_comps_i[j] : m_src_comps_i[j] - SrcData::NAI;
            int dst_comp = (m_dst_comps_i[j] < DstData::NAI) ? m_dst_comps_i[j] : m_dst_comps_i[j] - DstData::NAI;

            int* dst_data = (m_dst_comps_i[j] < DstData::NAI) ? dst.m_idata : dst.m_runtime_idata;
            int* src_data = (m_src_comps_i[j] < SrcData::NAI) ? src.m_idata : src.m_runtime_idata;

            dst_data[dst_comp] = src_data[src_comp];
        }
    }
};

class SmartCopyFactory
{
    SmartCopyTag m_tag_real;
    SmartCopyTag m_tag_int;
    bool m_defined;

    template <class PC1, class PC2>
    SmartCopyFactory (const PC1& pc1, const PC2& pc2) noexcept
        : m_defined(false)
    {
        m_tag_real = getSmartCopyTag(pc1.getParticleComps(), pc2.getParticleComps());
        m_tag_int  = getSmartCopyTag(pc1.getParticleiComps(), pc2.getParticleiComps());
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
                         m_tag_int. dst_comps.dataPtr()};
    }

    bool isDefined () const noexcept { return m_defined; }
};

#endif
