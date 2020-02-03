#ifndef SMART_COPY_H_
#define SMART_COPY_H_

#include <DefaultInitialization.H>

#include <AMReX_GpuContainers.H>

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
    void operator() (const DstData& dst, const SrcData& src, int i_src, int i_dst) const noexcept
    {
		// the particle struct is always copied over
        dst.m_aos[i_dst] = src.m_aos[i_src];

		// initialize the real components
		for (int j = 0; j < DstData::NAR; ++j)
			dst.m_rdata[j] = initializeRealValue(m_policy_real[j]);
		for (int j = 0; j < dst.m_num_runtime_real; ++j)
			dst.m_rdata[j] = initializeRealValue(m_policy_real[j+DstData::NAR]);

		// initialize the int components
		for (int j = 0; j < DstData::NAI; ++j)
			dst.m_idata[j] = initializeIntValue(m_policy_int[j]);
		for (int j = 0; j < dst.m_num_runtime_int; ++j)
			dst.m_idata[j] = initializeIntValue(m_policy_int[j+DstData::NAI]);
		
		// copy the shared real components
        for (int j = 0; j < m_num_copy_real; ++j)
        {
            int src_comp = (m_src_comps_r[j] < SrcData::NAR) ? m_src_comps_r[j] : m_src_comps_r[j] - SrcData::NAR;
            int dst_comp = (m_dst_comps_r[j] < DstData::NAR) ? m_dst_comps_r[j] : m_dst_comps_r[j] - DstData::NAR;

            amrex::ParticleReal* dst_data = (m_dst_comps_r[j] < DstData::NAR) ? dst.m_rdata : dst.m_runtime_rdata;
            amrex::ParticleReal* src_data = (m_src_comps_r[j] < SrcData::NAR) ? src.m_rdata : src.m_runtime_rdata;

            dst_data[dst_comp] = src_data[src_comp];
        }

		// copy the shared int components
        for (int j = 0; j < m_num_copy_int; ++j)
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

#endif
