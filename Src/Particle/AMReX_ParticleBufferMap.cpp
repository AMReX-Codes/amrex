#include <AMReX_ParticleBufferMap.H>

using namespace amrex;

ParticleBufferMap::ParticleBufferMap (const BoxArray& ba, const DistributionMapping& dm)
{
    define(ba, dm);
}

void ParticleBufferMap::define (const BoxArray& ba, const DistributionMapping& dm)
{
    m_defined = true;
    m_ba = ba;
    m_dm = dm;
    int num_boxes = ba.size();        
    m_bucket_to_gid.resize(0);
    m_bucket_to_gid.resize(num_boxes);
    m_gid_to_bucket.resize(0);
    m_gid_to_bucket.resize(num_boxes);
    Gpu::DeviceVector<int> box_proc_ids(num_boxes);
    auto& pmap = dm.ProcessorMap();
    Gpu::thrust_copy(pmap.begin(), pmap.end(), box_proc_ids.begin());
    thrust::sequence(thrust::device, m_bucket_to_gid.begin(), m_bucket_to_gid.end());	
    thrust::sort_by_key(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                        box_proc_ids.begin(),
                        box_proc_ids.end(),
                        m_bucket_to_gid.begin());
    
    auto b_to_gid_ptr = m_bucket_to_gid.dataPtr();
    auto gid_to_b_ptr = m_gid_to_bucket.dataPtr();
    AMREX_FOR_1D (num_boxes, i,
    {
        gid_to_b_ptr[b_to_gid_ptr[i]] = i;
    });	 

    m_proc_box_counts.resize(0);
    m_proc_box_counts.resize(ParallelDescriptor::NProcs(), 0);
    
    m_proc_box_offsets.resize(0);
    m_proc_box_offsets.resize(ParallelDescriptor::NProcs());
    for (auto& val : dm.ProcessorMap() )
    {
        m_proc_box_counts[val]++;
    } 

    amrex::Gpu::exclusive_scan(m_proc_box_counts.begin(), m_proc_box_counts.end(),
                               m_proc_box_offsets.begin());
 
    m_proc_box_offsets.push_back(num_boxes);        
}

bool ParticleBufferMap::isValid (const BoxArray& ba, const DistributionMapping& dm) const
{
    if (!m_defined) return false;
    bool same_ba = BoxArray::SameRefs(ba, m_ba);
    bool same_dm = DistributionMapping::SameRefs(dm, m_dm);
    return same_ba && same_dm;
}
