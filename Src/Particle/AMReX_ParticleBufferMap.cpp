#include <AMReX_ParticleBufferMap.H>

using namespace amrex;

ParticleBufferMap::ParticleBufferMap (const BoxArray& ba, const DistributionMapping& dm)
{
    define(ba, dm);
}

void ParticleBufferMap::define (const BoxArray& ba, const DistributionMapping& dm)
{
    BL_PROFILE("ParticleBufferMap::define");

    m_defined = true;
    m_ba = ba;
    m_dm = dm;
    int num_boxes = ba.size();        
    m_bucket_to_gid.resize(0);
    m_bucket_to_gid.resize(num_boxes);
    m_gid_to_bucket.resize(0);
    m_gid_to_bucket.resize(num_boxes);

    using int2 = std::pair<int, int>;
    std::vector<int2> box_proc_ids;
    for (int i = 0; i < num_boxes; ++i) 
        box_proc_ids.push_back(std::make_pair(i, dm[i]));
    
    std::sort(box_proc_ids.begin(), box_proc_ids.end(), 
              [](const int2& a, const int2& b) -> bool
              {
                  return a.second < b.second;
              });

    for (int i = 0; i < num_boxes; ++i) 
        m_bucket_to_gid[i] = box_proc_ids[i].first;
    
    for (int i = 0; i < num_boxes; ++i)
        m_gid_to_bucket[m_bucket_to_gid[i]] = i;

    m_proc_box_counts.resize(0);
    m_proc_box_counts.resize(ParallelDescriptor::NProcs(), 0);    
    for (auto& val : dm.ProcessorMap() )
    {
        m_proc_box_counts[val]++;
    } 

    m_proc_box_offsets.resize(0);
    m_proc_box_offsets.push_back(0);
    for (auto count : m_proc_box_counts)
        m_proc_box_offsets.push_back(m_proc_box_offsets.back() + count);

    d_bucket_to_gid.resize(0);
    d_bucket_to_gid.resize(num_boxes);
    d_gid_to_bucket.resize(0);
    d_gid_to_bucket.resize(num_boxes);    

    Gpu::thrust_copy(m_bucket_to_gid.begin(), m_bucket_to_gid.end(), d_bucket_to_gid.begin());
    Gpu::thrust_copy(m_gid_to_bucket.begin(), m_gid_to_bucket.end(), d_gid_to_bucket.begin());
}

bool ParticleBufferMap::isValid (const BoxArray& ba, const DistributionMapping& dm) const
{
    if (!m_defined) return false;
    bool same_ba = BoxArray::SameRefs(ba, m_ba);
    bool same_dm = DistributionMapping::SameRefs(dm, m_dm);
    return same_ba && same_dm;
}
