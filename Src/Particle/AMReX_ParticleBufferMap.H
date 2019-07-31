#ifndef AMREX_PARTICLEBUFFERMAP_H_
#define AMREX_PARTICLEBUFFERMAP_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_CudaContainers.H>

namespace amrex {

class ParticleBufferMap {
    
    bool m_defined;
    BoxArray m_ba;
    DistributionMapping m_dm;
    
    Vector<int> m_bucket_to_gid;
    Vector<int> m_gid_to_bucket;

    Vector<int> m_proc_box_counts;
    Vector<int> m_proc_box_offsets;

    Gpu::DeviceVector<int> d_bucket_to_gid;
    Gpu::DeviceVector<int> d_gid_to_bucket;    
    
public:
    ParticleBufferMap ()
        : m_defined(false), m_ba(), m_dm()
        {}

    ParticleBufferMap (const BoxArray& ba, const DistributionMapping& dm);

    void define (const BoxArray& ba, const DistributionMapping& dm);

    bool isValid (const BoxArray& ba, const DistributionMapping& dm) const;
    
    AMREX_INLINE
    int bucketToGrid (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_gid[bid];
    }

    AMREX_INLINE
    int gridToBucket (int gid) const
    {
        AMREX_ASSERT(m_defined);
        return m_gid_to_bucket[gid];
    }

    AMREX_INLINE
    const int* bucketToGridPtr () const
    {
        AMREX_ASSERT(m_defined);       
        return d_bucket_to_gid.dataPtr();
    }

    AMREX_INLINE
    const int* gridToBucketPtr () const
    {
        AMREX_ASSERT(m_defined);
        return d_gid_to_bucket.dataPtr();
    }

    AMREX_INLINE
    int firstBucketOnProc (int pid) const
    {
        AMREX_ASSERT(m_defined);
        return m_proc_box_offsets[pid];
    }

    AMREX_INLINE
    int numBoxesOnProc(int pid) const
    {
        AMREX_ASSERT(m_defined);
        return m_proc_box_counts[pid];
    }

    AMREX_INLINE
    Vector<int> allBucketsOnProc(int pid) const
    {
        AMREX_ASSERT(m_defined);
        Vector<int> buckets;        
        for (int i = 0; i < numBoxesOnProc(pid); ++i)
        {
            buckets.push_back(i+firstBucketOnProc(pid));
        }
        return buckets;
    }

    AMREX_INLINE
    int procID(int gid) const 
    { 
        AMREX_ASSERT(m_defined);
        return m_dm[gid]; 
    }
};

} // namespace amrex

#endif // AMREX_PARTICLEBUFFERMAP_H_
