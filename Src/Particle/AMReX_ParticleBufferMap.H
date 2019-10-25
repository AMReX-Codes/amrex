#ifndef AMREX_PARTICLEBUFFERMAP_H_
#define AMREX_PARTICLEBUFFERMAP_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_ParGDB.H>

namespace amrex {

class ParticleBufferMap {
    
    bool m_defined;
    Vector<BoxArray> m_ba;
    Vector<DistributionMapping> m_dm;

    Vector<int> m_bucket_to_gid;
    Vector<int> m_bucket_to_lev;
    Vector<int> m_bucket_to_pid;

    Vector<int> m_lev_gid_to_bucket;
    Vector<int> m_lev_offsets;

    Vector<int> m_proc_box_counts;
    Vector<int> m_proc_box_offsets;

    Gpu::DeviceVector<int> d_bucket_to_pid;
    Gpu::DeviceVector<int> d_lev_gid_to_bucket;    
    Gpu::DeviceVector<int> d_lev_offsets;
    
public:
    ParticleBufferMap ()
        : m_defined(false), m_ba(), m_dm()
        {}

    ParticleBufferMap (const ParGDBBase* a_gdb);
    
    void define (const ParGDBBase* a_gdb);

    bool isValid (const ParGDBBase* a_gdb) const;

    AMREX_INLINE
    int numLevels () const
    {
        AMREX_ASSERT(m_defined);
        return m_lev_offsets.size()-1;
    }

    AMREX_INLINE
    int numBuckets () const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_gid.size();
    }
        
    AMREX_INLINE
    int bucketToGrid (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_gid[bid];
    }

    AMREX_INLINE
    int bucketToLevel (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_lev[bid];
    }

    AMREX_INLINE
    int bucketToProc (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_pid[bid];
    }

    AMREX_INLINE
    int gridAndLevToBucket (int gid, int lev) const
    {
        AMREX_ASSERT(m_defined);        
        return m_lev_gid_to_bucket[m_lev_offsets[lev] + gid];
    }

    AMREX_INLINE
    const int* levGridToBucketPtr () const
    {
        AMREX_ASSERT(m_defined);
        return d_lev_gid_to_bucket.dataPtr();
    }

    AMREX_INLINE
    const int* levelOffsetsPtr () const
    {
        AMREX_ASSERT(m_defined);
        return d_lev_offsets.dataPtr();
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
    int procID(int gid, int lev) const 
    { 
        AMREX_ASSERT(m_defined);
        return m_dm[lev][gid]; 
    }

    AMREX_INLINE
    const int* bucketToPIDPtr () const
    {
        AMREX_ASSERT(m_defined);
        return d_bucket_to_pid.dataPtr();
    }
};

} // namespace amrex

#endif // AMREX_PARTICLEBUFFERMAP_H_
