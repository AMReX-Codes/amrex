#ifndef AMREX_PARTICLEBUFFERMAP_H_
#define AMREX_PARTICLEBUFFERMAP_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParGDB.H>

namespace amrex {

struct GetPID
{
    const int* m_bucket_to_pid;
    const int* m_lev_gid_to_bucket;
    const int* m_lev_offsets;

    GetPID (const Gpu::DeviceVector<int>& bucket_to_pid,
            const Gpu::DeviceVector<int>& lev_gid_to_bucket,
            const Gpu::DeviceVector<int>& lev_offsets)
        : m_bucket_to_pid(bucket_to_pid.dataPtr()),
          m_lev_gid_to_bucket(lev_gid_to_bucket.dataPtr()),
          m_lev_offsets(lev_offsets.dataPtr())
        {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const int lev, const int gid) const noexcept
    {
        return m_bucket_to_pid[m_lev_gid_to_bucket[m_lev_offsets[lev]+gid]];
    }
};

struct GetBucket
{
    const int* m_lev_gid_to_bucket;
    const int* m_lev_offsets;

    GetBucket (const Gpu::DeviceVector<int>& lev_gid_to_bucket,
               const Gpu::DeviceVector<int>& lev_offsets)
        : m_lev_gid_to_bucket(lev_gid_to_bucket.dataPtr()),
          m_lev_offsets(lev_offsets.dataPtr())
        {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const int lev, const int gid) const noexcept
    {
        return m_lev_gid_to_bucket[m_lev_offsets[lev]+gid];
    }
};

class ParticleBufferMap
{
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

    AMREX_FORCE_INLINE
    int numLevels () const
    {
        AMREX_ASSERT(m_defined);
        return m_lev_offsets.size()-1;
    }

    AMREX_FORCE_INLINE
    int numBuckets () const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_gid.size();
    }

    AMREX_FORCE_INLINE
    int bucketToGrid (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_gid[bid];
    }

    AMREX_FORCE_INLINE
    int bucketToLevel (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_lev[bid];
    }

    AMREX_FORCE_INLINE
    int bucketToProc (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_pid[bid];
    }

    AMREX_FORCE_INLINE
    int gridAndLevToBucket (int gid, int lev) const
    {
        AMREX_ASSERT(m_defined);
        return m_lev_gid_to_bucket[m_lev_offsets[lev] + gid];
    }

    AMREX_FORCE_INLINE
    int firstBucketOnProc (int pid) const
    {
        AMREX_ASSERT(m_defined);
        return m_proc_box_offsets[pid];
    }

    AMREX_FORCE_INLINE
    int numBoxesOnProc (int pid) const
    {
        AMREX_ASSERT(m_defined);
        return m_proc_box_counts[pid];
    }

    AMREX_FORCE_INLINE
    Vector<int> allBucketsOnProc (int pid) const
    {
        AMREX_ASSERT(m_defined);
        Vector<int> buckets;
        for (int i = 0; i < numBoxesOnProc(pid); ++i)
        {
            buckets.push_back(i+firstBucketOnProc(pid));
        }
        return buckets;
    }

    AMREX_FORCE_INLINE
    int procID (int gid, int lev) const
    {
        AMREX_ASSERT(m_defined);
        return m_dm[lev][gid];
    }

    GetPID getPIDFunctor () const noexcept { return GetPID(d_bucket_to_pid, d_lev_gid_to_bucket, d_lev_offsets);}
    GetBucket getBucketFunctor () const noexcept { return GetBucket(d_lev_gid_to_bucket, d_lev_offsets);}
};

} // namespace amrex

#endif // AMREX_PARTICLEBUFFERMAP_H_
