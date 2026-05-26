#ifndef AMREX_PARTICLEBUFFERMAP_H_
#define AMREX_PARTICLEBUFFERMAP_H_
#include <AMReX_Config.H>

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParGDB.H>

namespace amrex {

struct GetPID
{
    const int* m_bucket_to_pid;
    const int* m_lev_gid_tid_to_bucket;
    const int* m_lev_offsets;
    const int* m_gid_offsets;

    GetPID (const Gpu::DeviceVector<int>& bucket_to_pid,
            const Gpu::DeviceVector<int>& lev_gid_tid_to_bucket,
            const Gpu::DeviceVector<int>& lev_offsets,
            const Gpu::DeviceVector<int>& gid_offsets)
        : m_bucket_to_pid(bucket_to_pid.dataPtr()),
          m_lev_gid_tid_to_bucket(lev_gid_tid_to_bucket.dataPtr()),
          m_lev_offsets(lev_offsets.dataPtr()),
          m_gid_offsets(gid_offsets.dataPtr())
        {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const int lev, const int gid, const int tid) const noexcept
    {
        return m_bucket_to_pid[m_lev_gid_tid_to_bucket[m_gid_offsets[m_lev_offsets[lev]+gid] + tid]];
    }
};

struct GetBucket
{
    const int* m_lev_gid_tid_to_bucket;
    const int* m_lev_offsets;
    const int* m_gid_offsets;

    GetBucket (const int* lev_gid_tid_to_bucket_ptr,
               const int* lev_offsets_ptr,
               const int* gid_offsets_ptr)
        : m_lev_gid_tid_to_bucket(lev_gid_tid_to_bucket_ptr),
          m_lev_offsets(lev_offsets_ptr),
          m_gid_offsets(gid_offsets_ptr)
        {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int operator() (const int lev, const int gid, const int tid) const noexcept
    {
        return m_lev_gid_tid_to_bucket[m_gid_offsets[m_lev_offsets[lev]+gid] + tid];
    }
};

class ParticleBufferMap
{
    bool m_defined{false};
    Vector<BoxArray> m_ba;
    Vector<DistributionMapping> m_dm;

    Vector<int> m_bucket_to_gid;
    Vector<int> m_bucket_to_tid;
    Vector<int> m_bucket_to_lev;
    Vector<int> m_bucket_to_pid;

    Vector<int> m_lev_gid_tid_to_bucket;
    Vector<int> m_lev_offsets;
    Vector<int> m_gid_offsets;

    Vector<int> m_proc_box_counts;
    Vector<int> m_proc_box_offsets;

    bool m_do_tiling{false};
    IntVect m_tile_size{AMREX_D_DECL(1024000, 1024000, 1024000)};

    Gpu::DeviceVector<int> d_bucket_to_pid;
    Gpu::DeviceVector<int> d_lev_gid_tid_to_bucket;
    Gpu::DeviceVector<int> d_lev_offsets;
    Gpu::DeviceVector<int> d_gid_offsets;

public:
    ParticleBufferMap () = default;

    ParticleBufferMap (const ParGDBBase* a_gdb);

    ParticleBufferMap (const ParGDBBase* a_gdb, bool a_do_tiling, const IntVect& a_tile_size);

    void define (const ParGDBBase* a_gdb);

    void define (const ParGDBBase* a_gdb, bool a_do_tiling, const IntVect& a_tile_size);

    bool isValid (const ParGDBBase* a_gdb) const;

    bool isValid (const ParGDBBase* a_gdb, bool a_do_tiling, const IntVect& a_tile_size) const;

    [[nodiscard]] AMREX_FORCE_INLINE
    int numLevels () const
    {
        AMREX_ASSERT(m_defined);
        return static_cast<int>(m_lev_offsets.size()-1);
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int numBuckets () const
    {
        AMREX_ASSERT(m_defined);
        return static_cast<int>(m_bucket_to_gid.size());
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int bucketToGrid (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_gid[bid];
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int bucketToTile (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_tid[bid];
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int bucketToLevel (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_lev[bid];
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int bucketToProc (int bid) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_pid[bid];
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int gridAndLevToBucket (int gid, int lev) const
    {
        AMREX_ASSERT(!m_do_tiling);
        return gridAndTileAndLevToBucket(gid, 0, lev);
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int gridAndTileAndLevToBucket (int gid, int tid, int lev) const
    {
        AMREX_ASSERT(m_defined);
        return m_lev_gid_tid_to_bucket[m_gid_offsets[m_lev_offsets[lev] + gid] + tid];
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int firstBucketOnProc (int pid) const
    {
        AMREX_ASSERT(m_defined);
        return m_proc_box_offsets[pid];
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int numBoxesOnProc (int pid) const
    {
        AMREX_ASSERT(m_defined);
        return m_proc_box_counts[pid];
    }

    [[nodiscard]] AMREX_FORCE_INLINE
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

    [[nodiscard]] AMREX_FORCE_INLINE
    int procID (int gid, int lev) const
    {
        AMREX_ASSERT(!m_do_tiling);
        return procID(gid, 0, lev);
    }

    [[nodiscard]] AMREX_FORCE_INLINE
    int procID (int gid, int tid, int lev) const
    {
        AMREX_ASSERT(m_defined);
        return m_bucket_to_pid[gridAndTileAndLevToBucket(gid, tid, lev)];
    }

    [[nodiscard]] GetPID getPIDFunctor () const noexcept { return GetPID(d_bucket_to_pid, d_lev_gid_tid_to_bucket, d_lev_offsets, d_gid_offsets);}
    [[nodiscard]] GetBucket getBucketFunctor () const noexcept { return GetBucket(d_lev_gid_tid_to_bucket.data(), d_lev_offsets.data(), d_gid_offsets.data());}
    [[nodiscard]] GetBucket getHostBucketFunctor () const noexcept { return GetBucket(m_lev_gid_tid_to_bucket.data(), m_lev_offsets.data(), m_gid_offsets.data());}

};

} // namespace amrex

#endif // AMREX_PARTICLEBUFFERMAP_H_
