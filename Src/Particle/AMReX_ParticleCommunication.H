#ifndef AMREX_PARTICLECOMMUNICATION_H_
#define AMREX_PARTICLECOMMUNICATION_H_
#include <AMReX_Config.H>

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_MFIter.H>
#include <AMReX_OpenMP.H>
#include <AMReX_Scan.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_MakeParticle.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <iterator>
#include <map>
#include <numeric>
#include <utility>

namespace amrex {

class ParticleContainerBase;

struct NeighborUnpackPolicy
{
    template <class PTile>
    void resizeTiles (std::vector<PTile*>& tiles, const std::vector<int>& sizes, std::vector<int>& offsets) const
    {
        for(int i = 0; i < std::ssize(sizes); ++i)
        {
            int offset = tiles[i]->numTotalParticles();
            int nn = tiles[i]->getNumNeighbors();
            tiles[i]->setNumNeighbors(nn + sizes[i]);
            offsets.push_back(offset);
        }
    }
};

struct RedistributeUnpackPolicy
{
    template <class PTile>
    void resizeTiles (std::vector<PTile*>& tiles, const std::vector<int>& sizes, std::vector<int>& offsets) const
    {
        int N = static_cast<int>(sizes.size());

        std::map<PTile*, int> tile_sizes;
        for(int i = 0; i < N; ++i) {
            tile_sizes[tiles[i]] = tiles[i]->numParticles();
        }

        for(int i = 0; i < N; ++i)
        {
            offsets.push_back(tile_sizes[tiles[i]]);
            tile_sizes[tiles[i]] += sizes[i];
        }

        for (auto& kv : tile_sizes) {
            kv.first->resize(kv.second);
        }
    }
};

struct ParticleCopyOp
{
    using TileKey = std::pair<int, int>;

    Vector<std::map<TileKey, Gpu::DeviceVector<int> > > m_boxes;
    Vector<std::map<TileKey, Gpu::DeviceVector<int> > > m_levels;
    Vector<std::map<TileKey, Gpu::DeviceVector<int> > > m_tiles;
    Vector<std::map<TileKey, Gpu::DeviceVector<int> > > m_src_indices;
    Vector<std::map<TileKey, Gpu::DeviceVector<IntVect> > > m_periodic_shift;

    void clear ();

    void setNumLevels (int num_levels);

    void resize (int gid, int tid, int lev, int size);

    [[nodiscard]] int numCopies (TileKey const& index, int lev) const
    {
        if (m_boxes.size() <= lev) { return 0; }
        auto mit = m_boxes[lev].find(index);
        return mit == m_boxes[lev].end() ? 0 : int(mit->second.size());
    }

    [[nodiscard]] int numLevels () const { return int(m_boxes.size()); }
};

struct ParticleCopyPlan
{
    using TileKey = std::pair<int, int>;

    struct StableOrderedAlgorithm {};
    struct TwoPassHostAlgorithm {};
    struct AtomicScatterAlgorithm {};

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    struct OmpCopyWork
    {
        const int* boxes = nullptr;
        const int* levs = nullptr;
        const int* tiles = nullptr;
        int* dst_indices = nullptr;
        int num_copies = 0;
    };
#endif

    struct BuildWorkspace
    {
        explicit BuildWorkspace (int a_num_buckets)
            : num_buckets(a_num_buckets)
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
            , omp_copy_offsets(1, 0)
#endif
        {}

        int num_buckets = 0;
        Gpu::HostVector<unsigned int> h_box_counts;

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
        Vector<unsigned int> omp_thread_box_counts;
        Vector<OmpCopyWork> omp_copy_work;
        Vector<Long> omp_copy_offsets;
#endif
    };

    template <class PC, class F>
    void forEachCopyBatch (const PC& pc, const ParticleCopyOp& op, F&& f)
    {
        const int num_levels = op.numLevels();
        m_dst_indices.resize(num_levels);
        auto&& batch_fn = std::forward<F>(f);

        for (int lev = 0; lev < num_levels; ++lev)
        {
            for (const auto& kv : pc.GetParticles(lev))
            {
                auto index = kv.first;
                int num_copies = op.numCopies(index, lev);
                if (num_copies == 0) { continue; }

                auto& dst_indices = m_dst_indices[lev][index];
                dst_indices.resize(num_copies);

                batch_fn(lev, index, num_copies, dst_indices);
            }
        }
    }

    template <class PC, class GetBucket>
    void buildCopies (const PC& pc, const ParticleCopyOp& op,
                      StableOrderedAlgorithm, BuildWorkspace& workspace, GetBucket const& getBucket)
    {
        BL_PROFILE("ParticleCopyPlan::buildCopiesStableOrdered");

        forEachCopyBatch(pc, op,
        [&] (int lev, TileKey const& index, int num_copies, Gpu::DeviceVector<int>& dst_indices)
        {
#ifdef AMREX_USE_GPU
            const Gpu::DeviceVector<int>& d_boxes = op.m_boxes[lev].at(index);
            Gpu::HostVector<int> h_boxes(d_boxes.size());
            Gpu::copy(Gpu::deviceToHost, d_boxes.begin(), d_boxes.end(), h_boxes.begin());

            const Gpu::DeviceVector<int>& d_levs = op.m_levels[lev].at(index);
            Gpu::HostVector<int> h_levs(d_levs.size());
            Gpu::copy(Gpu::deviceToHost, d_levs.begin(), d_levs.end(), h_levs.begin());

            const Gpu::DeviceVector<int>& d_tiles = op.m_tiles[lev].at(index);
            Gpu::HostVector<int> h_tiles(d_tiles.size());
            Gpu::copy(Gpu::deviceToHost, d_tiles.begin(), d_tiles.end(), h_tiles.begin());

        Gpu::HostVector<int> h_dst_indices(num_copies);
#else
            const Gpu::DeviceVector<int>& h_boxes = op.m_boxes[lev].at(index);
            const Gpu::DeviceVector<int>& h_levs = op.m_levels[lev].at(index);
            const Gpu::DeviceVector<int>& h_tiles = op.m_tiles[lev].at(index);

        Gpu::DeviceVector<int>& h_dst_indices = dst_indices;
#endif
            for (int i = 0; i < num_copies; ++i) {
                int dst_box = h_boxes[i];
                if (dst_box >= 0) {
                    int dst_tile = h_tiles[i];
                    int dst_lev = h_levs[i];
                    int dst_index = static_cast<int>(workspace.h_box_counts[getBucket(dst_lev, dst_box, dst_tile)]++);
                    h_dst_indices[i] = dst_index;
                }
            }

#ifdef AMREX_USE_GPU
            Gpu::copy(Gpu::hostToDevice,
                      h_dst_indices.begin(), h_dst_indices.end(),
                      dst_indices.begin());
#endif
        });
    }

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    template <class PC, class GetBucket>
    void buildCopies (const PC& pc, const ParticleCopyOp& op,
                      TwoPassHostAlgorithm, BuildWorkspace& workspace, GetBucket const& getBucket)
    {
        BL_PROFILE("ParticleCopyPlan::buildCopiesTwoPassHost");

        workspace.omp_thread_box_counts.resize(OpenMP::get_max_threads()*workspace.num_buckets, 0U);

        forEachCopyBatch(pc, op,
        [&op, &workspace] (int lev, TileKey const& index, int num_copies, Gpu::DeviceVector<int>& dst_indices)
        {
            const auto* p_boxes = op.m_boxes[lev].at(index).dataPtr();
            const auto* p_levs = op.m_levels[lev].at(index).dataPtr();
            const auto* p_tiles = op.m_tiles[lev].at(index).dataPtr();
            auto* p_dst_indices = dst_indices.dataPtr();

            workspace.omp_copy_work.push_back({p_boxes, p_levs, p_tiles, p_dst_indices, num_copies});
            workspace.omp_copy_offsets.push_back(workspace.omp_copy_offsets.back() + num_copies);
        });

        if (workspace.omp_copy_work.empty()) { return; }

        std::fill(workspace.omp_thread_box_counts.begin(), workspace.omp_thread_box_counts.end(), 0U);
        auto* p_omp_thread_box_counts = workspace.omp_thread_box_counts.data();
        const auto* p_omp_copy_work = workspace.omp_copy_work.data();
        const auto* p_omp_copy_offsets = workspace.omp_copy_offsets.data();
        const Long total_num_copies = workspace.omp_copy_offsets.back();

#pragma omp parallel
        {
            int thread_num = OpenMP::get_thread_num();
            int num_threads = OpenMP::get_num_threads();
            Long ibegin = thread_num*total_num_copies/num_threads;
            Long iend = (thread_num+1)*total_num_copies/num_threads;
            auto* p_thread_box_counts = p_omp_thread_box_counts + thread_num*workspace.num_buckets;

            if (ibegin < iend)
            {
                int iwork = static_cast<int>(std::upper_bound(workspace.omp_copy_offsets.begin(),
                                                              workspace.omp_copy_offsets.end(),
                                                              ibegin)
                                             - workspace.omp_copy_offsets.begin()) - 1;
                while (iwork < static_cast<int>(workspace.omp_copy_work.size()) &&
                       p_omp_copy_offsets[iwork] < iend)
                {
                    auto const& work = p_omp_copy_work[iwork];
                    Long global_begin = std::max(ibegin, p_omp_copy_offsets[iwork]);
                    Long global_end = std::min(iend, p_omp_copy_offsets[iwork+1]);
                    int local_begin = static_cast<int>(global_begin - p_omp_copy_offsets[iwork]);
                    int local_end = static_cast<int>(global_end - p_omp_copy_offsets[iwork]);
                    for (int i = local_begin; i < local_end; ++i)
                    {
                        int dst_box = work.boxes[i];
                        if (dst_box >= 0)
                        {
                            int dst_tile = work.tiles[i];
                            int dst_lev = work.levs[i];
                            ++p_thread_box_counts[getBucket(dst_lev, dst_box, dst_tile)];
                        }
                    }
                    ++iwork;
                }
            }

#pragma omp barrier
#pragma omp for
            for (int ibucket = 0; ibucket < workspace.num_buckets; ++ibucket)
            {
                unsigned int offset = workspace.h_box_counts[ibucket];
                for (int tid = 0; tid < num_threads; ++tid)
                {
                    auto& count = p_omp_thread_box_counts[tid*workspace.num_buckets + ibucket];
                    unsigned int total = count;
                    count = offset;
                    offset += total;
                }
                workspace.h_box_counts[ibucket] = offset;
            }

            if (ibegin < iend)
            {
                int iwork = static_cast<int>(std::upper_bound(workspace.omp_copy_offsets.begin(),
                                                              workspace.omp_copy_offsets.end(),
                                                              ibegin)
                                             - workspace.omp_copy_offsets.begin()) - 1;
                while (iwork < static_cast<int>(workspace.omp_copy_work.size()) &&
                       p_omp_copy_offsets[iwork] < iend)
                {
                    auto const& work = p_omp_copy_work[iwork];
                    Long global_begin = std::max(ibegin, p_omp_copy_offsets[iwork]);
                    Long global_end = std::min(iend, p_omp_copy_offsets[iwork+1]);
                    int local_begin = static_cast<int>(global_begin - p_omp_copy_offsets[iwork]);
                    int local_end = static_cast<int>(global_end - p_omp_copy_offsets[iwork]);
                    for (int i = local_begin; i < local_end; ++i)
                    {
                        int dst_box = work.boxes[i];
                        if (dst_box >= 0)
                        {
                            int dst_tile = work.tiles[i];
                            int dst_lev = work.levs[i];
                            int bucket = getBucket(dst_lev, dst_box, dst_tile);
                            work.dst_indices[i] = static_cast<int>(p_thread_box_counts[bucket]++);
                        }
                    }
                    ++iwork;
                }
            }
        }
    }
#endif

    template <class PC, class GetBucket>
    void buildCopies (const PC& pc, const ParticleCopyOp& op,
                      AtomicScatterAlgorithm, BuildWorkspace&, GetBucket const& getBucket)
    {
        BL_PROFILE("ParticleCopyPlan::buildCopiesAtomicScatter");

        auto* p_dst_box_counts = m_box_counts_d.dataPtr();

        forEachCopyBatch(pc, op,
        [&op, &getBucket, p_dst_box_counts] (int lev, TileKey const& index,
                                             int num_copies, Gpu::DeviceVector<int>& dst_indices)
        {
            const auto* p_boxes = op.m_boxes[lev].at(index).dataPtr();
            const auto* p_levs = op.m_levels[lev].at(index).dataPtr();
            const auto* p_tiles = op.m_tiles[lev].at(index).dataPtr();
            auto* p_dst_indices = dst_indices.dataPtr();

            amrex::ParallelFor(num_copies, [=] AMREX_GPU_DEVICE (int i)
            {
                int dst_box = p_boxes[i];
                if (dst_box >= 0)
                {
                    int dst_tile = p_tiles[i];
                    int dst_lev = p_levs[i];
                    int dst_index = static_cast<int>(HostDevice::Atomic::FetchAdd(
                        &p_dst_box_counts[getBucket(dst_lev, dst_box, dst_tile)], 1U));
                    p_dst_indices[i] = dst_index;
                }
            });
        });
    }

    void finalizeBuildBoxCounts (BuildWorkspace const& workspace, bool use_host_box_counters)
    {
        if (use_host_box_counters) {
            Gpu::copy(Gpu::hostToDevice,
                      workspace.h_box_counts.begin(), workspace.h_box_counts.end(),
                      m_box_counts_d.begin());
        }

        amrex::Gpu::exclusive_scan(m_box_counts_d.begin(), m_box_counts_d.end(),
                                   m_box_offsets.begin());
    }

public:

    Vector<std::map<TileKey, Gpu::DeviceVector<int> > > m_dst_indices;

    Gpu::DeviceVector<unsigned int> m_box_counts_d;
    Gpu::HostVector<unsigned int>   m_box_counts_h;
    Gpu::DeviceVector<unsigned int> m_box_offsets;

    Vector<int> m_rcv_box_counts;
    Vector<int> m_rcv_box_offsets;
    Vector<int> m_rcv_box_ids;
    Vector<int> m_rcv_box_tids;
    Vector<int> m_rcv_box_pids;
    Vector<int> m_rcv_box_levs;

    Long m_NumSnds = 0;
    int m_nrcvs = 0;
    mutable Vector<MPI_Status> m_build_stats;
    mutable Vector<MPI_Request> m_build_rreqs;

    mutable Vector<MPI_Status> m_particle_rstats;
    mutable Vector<MPI_Request> m_particle_rreqs;

    mutable Vector<MPI_Status> m_particle_sstats;
    mutable Vector<MPI_Request> m_particle_sreqs;

    Vector<Long> m_snd_num_particles;
    Vector<Long> m_rcv_num_particles;

    Vector<int> m_neighbor_procs;

    Vector<Long> m_Snds;
    Vector<Long> m_Rcvs;
    Vector<int> m_RcvProc;
    Vector<std::size_t> m_rOffset;
    Gpu::HostVector<int> m_rcv_data;

    Vector<std::size_t> m_snd_offsets;
    Vector<std::size_t> m_snd_counts;

    Vector<std::size_t> m_snd_pad_correction_h;
    Gpu::DeviceVector<std::size_t> m_snd_pad_correction_d;

    Vector<std::size_t> m_rcv_pad_correction_h;
    Gpu::DeviceVector<std::size_t> m_rcv_pad_correction_d;

    Gpu::DeviceVector<int> d_int_comp_mask, d_real_comp_mask;
    Long m_superparticle_size;

    Long superParticleSize() const { return m_superparticle_size; }

    template <class PC>
    requires (IsParticleContainer<PC>::value)
    void build (const PC& pc,
                const ParticleCopyOp& op,
                const Vector<int>& int_comp_mask,
                const Vector<int>& real_comp_mask,
                int local)
    {
        BL_PROFILE("ParticleCopyPlan::build");

        ParmParse pp("particles");
        pp.query("do_one_sided_comms", m_do_one_sided_comms);
          const int num_buckets = pc.BufferMap().numBuckets();

        m_local = local;
        if (local)
        {
            m_neighbor_procs = pc.NeighborProcs(local);
        }
        else
        {
            m_neighbor_procs.resize(ParallelContext::NProcsSub());
            std::iota(m_neighbor_procs.begin(), m_neighbor_procs.end(), 0);
        }

        m_box_counts_d.resize(0);
        m_box_counts_d.resize(num_buckets+1, 0);
        m_box_offsets.resize(num_buckets+1);

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
        constexpr bool use_host_bucket_counters = true;
#else
        constexpr bool use_host_bucket_counters = false;
#endif
        BuildWorkspace workspace(num_buckets);
        bool use_host_box_counters = pc.stableRedistribute() || use_host_bucket_counters;
        if (use_host_box_counters) {
            workspace.h_box_counts.resize(m_box_counts_d.size(), 0);
        }

        if (pc.stableRedistribute())
        {
            buildCopies(pc, op, StableOrderedAlgorithm{}, workspace, pc.BufferMap().getHostBucketFunctor());
        }
        else
        {
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
            buildCopies(pc, op, TwoPassHostAlgorithm{}, workspace, pc.BufferMap().getBucketFunctor());
#else
            buildCopies(pc, op, AtomicScatterAlgorithm{}, workspace, pc.BufferMap().getBucketFunctor());
#endif
        }

        finalizeBuildBoxCounts(workspace, use_host_box_counters);

        m_box_counts_h.resize(m_box_counts_d.size());
        Gpu::copyAsync(Gpu::deviceToHost, m_box_counts_d.begin(), m_box_counts_d.end(),
                       m_box_counts_h.begin());

        m_snd_pad_correction_h.resize(0);
        m_snd_pad_correction_h.resize(ParallelContext::NProcsSub()+1, 0);

        m_snd_pad_correction_d.resize(m_snd_pad_correction_h.size());
        Gpu::copyAsync(Gpu::hostToDevice, m_snd_pad_correction_h.begin(), m_snd_pad_correction_h.end(),
                       m_snd_pad_correction_d.begin());

        d_int_comp_mask.resize(int_comp_mask.size());
        Gpu::copyAsync(Gpu::hostToDevice,  int_comp_mask.begin(),  int_comp_mask.end(),
                       d_int_comp_mask.begin());
        d_real_comp_mask.resize(real_comp_mask.size());
        Gpu::copyAsync(Gpu::hostToDevice, real_comp_mask.begin(), real_comp_mask.end(),
                       d_real_comp_mask.begin());

        Gpu::streamSynchronize();

        int NStructReal = PC::ParticleContainerType::NStructReal;
        int NStructInt  = PC::ParticleContainerType::NStructInt;

        int num_real_comm_comp = 0;
        int comm_comps_start = 0;
        if constexpr (!PC::ParticleType::is_soa_particle) {
            comm_comps_start += AMREX_SPACEDIM + NStructReal;
        }
        for (int i = comm_comps_start; i < std::ssize(real_comp_mask); ++i) {
            if (real_comp_mask[i]) {++num_real_comm_comp;}
        }

        int num_int_comm_comp = 0;
        for (int i = 2 + NStructInt; i < std::ssize(int_comp_mask); ++i) {
            if (int_comp_mask[i])  {++num_int_comm_comp;}
        }

        if constexpr (PC::ParticleType::is_soa_particle) {
            m_superparticle_size = sizeof(uint64_t);  // idcpu
        } else {
            m_superparticle_size = sizeof(typename PC::ParticleType);
        }
        m_superparticle_size += num_real_comm_comp * sizeof(typename PC::ParticleType::RealType)
                              + num_int_comm_comp  * sizeof(int);

        buildMPIStart(pc, pc.BufferMap(), m_superparticle_size);
    }

    void clear ();

    void buildMPIFinish (const ParticleBufferMap& map);

private:

    void buildMPIStart (const ParticleContainerBase& pc, const ParticleBufferMap& map, Long psize);

    //
    // Snds - a Vector with the number of bytes that is process will send to each proc.
    // Rcvs - a Vector that, after calling this method, will contain the
    //        number of bytes this process will receive from each proc.
    //
    void doHandShake (const ParticleContainerBase& pc, const Vector<Long>& Snds, Vector<Long>& Rcvs) const;

    //
    // In the local version of this method, each proc knows which other
    // procs it could possibly receive messages from, meaning we can do
    // this purely with point-to-point communication.
    //
    void doHandShakeLocal (const Vector<Long>& Snds, Vector<Long>& Rcvs) const;

    //
    // In the global version, we don't know who we'll receive from, so we
    // need to do some collective communication first.
    //
    static void doHandShakeReduceScatter (const Vector<Long>& Snds, Vector<Long>& Rcvs);

    //
    // Another version of the global handshake implemented with MPI-3
    // one-sided communication.
    //
    static void doHandShakeOneSided (const ParticleContainerBase& pc,
                                     const Vector<Long>& Snds, Vector<Long>& Rcvs);

    //
    // Another version of the above that is implemented using MPI All-to-All
    //
    static void doHandShakeAllToAll (const Vector<Long>& Snds, Vector<Long>& Rcvs);

    bool m_local = false;
    int m_do_one_sided_comms = 0;
};

struct GetSendBufferOffset
{
    const unsigned int* m_box_offsets;
    const std::size_t* m_pad_correction;

    GetPID m_get_pid;
    GetBucket m_get_bucket;

    GetSendBufferOffset (const ParticleCopyPlan& plan, const ParticleBufferMap& map)
        : m_box_offsets(plan.m_box_offsets.dataPtr()),
          m_pad_correction(plan.m_snd_pad_correction_d.dataPtr()),
          m_get_pid(map.getPIDFunctor()),
          m_get_bucket(map.getBucketFunctor())
    {}

    AMREX_FORCE_INLINE AMREX_GPU_DEVICE
    Long operator() (int dst_box, int dst_tile, int dst_lev, std::size_t psize, int i) const
    {
        int dst_pid = m_get_pid(dst_lev, dst_box, dst_tile);
        Long dst_offset = Long(psize)*(m_box_offsets[m_get_bucket(dst_lev, dst_box, dst_tile)] + i);
        dst_offset += Long(m_pad_correction[dst_pid]);
        return dst_offset;
    }
};

template <class PC, class Buffer>
requires (IsParticleContainer<PC>::value &&
          std::derived_from<Buffer,PolymorphicArenaAllocator<typename Buffer::value_type>>)
void packBuffer (const PC& pc, const ParticleCopyOp& op, const ParticleCopyPlan& plan,
                 Buffer& snd_buffer)
{
    BL_PROFILE("amrex::packBuffer");

    Long psize = plan.superParticleSize();

    int num_levels = op.numLevels();
    int num_buckets = pc.BufferMap().numBuckets();

    std::size_t total_buffer_size = 0;
    if (plan.m_snd_offsets.empty())
    {
        unsigned int np = 0;
        Gpu::copy(Gpu::deviceToHost, plan.m_box_offsets.begin() + num_buckets,
                  plan.m_box_offsets.begin() + num_buckets + 1, &np);
        total_buffer_size = np*psize;
    }
    else
    {
        total_buffer_size = plan.m_snd_offsets.back();
    }

    if (! snd_buffer.arena()->hasFreeDeviceMemory(total_buffer_size)) {
        snd_buffer.clear();
        snd_buffer.setArena(The_Pinned_Arena());
    }
    snd_buffer.resize(total_buffer_size);

    const auto* p_comm_real = plan.d_real_comp_mask.dataPtr();
    const auto* p_comm_int  = plan.d_int_comp_mask.dataPtr();

    const auto plo = pc.Geom(0).ProbLoArray();
    const auto phi = pc.Geom(0).ProbHiArray();
    const auto is_per = pc.Geom(0).isPeriodicArray();
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    struct OmpPackWork
    {
        typename PC::ParticleTileType const* src_tile = nullptr;
        const int* boxes = nullptr;
        const int* levels = nullptr;
        const int* tiles = nullptr;
        const int* src_indices = nullptr;
        const IntVect* periodic_shift = nullptr;
        const int* dst_indices = nullptr;
        int num_copies = 0;
    };
    Vector<OmpPackWork> omp_pack_work;
    Vector<Long> omp_pack_offsets(1, 0);
#endif
    for (int lev = 0; lev < num_levels; ++lev)
    {
        auto& plev = pc.GetParticles(lev);
        for (auto& kv : plev)
        {
            auto index = kv.first;
            auto& src_tile = plev.at(index);
            int num_copies = op.numCopies(index, lev);
            if (num_copies == 0) { continue; }

            const auto* p_boxes = op.m_boxes[lev].at(index).dataPtr();
            const auto* p_levels = op.m_levels[lev].at(index).dataPtr();
            const auto* p_tiles = op.m_tiles[lev].at(index).dataPtr();
            const auto* p_src_indices = op.m_src_indices[lev].at(index).dataPtr();
            const auto* p_periodic_shift = op.m_periodic_shift[lev].at(index).dataPtr();
            const auto* p_dst_indices = plan.m_dst_indices[lev].at(index).dataPtr();
#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
            omp_pack_work.push_back({&src_tile, p_boxes, p_levels, p_tiles,
                                     p_src_indices, p_periodic_shift, p_dst_indices, num_copies});
            omp_pack_offsets.push_back(omp_pack_offsets.back() + num_copies);
#else
            const auto& ptd = src_tile.getConstParticleTileData();
            auto* p_snd_buffer = snd_buffer.dataPtr();
            GetSendBufferOffset get_offset(plan, pc.BufferMap());
            amrex::ParallelFor(num_copies, [=] AMREX_GPU_DEVICE (int i)
            {
                int dst_box = p_boxes[i];
                if (dst_box >= 0)
                {
                    int dst_tile = p_tiles[i];
                    int dst_lev = p_levels[i];
                    auto dst_offset = get_offset(dst_box, dst_tile, dst_lev, psize, p_dst_indices[i]);
                    int src_index = p_src_indices[i];
                    ptd.packParticleData(p_snd_buffer, src_index, dst_offset, p_comm_real, p_comm_int);

                    const IntVect& pshift = p_periodic_shift[i];
                    bool do_periodic_shift =
                        AMREX_D_TERM( (is_per[0] && pshift[0] != 0),
                                   || (is_per[1] && pshift[1] != 0),
                                   || (is_per[2] && pshift[2] != 0) );

                    if (do_periodic_shift)
                    {
                        ParticleReal pos[AMREX_SPACEDIM];
                        Long pos_offset = dst_offset;
                        // for pure SoA positions come after idcpu
                        if constexpr (PC::ParticleType::is_soa_particle) {
                            pos_offset += sizeof(uint64_t);
                        }
                        amrex::Gpu::memcpy(&pos[0], &p_snd_buffer[pos_offset],
                                           AMREX_SPACEDIM*sizeof(ParticleReal));
                        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                        {
                            if (! is_per[idim]) { continue; }
                            if (pshift[idim] > 0) {
                                pos[idim] += phi[idim] - plo[idim];
                            } else if (pshift[idim] < 0) {
                                pos[idim] -= phi[idim] - plo[idim];
                            }
                        }
                        amrex::Gpu::memcpy(&p_snd_buffer[pos_offset], &pos[0],
                                           AMREX_SPACEDIM*sizeof(ParticleReal));
                    }
                }
            });
#endif
        }
    }

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    if (!omp_pack_work.empty())
    {
        auto* p_snd_buffer = snd_buffer.dataPtr();
        GetSendBufferOffset get_offset(plan, pc.BufferMap());
        const Long total_num_copies = omp_pack_offsets.back();

#pragma omp parallel
        {
            int thread_num = OpenMP::get_thread_num();
            int num_threads = OpenMP::get_num_threads();
            Long ibegin = thread_num*total_num_copies/num_threads;
            Long iend = (thread_num+1)*total_num_copies/num_threads;

            if (ibegin < iend)
            {
                int iwork = static_cast<int>(std::upper_bound(omp_pack_offsets.begin(),
                                                              omp_pack_offsets.end(),
                                                              ibegin)
                                             - omp_pack_offsets.begin()) - 1;
                while (iwork < static_cast<int>(omp_pack_work.size()) &&
                       omp_pack_offsets[iwork] < iend)
                {
                    auto const& work = omp_pack_work[iwork];
                    auto const& ptd = work.src_tile->getConstParticleTileData();
                    Long global_begin = std::max(ibegin, omp_pack_offsets[iwork]);
                    Long global_end = std::min(iend, omp_pack_offsets[iwork+1]);
                    int local_begin = static_cast<int>(global_begin - omp_pack_offsets[iwork]);
                    int local_end = static_cast<int>(global_end - omp_pack_offsets[iwork]);
                    for (int i = local_begin; i < local_end; ++i)
                    {
                        int dst_box = work.boxes[i];
                        if (dst_box >= 0)
                        {
                            int dst_tile = work.tiles[i];
                            int dst_lev = work.levels[i];
                            auto dst_offset = get_offset(dst_box, dst_tile, dst_lev, psize,
                                                         work.dst_indices[i]);
                            int src_index = work.src_indices[i];
                            ptd.packParticleData(p_snd_buffer, src_index, dst_offset,
                                                 p_comm_real, p_comm_int);

                            const IntVect& pshift = work.periodic_shift[i];
                            bool do_periodic_shift =
                                AMREX_D_TERM( (is_per[0] && pshift[0] != 0),
                                           || (is_per[1] && pshift[1] != 0),
                                           || (is_per[2] && pshift[2] != 0) );

                            if (do_periodic_shift)
                            {
                                ParticleReal pos[AMREX_SPACEDIM];
                                Long pos_offset = dst_offset;
                                if constexpr (PC::ParticleType::is_soa_particle) {
                                    pos_offset += sizeof(uint64_t);
                                }
                                amrex::Gpu::memcpy(&pos[0], &p_snd_buffer[pos_offset],
                                                   AMREX_SPACEDIM*sizeof(ParticleReal));
                                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                                {
                                    if (! is_per[idim]) { continue; }
                                    if (pshift[idim] > 0) {
                                        pos[idim] += phi[idim] - plo[idim];
                                    } else if (pshift[idim] < 0) {
                                        pos[idim] -= phi[idim] - plo[idim];
                                    }
                                }
                                amrex::Gpu::memcpy(&p_snd_buffer[pos_offset], &pos[0],
                                                   AMREX_SPACEDIM*sizeof(ParticleReal));
                            }
                        }
                    }
                    ++iwork;
                }
            }
        }
    }
#endif
}

template <class PC, class Buffer, class UnpackPolicy>
requires (IsParticleContainer<PC>::value)
void unpackBuffer (PC& pc, const ParticleCopyPlan& plan, const Buffer& snd_buffer, UnpackPolicy const& policy)
{
    BL_PROFILE("amrex::unpackBuffer");

    using PTile = typename PC::ParticleTileType;

    int num_levels = pc.BufferMap().numLevels();
    Long psize = plan.superParticleSize();

    // count how many particles we have to add to each tile
    std::vector<int> sizes;
    std::vector<PTile*> tiles;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            int tid = mfi.LocalTileIndex();
            auto& tile = pc.DefineAndReturnParticleTile(lev, gid, tid);
            int num_copies = plan.m_box_counts_h[pc.BufferMap().gridAndTileAndLevToBucket(gid, tid, lev)];
            sizes.push_back(num_copies);
            tiles.push_back(&tile);
        }
    }

    // resize the tiles and compute offsets
    std::vector<int> offsets;
    policy.resizeTiles(tiles, sizes, offsets);

    const auto* p_comm_real = plan.d_real_comp_mask.dataPtr();
    const auto* p_comm_int  = plan.d_int_comp_mask.dataPtr();

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    struct OmpUnpackWork
    {
        PTile* tile = nullptr;
        int gid = 0;
        int tid = 0;
        int lev = 0;
        int offset = 0;
        int size = 0;
    };
    Vector<OmpUnpackWork> omp_unpack_work;
    Vector<Long> omp_unpack_offsets(1, 0);
#endif

    // local unpack
    int uindex = 0;
    for (int lev = 0; lev < num_levels; ++lev)
    {
        auto& plev  = pc.GetParticles(lev);
        for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            int tid = mfi.LocalTileIndex();
            auto index = std::make_pair(gid, tid);

            auto& tile = plev[index];

            GetSendBufferOffset get_offset(plan, pc.BufferMap());

            int offset = offsets[uindex];
            int size = sizes[uindex];
            ++uindex;

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
            omp_unpack_work.push_back({&tile, gid, tid, lev, offset, size});
            omp_unpack_offsets.push_back(omp_unpack_offsets.back() + size);
#else
            auto p_snd_buffer = snd_buffer.dataPtr();
            auto ptd = tile.getParticleTileData();
            amrex::ParallelFor(size, [=] AMREX_GPU_DEVICE (int i)
            {
                auto src_offset = get_offset(gid, tid, lev, psize, i);
                int dst_index = offset + i;
                ptd.unpackParticleData(p_snd_buffer, src_offset, dst_index, p_comm_real, p_comm_int);
            });
#endif
        }
    }

#if defined(AMREX_USE_OMP) && !defined(AMREX_USE_GPU)
    if (!omp_unpack_work.empty())
    {
        GetSendBufferOffset get_offset(plan, pc.BufferMap());
        auto p_snd_buffer = snd_buffer.dataPtr();
        const Long total_num_copies = omp_unpack_offsets.back();

#pragma omp parallel
        {
            int thread_num = OpenMP::get_thread_num();
            int num_threads = OpenMP::get_num_threads();
            Long ibegin = thread_num*total_num_copies/num_threads;
            Long iend = (thread_num+1)*total_num_copies/num_threads;

            if (ibegin < iend)
            {
                int iwork = static_cast<int>(std::upper_bound(omp_unpack_offsets.begin(),
                                                              omp_unpack_offsets.end(),
                                                              ibegin)
                                             - omp_unpack_offsets.begin()) - 1;
                while (iwork < static_cast<int>(omp_unpack_work.size()) &&
                       omp_unpack_offsets[iwork] < iend)
                {
                    auto const& work = omp_unpack_work[iwork];
                    auto ptd = work.tile->getParticleTileData();
                    Long global_begin = std::max(ibegin, omp_unpack_offsets[iwork]);
                    Long global_end = std::min(iend, omp_unpack_offsets[iwork+1]);
                    int local_begin = static_cast<int>(global_begin - omp_unpack_offsets[iwork]);
                    int local_end = static_cast<int>(global_end - omp_unpack_offsets[iwork]);
                    for (int i = local_begin; i < local_end; ++i)
                    {
                        auto src_offset = get_offset(work.gid, work.tid, work.lev, psize, i);
                        int dst_index = work.offset + i;
                        ptd.unpackParticleData(p_snd_buffer, src_offset, dst_index,
                                               p_comm_real, p_comm_int);
                    }
                    ++iwork;
                }
            }
        }
    }
#endif
}

template <class PC, class SndBuffer, class RcvBuffer>
requires (IsParticleContainer<PC>::value)
void communicateParticlesStart (const PC& pc, ParticleCopyPlan& plan, const SndBuffer& snd_buffer, RcvBuffer& rcv_buffer)
{
    BL_PROFILE("amrex::communicateParticlesStart");

#ifdef AMREX_USE_MPI
    Long psize = plan.superParticleSize();
    const int NProcs = ParallelContext::NProcsSub();
    const int MyProc = ParallelContext::MyProcSub();

    if (NProcs == 1) { return; }

    Vector<int> RcvProc;
    Vector<Long> rOffset;

    plan.m_rcv_pad_correction_h.resize(0);
    plan.m_rcv_pad_correction_h.push_back(0);

    Long TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (plan.m_rcv_num_particles[i] > 0) {
            RcvProc.push_back(i);
            Long nbytes = plan.m_rcv_num_particles[i]*psize;
            std::size_t acd = ParallelDescriptor::sizeof_selected_comm_data_type(nbytes);
            TotRcvBytes = Long(amrex::aligned_size(acd, TotRcvBytes));
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += Long(amrex::aligned_size(acd, nbytes));
            plan.m_rcv_pad_correction_h.push_back(plan.m_rcv_pad_correction_h.back() + nbytes);
        }
    }

    for (int i = 0; i < plan.m_nrcvs; ++i)
    {
        plan.m_rcv_pad_correction_h[i] = rOffset[i] - plan.m_rcv_pad_correction_h[i];
    }

    plan.m_rcv_pad_correction_d.resize(plan.m_rcv_pad_correction_h.size());
    Gpu::copy(Gpu::hostToDevice, plan.m_rcv_pad_correction_h.begin(), plan.m_rcv_pad_correction_h.end(),
              plan.m_rcv_pad_correction_d.begin());

    rcv_buffer.resize(TotRcvBytes);

    plan.m_nrcvs = int(RcvProc.size());

    plan.m_particle_rstats.resize(0);
    plan.m_particle_rstats.resize(plan.m_nrcvs);

    plan.m_particle_rreqs.resize(0);
    plan.m_particle_rreqs.resize(plan.m_nrcvs);

    plan.m_particle_sstats.resize(0);
    plan.m_particle_sreqs.resize(0);

    const int SeqNum = ParallelDescriptor::SeqNum();

    // Post receives.
    for (int i = 0; i < plan.m_nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        Long nbytes       = plan.m_rcv_num_particles[Who]*psize;
        std::size_t acd   = ParallelDescriptor::sizeof_selected_comm_data_type(nbytes);
        const auto Cnt    = amrex::aligned_size(acd, nbytes);

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(amrex::aligned_size(acd, nbytes) % acd == 0);

        plan.m_particle_rreqs[i] =
            ParallelDescriptor::Arecv((char*) (rcv_buffer.dataPtr() + offset), Cnt, Who, SeqNum, ParallelContext::CommunicatorSub()).req();
    }

    if (plan.m_NumSnds == 0) { return; }

    // Send.
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) { continue; }
        const auto Who  = i;
        const auto Cnt  = plan.m_snd_counts[i];
        if (Cnt == 0) { continue; }

        auto snd_offset = plan.m_snd_offsets[i];
        AMREX_ASSERT(plan.m_snd_counts[i] % ParallelDescriptor::sizeof_selected_comm_data_type(plan.m_snd_num_particles[i]*psize) == 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        plan.m_particle_sreqs.push_back(ParallelDescriptor::Asend((char const*)(snd_buffer.dataPtr()+snd_offset), Cnt, Who, SeqNum,
                                                                  ParallelContext::CommunicatorSub()).req());
    }

    plan.m_particle_sstats.resize(plan.m_particle_sreqs.size());

    amrex::ignore_unused(pc);
#else
    amrex::ignore_unused(pc,plan,snd_buffer,rcv_buffer);
#endif // MPI
}

void communicateParticlesFinish (const ParticleCopyPlan& plan);

template <class PC, class Buffer, class UnpackPolicy>
requires (IsParticleContainer<PC>::value)
void unpackRemotes (PC& pc, const ParticleCopyPlan& plan, Buffer& rcv_buffer, UnpackPolicy const& policy)
{
    BL_PROFILE("amrex::unpackRemotes");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs == 1) { return; }

    const int MyProc = ParallelContext::MyProcSub();
    amrex::ignore_unused(MyProc);
    using PTile = typename PC::ParticleTileType;

    if (plan.m_nrcvs > 0)
    {
        const auto* p_comm_real = plan.d_real_comp_mask.dataPtr();
        const auto* p_comm_int  = plan.d_int_comp_mask.dataPtr();
        auto* p_rcv_buffer = rcv_buffer.dataPtr();

        std::vector<int> sizes;
        std::vector<PTile*> tiles;
        for (int i = 0; i < std::ssize(plan.m_rcv_box_counts); ++i)
        {
            int copy_size = plan.m_rcv_box_counts[i];
            int lev = plan.m_rcv_box_levs[i];
            int gid = plan.m_rcv_box_ids[i];
            int tid = plan.m_rcv_box_tids[i];
            auto& tile = pc.DefineAndReturnParticleTile(lev, gid, tid);
            sizes.push_back(copy_size);
            tiles.push_back(&tile);
        }

        Vector<int> offsets;
        policy.resizeTiles(tiles, sizes, offsets);
        Gpu::streamSynchronize();
        int uindex = 0;
        int procindex = 0, rproc = plan.m_rcv_box_pids[0];
        for (int i = 0; i < std::ssize(plan.m_rcv_box_counts); ++i)
        {
            int lev = plan.m_rcv_box_levs[i];
            int gid = plan.m_rcv_box_ids[i];
            int tid = plan.m_rcv_box_tids[i];
            auto offset = plan.m_rcv_box_offsets[i];
            procindex = (rproc == plan.m_rcv_box_pids[i]) ? procindex : procindex+1;
            rproc = plan.m_rcv_box_pids[i];

            auto& tile = pc.DefineAndReturnParticleTile(lev, gid, tid);
            auto ptd = tile.getParticleTileData();

            AMREX_ASSERT(MyProc ==
                ParallelContext::global_to_local_rank(pc.ParticleDistributionMap(lev)[gid]));

            int dst_offset = offsets[uindex];
            int size = sizes[uindex];
            ++uindex;

            Long psize = plan.superParticleSize();
            const auto* p_pad_adjust = plan.m_rcv_pad_correction_d.dataPtr();

            amrex::ParallelForOMP(size, [=] AMREX_GPU_DEVICE (int ip) {
                Long src_offset = psize * static_cast<Long>(offset + ip)
                    + static_cast<Long>(p_pad_adjust[procindex]);
                int dst_index = dst_offset + ip;
                ptd.unpackParticleData(p_rcv_buffer, src_offset, dst_index,
                                       p_comm_real, p_comm_int);
              });
        }
    }
#else
    amrex::ignore_unused(pc,plan,rcv_buffer,policy);
#endif // MPI
}

} // namespace amrex

#endif // AMREX_PARTICLECOMMUNICATION_H_
