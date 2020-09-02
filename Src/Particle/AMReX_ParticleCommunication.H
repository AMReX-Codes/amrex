#ifndef AMREX_PARTICLECOMMUNICATION_H_
#define AMREX_PARTICLECOMMUNICATION_H_

#include <AMReX_Gpu.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_MFIter.H>
#include <AMReX_TypeTraits.H>

#include <map>

namespace amrex {

struct NeighborUnpackPolicy
{
    template <class PTile>
    void resizeTiles (std::vector<PTile*>& tiles, const std::vector<int>& sizes, std::vector<int>& offsets) const
    {
        for(int i = 0; i < static_cast<int>(sizes.size()); ++i)
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
        for(int i = 0; i < N; ++i)
            tile_sizes[tiles[i]] = tiles[i]->numParticles();

        for(int i = 0; i < N; ++i)
        {
            offsets.push_back(tile_sizes[tiles[i]]);
            tile_sizes[tiles[i]] += sizes[i];
        }

        for (auto& kv : tile_sizes)
            kv.first->resize(kv.second);
    }
};

struct ParticleCopyOp
{
    Vector<std::map<int, Gpu::DeviceVector<int> > > m_boxes;
    Vector<std::map<int, Gpu::DeviceVector<int> > > m_levels;
    Vector<std::map<int, Gpu::DeviceVector<int> > > m_src_indices;
    Vector<std::map<int, Gpu::DeviceVector<IntVect> > > m_periodic_shift;

    void clear ();

    void setNumLevels (const int num_levels);

    void resize (const int gid, const int lev, const int size);

    int numCopies (const int gid, const int lev) const
    {
        if (m_boxes.size() <= lev) return 0;
        auto mit = m_boxes[lev].find(gid);
        return mit == m_boxes[lev].end() ? 0 : mit->second.size();
    }
};

struct ParticleCopyPlan
{
    Vector<std::map<int, Gpu::DeviceVector<int> > > m_dst_indices;

    Gpu::DeviceVector<unsigned int> m_box_counts;
    Gpu::DeviceVector<unsigned int> m_box_offsets;

    Gpu::DeviceVector<int> m_rcv_box_counts;
    Gpu::DeviceVector<int> m_rcv_box_offsets;
    Gpu::DeviceVector<int> m_rcv_box_ids;
    Gpu::DeviceVector<int> m_rcv_box_pids;
    Gpu::DeviceVector<int> m_rcv_box_levs;

    Long m_NumSnds;
    int m_nrcvs;
    mutable Vector<MPI_Status> m_build_stats;
    mutable Vector<MPI_Request> m_build_rreqs;

    mutable Vector<MPI_Status> m_particle_stats;
    mutable Vector<MPI_Request> m_particle_rreqs;

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

    template <class PC, EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
    void build (const PC& pc, const ParticleCopyOp& op, bool local)
    {
        BL_PROFILE("ParticleCopyPlan::build");

        m_local = local;

        const int ngrow = 1;  // note - fix

        const int num_levels = pc.BufferMap().numLevels();
        const int num_buckets = pc.BufferMap().numBuckets();

        if (m_local)
        {
            m_neighbor_procs = pc.NeighborProcs(ngrow);
        }
        else
        {
            m_neighbor_procs.resize(ParallelContext::NProcsSub());
            std::iota(m_neighbor_procs.begin(), m_neighbor_procs.end(), 0);
        }

        m_box_counts.resize(0);
        m_box_counts.resize(num_buckets+1, 0);
        m_box_offsets.resize(num_buckets+1);
        auto p_dst_box_counts = m_box_counts.dataPtr();
        auto getBucket = pc.BufferMap().getBucketFunctor();

        constexpr unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();

        m_dst_indices.resize(num_levels);
        for (int lev = 0; lev < num_levels; ++lev)
        {
            for (const auto& kv : pc.GetParticles(lev))
            {
                int gid = kv.first.first;
                int num_copies = op.numCopies(gid, lev);
                if (num_copies == 0) continue;
                m_dst_indices[lev][gid].resize(num_copies);

                auto p_boxes = op.m_boxes[lev].at(gid).dataPtr();
                auto p_levs = op.m_levels[lev].at(gid).dataPtr();
                auto p_dst_indices = m_dst_indices[lev][gid].dataPtr();

                AMREX_FOR_1D ( num_copies, i,
                {
                    int dst_box = p_boxes[i];
                    if (dst_box >= 0)
                    {
                        int dst_lev = p_levs[i];
                        int index = Gpu::Atomic::Inc(
                            &p_dst_box_counts[getBucket(dst_lev, dst_box)], max_unsigned_int);
                        p_dst_indices[i] = index;
                    }
                });
            }
        }

        amrex::Gpu::exclusive_scan(m_box_counts.begin(), m_box_counts.end(), m_box_offsets.begin());

        m_snd_pad_correction_h.resize(0);
        m_snd_pad_correction_h.resize(ParallelContext::NProcsSub()+1, 0);

        m_snd_pad_correction_d.resize(m_snd_pad_correction_h.size());
        Gpu::copy(Gpu::hostToDevice, m_snd_pad_correction_h.begin(), m_snd_pad_correction_h.end(),
                  m_snd_pad_correction_d.begin());

        buildMPIStart(pc.BufferMap(), pc.superParticleSize());
    }

    void clear ();

    void buildMPIFinish (const ParticleBufferMap& map);

private:

    void buildMPIStart (const ParticleBufferMap& map, Long psize);

    //
    // Snds - a Vector with the number of bytes that is process will send to each proc.
    // Rcvs - a Vector that, after calling this method, will contain the
    //        number of bytes this process will reveive from each proc.
    //
    void doHandShake (const Vector<Long>& Snds, Vector<Long>& Rcvs) const;

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
    void doHandShakeGlobal (const Vector<Long>& Snds, Vector<Long>& Rcvs) const;

    //
    // Another version of the above that is implemented using MPI All-to-All
    //
    void doHandShakeAllToAll (const Vector<Long>& Snds, Vector<Long>& Rcvs) const;

    bool m_local;
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

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Long operator() (int dst_box, int dst_lev, std::size_t psize, int i) const
    {
        int dst_pid = m_get_pid(dst_lev, dst_box);
        Long dst_offset = psize*(m_box_offsets[m_get_bucket(dst_lev, dst_box)] + i);
        dst_offset += m_pad_correction[dst_pid];
        return dst_offset;
    }
};

template <class PC, class Buffer, EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
void packBuffer (const PC& pc, const ParticleCopyOp& op, const ParticleCopyPlan& plan,
                 Buffer& snd_buffer)
{
    BL_PROFILE("amrex::packBuffer");

    using ParticleType = typename PC::ParticleType;

    Long psize = pc.superParticleSize();

    int num_levels = pc.BufferMap().numLevels();
    int num_buckets = pc.BufferMap().numBuckets();
    Long total_buffer_size = (plan.m_snd_offsets.size() == 0) ? plan.m_box_offsets[num_buckets]*psize : plan.m_snd_offsets.back();
    snd_buffer.resize(total_buffer_size);

    auto p_comm_real = pc.d_communicate_real_comp.dataPtr();
    auto p_comm_int  = pc.d_communicate_int_comp.dataPtr();

    for (int lev = 0; lev < num_levels; ++lev)
    {
        const auto& geom = pc.Geom(lev);
        auto& plev = pc.GetParticles(lev);
        auto& ba = pc.ParticleBoxArray(lev);
        const auto plo = geom.ProbLoArray();
        const auto phi = geom.ProbHiArray();
        const auto is_per = geom.isPeriodicArray();

        for (auto& kv : plev)
        {
            int gid = kv.first.first;
            int tid = kv.first.second;
            auto index = std::make_pair(gid, tid);

            auto& src_tile = plev.at(index);
            auto& aos   = src_tile.GetArrayOfStructs();
            const auto ptd = src_tile.getConstParticleTileData();

            int num_copies = op.numCopies(gid, lev);
            if (num_copies == 0) continue;

            auto p_boxes = op.m_boxes[lev].at(gid).dataPtr();
            auto p_levels = op.m_levels[lev].at(gid).dataPtr();
            auto p_src_indices = op.m_src_indices[lev].at(gid).dataPtr();
            auto p_periodic_shift = op.m_periodic_shift[lev].at(gid).dataPtr();
            auto p_dst_indices = plan.m_dst_indices[lev].at(gid).dataPtr();
            auto p_snd_buffer = snd_buffer.dataPtr();
            GetSendBufferOffset get_offset(plan, pc.BufferMap());

            AMREX_FOR_1D ( num_copies, i,
            {
                int dst_box = p_boxes[i];
                if (dst_box >= 0)
                {
                    int dst_lev = p_levels[i];
                    auto dst_offset = get_offset(dst_box, dst_lev, psize, p_dst_indices[i]);
                    int src_index = p_src_indices[i];
                    ptd.packParticleData(p_snd_buffer, src_index, dst_offset, p_comm_real, p_comm_int);

                    ParticleType* p = (ParticleType*) &p_snd_buffer[dst_offset];
                    const IntVect& pshift = p_periodic_shift[i];
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
                    {
                        if (not is_per[idim]) continue;
                        if (pshift[idim] > 0)
                            p->pos(idim) += phi[idim] - plo[idim];
                        else if (pshift[idim] < 0)
                            p->pos(idim) -= phi[idim] - plo[idim];
                    }
                }
            });
        }
    }
}

template <class PC, class Buffer, class UnpackPolicy,
          EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
void unpackBuffer (PC& pc, const ParticleCopyPlan& plan, const Buffer& snd_buffer, const UnpackPolicy&& policy)
{
    BL_PROFILE("amrex::unpackBuffer");

    using PTile = typename PC::ParticleTileType;

    int num_levels = pc.BufferMap().numLevels();
    Long psize = pc.superParticleSize();

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
            int num_copies = plan.m_box_counts[pc.BufferMap().gridAndLevToBucket(gid, lev)];
            sizes.push_back(num_copies);
            tiles.push_back(&tile);
        }
    }

    // resize the tiles and compute offsets
    std::vector<int> offsets;
    policy.resizeTiles(tiles, sizes, offsets);

    auto p_comm_real = pc.d_communicate_real_comp.dataPtr();
    auto p_comm_int  = pc.d_communicate_int_comp.dataPtr();

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
            auto& aos   = tile.GetArrayOfStructs();

            GetSendBufferOffset get_offset(plan, pc.BufferMap());
            auto p_snd_buffer = snd_buffer.dataPtr();

            int offset = offsets[uindex];
            int size = sizes[uindex];
            ++uindex;

            auto ptd = tile.getParticleTileData();
            AMREX_FOR_1D ( size, i,
            {
                auto src_offset = get_offset(gid, lev, psize, i);
                int dst_index = offset + i;
                ptd.unpackParticleData(p_snd_buffer, src_offset, dst_index, p_comm_real, p_comm_int);
            });
        }
    }
}

template <class PC, class Buffer,
          EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
void communicateParticlesStart (const PC& pc, ParticleCopyPlan& plan, const Buffer& snd_buffer, Buffer& rcv_buffer)
{
    BL_PROFILE("amrex::communicateParticlesStart");

#ifdef AMREX_USE_MPI
    Long psize = pc.superParticleSize();
    const int NProcs = ParallelContext::NProcsSub();
    const int MyProc = ParallelContext::MyProcSub();

    if (NProcs == 1) return;

    Vector<int> RcvProc;
    Vector<Long> rOffset;

    plan.m_rcv_pad_correction_h.resize(0);
    plan.m_rcv_pad_correction_h.push_back(0);

    Long TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (plan.m_rcv_num_particles[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            Long nbytes = plan.m_rcv_num_particles[i]*psize;
            std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
            TotRcvBytes = amrex::aligned_size(acd, TotRcvBytes);
            TotRcvBytes += amrex::aligned_size(acd, nbytes);
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

    plan.m_nrcvs = RcvProc.size();

    plan.m_particle_stats.resize(0);
    plan.m_particle_stats.resize(plan.m_nrcvs);

    plan.m_particle_rreqs.resize(0);
    plan.m_particle_rreqs.resize(plan.m_nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();

    // Post receives.
    for (int i = 0; i < plan.m_nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        Long nbytes       = plan.m_rcv_num_particles[Who]*psize;
        std::size_t acd   = ParallelDescriptor::alignof_comm_data(nbytes);
        const auto Cnt    = amrex::aligned_size(acd, nbytes) / acd;

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(amrex::aligned_size(acd, nbytes) % acd == 0);

        const int comm_data_type = ParallelDescriptor::select_comm_data_type(nbytes);
        if (comm_data_type == 1) {
            plan.m_particle_rreqs[i] =
                ParallelDescriptor::Arecv((char*) (rcv_buffer.dataPtr() + offset), Cnt, Who, SeqNum, ParallelContext::CommunicatorSub()).req();
        } else if (comm_data_type == 2) {
            plan.m_particle_rreqs[i] =
                ParallelDescriptor::Arecv((unsigned long long*) (rcv_buffer.dataPtr() + offset), Cnt, Who, SeqNum, ParallelContext::CommunicatorSub()).req();
        } else if (comm_data_type == 3) {
            plan.m_particle_rreqs[i] =
                ParallelDescriptor::Arecv((ParallelDescriptor::lull_t *) (rcv_buffer.dataPtr() + offset), Cnt, Who, SeqNum, ParallelContext::CommunicatorSub()).req();
        } else {
            amrex::Abort("TODO: message size is too big");
        }
    }

    if (plan.m_NumSnds == 0) return;

    // Send.
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        const auto Who  = i;
        Long nbytes     = plan.m_snd_num_particles[i]*psize;
        std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
        const auto Cnt  = plan.m_snd_counts[i] / acd;
        if (Cnt == 0) continue;

        auto snd_offset = plan.m_snd_offsets[i];
        AMREX_ASSERT(plan.m_snd_counts[i] % acd == 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(snd_offset % acd == 0);

        const int comm_data_type = ParallelDescriptor::select_comm_data_type(nbytes);
        if (comm_data_type == 1) {
            ParallelDescriptor::Send((char*)(snd_buffer.dataPtr()+snd_offset), Cnt, Who, SeqNum,
                                     ParallelContext::CommunicatorSub());
        } else if (comm_data_type == 2) {
            ParallelDescriptor::Send((unsigned long long*)(snd_buffer.dataPtr()+snd_offset), Cnt, Who, SeqNum, ParallelContext::CommunicatorSub());
        } else if (comm_data_type == 3) {
            ParallelDescriptor::Send((ParallelDescriptor::lull_t *)(snd_buffer.dataPtr()+snd_offset), Cnt, Who, SeqNum, ParallelContext::CommunicatorSub());
        } else {
            amrex::Abort("TODO: message size is too big");
        }
    }
#else
    amrex::ignore_unused(pc,plan,snd_buffer,rcv_buffer);
#endif // MPI
}

void communicateParticlesFinish (const ParticleCopyPlan& plan);

template <class PC, class Buffer, class UnpackPolicy,
          EnableIf_t<IsParticleContainer<PC>::value, int> foo = 0>
void unpackRemotes (PC& pc, const ParticleCopyPlan& plan, Buffer& rcv_buffer, UnpackPolicy&& policy)
{
    BL_PROFILE("amrex::unpackRemotes");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs == 1) return;

    const int MyProc = ParallelContext::MyProcSub();
    amrex::ignore_unused(MyProc);
    using PTile = typename PC::ParticleTileType;

    if (plan.m_nrcvs > 0)
    {
        auto p_comm_real = pc.d_communicate_real_comp.dataPtr();
        auto p_comm_int  = pc.d_communicate_int_comp.dataPtr();
        auto p_rcv_buffer = rcv_buffer.dataPtr();

        std::vector<int> sizes;
        std::vector<PTile*> tiles;
        for (int i = 0; i < plan.m_rcv_box_counts.size(); ++i)
          {
            int copy_size = plan.m_rcv_box_counts[i];
            int lev = plan.m_rcv_box_levs[i];
            int gid = plan.m_rcv_box_ids[i];
            int tid = 0;
            auto& tile = pc.DefineAndReturnParticleTile(lev, gid, tid);
            sizes.push_back(copy_size);
            tiles.push_back(&tile);
        }

        Vector<int> offsets;
        policy.resizeTiles(tiles, sizes, offsets);
        Gpu::Device::synchronize();
        int uindex = 0;
        int procindex = 0, rproc = plan.m_rcv_box_pids[0];
        for (int i = 0; i < plan.m_rcv_box_counts.size(); ++i)
          {
            int lev = plan.m_rcv_box_levs[i];
            int gid = plan.m_rcv_box_ids[i];
            int tid = 0;
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

            Long psize = pc.superParticleSize();
            auto p_pad_adjust = plan.m_rcv_pad_correction_d.dataPtr();

            AMREX_FOR_1D ( size, ip, {
                Long src_offset = psize*(offset + ip) + p_pad_adjust[procindex];
                int dst_index = dst_offset + ip;
                ptd.unpackParticleData(p_rcv_buffer, src_offset, dst_index,
                                       p_comm_real, p_comm_int);
              });

            Gpu::synchronize();
          }
    }
#else
    amrex::ignore_unused(pc,plan,rcv_buffer,policy);
#endif // MPI
}

} // namespace amrex

#endif // AMREX_PARTICLECOMMUNICATION_H_
