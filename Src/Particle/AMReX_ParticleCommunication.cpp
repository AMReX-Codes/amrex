#include <AMReX_ParticleCommunication.H>
#include <AMReX_ParticleContainerBase.H>
#include <AMReX_ParallelDescriptor.H>

#include <iterator>

namespace amrex {

void ParticleCopyOp::clear ()
{
    m_boxes.resize(0);
    m_levels.resize(0);
    m_tiles.resize(0);
    m_src_indices.resize(0);
    m_periodic_shift.resize(0);
}

void ParticleCopyOp::setNumLevels (int num_levels)
{
    m_boxes.resize(num_levels);
    m_levels.resize(num_levels);
    m_tiles.resize(num_levels);
    m_src_indices.resize(num_levels);
    m_periodic_shift.resize(num_levels);
}

void ParticleCopyOp::resize (int gid, int tid, int lev, int size)
{
    if (lev >= std::ssize(m_boxes))
    {
        setNumLevels(lev+1);
    }
    const auto index = std::make_pair(gid, tid);
    m_boxes[lev][index].resize(size);
    m_levels[lev][index].resize(size);
    m_tiles[lev][index].resize(size);
    m_src_indices[lev][index].resize(size);
    m_periodic_shift[lev][index].resize(size, IntVect(AMREX_D_DECL(0,0,0)));
}

void ParticleCopyPlan::clear ()
{
    m_dst_indices.clear();
    m_box_counts_d.clear();
    m_box_counts_h.clear();
    m_box_offsets.clear();

    m_rcv_box_counts.clear();
    m_rcv_box_offsets.clear();
    m_rcv_box_ids.clear();
    m_rcv_box_tids.clear();
    m_rcv_box_pids.clear();
    m_rcv_box_levs.clear();
}

void ParticleCopyPlan::buildMPIStartHandshake (const ParticleContainerBase& pc,
                                               const ParticleBufferMap& map)
{
    BL_PROFILE("ParticleCopyPlan::buildMPIStartHandshake");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    const int MyProc = ParallelContext::MyProcSub();

    if (NProcs == 1) { return; }

    m_Snds.assign(NProcs, 0);
    m_Rcvs.assign(NProcs, 0);
    m_snd_num_particles.assign(NProcs, 0);
    m_rcv_num_particles.assign(NProcs, 0);
    m_snd_data.clear();
    m_NumSnds = 0;

    for (auto i : m_neighbor_procs)
    {
        auto box_buffer_indices = map.allBucketsOnProc(i);
        Long nbytes = 0;
        for (auto bucket : box_buffer_indices)
        {
            int dst   = map.bucketToGrid(bucket);
            int tid   = map.bucketToTile(bucket);
            int lev   = map.bucketToLevel(bucket);
            AMREX_ASSERT(m_box_counts_h[bucket] <= static_cast<unsigned int>(std::numeric_limits<int>::max()));
            int npart = static_cast<int>(m_box_counts_h[bucket]);
            if (npart == 0) { continue; }
            m_snd_num_particles[i] += npart;
            if (i == MyProc) { continue; }
            m_snd_data[i].push_back(npart);
            m_snd_data[i].push_back(dst);
            m_snd_data[i].push_back(tid);
            m_snd_data[i].push_back(lev);
            m_snd_data[i].push_back(MyProc);
            nbytes += 5*sizeof(int);
        }
        m_Snds[i] = nbytes;
        m_NumSnds += nbytes;
    }

    doHandShakeStart(pc, m_Snds);
#else
    amrex::ignore_unused(pc, map);
#endif
}

void ParticleCopyPlan::buildMPIStartRest (const ParticleContainerBase& pc,
                                          const ParticleBufferMap& /*map*/, Long psize)
{
    BL_PROFILE("ParticleCopyPlan::buildMPIStartRest");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    const int MyProc = ParallelContext::MyProcSub();
    const auto NNeighborProcs = static_cast<int>(m_neighbor_procs.size());

    if (NProcs == 1) { return; }

    doHandShakeFinish(pc, m_Snds, m_Rcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();
    Long tot_snds_this_proc = 0;
    Long tot_rcvs_this_proc = 0;

    if (m_local)
    {
        for (int i = 0; i < NNeighborProcs; ++i)
        {
            tot_snds_this_proc += m_Snds[m_neighbor_procs[i]];
            tot_rcvs_this_proc += m_Rcvs[m_neighbor_procs[i]];
        }
    } else {
        for (int i = 0; i < NProcs; ++i)
        {
            tot_snds_this_proc += m_Snds[i];
            tot_rcvs_this_proc += m_Rcvs[i];
        }
    }

    if ( (tot_snds_this_proc == 0) && (tot_rcvs_this_proc == 0) )
    {
        m_nrcvs = 0;
        m_NumSnds = 0;
        return;
    }

    m_RcvProc.resize(0);
    m_rOffset.resize(0);
    std::size_t TotRcvBytes = 0;
    for (auto i : m_neighbor_procs)
    {
        if (m_Rcvs[i] > 0)
        {
            m_RcvProc.push_back(i);
            m_rOffset.push_back(TotRcvBytes/sizeof(int));
            TotRcvBytes += m_Rcvs[i];
        }
    }

    m_nrcvs = static_cast<int>(m_RcvProc.size());

    m_build_stats.assign(m_nrcvs, MPI_Status{});
    m_build_rreqs.assign(m_nrcvs, MPI_REQUEST_NULL);
    m_rcv_data.resize(TotRcvBytes/sizeof(int));

    for (int i = 0; i < m_nrcvs; ++i)
    {
        const auto Who    = m_RcvProc[i];
        const auto offset = m_rOffset[i];
        const auto Cnt    = m_Rcvs[Who];

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        m_build_rreqs[i] = ParallelDescriptor::Arecv(
            (char*)(m_rcv_data.dataPtr() + offset), Cnt, Who, SeqNum,
            ParallelContext::CommunicatorSub()).req();
    }

    Vector<MPI_Request> snd_reqs;
    Vector<MPI_Status>  snd_stats;
    for (auto i : m_neighbor_procs)
    {
        if (i == MyProc) { continue; }
        const auto Cnt = m_Snds[i];
        if (Cnt == 0) { continue; }

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(i >= 0 && i < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        snd_reqs.push_back(ParallelDescriptor::Asend(
            (char*) m_snd_data[i].data(), Cnt, i, SeqNum,
            ParallelContext::CommunicatorSub()).req());
    }

    m_snd_counts.resize(0);
    m_snd_offsets.resize(0);
    m_snd_pad_correction_h.resize(0);

    m_snd_offsets.push_back(0);
    m_snd_pad_correction_h.push_back(0);
    for (int i = 0; i < NProcs; ++i)
    {
        Long nbytes = m_snd_num_particles[i]*psize;
        std::size_t acd = ParallelDescriptor::sizeof_selected_comm_data_type(nbytes);
        auto Cnt = static_cast<Long>(amrex::aligned_size(acd, nbytes));
        Long bytes_to_send = (i == MyProc) ? 0 : Cnt;
        m_snd_counts.push_back(bytes_to_send);
        m_snd_offsets.push_back(amrex::aligned_size(Arena::align_size,
                                                    m_snd_offsets.back() + Cnt));
        m_snd_pad_correction_h.push_back(m_snd_pad_correction_h.back() + nbytes);
    }

    for (int i = 0; i < NProcs; ++i)
    {
        m_snd_pad_correction_h[i] = m_snd_offsets[i] - m_snd_pad_correction_h[i];
    }

    m_snd_pad_correction_d.resize(m_snd_pad_correction_h.size());
    Gpu::copy(Gpu::hostToDevice, m_snd_pad_correction_h.begin(), m_snd_pad_correction_h.end(),
              m_snd_pad_correction_d.begin());

    snd_stats.resize(snd_reqs.size());
    ParallelDescriptor::Waitall(snd_reqs, snd_stats);
    m_snd_data.clear();
#else
    amrex::ignore_unused(pc, psize);
#endif
}

void ParticleCopyPlan::buildMPIStart (const ParticleContainerBase& pc, const ParticleBufferMap& map, Long psize) // NOLINT(readability-convert-member-functions-to-static)
{
    BL_PROFILE("ParticleCopyPlan::buildMPIStart");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    const int MyProc = ParallelContext::MyProcSub();
    const auto NNeighborProcs = static_cast<int>(m_neighbor_procs.size());

    if (NProcs == 1) { return; }

    m_Snds.resize(0);
    m_Snds.resize(NProcs, 0);

    m_Rcvs.resize(0);
    m_Rcvs.resize(NProcs, 0);

    m_snd_num_particles.resize(0);
    m_snd_num_particles.resize(NProcs, 0);

    m_rcv_num_particles.resize(0);
    m_rcv_num_particles.resize(NProcs, 0);

    std::map<int, Vector<int> > snd_data;

    m_NumSnds = 0;
    for (auto i : m_neighbor_procs)
    {
        auto box_buffer_indices = map.allBucketsOnProc(i);
        Long nbytes = 0;
        for (auto bucket : box_buffer_indices)
        {
            int dst = map.bucketToGrid(bucket);
            int tid = map.bucketToTile(bucket);
            int lev = map.bucketToLevel(bucket);
            AMREX_ASSERT(m_box_counts_h[bucket] <= static_cast<unsigned int>(std::numeric_limits<int>::max()));
            int npart = static_cast<int>(m_box_counts_h[bucket]);
            if (npart == 0) { continue; }
            m_snd_num_particles[i] += npart;
            if (i == MyProc) { continue; }
            snd_data[i].push_back(npart);
            snd_data[i].push_back(dst);
            snd_data[i].push_back(tid);
            snd_data[i].push_back(lev);
            snd_data[i].push_back(MyProc);
            nbytes += 5*sizeof(int);
        }
        m_Snds[i] = nbytes;
        m_NumSnds += nbytes;
    }

    doHandShake(pc, m_Snds, m_Rcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();
    Long tot_snds_this_proc = 0;
    Long tot_rcvs_this_proc = 0;

    if (m_local)
    {
        for (int i = 0; i < NNeighborProcs; ++i)
        {
            tot_snds_this_proc += m_Snds[m_neighbor_procs[i]];
            tot_rcvs_this_proc += m_Rcvs[m_neighbor_procs[i]];
        }
    } else {
        for (int i = 0; i < NProcs; ++i)
        {
            tot_snds_this_proc += m_Snds[i];
            tot_rcvs_this_proc += m_Rcvs[i];
        }
    }

    if ( (tot_snds_this_proc == 0) && (tot_rcvs_this_proc == 0) )
    {
        m_nrcvs = 0;
        m_NumSnds = 0;
        return;
    }

    m_RcvProc.resize(0);
    m_rOffset.resize(0);
    std::size_t TotRcvBytes = 0;
    for (auto i : m_neighbor_procs)
    {
        if (m_Rcvs[i] > 0)
        {
            m_RcvProc.push_back(i);
            m_rOffset.push_back(TotRcvBytes/sizeof(int));
            TotRcvBytes += m_Rcvs[i];
        }
    }

    m_nrcvs = static_cast<int>(m_RcvProc.size());

    m_build_stats.resize(0);
    m_build_stats.resize(m_nrcvs);

    m_build_rreqs.resize(0);
    m_build_rreqs.resize(m_nrcvs);

    m_rcv_data.resize(TotRcvBytes/sizeof(int));

    for (int i = 0; i < m_nrcvs; ++i)
    {
        const auto Who    = m_RcvProc[i];
        const auto offset = m_rOffset[i];
        const auto Cnt    = m_Rcvs[Who];

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        m_build_rreqs[i] = ParallelDescriptor::Arecv((char*) (m_rcv_data.dataPtr() + offset), Cnt, Who, SeqNum, ParallelContext::CommunicatorSub()).req();
    }

    Vector<MPI_Request> snd_reqs;
    Vector<MPI_Status>  snd_stats;
    for (auto i : m_neighbor_procs)
    {
        if (i == MyProc) { continue; }
        const auto Who = i;
        const auto Cnt = m_Snds[i];
        if (Cnt == 0) { continue; }

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        snd_reqs.push_back(ParallelDescriptor::Asend((char*) snd_data[i].data(), Cnt, Who, SeqNum,
                                                      ParallelContext::CommunicatorSub()).req());
    }

    m_snd_counts.resize(0);
    m_snd_offsets.resize(0);
    m_snd_pad_correction_h.resize(0);

    m_snd_offsets.push_back(0);
    m_snd_pad_correction_h.push_back(0);
    for (int i = 0; i < NProcs; ++i)
    {
        Long nbytes = m_snd_num_particles[i]*psize;
        std::size_t acd = ParallelDescriptor::sizeof_selected_comm_data_type(nbytes);
        auto Cnt = static_cast<Long>(amrex::aligned_size(acd, nbytes));
        Long bytes_to_send = (i == MyProc) ? 0 : Cnt;
        m_snd_counts.push_back(bytes_to_send);
        m_snd_offsets.push_back(amrex::aligned_size(Arena::align_size,
                                                    m_snd_offsets.back() + Cnt));
        m_snd_pad_correction_h.push_back(m_snd_pad_correction_h.back() + nbytes);
    }

    for (int i = 0; i < NProcs; ++i)
    {
        m_snd_pad_correction_h[i] = m_snd_offsets[i] - m_snd_pad_correction_h[i];
    }

    m_snd_pad_correction_d.resize(m_snd_pad_correction_h.size());
    Gpu::copy(Gpu::hostToDevice, m_snd_pad_correction_h.begin(), m_snd_pad_correction_h.end(),
              m_snd_pad_correction_d.begin());

    snd_stats.resize(0);
    snd_stats.resize(snd_reqs.size());
    ParallelDescriptor::Waitall(snd_reqs, snd_stats);
#else
    amrex::ignore_unused(pc,map,psize);
#endif
}

void ParticleCopyPlan::buildMPIFinish (const ParticleBufferMap& map) // NOLINT(readability-convert-member-functions-to-static)
{
    amrex::ignore_unused(map);

    BL_PROFILE("ParticleCopyPlan::buildMPIFinish");

#ifdef AMREX_USE_MPI

    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs == 1) { return; }

    if (m_nrcvs > 0)
    {
        ParallelDescriptor::Waitall(m_build_rreqs, m_build_stats);

        m_rcv_box_offsets.resize(0);
        m_rcv_box_counts.resize(0);
        m_rcv_box_ids.resize(0);
        m_rcv_box_tids.resize(0);
        m_rcv_box_levs.resize(0);
        m_rcv_box_pids.resize(0);

        m_rcv_box_offsets.push_back(0);
        for (int i = 0, N = static_cast<int>(m_rcv_data.size()); i < N; i+=5)
        {
            m_rcv_box_counts.push_back(m_rcv_data[i]);
            AMREX_ASSERT(ParallelContext::MyProcSub() == map.procID(m_rcv_data[i+1], m_rcv_data[i+2], m_rcv_data[i+3]));
            m_rcv_box_ids.push_back(m_rcv_data[i+1]);
            m_rcv_box_tids.push_back(m_rcv_data[i+2]);
            m_rcv_box_levs.push_back(m_rcv_data[i+3]);
            m_rcv_box_pids.push_back(m_rcv_data[i+4]);
            m_rcv_box_offsets.push_back(m_rcv_box_offsets.back() + m_rcv_box_counts.back());
        }
    }

    for (int j = 0; j < m_nrcvs; ++j)
    {
        const auto Who    = m_RcvProc[j];
        const auto offset = m_rOffset[j];
        const auto Cnt    = m_Rcvs[Who]/sizeof(int);

        Long nparticles = 0;
        for (auto i = offset; i < offset + Cnt; i +=5)
        {
            nparticles += m_rcv_data[i];
        }
        m_rcv_num_particles[Who] = nparticles;
    }
#endif // MPI
}

void ParticleCopyPlan::doHandShake (const ParticleContainerBase& pc,
                                    const Vector<Long>& Snds,
                                    Vector<Long>& Rcvs) const // NOLINT(readability-convert-member-functions-to-static)
{
    BL_PROFILE("ParticleCopyPlan::doHandShake");
    if (m_local) { doHandShakeLocal(Snds, Rcvs); }
    else if (m_do_one_sided_comms) {
        doHandShakeOneSided(pc, Snds, Rcvs);
    }
    else {
        doHandShakeReduceScatter(Snds, Rcvs);
    }
}

void ParticleCopyPlan::doHandShakeStart (const ParticleContainerBase& pc,
                                         const Vector<Long>& Snds)
{
    BL_PROFILE("ParticleCopyPlan::doHandShakeStart");
    if (m_local) { doHandShakeLocalStart(Snds); }
    else if (m_do_one_sided_comms) { doHandShakeOneSidedStart(pc, Snds); }
    else { doHandShakeReduceScatterStart(Snds); }
}

void ParticleCopyPlan::doHandShakeFinish (const ParticleContainerBase& pc,
                                          const Vector<Long>& Snds,
                                          Vector<Long>& Rcvs)
{
    BL_PROFILE("ParticleCopyPlan::doHandShakeFinish");
    if (m_local) { doHandShakeLocalFinish(Rcvs); }
    else if (m_do_one_sided_comms) { doHandShakeOneSidedFinish(pc, Rcvs); }
    else { doHandShakeReduceScatterFinish(Snds, Rcvs); }
}

void ParticleCopyPlan::doHandShakeLocal (const Vector<Long>& Snds, Vector<Long>& Rcvs) const // NOLINT(readability-convert-member-functions-to-static)
{
#ifdef AMREX_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const auto num_rcvs = static_cast<int>(m_neighbor_procs.size());
    Vector<MPI_Status>  rstats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);
    Vector<MPI_Status>  sstats(num_rcvs);
    Vector<MPI_Request> sreqs(num_rcvs);

    // Post receives
    for (int i = 0; i < num_rcvs; ++i)
    {
        const int Who = m_neighbor_procs[i];
        const Long Cnt = 1;

        AMREX_ASSERT(Who >= 0 && Who < ParallelContext::NProcsSub());

        rreqs[i] = ParallelDescriptor::Arecv(&Rcvs[Who], Cnt, Who, SeqNum,
                                             ParallelContext::CommunicatorSub()).req();
    }

    // Send.
    for (int i = 0; i < num_rcvs; ++i)
    {
        const int Who = m_neighbor_procs[i];
        const Long Cnt = 1;

        AMREX_ASSERT(Who >= 0 && Who < ParallelContext::NProcsSub());

        sreqs[i] = ParallelDescriptor::Asend(&Snds[Who], Cnt, Who, SeqNum,
                                             ParallelContext::CommunicatorSub()).req();
    }

    if (num_rcvs > 0)
    {
        ParallelDescriptor::Waitall(sreqs, sstats);
        ParallelDescriptor::Waitall(rreqs, rstats);
    }
#else
    amrex::ignore_unused(Snds,Rcvs);
#endif
}

void ParticleCopyPlan::doHandShakeLocalStart (const Vector<Long>& Snds)
{
#ifdef AMREX_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const auto num_neighbors = static_cast<int>(m_neighbor_procs.size());

    m_hs_rreqs.resize(num_neighbors);
    m_hs_sreqs.resize(num_neighbors);
    m_hs_rstats.resize(num_neighbors);
    m_hs_sstats.resize(num_neighbors);

    for (int i = 0; i < num_neighbors; ++i)
    {
        const int Who = m_neighbor_procs[i];
        AMREX_ASSERT(Who >= 0 && Who < ParallelContext::NProcsSub());
        // receives write directly into m_Rcvs, which is a member
        m_hs_rreqs[i] = ParallelDescriptor::Arecv(&m_Rcvs[Who], 1, Who, SeqNum,
                                                    ParallelContext::CommunicatorSub()).req();
    }

    for (int i = 0; i < num_neighbors; ++i)
    {
        const int Who = m_neighbor_procs[i];
        AMREX_ASSERT(Who >= 0 && Who < ParallelContext::NProcsSub());
        m_hs_sreqs[i] = ParallelDescriptor::Asend(&Snds[Who], 1, Who, SeqNum,
                                                    ParallelContext::CommunicatorSub()).req();
    }
#else
    amrex::ignore_unused(Snds);
#endif
}

void ParticleCopyPlan::doHandShakeLocalFinish (Vector<Long>& /*Rcvs*/)
{
#ifdef AMREX_USE_MPI
    // Results land directly in m_Rcvs via the Arecv in doHandShakeLocalStart;
    // the Rcvs argument is the same vector and needs no separate copy.
    if (!m_hs_rreqs.empty())
    {
        ParallelDescriptor::Waitall(m_hs_sreqs, m_hs_sstats);
        ParallelDescriptor::Waitall(m_hs_rreqs, m_hs_rstats);
    }
    m_hs_rreqs.clear();
    m_hs_sreqs.clear();
    m_hs_rstats.clear();
    m_hs_sstats.clear();
#endif
}

void ParticleCopyPlan::doHandShakeAllToAll (const Vector<Long>& Snds, Vector<Long>& Rcvs)
{
#ifdef AMREX_USE_MPI
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(Long),
                    ParallelContext::MyProcSub(), BLProfiler::BeforeCall());

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<Long>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<Long>::type(),
                                 ParallelContext::CommunicatorSub()) );

    AMREX_ASSERT(Rcvs[ParallelContext::MyProcSub()] == 0);

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(Long),
                    ParallelContext::MyProcSub(), BLProfiler::AfterCall());
#else
    amrex::ignore_unused(Snds,Rcvs);
#endif
}

void ParticleCopyPlan::doHandShakeReduceScatter (const Vector<Long>& Snds, Vector<Long>& Rcvs)
{
#ifdef AMREX_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const int NProcs = ParallelContext::NProcsSub();

    Vector<Long> snd_connectivity(NProcs, 0);
    Vector<int > rcv_connectivity(NProcs, 1);
    for (int i = 0; i < NProcs; ++i) { if (Snds[i] > 0) { snd_connectivity[i] = 1; } }

    Long num_rcvs = 0;
    MPI_Reduce_scatter(snd_connectivity.data(), &num_rcvs, rcv_connectivity.data(),
                       ParallelDescriptor::Mpi_typemap<Long>::type(), MPI_SUM,
                       ParallelContext::CommunicatorSub());

    Vector<MPI_Status>  rstats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);
    Vector<MPI_Status>  sstats;
    Vector<MPI_Request> sreqs;

    Vector<Long> num_bytes_rcv(num_rcvs);
    for (int i = 0; i < static_cast<int>(num_rcvs); ++i)
    {
        BL_MPI_REQUIRE(MPI_Irecv( &num_bytes_rcv[i], 1, ParallelDescriptor::Mpi_typemap<Long>::type(),
                                  MPI_ANY_SOURCE, SeqNum, ParallelContext::CommunicatorSub(), &rreqs[i] ));
    }
    for (int i = 0; i < NProcs; ++i)
    {
        if (Snds[i] == 0) { continue; }
        const Long Cnt = 1;
        sreqs.push_back(ParallelDescriptor::Asend( &Snds[i], Cnt, i, SeqNum, ParallelContext::CommunicatorSub()).req());
    }

    sstats.resize(0);
    sstats.resize(sreqs.size());
    ParallelDescriptor::Waitall(sreqs, sstats);
    ParallelDescriptor::Waitall(rreqs, rstats);

    for (int i = 0; i < num_rcvs; ++i)
    {
        const auto Who = rstats[i].MPI_SOURCE;
        Rcvs[Who] = num_bytes_rcv[i];
    }
#else
    amrex::ignore_unused(Snds,Rcvs);
#endif
}

void ParticleCopyPlan::doHandShakeReduceScatterStart (const Vector<Long>& Snds)
{
#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();

    m_snd_connectivity.assign(NProcs, 0);
    m_rcv_connectivity.assign(NProcs, 1);
    for (int i = 0; i < NProcs; ++i) { if (Snds[i] > 0) { m_snd_connectivity[i] = 1; } }

    m_num_rcvs_rs = 0;
    BL_MPI_REQUIRE(MPI_Ireduce_scatter(
        m_snd_connectivity.data(), &m_num_rcvs_rs,
        m_rcv_connectivity.data(),
        ParallelDescriptor::Mpi_typemap<Long>::type(), MPI_SUM,
        ParallelContext::CommunicatorSub(), &m_reduce_scatter_req));
#else
    amrex::ignore_unused(Snds);
#endif
}

void ParticleCopyPlan::doHandShakeReduceScatterFinish (const Vector<Long>& Snds,
                                                       Vector<Long>& Rcvs)
{
#ifdef AMREX_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const int NProcs = ParallelContext::NProcsSub();

    BL_MPI_REQUIRE(MPI_Wait(&m_reduce_scatter_req, MPI_STATUS_IGNORE));
    m_snd_connectivity.clear();
    m_rcv_connectivity.clear();

    Long num_rcvs = m_num_rcvs_rs;
    Vector<Long>        num_bytes_rcv(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);
    Vector<MPI_Status>  rstats(num_rcvs);
    Vector<MPI_Request> sreqs;
    Vector<MPI_Status>  sstats;

    for (int i = 0; i < static_cast<int>(num_rcvs); ++i)
    {
        BL_MPI_REQUIRE(MPI_Irecv(&num_bytes_rcv[i], 1,
                                  ParallelDescriptor::Mpi_typemap<Long>::type(),
                                  MPI_ANY_SOURCE, SeqNum,
                                  ParallelContext::CommunicatorSub(), &rreqs[i]));
    }
    for (int i = 0; i < NProcs; ++i)
    {
        if (Snds[i] == 0) { continue; }
        sreqs.push_back(ParallelDescriptor::Asend(&Snds[i], 1, i, SeqNum,
                                                   ParallelContext::CommunicatorSub()).req());
    }

    sstats.resize(sreqs.size());
    ParallelDescriptor::Waitall(sreqs, sstats);
    ParallelDescriptor::Waitall(rreqs, rstats);

    for (int i = 0; i < num_rcvs; ++i)
    {
        Rcvs[rstats[i].MPI_SOURCE] = num_bytes_rcv[i];
    }
#else
    amrex::ignore_unused(Snds,Rcvs);
#endif
}

void ParticleCopyPlan::doHandShakeOneSided (const ParticleContainerBase& pc,
                                            const Vector<Long>& Snds,
                                            Vector<Long>& Rcvs)
{
#if defined(AMREX_USE_MPI)
    const int MyProc = ParallelContext::MyProcSub();
    const int NProcs = ParallelContext::NProcsSub();

    AMREX_ALWAYS_ASSERT(std::ssize(Snds) == NProcs);
    AMREX_ALWAYS_ASSERT(std::ssize(Rcvs) == NProcs);

    pc.ensureParticleHandshakeWindow();
    auto* handshake_buffer = pc.particleHandshakeBuffer();
    AMREX_ALWAYS_ASSERT(handshake_buffer != nullptr);
    std::fill_n(handshake_buffer, NProcs, Long(0));

    MPI_Win win = pc.particleHandshakeWindow();
    BL_MPI_REQUIRE(MPI_Win_fence(0, win));

    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc || Snds[i] == 0) { continue; }

        BL_MPI_REQUIRE(MPI_Put(&Snds[i],
                               1,
                               ParallelDescriptor::Mpi_typemap<Long>::type(),
                               i,
                               MyProc,
                               1,
                               ParallelDescriptor::Mpi_typemap<Long>::type(),
                               win));
    }

    BL_MPI_REQUIRE(MPI_Win_fence(0, win));
    std::copy_n(handshake_buffer, NProcs, Rcvs.begin());

    AMREX_ASSERT(Rcvs[MyProc] == 0);
#else
    amrex::ignore_unused(pc,Snds,Rcvs);
#endif
}

void ParticleCopyPlan::doHandShakeOneSidedStart (const ParticleContainerBase& pc,
                                                  const Vector<Long>& Snds)
{
#if defined(AMREX_USE_MPI)
    const int MyProc = ParallelContext::MyProcSub();
    const int NProcs = ParallelContext::NProcsSub();

    pc.ensureParticleHandshakeWindow();
    auto* handshake_buffer = pc.particleHandshakeBuffer();
    AMREX_ALWAYS_ASSERT(handshake_buffer != nullptr);
    std::fill_n(handshake_buffer, NProcs, Long(0));

    MPI_Win win = pc.particleHandshakeWindow();
    BL_MPI_REQUIRE(MPI_Win_fence(0, win));  // open RMA epoch

    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc || Snds[i] == 0) { continue; }
        BL_MPI_REQUIRE(MPI_Put(&Snds[i], 1,
                                ParallelDescriptor::Mpi_typemap<Long>::type(),
                                i, MyProc, 1,
                                ParallelDescriptor::Mpi_typemap<Long>::type(),
                                win));
    }
    // epoch left open; puts are in flight until doHandShakeOneSidedFinish closes it
#else
    amrex::ignore_unused(pc, Snds);
#endif
}

void ParticleCopyPlan::doHandShakeOneSidedFinish (const ParticleContainerBase& pc,
                                                   Vector<Long>& Rcvs)
{
#if defined(AMREX_USE_MPI)
    const int NProcs = ParallelContext::NProcsSub();
    MPI_Win win = pc.particleHandshakeWindow();
    BL_MPI_REQUIRE(MPI_Win_fence(0, win));  // close epoch; all puts complete
    std::copy_n(pc.particleHandshakeBuffer(), NProcs, Rcvs.begin());
    AMREX_ASSERT(Rcvs[ParallelContext::MyProcSub()] == 0);
#else
    amrex::ignore_unused(pc, Rcvs);
#endif
}

void communicateParticlesFinish (const ParticleCopyPlan& plan)
{
    BL_PROFILE("amrex::communicateParticlesFinish");
#ifdef AMREX_USE_MPI
    if (plan.m_NumSnds > 0)
    {
        ParallelDescriptor::Waitall(plan.m_particle_sreqs, plan.m_particle_sstats);
    }
    if (plan.m_nrcvs > 0)
    {
        ParallelDescriptor::Waitall(plan.m_particle_rreqs, plan.m_particle_rstats);
    }
#else
    amrex::ignore_unused(plan);
#endif
}

}
