#include <AMReX_ParticleCommunication.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

void ParticleCopyOp::clear ()
{
    m_boxes.resize(0);
    m_levels.resize(0);
    m_src_indices.resize(0);
    m_periodic_shift.resize(0);
}

void ParticleCopyOp::setNumLevels (const int num_levels)
{
    m_boxes.resize(num_levels);
    m_levels.resize(num_levels);
    m_src_indices.resize(num_levels);
    m_periodic_shift.resize(num_levels);
}

void ParticleCopyOp::resize (const int gid, const int lev, const int size)
{
    if (lev >= m_boxes.size())
    {
        setNumLevels(lev+1);
    }
    m_boxes[lev][gid].resize(size);
    m_levels[lev][gid].resize(size);
    m_src_indices[lev][gid].resize(size);
    m_periodic_shift[lev][gid].resize(size);
}

void ParticleCopyPlan::clear ()
{
    m_dst_indices.clear();
    m_box_counts.clear();
    m_box_offsets.clear();

    m_rcv_box_counts.clear();
    m_rcv_box_offsets.clear();
    m_rcv_box_ids.clear();
}

void ParticleCopyPlan::buildMPIStart (const ParticleBufferMap& map, Long psize)
{
    BL_PROFILE("ParticleCopyPlan::buildMPIStart");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();
    const int MyProc = ParallelContext::MyProcSub();
    const int NNeighborProcs = m_neighbor_procs.size();

    if (NProcs == 1) return;

    m_Snds.resize(0);
    m_Snds.resize(NProcs, 0);

    m_Rcvs.resize(0);
    m_Rcvs.resize(NProcs, 0);

    m_snd_num_particles.resize(0);
    m_snd_num_particles.resize(NProcs, 0);

    m_rcv_num_particles.resize(0);
    m_rcv_num_particles.resize(NProcs, 0);

    Gpu::HostVector<int> box_counts(m_box_counts.size());
    Gpu::copy(Gpu::deviceToHost, m_box_counts.begin(), m_box_counts.end(), box_counts.begin());
    std::map<int, Vector<int> > snd_data;

    m_NumSnds = 0;
    for (auto i : m_neighbor_procs)
    {
        auto box_buffer_indices = map.allBucketsOnProc(i);
        Long nbytes = 0;
        for (auto bucket : box_buffer_indices)
	{
            int dst = map.bucketToGrid(bucket);
            int lev = map.bucketToLevel(bucket);
            int npart = box_counts[bucket];
            if (npart == 0) continue;
            m_snd_num_particles[i] += npart;
            if (i == MyProc) continue;
            snd_data[i].push_back(npart);
            snd_data[i].push_back(dst);
            snd_data[i].push_back(lev);
            snd_data[i].push_back(MyProc);
            nbytes += 4*sizeof(int);
	}
	m_Snds[i] = nbytes;
	m_NumSnds += nbytes;
    }

    doHandShake(m_Snds, m_Rcvs);

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

    if ( (tot_snds_this_proc == 0) and (tot_rcvs_this_proc == 0) )
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
    
    m_nrcvs = m_RcvProc.size();

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
    
    for (auto i : m_neighbor_procs)
    {
        if (i == MyProc) continue;
        const auto Who = i;
        const auto Cnt = m_Snds[i];
        if (Cnt == 0) continue;

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send((char*) snd_data[i].data(), Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
    }

    m_snd_counts.resize(0);
    m_snd_offsets.resize(0);
    m_snd_pad_correction_h.resize(0);
    
    m_snd_offsets.push_back(0);
    m_snd_pad_correction_h.push_back(0);
    for (int i = 0; i < NProcs; ++i)
    {
        Long nbytes = m_snd_num_particles[i]*psize;
        std::size_t acd = ParallelDescriptor::alignof_comm_data(nbytes);
        Long Cnt = amrex::aligned_size(acd, nbytes);
        Long bytes_to_send = (i == MyProc) ? 0 : Cnt;
        m_snd_counts.push_back(bytes_to_send);
        m_snd_offsets.push_back(amrex::aligned_size(acd, m_snd_offsets.back()) + Cnt);
        m_snd_pad_correction_h.push_back(m_snd_pad_correction_h.back() + nbytes);
    }

    for (int i = 0; i < NProcs; ++i)
    {
        m_snd_pad_correction_h[i] = m_snd_offsets[i] - m_snd_pad_correction_h[i];
    }
    
    m_snd_pad_correction_d.resize(m_snd_pad_correction_h.size());
    Gpu::copy(Gpu::hostToDevice, m_snd_pad_correction_h.begin(), m_snd_pad_correction_h.end(),
              m_snd_pad_correction_d.begin());
#else
    amrex::ignore_unused(map,psize);
#endif
}

void ParticleCopyPlan::buildMPIFinish (const ParticleBufferMap& map)
{
    amrex::ignore_unused(map);

    BL_PROFILE("ParticleCopyPlan::buildMPIFinish");

#ifdef AMREX_USE_MPI

    const int NProcs = ParallelContext::NProcsSub();
    if (NProcs == 1) return;

    if (m_nrcvs > 0)
    {
        ParallelDescriptor::Waitall(m_build_rreqs, m_build_stats);

        Gpu::HostVector<int> rcv_box_offsets;
        Gpu::HostVector<int> rcv_box_counts;
        Gpu::HostVector<int> rcv_box_ids;
        Gpu::HostVector<int> rcv_box_levs;
        Gpu::HostVector<int> rcv_box_pids;

        rcv_box_offsets.push_back(0);
        for (int i = 0; i < m_rcv_data.size(); i+=4)
        {
            rcv_box_counts.push_back(m_rcv_data[i]);
            AMREX_ASSERT(ParallelContext::MyProcSub() == map.procID(m_rcv_data[i+1], m_rcv_data[i+2]));
            rcv_box_ids.push_back(m_rcv_data[i+1]);
            rcv_box_levs.push_back(m_rcv_data[i+2]);
            rcv_box_pids.push_back(m_rcv_data[i+3]);
            rcv_box_offsets.push_back(rcv_box_offsets.back() + rcv_box_counts.back());
        }

        m_rcv_box_counts.resize(rcv_box_counts.size());
        Gpu::copy(Gpu::hostToDevice, rcv_box_counts.begin(), rcv_box_counts.end(), m_rcv_box_counts.begin());

        m_rcv_box_offsets.resize(rcv_box_offsets.size());
        Gpu::copy(Gpu::hostToDevice, rcv_box_offsets.begin(), rcv_box_offsets.end(), m_rcv_box_offsets.begin());

        m_rcv_box_ids.resize(rcv_box_ids.size());
        Gpu::copy(Gpu::hostToDevice, rcv_box_ids.begin(), rcv_box_ids.end(), m_rcv_box_ids.begin());

        m_rcv_box_levs.resize(rcv_box_levs.size());
        Gpu::copy(Gpu::hostToDevice, rcv_box_levs.begin(), rcv_box_levs.end(), m_rcv_box_levs.begin());

        m_rcv_box_pids.resize(rcv_box_pids.size());
        Gpu::copy(Gpu::hostToDevice, rcv_box_pids.begin(), rcv_box_pids.end(), m_rcv_box_pids.begin());
    }

    for (int j = 0; j < m_nrcvs; ++j)
    {
        const auto Who    = m_RcvProc[j];
        const auto offset = m_rOffset[j];
        const auto Cnt    = m_Rcvs[Who]/sizeof(int);

        Long nparticles = 0;
        for (int i = offset; i < offset + Cnt; i +=4)
        {
            nparticles += m_rcv_data[i];
        }
        m_rcv_num_particles[Who] = nparticles;
    }
#endif // MPI
}

void ParticleCopyPlan::doHandShake (const Vector<Long>& Snds, Vector<Long>& Rcvs) const
{
    BL_PROFILE("ParticleCopyPlan::doHandShake");
    if (m_local) doHandShakeLocal(Snds, Rcvs);
    else doHandShakeGlobal(Snds, Rcvs);
}

void ParticleCopyPlan::doHandShakeLocal (const Vector<Long>& Snds, Vector<Long>& Rcvs) const
{
#ifdef AMREX_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const int num_rcvs = m_neighbor_procs.size();
    Vector<MPI_Status>  stats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);

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

        ParallelDescriptor::Send(&Snds[Who], Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
    }

    if (num_rcvs > 0)
    {
        ParallelDescriptor::Waitall(rreqs, stats);
    }
#else
    amrex::ignore_unused(Snds,Rcvs);
#endif
}

void ParticleCopyPlan::doHandShakeAllToAll (const Vector<Long>& Snds, Vector<Long>& Rcvs) const
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

void ParticleCopyPlan::doHandShakeGlobal (const Vector<Long>& Snds, Vector<Long>& Rcvs) const
{
#ifdef AMREX_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const int NProcs = ParallelContext::NProcsSub();

    Vector<Long> snd_connectivity(NProcs, 0);
    Vector<int > rcv_connectivity(NProcs, 1);
    for (int i = 0; i < NProcs; ++i) { if (Snds[i] > 0) snd_connectivity[i] = 1; }

    Long num_rcvs = 0;
    MPI_Reduce_scatter(snd_connectivity.data(), &num_rcvs, rcv_connectivity.data(),
                       ParallelDescriptor::Mpi_typemap<Long>::type(), MPI_SUM,
                       ParallelContext::CommunicatorSub());

    Vector<MPI_Status>  stats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);

    Vector<Long> num_bytes_rcv(num_rcvs);
    for (int i = 0; i < num_rcvs; ++i)
    {
        MPI_Irecv( &num_bytes_rcv[i], 1, ParallelDescriptor::Mpi_typemap<Long>::type(),
                   MPI_ANY_SOURCE, SeqNum, ParallelContext::CommunicatorSub(), &rreqs[i] );
    }
    for (int i = 0; i < NProcs; ++i)
    {
        if (Snds[i] == 0) continue;
        const Long Cnt = 1;
        MPI_Send( &Snds[i], Cnt, ParallelDescriptor::Mpi_typemap<Long>::type(), i, SeqNum,
                  ParallelContext::CommunicatorSub());
    }

    MPI_Waitall(num_rcvs, rreqs.data(), stats.data());

    for (int i = 0; i < num_rcvs; ++i)
    {
        const auto Who = stats[i].MPI_SOURCE;
        Rcvs[Who] = num_bytes_rcv[i];
    }
#else
    amrex::ignore_unused(Snds,Rcvs);
#endif
}

void amrex::communicateParticlesFinish (const ParticleCopyPlan& plan)
{
    BL_PROFILE("amrex::communicateParticlesFinish");
#ifdef AMREX_USE_MPI
    if (plan.m_nrcvs > 0)
    {
        ParallelDescriptor::Waitall(plan.m_particle_rreqs, plan.m_particle_stats);
    }
#else
    amrex::ignore_unused(plan);
#endif
}
