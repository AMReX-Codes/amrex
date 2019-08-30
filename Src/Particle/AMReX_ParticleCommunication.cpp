#include <AMReX_ParticleCommunication.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

void ParticleCopyOp::clear ()
{
    m_boxes.clear();
    m_src_indices.clear();
    m_periodic_shift.clear();
}

void ParticleCopyOp::resize (const int gid, const int size)
{
    m_boxes[gid].resize(size);
    m_src_indices[gid].resize(size);
    m_periodic_shift[gid].resize(size);
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

void ParticleCopyPlan::buildMPIStart (const ParticleBufferMap& map, bool do_handshake)
{
    BL_PROFILE("ParticleCopyPlan::buildMPIStart");

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    const int MyProc = ParallelDescriptor::MyProc();
    const int NNeighborProcs = m_neighbor_procs.size();

    if (NProcs == 1) return;

    if (do_handshake)
    {
        m_Snds.resize(0);
        m_Snds.resize(NProcs, 0);

        m_Rcvs.resize(0);
        m_Rcvs.resize(NProcs, 0);
    }
    
    m_snd_num_particles.resize(0);
    m_snd_num_particles.resize(NProcs, 0);

    m_rcv_num_particles.resize(0);
    m_rcv_num_particles.resize(NProcs, 0);

    Gpu::HostVector<int> box_counts(m_box_counts.size());
    Gpu::thrust_copy(m_box_counts.begin(), m_box_counts.end(), box_counts.begin());
    std::map<int, Vector<int> > snd_data;

    m_NumSnds = 0;
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        auto box_buffer_indices = map.allBucketsOnProc(i);
        int nboxes = box_buffer_indices.size();
        for (auto bucket : box_buffer_indices)
	{
            int dst = map.bucketToGrid(bucket);
            int npart = box_counts[bucket];
            m_snd_num_particles[i] += npart;
            snd_data[i].push_back(npart);
            snd_data[i].push_back(dst);
	}
	long nbytes = 2*nboxes*sizeof(int);
	m_Snds[i] = nbytes;
	m_NumSnds += nbytes;
    }

    if (do_handshake) doHandShake(m_Snds, m_Rcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();
    long tot_snds_this_proc = 0;
    long tot_rcvs_this_proc = 0;
    for (int i = 0; i < NNeighborProcs; ++i)
    {
        tot_snds_this_proc += m_Snds[m_neighbor_procs[i]];
        tot_rcvs_this_proc += m_Rcvs[m_neighbor_procs[i]];
    }
    if ( (tot_snds_this_proc == 0) and (tot_rcvs_this_proc == 0) )
    {
        m_nrcvs = 0;
        return;
    } 

    m_RcvProc.resize(0);
    m_rOffset.resize(0);    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i)
    {
        if (m_Rcvs[i] > 0)
        {
            m_RcvProc.push_back(i);
            m_rOffset.push_back(TotRcvBytes/sizeof(int));
            TotRcvBytes += m_Rcvs[i];
        }
    }
    
    m_nrcvs = m_RcvProc.size();

    m_stats.resize(0);
    m_stats.resize(m_nrcvs);

    m_rreqs.resize(0);
    m_rreqs.resize(m_nrcvs);

    m_rcv_data.resize(TotRcvBytes/sizeof(int));
 
    for (int i = 0; i < m_nrcvs; ++i)
    {
        const auto Who    = m_RcvProc[i];
        const auto offset = m_rOffset[i];
        const auto Cnt    = m_Rcvs[Who];
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        m_rreqs[i] = ParallelDescriptor::Arecv((char*) (m_rcv_data.dataPtr() + offset), Cnt, Who, SeqNum).req();
    }
    
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        const auto Who = i;
        const auto Cnt = m_Snds[i];
        if (Cnt == 0) continue;

        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Asend((char*) snd_data[i].data(), Cnt, Who, SeqNum);
    }
#endif
}

void ParticleCopyPlan::buildMPIFinish (const ParticleBufferMap& map)
{
    BL_PROFILE("ParticleCopyPlan::buildMPIFinish");

#ifdef BL_USE_MPI

    const int NProcs = ParallelDescriptor::NProcs();
    if (NProcs == 1) return;

    if (m_nrcvs > 0)
    {
        ParallelDescriptor::Waitall(m_rreqs, m_stats);

        Gpu::HostVector<int> rcv_box_offsets;
        Gpu::HostVector<int> rcv_box_counts;
        Gpu::HostVector<int> rcv_box_ids;        
        rcv_box_offsets.push_back(0);
        for (int i = 0; i < m_rcv_data.size(); i+=2)
        {
            rcv_box_counts.push_back(m_rcv_data[i]);
            AMREX_ASSERT(ParallelDescriptor::MyProc() == map.procID(m_rcv_data[i+1]));
            rcv_box_ids.push_back(m_rcv_data[i+1]);
            rcv_box_offsets.push_back(rcv_box_offsets.back() + rcv_box_counts.back());
        }
        
        m_rcv_box_counts.resize(rcv_box_counts.size());
        Gpu::thrust_copy(rcv_box_counts.begin(), rcv_box_counts.end(), m_rcv_box_counts.begin());
        
        m_rcv_box_offsets.resize(rcv_box_offsets.size());
        Gpu::thrust_copy(rcv_box_offsets.begin(), rcv_box_offsets.end(), m_rcv_box_offsets.begin());
        
        m_rcv_box_ids.resize(rcv_box_ids.size());
        Gpu::thrust_copy(rcv_box_ids.begin(), rcv_box_ids.end(), m_rcv_box_ids.begin());
    }
    
    for (int j = 0; j < m_nrcvs; ++j)
    {
        const auto Who    = m_RcvProc[j];
        const auto offset = m_rOffset[j];
        const auto Cnt    = m_Rcvs[Who]/sizeof(int);
        
        long nparticles = 0;
        for (int i = offset; i < offset + Cnt; i +=2)
        {
            nparticles += m_rcv_data[i];
        }
        m_rcv_num_particles[Who] = nparticles;
    }
#endif // MPI
}

void ParticleCopyPlan::doHandShake (const Vector<long>& Snds, Vector<long>& Rcvs) const
{
    BL_PROFILE("ParticleCopyPlan::doHandShake");
    if (m_local) doHandShakeLocal(Snds, Rcvs);
    else doHandShakeGlobal(Snds, Rcvs);
}

void ParticleCopyPlan::doHandShakeLocal (const Vector<long>& Snds, Vector<long>& Rcvs) const
{
#ifdef BL_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const int num_rcvs = m_neighbor_procs.size();
    Vector<MPI_Status>  stats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);
    
    // Post receives
    for (int i = 0; i < num_rcvs; ++i)
    {
        const int Who = m_neighbor_procs[i];
        const long Cnt = 1;
        
        BL_ASSERT(Who >= 0 && Who < ParallelDescriptor::NProcs());
        
        rreqs[i] = ParallelDescriptor::Arecv(&Rcvs[Who], Cnt, Who, SeqNum).req();
    }
        
    // Send.
    for (int i = 0; i < num_rcvs; ++i)
    {
        const int Who = m_neighbor_procs[i];
        const long Cnt = 1;
        
        BL_ASSERT(Who >= 0 && Who < ParallelDescriptor::NProcs());
        
        ParallelDescriptor::Send(&Snds[Who], Cnt, Who, SeqNum);        
    }
        
    if (num_rcvs > 0)
    {
        ParallelDescriptor::Waitall(rreqs, stats);
    }
#endif
}

void ParticleCopyPlan::doHandShakeGlobal (const Vector<long>& Snds, Vector<long>& Rcvs) const
{
#ifdef BL_USE_MPI
    const int SeqNum = ParallelDescriptor::SeqNum();
    const int NProcs = ParallelDescriptor::NProcs();

    Vector<long> snd_connectivity(NProcs, 0);
    Vector<int > rcv_connectivity(NProcs, 1);
    for (int i = 0; i < NProcs; ++i) { if (Snds[i] > 0) snd_connectivity[i] = 1; }

    long num_rcvs = 0;
    MPI_Reduce_scatter(snd_connectivity.data(), &num_rcvs, rcv_connectivity.data(), 
                       MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

    Vector<MPI_Status>  stats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);

    Vector<std::size_t> num_bytes_rcv(num_rcvs);
    for (int i = 0; i < num_rcvs; ++i)
        MPI_Irecv( &num_bytes_rcv[i], 1, MPI_UNSIGNED_LONG, MPI_ANY_SOURCE,
                   SeqNum, MPI_COMM_WORLD, &rreqs[i] );
        
    for (int i = 0; i < NProcs; ++i)
    {
        if (Snds[i] == 0) continue;
        const long Cnt = 1;
        ParallelDescriptor::Send(&Snds[i], Cnt, i, SeqNum);
    }

    MPI_Waitall(num_rcvs, rreqs.data(), stats.data());

    for (int i = 0; i < num_rcvs; ++i)
    {
        const auto Who = stats[i].MPI_SOURCE;
        Rcvs[Who] = num_bytes_rcv[i];
    }
#endif
}

void amrex::communicateParticlesFinish (const ParticleCopyPlan& plan)
{
    BL_PROFILE("amrex::communicateParticlesFinish");
#ifdef BL_USE_MPI
    if (plan.m_nrcvs > 0)
    {
        ParallelDescriptor::Waitall(plan.m_rreqs, plan.m_stats);
    }
#endif
}
