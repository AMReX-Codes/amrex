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

void ParticleCopyPlan::buildMPI (const ParticleBufferMap& map)
{
    BL_PROFILE("ParticleCopyPlan::buildMPI");

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    const int MyProc = ParallelDescriptor::MyProc();
    
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0); // number of bytes per snd / receive
    
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
	Snds[i] = nbytes;
	m_NumSnds += nbytes;
    }

    ParallelDescriptor::ReduceLongMax(m_NumSnds);

    if (m_NumSnds == 0) return;

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());
    
    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 ParallelDescriptor::Communicator()) );
    
    BL_ASSERT(Rcvs[ParallelDescriptor::MyProc()] == 0);
    
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset;
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i)
    {
        if (Rcvs[i] > 0)
        {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes/sizeof(int));
            TotRcvBytes += Rcvs[i];
        }
    }
    
    m_nrcvs = RcvProc.size();
    m_stats.resize(m_nrcvs);
    m_rreqs.resize(m_nrcvs);
    const int SeqNum = ParallelDescriptor::SeqNum();
    Gpu::HostVector<int> rcv_data(TotRcvBytes/sizeof(int));

    for (int i = 0; i < m_nrcvs; ++i)
    {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = Rcvs[Who];
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        m_rreqs[i] = ParallelDescriptor::Arecv((char*) (rcv_data.dataPtr() + offset), Cnt, Who, SeqNum).req();
    }
    
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        const auto Who = i;
        const auto Cnt = Snds[i];

        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send((char*) snd_data[i].data(), Cnt, Who, SeqNum);
    }

    if (m_nrcvs > 0)
    {
        ParallelDescriptor::Waitall(m_rreqs, m_stats);

        Gpu::HostVector<int> rcv_box_offsets;
        Gpu::HostVector<int> rcv_box_counts;
        Gpu::HostVector<int> rcv_box_ids;        
        rcv_box_offsets.push_back(0);
        for (int i = 0; i < rcv_data.size(); i +=2)
        {
            rcv_box_counts.push_back(rcv_data[i]);
            AMREX_ASSERT(MyProc == map.dm[rcv_data[i+1]]);
            rcv_box_ids.push_back(rcv_data[i+1]);
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
        const auto Who    = RcvProc[j];
        const auto offset = rOffset[j];
        const auto Cnt    = Rcvs[Who]/sizeof(int);
        
        long nparticles = 0;
        for (int i = offset; i < offset + Cnt; i +=2)
        {
            nparticles += rcv_data[i];
        }
        m_rcv_num_particles[Who] = nparticles;
    }
#endif // MPI
}
