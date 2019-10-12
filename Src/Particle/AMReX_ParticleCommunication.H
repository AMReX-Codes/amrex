#ifndef AMREX_PARTICLECOMMUNICATION_H_
#define AMREX_PARTICLECOMMUNICATION_H_

#include <AMReX_Gpu.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_MFIter.H>

#include <map>

namespace amrex {

struct NeighborUnpackPolicy
{
    template <class PTile> 
    void resizeTiles(std::vector<PTile*>& tiles, const std::vector<int>& sizes, std::vector<int>& offsets) const
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
    void resizeTiles(std::vector<PTile*>& tiles, const std::vector<int>& sizes, std::vector<int>& offsets) const
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
    std::map<int, Gpu::ManagedDeviceVector<int> > m_boxes;
    std::map<int, Gpu::ManagedDeviceVector<int> > m_src_indices;
    std::map<int, Gpu::ManagedDeviceVector<IntVect> > m_periodic_shift;

    void clear ();

    void resize (const int gid, const int size);

    int numCopies (const int gid) const { return m_boxes.at(gid).size(); }
};

struct ParticleCopyPlan 
{
    std::map<int, Gpu::ManagedDeviceVector<int> > m_dst_indices;
    
    Gpu::ManagedDeviceVector<unsigned int> m_box_counts;
    Gpu::ManagedDeviceVector<unsigned int> m_box_offsets;

    Gpu::ManagedDeviceVector<int> m_rcv_box_counts;
    Gpu::ManagedDeviceVector<int> m_rcv_box_offsets;
    Gpu::ManagedDeviceVector<int> m_rcv_box_ids;

    long m_NumSnds;
    int m_nrcvs;
    mutable Vector<MPI_Status> m_build_stats;
    mutable Vector<MPI_Request> m_build_rreqs;

    mutable Vector<MPI_Status> m_particle_stats;
    mutable Vector<MPI_Request> m_particle_rreqs;

    Vector<long> m_snd_num_particles;
    Vector<long> m_rcv_num_particles;

    Vector<int> m_neighbor_procs;
        
    Vector<long> m_Snds;
    Vector<long> m_Rcvs;
    Vector<int> m_RcvProc;
    Vector<std::size_t> m_rOffset;
    Gpu::HostVector<int> m_rcv_data;

    template <class PC>
    void build (const PC& pc, const ParticleCopyOp& op, bool do_handshake)
    {
        BL_PROFILE("ParticleCopyPlan::build");
        
        AMREX_ASSERT(pc.numLevels() == 1);
        
        const int lev = 0;
        const int ngrow = 1;
        const auto& geom = pc.Geom(lev);
                
        m_neighbor_procs.assign(pc.NeighborProcs(lev, ngrow).begin(), pc.NeighborProcs(lev, ngrow).end());

        int num_boxes = pc.ParticleBoxArray(lev).size();
        if (num_boxes == 1 and (not geom.isAnyPeriodic()) ) return;
        
        m_box_counts.resize(0);
        m_box_counts.resize(num_boxes, 0);
        m_box_offsets.resize(num_boxes+1);
        auto p_dst_box_counts = m_box_counts.dataPtr();
        auto p_box_perm = pc.BufferMap().gridToBucketPtr();
        
        constexpr unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();

        for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            int num_copies = op.numCopies(gid);
            m_dst_indices[gid].resize(num_copies);
        }
        
        for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int gid = mfi.index();
            int num_copies = op.numCopies(gid);
            
            auto p_boxes = op.m_boxes.at(gid).dataPtr();
            auto p_dst_indices = m_dst_indices[gid].dataPtr();
            
            AMREX_FOR_1D ( num_copies, i,
            {
                int dst_box = p_boxes[i];
                if (dst_box < 0) continue;
                int index = Gpu::Atomic::Inc(&p_dst_box_counts[p_box_perm[dst_box]], max_unsigned_int);
                p_dst_indices[i] = index;
            });
        }

        amrex::Gpu::exclusive_scan(m_box_counts.begin(), m_box_counts.end()+1, m_box_offsets.begin());

        buildMPIStart(pc.BufferMap(), do_handshake);
    }

    void clear();

    void buildMPIFinish (const ParticleBufferMap& map);

private:
    
    void buildMPIStart (const ParticleBufferMap& map, bool do_handshake=true);

    //
    // Snds - a Vector with the number of bytes that is process will send to each proc.
    // Rcvs - a Vector that, after calling this method, will contain the
    //        number of bytes this process will reveive from each proc.
    //
    void doHandShake (const Vector<long>& Snds, Vector<long>& Rcvs) const;

    //
    // In the local version of this method, each proc knows which other
    // procs it could possibly receive messages from, meaning we can do
    // this purely with point-to-point communication.
    //
    void doHandShakeLocal (const Vector<long>& Snds, Vector<long>& Rcvs) const;
    
    //
    // In the global version, we don't know who we'll receive from, so we
    // need to do some collective communication first.
    //
    void doHandShakeGlobal (const Vector<long>& Snds, Vector<long>& Rcvs) const;

    bool m_local = true;
};

template <class PC, class Buffer>
void packBuffer (const PC& pc, const ParticleCopyOp& op, const ParticleCopyPlan& plan, Buffer& snd_buffer)
{
    BL_PROFILE("amrex::packBuffer");

    AMREX_ASSERT(pc.numLevels() == 1);

    const int lev = 0;
    const auto& geom = pc.Geom(lev);
    auto& plev = pc.GetParticles(lev);
    auto& ba = pc.ParticleBoxArray(lev);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto is_per = geom.isPeriodicArray();

    if (ba.size() == 1 and (not geom.isAnyPeriodic()) ) return;
    
    int total_buffer_size = plan.m_box_offsets[ba.size()];
    snd_buffer.resize(total_buffer_size);

    for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& src_tile = plev.at(index);
        auto& aos   = src_tile.GetArrayOfStructs();
        const auto ptd = src_tile.getConstParticleTileData();

	int num_copies = op.numCopies(gid);

	auto p_boxes = op.m_boxes.at(gid).dataPtr();
	auto p_src_indices = op.m_src_indices.at(gid).dataPtr();
        auto p_dst_indices = plan.m_dst_indices.at(gid).dataPtr();
	auto p_periodic_shift = op.m_periodic_shift.at(gid).dataPtr();
	auto p_box_offsets = plan.m_box_offsets.dataPtr();
        auto p_box_perm = pc.BufferMap().gridToBucketPtr();
	auto p_snd_buffer = snd_buffer.dataPtr();

	AMREX_FOR_1D ( num_copies, i,
        {
            int dst_box = p_boxes[i];
            if (dst_box < 0) continue;
            int box_offset = p_box_offsets[p_box_perm[dst_box]];
            int dst_index = box_offset + p_dst_indices[i];
            int src_index = p_src_indices[i];
            p_snd_buffer[dst_index] = ptd.getSuperParticle(src_index);
            
            auto& p = p_snd_buffer[dst_index];
            const IntVect& pshift = p_periodic_shift[i];
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
                if (not is_per[idim]) continue;
                if (pshift[idim] > 0) 
                    p.pos(idim) += phi[idim] - plo[idim]; 
                else if (pshift[idim] < 0)
                    p.pos(idim) -= phi[idim] - plo[idim];
            }
	});
    }
}

template <class PC, class Buffer, class UnpackPolicy>
void unpackBuffer (PC& pc, const ParticleCopyPlan& plan, const Buffer& snd_buffer, const UnpackPolicy policy)
{
    BL_PROFILE("amrex::unpackBuffer");

    AMREX_ASSERT(pc.numLevels() == 1);

    using PTile = typename PC::ParticleTileType;

    const int lev = 0;
    auto& plev  = pc.GetParticles(lev);

    if (pc.ParticleBoxArray(lev).size() == 1 and (not pc.Geom(lev).isAnyPeriodic()) ) return;
    
    // count how many particles we have to add to each tile
    std::vector<int> sizes;
    std::vector<PTile*> tiles;
    for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);
        auto& tile = plev[index];
	int num_copies = plan.m_box_counts[pc.BufferMap().gridToBucket(gid)];
        sizes.push_back(num_copies);
        tiles.push_back(&tile);
    }

    // resize the tiles and compute offsets
    std::vector<int> offsets;
    policy.resizeTiles(tiles, sizes, offsets);

    // local unpack
    int uindex = 0;
    auto p_dst_box_counts = plan.m_box_counts.dataPtr();
    for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& tile = plev[index];
        auto& aos   = tile.GetArrayOfStructs();

	auto p_box_offsets = plan.m_box_offsets.dataPtr();
	auto p_box_perm = pc.BufferMap().gridToBucketPtr();
	auto p_snd_buffer = snd_buffer.dataPtr();

        int offset = offsets[uindex];
        int size = sizes[uindex];
        ++uindex;
        
        auto ptd = tile.getParticleTileData();
        AMREX_FOR_1D ( size, i,
        {
            int box_offset = p_box_offsets[p_box_perm[gid]];
            int src_index = box_offset + i;
            int dst_index = offset + i;
            ptd.setSuperParticle(p_snd_buffer[src_index], dst_index);
        });
    }    
}

template <class PC, class Buffer>
void communicateParticlesStart (const PC& pc, ParticleCopyPlan& plan, const Buffer& snd_buffer, Buffer& rcv_buffer)
{
    BL_PROFILE("amrex::communicateParticlesStart");
    
    AMREX_ASSERT(pc.numLevels() == 1);

    using ParticleType = typename PC::SuperParticleType;

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    const int MyProc = ParallelDescriptor::MyProc();
    
    if (NProcs == 1) return;

    Vector<int> RcvProc;
    Vector<long> rOffset;    
    long TotRcvParticles = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (plan.m_rcv_num_particles[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvParticles);
            TotRcvParticles += plan.m_rcv_num_particles[i];
        }
    }
    rcv_buffer.resize(TotRcvParticles);

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
        const auto Cnt    = plan.m_rcv_num_particles[Who]*sizeof(ParticleType);
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        plan.m_particle_rreqs[i] = 
            ParallelDescriptor::Arecv((char*) (rcv_buffer.dataPtr() + offset), Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        const auto Who = i;
        const auto Cnt = plan.m_snd_num_particles[i]*sizeof(ParticleType);
        if (Cnt == 0) continue;

        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        int bucket = pc.BufferMap().firstBucketOnProc(i);
        const auto snd_offset = plan.m_box_offsets[bucket];

        ParallelDescriptor::Asend( (char*) (snd_buffer.dataPtr() + snd_offset), Cnt, Who, SeqNum);
    }
#endif // MPI
}

void communicateParticlesFinish (const ParticleCopyPlan& plan);

template <class PC, class Buffer, class UnpackPolicy>
void unpackRemotes (PC& pc, const ParticleCopyPlan& plan, Buffer& rcv_buffer, const UnpackPolicy& policy)
{
    BL_PROFILE("amrex::unpackRemotes");
    
    AMREX_ASSERT(pc.numLevels() == 1);

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    if (NProcs == 1) return;

    const int MyProc = ParallelDescriptor::MyProc();
    using PTile = typename PC::ParticleTileType;

    if (plan.m_nrcvs > 0)
    {
 	const int lev = 0;
	auto& plev  = pc.GetParticles(lev);

	auto p_rcv_box_counts = plan.m_rcv_box_counts.dataPtr();
	auto p_rcv_box_offsets = plan.m_rcv_box_offsets.dataPtr();
	auto p_rcv_box_ids = plan.m_rcv_box_ids.dataPtr();
	auto p_rcv_buffer = rcv_buffer.dataPtr();

        std::vector<int> sizes;
        std::vector<PTile*> tiles;
	for (int i = 0; i < plan.m_rcv_box_counts.size(); ++i)
	{
            int copy_size = plan.m_rcv_box_counts[i];
            int gid = plan.m_rcv_box_ids[i];
            int tid = 0;
            auto index = std::make_pair(gid, tid);         
            auto& tile = plev[index];
            sizes.push_back(copy_size);
            tiles.push_back(&tile);
        }

        Vector<int> offsets;
        policy.resizeTiles(tiles, sizes, offsets);
        Gpu::Device::synchronize();
        int uindex = 0;
	for (int i = 0; i < plan.m_rcv_box_counts.size(); ++i)
	{
	  int gid = plan.m_rcv_box_ids[i];
	  int tid = 0;
	  int offset = plan.m_rcv_box_offsets[i];

	  auto index = std::make_pair(gid, tid);
	  
	  auto& tile = plev[index];
          auto ptd = tile.getParticleTileData();

	  AMREX_ASSERT(MyProc == pc.ParticleDistributionMap(0)[gid]);

          int dst_offset = offsets[uindex];
          int size = sizes[uindex];
          ++uindex;

	  AMREX_FOR_1D ( size, ip,
	  {
              int src_index = offset + ip;
              int dst_index = dst_offset + ip;
              ptd.setSuperParticle(p_rcv_buffer[src_index], dst_index);
	  });	 
	}	
    }    
#endif // MPI
}

} // namespace amrex

#endif // AMREX_PARTICLECOMMUNICATION_H_
