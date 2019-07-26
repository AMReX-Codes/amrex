#ifndef AMREX_PARTICLECOMMUNICATION_H_
#define AMREX_PARTICLECOMMUNICATION_H_

#include <AMReX_Gpu.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_IntVect.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_MFIter.H>

#include <map>

namespace amrex {

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

    template <class PC>
    void build (const PC& pc, const ParticleCopyOp& op)
    {
        BL_PROFILE("ParticleCopyPlan::build");
        
        AMREX_ASSERT(pc.numLevels() == 1);
        AMREX_ASSERT(pc.hasNeighbors() == false);
        
        const int lev = 0;
        const auto& geom = pc.Geom(lev);
        
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
            
            auto p_boxes = op.m_boxes.at(gid).dataPtr();
            auto p_dst_indices = m_dst_indices[gid].dataPtr();
            
            AMREX_FOR_1D ( num_copies, i,
            {
                int dst_box = p_boxes[i];
                int index = Gpu::Atomic::Inc(&p_dst_box_counts[p_box_perm[dst_box]], max_unsigned_int);
                p_dst_indices[i] = index;
            });
        }

        amrex::Gpu::exclusive_scan(m_box_counts.begin(), m_box_counts.end()+1, m_box_offsets.begin());

        buildMPI(pc.BufferMap());
    }

    void clear();

private:
    
    void buildMPI (const ParticleBufferMap& map);
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
        const auto p_ptr = &(aos[0]);

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
            int box_offset = p_box_offsets[p_box_perm[p_boxes[i]]];
            int dst_index = box_offset + p_dst_indices[i];
            int src_index = p_src_indices[i];
            p_snd_buffer[dst_index] = p_ptr[src_index];
            
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

template <class PC, class Buffer>
void unpackBuffer (PC& pc, const ParticleCopyOp& op, const ParticleCopyPlan& plan, const Buffer& snd_buffer)
{
    BL_PROFILE("amrex::packBuffer");

    AMREX_ASSERT(pc.numLevels() == 1);

    const int lev = 0;
    auto& plev  = pc.GetParticles(lev);

    if (pc.ParticleBoxArray(lev).size() == 1 and (not pc.Geom(lev).isAnyPeriodic()) ) return;
    
    // local unpack
    auto p_dst_box_counts = plan.m_box_counts.dataPtr();
    for(MFIter mfi = pc.MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& tile = plev[index];
        auto& aos   = tile.GetArrayOfStructs();
	auto np = tile.numTotalParticles();

	int num_copies = plan.m_box_counts[pc.BufferMap().gridToBucket(gid)];

	auto p_box_offsets = plan.m_box_offsets.dataPtr();
	auto p_box_perm = pc.BufferMap().gridToBucketPtr();
	auto p_snd_buffer = snd_buffer.dataPtr();

        int nn = tile.getNumNeighbors();
        tile.setNumNeighbors(nn + num_copies);	

        auto p_ptr = aos().dataPtr();
        AMREX_FOR_1D ( num_copies, i,
        {
            int box_offset = p_box_offsets[p_box_perm[gid]];
            int src_index = box_offset + i;
            int dst_index = np + i;
            p_ptr[dst_index] = p_snd_buffer[src_index];
        });
    }    
}

} // namespace amrex

#endif // AMREX_PARTICLECOMMUNICATION_H_
