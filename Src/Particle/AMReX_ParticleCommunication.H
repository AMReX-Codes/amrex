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

} // namespace amrex

#endif // AMREX_PARTICLECOMMUNICATION_H_
