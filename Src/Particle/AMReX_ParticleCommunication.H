#ifndef AMREX_PARTICLECOMMUNICATION_H_
#define AMREX_PARTICLECOMMUNICATION_H_

#include <AMReX_Gpu.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_IntVect.H>

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

    void build ();

    void clear();
};

} // namespace amrex

#endif // AMREX_PARTICLECOMMUNICATION_H_
