#include "NeighborList.H"

using namespace amrex;

void NeighborList::print ()
{
    BL_PROFILE("NeighborList::print");
    
    Gpu::HostVector<unsigned int> host_nbor_offsets(m_nbor_offsets.size());
    Gpu::HostVector<unsigned int> host_nbor_list(m_nbor_list.size());
    
    Cuda::thrust_copy(m_nbor_offsets.begin(),
                      m_nbor_offsets.end(),
                      host_nbor_offsets.begin());
    
    Cuda::thrust_copy(m_nbor_list.begin(),
                      m_nbor_list.end(),
                      host_nbor_list.begin());
    
    for (int i = 0; i < numParticles(); ++i) {
        amrex::Print() << "Particle " << i << " could collide with: ";
        for (int j = host_nbor_offsets[i]; j < host_nbor_offsets[i+1]; ++j) {
            amrex::Print() << host_nbor_list[j] << " ";
        }
        amrex::Print() << "\n";
    }
}
