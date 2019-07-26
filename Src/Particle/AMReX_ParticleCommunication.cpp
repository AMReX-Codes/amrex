#include <AMReX_ParticleCommunication.H>

using namespace amrex;

void ParticleCopyOp::clear ()
{
    m_boxes.clear();
    m_src_indices.clear();
    m_periodic_shift.clear();
}

void ParticleCopyPlan::build ()
{
    
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
