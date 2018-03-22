#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex {
namespace ParallelContext {

Vector<Frame> frames; // stack of communicator frames

Frame::Frame (MPI_Comm c)
    : comm(c),
      m_mpi_tag(ParallelDescriptor::MinTag())
{
#ifdef BL_USE_MPI
    MPI_Comm_rank(comm, &m_rank_me);
    MPI_Comm_size(comm, &m_nranks);
#else
    m_rank_me = 0;
    m_nranks = 1;
#endif
}

int
Frame::local_to_global_rank (int lrank) const
{
#ifdef BL_USE_MPI
    int r;
    MPI_Group ggroup, lgroup;
    MPI_Comm_group(ParallelContext::CommunicatorAll(), &ggroup);
    MPI_Comm_group(ParallelContext::Communicator(), &lgroup);
    MPI_Group_translate_ranks(lgroup, 1, &lrank, ggroup, &r);
    return r;
#else
    return 0;
#endif
}

int
Frame::global_to_local_rank (int grank) const
{
#ifdef BL_USE_MPI
    int r;
    MPI_Group ggroup, lgroup;
    MPI_Comm_group(ParallelContext::CommunicatorAll(), &ggroup);
    MPI_Comm_group(ParallelContext::Communicator(), &lgroup);
    MPI_Group_translate_ranks(ggroup, 1, &grank, lgroup, &r);
    return r;
#else
    return 0;
#endif    
}

int
Frame::get_inc_mpi_tag ()
{
    int cur_tag = m_mpi_tag;
    m_mpi_tag = (m_mpi_tag < ParallelDescriptor::MaxTag()) ?
        m_mpi_tag + 1 : ParallelDescriptor::MinTag();
    return cur_tag;
}

}}
