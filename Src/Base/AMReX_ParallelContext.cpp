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
    MPI_Comm_group(comm, &group);
    MPI_Comm_rank(comm, &m_rank_me);
    MPI_Comm_size(comm, &m_nranks);
#else
    m_rank_me = 0;
    m_nranks = 1;
#endif
}

Frame::Frame (Frame && rhs)
    : comm     (rhs.comm),
      group    (rhs.group),
      m_rank_me(rhs.m_rank_me),
      m_nranks (rhs.m_nranks),
      m_mpi_tag(rhs.m_mpi_tag)
{
    rhs.group = MPI_GROUP_NULL;
}

Frame::~Frame ()
{
#ifdef BL_USE_MPI
    if (group != MPI_GROUP_NULL) {
        MPI_Group_free(&group);
    }
#endif
}

int
Frame::local_to_global_rank (int lrank) const
{
    int r;
    local_to_global_rank(&r, &lrank, 1);
    return r;
}

void
Frame::local_to_global_rank (int* global, const int* local, int n) const
{
#ifdef BL_USE_MPI
    if (frames.size() > 1)
    {
        MPI_Group_translate_ranks(GroupTop(), n, local, GroupAll(), global);
    }
    else
    {
        for (int i = 0; i < n; ++i) global[i] = local[i];
    }
#else
    for (int i = 0; i < n; ++i) global[i] = 0;
#endif
}

int
Frame::global_to_local_rank (int grank) const
{
    int r;
    global_to_local_rank(&r, &grank, 1);
    return r;
}

void
Frame::global_to_local_rank (int* local, const int* global, int n) const
{
#ifdef BL_USE_MPI
    if (frames.size() > 1)
    {
        MPI_Group_translate_ranks(GroupAll(), n, global, GroupTop(), local);
    }
    else
    {
        for (int i = 0; i < n; ++i) local[i] = global[i];
    }
#else
    for (int i = 0; i < n; ++i) local[i] = 0;
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
