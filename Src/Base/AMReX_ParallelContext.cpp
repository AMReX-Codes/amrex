#include <sstream>

#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex {
namespace ParallelContext {

Vector<Frame> frames; // stack of communicator frames

Frame::Frame (MPI_Comm c, int id, int io_rank)
    : comm(c),
      m_id(id),
      m_mpi_tag(ParallelDescriptor::MinTag()),
      m_io_rank(io_rank),
      m_out_filename(""),
      m_out(nullptr)
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

Frame::Frame (MPI_Comm c)
    : Frame(c, -1, ParallelDescriptor::IOProcessorNumber()) {}

Frame::Frame (Frame && rhs) noexcept
    : comm     (rhs.comm),
      group    (rhs.group),
      m_id     (rhs.m_id),
      m_rank_me(rhs.m_rank_me),
      m_nranks (rhs.m_nranks),
      m_mpi_tag(rhs.m_mpi_tag),
      m_io_rank(rhs.m_io_rank),
      m_out_filename(std::move(rhs.m_out_filename)),
      m_out    (std::move(rhs.m_out))
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
Frame::local_to_global_rank (int* global, const int* local, std::size_t n) const
{
#ifdef BL_USE_MPI
    if (frames.size() > 1)
    {
      MPI_Group_translate_ranks(GroupSub(), n, const_cast<int*>(local), GroupAll(), global);
    }
    else
    {
        for (std::size_t i = 0; i < n; ++i) global[i] = local[i];
    }
#else
    amrex::ignore_unused(local);
    for (std::size_t i = 0; i < n; ++i) global[i] = 0;
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
Frame::global_to_local_rank (int* local, const int* global, std::size_t n) const
{
#ifdef BL_USE_MPI
    if (frames.size() > 1)
    {
      MPI_Group_translate_ranks(GroupAll(), n, const_cast<int*>(global), GroupSub(), local);
    }
    else
    {
        for (std::size_t i = 0; i < n; ++i) local[i] = global[i];
    }
#else
    amrex::ignore_unused(global);
    for (std::size_t i = 0; i < n; ++i) local[i] = 0;
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

void
Frame::set_ofs_name (std::string filename)
{
    m_out.reset(); // in case changing name, close previous file
    m_out_filename = std::move(filename);
}

std::ofstream *
Frame::get_ofs_ptr ()
{
    if (m_out_filename.empty()) {
        return nullptr;
    } else {
        if (!m_out) {
            m_out.reset(new std::ofstream(m_out_filename, std::ios_base::app));
        }
        return m_out.get();
    }
}

}}
