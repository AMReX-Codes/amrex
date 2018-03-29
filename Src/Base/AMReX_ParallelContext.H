#ifndef AMREX_PARALLELCONTEXT_H
#define AMREX_PARALLELCONTEXT_H

#include <AMReX_Vector.H>
#include <AMReX_ccse-mpi.H>

namespace amrex {
namespace ParallelContext {

class Frame
{
public:

    explicit Frame (MPI_Comm c);
    Frame (Frame const&) = delete;
    Frame (Frame && rhs);
    ~Frame ();

    int MyProc () const { return m_rank_me; }
    int NProcs () const { return m_nranks; }
    int local_to_global_rank (int lrank) const;
    void local_to_global_rank (int* global, const int* local, int n) const;
    void global_to_local_rank (int* local, const int* global, int n) const;
    int global_to_local_rank (int grank) const;
    int get_inc_mpi_tag ();

    // sub-communicator associated with frame
    MPI_Comm comm = MPI_COMM_NULL;
    MPI_Group group = MPI_GROUP_NULL; // to avoid repeatedly creating groups in rank translation

private:
    int m_rank_me = -1; // local rank
    int m_nranks  =  0; // local # of ranks
    int m_mpi_tag = -1;
};

extern Vector<Frame> frames; // stack of communicator frames

// world communicator
inline MPI_Comm CommunicatorAll () { return frames[0].comm; }
// world group
inline MPI_Comm GroupAll () { return frames[0].group; }
// number of ranks in world communicator
inline int NProcsAll () { return frames[0].NProcs(); }
// my rank in world communicator
inline int MyProcAll () { return frames[0].MyProc(); }

// sub-communicator for current frame
inline MPI_Comm CommunicatorTop () { return frames.back().comm; }
// sub-group for current frame
inline MPI_Group GroupTop () { return frames.back().group; }
// number of ranks in current frame
inline int NProcsTop () { return frames.back().NProcs(); }
// my sub-rank in current frame
inline int MyProcTop () { return frames.back().MyProc(); }

// get and increment mpi tag in current frame
inline int get_inc_mpi_tag () { return frames.back().get_inc_mpi_tag(); }
// translate between local rank and global rank
inline int local_to_global_rank (int rank) { return frames.back().local_to_global_rank(rank); }
inline void local_to_global_rank (int* global, const int* local, int n)
    { frames.back().local_to_global_rank(global, local, n); }
inline int global_to_local_rank (int rank) { return frames.back().global_to_local_rank(rank); }
inline void global_to_local_rank (int* local, const int* global, int n)
    { frames.back().global_to_local_rank(local, global, n); }

inline void push (MPI_Comm c) { frames.emplace_back(c); }
// Note that it's the user's responsibility to free the MPI_Comm
inline void pop () { frames.pop_back(); }

}}

#endif // AMREX_PARALLELCONTEXT_H
