#ifndef AMREX_PARALLELCONTEXT_H
#define AMREX_PARALLELCONTEXT_H

#include <AMReX_Vector.H>
#include <AMReX_ccse-mpi.H>

namespace amrex {
namespace ParallelContext {

struct Frame
{
    // sub-communicator associated with frame
    MPI_Comm comm = MPI_COMM_NULL;
    int glo_rank_lo = -1, glo_rank_hi = -1; // members of task have contiguous global ranks
    int loc_rank_me = -1; // local rank relative to current sub-communicator
    int mpi_tag = ParallelDescriptor::MinTag();

    Frame (MPI_Comm c, int glo, int ghi, int lme) :
      comm(c), glo_rank_lo(glo), glo_rank_hi(ghi), loc_rank_me(lme) { }

    int rank_n () { return glo_rank_hi - glo_rank_lo; }
    int rank_me () { return loc_rank_me; }
    int local_to_global_rank (int rank) {
        AMREX_ASSERT(rank >= 0 && rank <= rank_n()); // allow exclusive upper bound
        return rank + glo_rank_lo;
    }
    int global_to_local_rank (int rank) {
        AMREX_ASSERT(rank >= glo_rank_lo && rank <= glo_rank_hi); // allow exclusive upper bound
        return rank - glo_rank_lo;
    }
    int get_inc_mpi_tag () {
        // get and increment the mpi_tag in this frame
        auto cur_tag = mpi_tag;
        mpi_tag = mpi_tag < ParallelDescriptor::MaxTag() ?
                  mpi_tag + 1 : ParallelDescriptor::MinTag();
        return cur_tag;
    }
};

extern Vector<Frame> frames; // stack of communicator frames

void init ();
void finalize ();

// world communicator
inline MPI_Comm CommunicatorAll () { return frames[0].comm; }
// number of ranks in world communicator
inline int NProcsAll () { return frames[0].rank_n(); }
// my rank in world communicator
inline int MyProcAll () { return frames[0].rank_me(); }

// sub-communicator for current frame
inline MPI_Comm Communicator () { return frames.back().comm; }
// number of ranks in current frame
inline int NProcs () { return frames.back().rank_n(); }
// my sub-rank in current frame
inline int MyProc () { return frames.back().rank_me(); }

// get and increment mpi tag in current frame
inline int get_inc_mpi_tag () { return frames.back().get_inc_mpi_tag(); }
// translate between local rank and global rank
inline int local_to_global_rank (int rank) { return frames.back().local_to_global_rank(rank); }
inline int global_to_local_rank (int rank) { return frames.back().global_to_local_rank(rank); }

// split ranks in current frame into contiguous chunks
// task i has ranks over the interval [result[i], result[i+1])
Vector<int> get_split_bounds (const Vector<int> &task_rank_n);
// split top frame of stack and push new frame on top
int split (const Vector<int> &task_rank_n);
// pop split frame off the top
void unsplit ();

}}

#endif // AMREX_PARALLELCONTEXT_H
