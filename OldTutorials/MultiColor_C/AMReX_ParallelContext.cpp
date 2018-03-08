#include <AMReX_ParallelContext.H>

namespace {

template <class T>
T sum(const amrex::Vector<T> &v) {
    T s = 0;
    for (int i = 0; i < v.size(); ++i) {
        s += v[i];
    }
    return s;
}

}

namespace amrex {
namespace ParallelContext {

Vector<Frame> frames; // stack of communicator frames

// call at beginning of program, after MPI_Init and ParallelDescriptor::StartParallel()
// probably somewhere inside amrex::Initialize()
void init() {
    // initialize "global" first frame in stack to ParallelDescriptor's communicator
    int glo_rank_n, glo_rank_me;
    MPI_Comm_size(ParallelDescriptor::Communicator(), &glo_rank_n);
    MPI_Comm_rank(ParallelDescriptor::Communicator(), &glo_rank_me);
    frames.emplace_back(ParallelDescriptor::Communicator(), 0, glo_rank_n, glo_rank_me);
}

// split ranks in current frame into contiguous chunks
// task i has ranks over the interval [result[i], result[i+1])
Vector<int> get_split_bounds(const Vector<int> &task_rank_n)
{
    AMREX_ASSERT(sum(task_rank_n) == rank_n());

    const auto task_n = task_rank_n.size();
    Vector<int> result(task_n + 1);
    result[0] = 0;
    for (int i = 0; i < task_n; ++i) {
        result[i + 1] = result[i] + task_rank_n[i];
    }
    return result;
}

// split top frame of stack and push new frame on top
// TODO: write version that takes cached comm object as argument in case of repeated identical split calls
int split(const Vector<int> &task_rank_n)
{
    // figure out what color (task_me) to pass into MPI_Comm_split
    const auto task_n = task_rank_n.size();
    AMREX_ASSERT(sum(task_rank_n) == rank_n());
    auto split_bounds = get_split_bounds(task_rank_n);
    int new_glo_rank_lo, new_glo_rank_hi, new_loc_rank_me;
    int task_me;
    for (task_me = 0; task_me < task_n; ++task_me) {
        int lo = split_bounds[task_me];
        int hi = split_bounds[task_me + 1];
        if (rank_me() >= lo && rank_me() < hi) {
            new_glo_rank_lo = local_to_global_rank(lo);
            new_glo_rank_hi = local_to_global_rank(hi);
            new_loc_rank_me = rank_me() - lo;
            break;
        }
    }
    AMREX_ASSERT(task_me < task_n);

    MPI_Comm new_comm;
    MPI_Comm_split(comm(), task_me, rank_me(), &new_comm);

    frames.emplace_back(new_comm, new_glo_rank_lo, new_glo_rank_hi, new_loc_rank_me);
    return task_me;
}

void unsplit() {
    MPI_Comm_free(&frames.back().comm);
    frames.pop_back();
}

}}
