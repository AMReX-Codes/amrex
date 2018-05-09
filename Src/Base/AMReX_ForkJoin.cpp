#include <AMReX_ForkJoin.H>
#include <AMReX_Print.H>

namespace amrex {

ForkJoin::ForkJoin (Vector<int> trn)
    : task_rank_n(std::move(trn))
{
    init();
}

ForkJoin::ForkJoin (const Vector<double> &task_rank_pct)
{
    auto rank_n = ParallelContext::NProcsSub(); // number of ranks in current frame
    auto ntasks = task_rank_pct.size();
    task_rank_n.resize(ntasks);
    int prev = 0;
    double accum = 0;
    for (int i = 0; i < ntasks; ++i) {
        accum += task_rank_pct[i];
        int cur = std::round(rank_n * accum);
        task_rank_n[i] = cur - prev;
        prev = cur;
    }

    init();
}

void
ForkJoin::init()
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(NTasks() > 0,
                                     "ForkJoin must have at least 1 task");
    int min_task_rank_n = task_rank_n[0];
    for (int i = 1; i < NTasks(); ++i) {
      min_task_rank_n = std::min(min_task_rank_n, task_rank_n[i]);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(min_task_rank_n > 0,
                                     "All tasks must have non-negative ranks");
    auto rank_n = ParallelContext::NProcsSub(); // number of ranks in current frame
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::accumulate(task_rank_n.begin(),task_rank_n.end(),0) == rank_n,
                                     "Sum of ranks assigned to tasks must sum to parent number of ranks");
    compute_split_bounds();
}

void
ForkJoin::reg_mf (MultiFab &mf, const std::string &name, int idx,
                  Strategy strategy, Intent intent, int owner)
{
    if (idx >= data[name].size()) {
        data[name].resize(idx + 1);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(data[name][idx].empty(),
                                     "Can only register to a (name, index) pair once");
    data[name][idx] = MFFork(&mf, strategy, intent, owner);
}

void
ForkJoin::copy_data_to_tasks (MPI_Comm /*task_comm*/)
{
    if (flag_verbose) {
        amrex::Print() << "Copying data into fork-join tasks ...\n";
    }
    for (auto &p : data) { // for each name
        const auto &mf_name = p.first;
        for (int idx = 0; idx < p.second.size(); ++idx) { // for each index
            auto &mff = p.second[idx];
            const MultiFab &orig = *mff.orig;
            const auto &ba = orig.boxArray();
            Vector<MultiFab> &forked = mff.forked;
            int comp_n = orig.nComp(); // number of components in original

            if (mff.strategy == Strategy::split) {
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(NTasks() <= comp_n,
                                                 "Number of tasks cannot be larger than number of components!");
            }

            forked.reserve(NTasks()); // does nothing if forked MFs already created
            for (int i = 0; i < NTasks(); ++i) {
                // check if this task needs this MF
                if (mff.strategy != Strategy::single || i == mff.owner_task) {

                    // compute task's component lower and upper bound
                    int comp_lo, comp_hi;
                    if (mff.strategy == Strategy::split) {
                        // split components across tasks
                        comp_lo = comp_n *  i    / NTasks();
                        comp_hi = comp_n * (i+1) / NTasks();
                    } else {
                        // copy all components to task
                        comp_lo = 0;
                        comp_hi = comp_n;
                    }

                    // create task's MF if first time through
                    if (forked.size() <= i) {
                        if (flag_verbose) {
                            amrex::Print() << "  Creating forked " << mf_name << "[" << idx << "] for task " << i
                                           << (mff.strategy == Strategy::split ? " (split)" : " (whole)") << std::endl;
                        }
                        // look up the distribution mapping for this (box array, task) pair
                        const DistributionMapping &dm = get_dm(ba, i, orig.DistributionMap());
                        forked.emplace_back(ba, dm, comp_hi - comp_lo, 0);
                    } else if (flag_verbose) {
                        amrex::Print() << "  Forked " << mf_name << "[" << idx << "] for task " << i
                                       << " already created" << std::endl;
                    }
                    AMREX_ASSERT(i < forked.size());

                    // copy data if needed
                    if (mff.intent == Intent::in || mff.intent == Intent::inout) {
                        if (flag_verbose) {
                            amrex::Print() << "    Copying " << mf_name << "[" << idx << "] into to task " << i << std::endl;
                        }
                        // parallel copy data into forked MF
                        forked[i].copy(orig, comp_lo, 0, comp_hi - comp_lo);
                    }

                } else {
                    // this task doesn't use the MultiFab
                    if (forked.size() <= i) {
                        // first time through, push empty placeholder (not used)
                        forked.push_back(MultiFab());
                    }
                }
            }
            AMREX_ASSERT(forked.size() == NTasks());
        }
    }
}

// this is called after ParallelContext::unsplit
// the parent task is the top frame in ParallelContext's stack
void
ForkJoin::copy_data_from_tasks ()
{
    if (flag_verbose) {
        amrex::Print() << "Copying data out of fork-join tasks ...\n";
    }
    for (auto &p : data) { // for each name
        const auto &mf_name = p.first;
        for (int idx = 0; idx < p.second.size(); ++idx) { // for each index
            auto &mff = p.second[idx];
            if (mff.intent == Intent::out || mff.intent == Intent::inout) {
                MultiFab &orig = *mff.orig;
                int comp_n = orig.nComp(); // number of components in original
                const Vector<MultiFab> &forked = mff.forked;
                if (mff.strategy == Strategy::split) {
                    // gather components from across tasks
                    for (int i = 0; i < NTasks(); ++i) {
                        if (flag_verbose) {
                            amrex::Print() << "  Copying " << mf_name << "[" << idx << "] out from task " << i << "  (unsplit)" << std::endl;
                        }
                        int comp_lo = comp_n *  i    / NTasks();
                        int comp_hi = comp_n * (i+1) / NTasks();
                        orig.copy(forked[i], 0, comp_lo, comp_hi - comp_lo);
                    }
                } else { // mff.strategy == single or duplicate
                    // copy all components from owner_task
                    if (flag_verbose) {
                        amrex::Print() << "Copying " << mf_name << " out from task " << mff.owner_task << "  (whole)" << std::endl;
                    }
                    orig.copy(forked[mff.owner_task], 0, 0, comp_n);
                }
            }
        }
    }
}

// multiple MultiFabs may share the same box array
// only compute the DM once per unique (box array, task) pair and cache it
// create map from box array RefID to vector of DistributionMapping indexed by task ID
const DistributionMapping &
ForkJoin::get_dm (const BoxArray& ba, int task_idx, const DistributionMapping& dm_orig)
{
    auto &dm_vec = dms[ba.getRefID()];

    if (dm_vec.size() == 0) {
        // new entry
        dm_vec.resize(NTasks());
    }
    AMREX_ASSERT(task_idx < dm_vec.size());

    if (dm_vec[task_idx] == nullptr) {
        // create DM of current box array over current task's ranks
        int rank_lo = split_bounds[task_idx];  // note that these ranks are not necessarily global
        int nprocs_task = task_rank_n[task_idx];

        Vector<int> pmap = dm_orig.ProcessorMap(); // DistributionMapping stores global ranks
        for (auto& r : pmap) {
            int lr = ParallelContext::global_to_local_rank(r);
            lr = lr%nprocs_task + rank_lo;
            r = ParallelContext::local_to_global_rank(lr);
        }

        dm_vec[task_idx].reset(new DistributionMapping(std::move(pmap)));

        if (flag_verbose) {
            amrex::Print() << "    Creating DM for (box array, task id) = ("
                      << ba.getRefID() << ", " << task_idx << ")" << std::endl;
        }

//        amrex::Print() << " xxxxx get_dm " << task_idx << ", " << *dm_vec[task_idx] << "\n";

    } else {
        // DM has already been created
        if (flag_verbose) {
            amrex::Print() << "    DM for (box array, task id) = (" << ba.getRefID() << ", " << task_idx
                           << ") already created" << std::endl;
        }
    }
    AMREX_ASSERT(dm_vec[task_idx] != nullptr);

    return *dm_vec[task_idx];
}

// split ranks in current frame into contiguous chunks
// task i has ranks over the interval [result[i], result[i+1])
void
ForkJoin::compute_split_bounds ()
{
    AMREX_ASSERT(std::accumulate(task_rank_n.begin(),task_rank_n.end(),0) == ParallelContext::NProcsSub());

    const auto ntasks = task_rank_n.size();
    split_bounds.resize(ntasks + 1);
    split_bounds[0] = 0;
    for (int i = 0; i < ntasks; ++i) {
        split_bounds[i + 1] = split_bounds[i] + task_rank_n[i];
    }
}

// split top frame of stack
// TODO: write version that takes cached comm object as argument in case of repeated identical split calls
MPI_Comm
ForkJoin::split_tasks ()
{
    const auto ntasks = task_rank_n.size();
    int myproc = ParallelContext::MyProcSub();
    for (task_me = 0; task_me < ntasks; ++task_me) {
        int lo = split_bounds[task_me];
        int hi = split_bounds[task_me + 1];
        if (myproc >= lo && myproc < hi) {
            break;
        }
    }
    AMREX_ASSERT(task_me < ntasks);

#ifdef BL_USE_MPI
    MPI_Comm new_comm;
    MPI_Comm_split(ParallelContext::CommunicatorSub(), task_me, myproc, &new_comm);
#else
    MPI_Comm new_comm = ParallelContext::CommunicatorSub();
#endif

    return new_comm;
}

}
