#include <AMReX_ForkJoin.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

namespace {

inline bool file_exists(std::string file_path) {
  std::ifstream ifs(file_path);
  return ifs.good();
}

template <class T>
std::string
str_join (Vector<T> xs, std::string sep)
{
    std::ostringstream ss;
    bool flag_first = true;
    for (int i = 0; i < xs.size(); ++i) {
        if (!flag_first) {
            ss << sep;
        }
        flag_first = false;
        ss << xs[i];
    }
    return ss.str();
}

Vector<int>
get_frame_id_vec ()
{
    const auto &frames = amrex::ParallelContext::frames;
    Vector<int> result;
    // ignore first (global) frame
    for (int i = 1; i < frames.size(); ++i) {
        result.push_back(frames[i].MyID());
    }
    return result;
}

}

namespace amrex {

ForkJoin::ForkJoin (const Vector<int> &task_rank_n)
{
    init(task_rank_n);
}

ForkJoin::ForkJoin (const Vector<double> &task_rank_pct)
{
    auto rank_n = ParallelContext::NProcsSub(); // number of ranks in current frame
    auto ntasks = task_rank_pct.size();
    Vector<int> task_rank_n(ntasks);
    int prev = 0;
    double accum = 0;
    for (int i = 0; i < ntasks; ++i) {
        accum += task_rank_pct[i];
        int cur = std::round(rank_n * accum);
        task_rank_n[i] = cur - prev;
        prev = cur;
    }

    init(task_rank_n);
}

void
ForkJoin::init(const Vector<int> &task_rank_n)
{
    ParmParse pp("forkjoin");
    pp.query("verbose", flag_verbose);

    const auto task_n = task_rank_n.size();
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(task_n > 0,
                                     "ForkJoin must have at least 1 task");
    int min_task_rank_n = task_rank_n[0];
    for (int i = 1; i < task_n; ++i) {
      min_task_rank_n = std::min(min_task_rank_n, task_rank_n[i]);
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(min_task_rank_n > 0,
                                     "All tasks must have at least one rank");
    auto rank_n = ParallelContext::NProcsSub(); // number of ranks in current frame
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::accumulate(task_rank_n.begin(),task_rank_n.end(),0) == rank_n,
                                     "Sum of ranks assigned to tasks must sum to parent number of ranks");

    // split ranks into contiguous chunks
    // task i has ranks over the interval [split_bounds[i], split_bounds[i+1])
    split_bounds.resize(task_n + 1);
    split_bounds[0] = 0;
    for (int i = 0; i < task_n; ++i) {
        split_bounds[i + 1] = split_bounds[i] + task_rank_n[i];
    }

    if (flag_verbose) {
        amrex::Print() << "Initialized ForkJoin:\n";
        for (int i = 0; i < task_n; ++i) {
            int glo_rank_lo = ParallelContext::local_to_global_rank(split_bounds[i]);
            int glo_rank_hi = ParallelContext::local_to_global_rank(split_bounds[i+1]-1);
            amrex::Print() << "  Task " << i << " has " << NProcsTask(i)
                           << " Ranks: [" << glo_rank_lo << ", " << glo_rank_hi << "]\n";
        }
    }
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

    IntVect ngrow = mf.nGrowVect(); // default use original MultiFab's grow cells

    // compute how components are copied to tasks
    // default evenly split components across tasks
    int comp_n = mf.nComp(); // number of components in original
    Vector<ComponentSet> comp_split;
    comp_split.resize(NTasks());
    for (int i = 0; i < NTasks(); ++i) {
        if (strategy == Strategy::split) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(NTasks() <= comp_n,
                                             "Number of tasks cannot be larger than number of components!");
            // split components across tasks
            comp_split[i].lo = comp_n *  i    / NTasks();
            comp_split[i].hi = comp_n * (i+1) / NTasks();
        } else {
            // copy all components to task
            comp_split[i].lo = 0;
            comp_split[i].hi = comp_n;
        }
    }

    data[name][idx] = MFFork(&mf, strategy, intent, owner, ngrow, std::move(comp_split));
}

void
ForkJoin::modify_ngrow (const std::string &name, int idx, IntVect ngrow)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(data.count(name) > 0 && data[name].size() > idx,
                                     "(name, index) pair doesn't exist");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!flag_invoked,
                                     "Can only specify grow cells before first forkjoin() invocation");
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ngrow[i] >= 0,
                                         "ngrow[i] must be non-negative");
    }
    data[name][idx].ngrow = ngrow;
}

void
ForkJoin::modify_split (const std::string &name, int idx, Vector<ComponentSet> comp_split)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(data.count(name) > 0 && data[name].size() > idx,
                                     "(name, index) pair doesn't exist");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!flag_invoked,
                                     "Can only specify custom split before first forkjoin() invocation");
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(comp_split.size() == NTasks(),
                                     "comp_split must be same length as number of tasks");
    for (int i = 0; i < NTasks(); ++i) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(comp_split[i].hi - comp_split[i].lo > 0,
                                         "comp_split[i] must have positive number of components");
    }
    data[name][idx].comp_split = std::move(comp_split);
}

ForkJoin::ComponentSet
ForkJoin::ComponentBounds(const std::string& name, int idx) const
{
    const auto & mffork_vec = data.at(name);
    BL_ASSERT(idx >= 0 && idx < mffork_vec.size());
    return mffork_vec[idx].comp_split[task_me];
}

void
ForkJoin::copy_data_to_tasks ()
{
    BL_PROFILE("ForkJoin::copy_data_to_tasks()");
    if (flag_verbose) {
        amrex::Print() << "Copying data into fork-join tasks ...\n";
    }
    for (auto &p : data) { // for each name
        const auto &mf_name = p.first;
        for (int idx = 0; idx < p.second.size(); ++idx) { // for each index
            auto &mff = p.second[idx];
            const auto &orig = *mff.orig;
            const auto &ba = orig.boxArray();
            const auto &comp_split = mff.comp_split;
            auto &forked = mff.forked;

            forked.reserve(NTasks()); // does nothing if forked MFs already created
            for (int i = 0; i < NTasks(); ++i) {
                // check if this task needs this MF
                if (mff.strategy != Strategy::single || i == mff.owner_task) {
                    int task_comp_n = comp_split[i].hi - comp_split[i].lo;

                    // create task's MF if first time through
                    if (forked.size() <= i) {
                        if (flag_verbose) {
                            amrex::Print() << "  Creating forked " << mf_name << "[" << idx << "] for task " << i
                                           << (mff.strategy == Strategy::split ? " (split)" : " (whole)") << std::endl;
                        }
                        // look up the distribution mapping for this (box array, task) pair
                        const DistributionMapping &dm = get_dm(ba, i, orig.DistributionMap());
                        forked.emplace_back(ba, dm, task_comp_n, mff.ngrow);
                    } else if (flag_verbose) {
                        amrex::Print() << "  Forked " << mf_name << "[" << idx << "] for task " << i
                                       << " already created" << std::endl;
                    }
                    AMREX_ASSERT(i < forked.size());

                    // copy data if needed
                    if (mff.intent == Intent::in || mff.intent == Intent::inout) {
                        if (flag_verbose) {
                            amrex::Print() << "    Copying " << mf_name << "[" << idx << "] components ["
                                           << comp_split[i].lo << ", " << comp_split[i].hi << ") into to task " << i << std::endl;
                        }
                        // parallel copy data into forked MF
                        forked[i].Redistribute(orig, comp_split[i].lo, 0, task_comp_n, mff.ngrow);
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
    BL_PROFILE("ForkJoin::copy_data_from_tasks()");
    if (flag_verbose) {
        amrex::Print() << "Copying data out of fork-join tasks ...\n";
    }
    for (auto &p : data) { // for each name
        const auto &mf_name = p.first;
        for (int idx = 0; idx < p.second.size(); ++idx) { // for each index
            auto &mff = p.second[idx];
            if (mff.intent == Intent::out || mff.intent == Intent::inout) {
                MultiFab &orig = *mff.orig;
                const auto &comp_split = mff.comp_split;
                const Vector<MultiFab> &forked = mff.forked;
                if (mff.strategy == Strategy::split) {
                    // gather components from across tasks
                    for (int i = 0; i < NTasks(); ++i) {
                        if (flag_verbose) {
                            amrex::Print() << "  Copying " << mf_name << "[" << idx << "] components ["
                                           << comp_split[i].lo << ", " << comp_split[i].hi << ") out from task " << i << " (unsplit)" << std::endl;
                        }
                        int task_comp_n = comp_split[i].hi - comp_split[i].lo;
                        AMREX_ASSERT(forked[i].nComp() == task_comp_n);
                        orig.Redistribute(forked[i], 0, comp_split[i].lo, task_comp_n, mff.ngrow);
                    }
                } else { // mff.strategy == single or duplicate
                    // copy all components from owner_task
                    if (flag_verbose) {
                        amrex::Print() << "Copying " << mf_name << " out from task " << mff.owner_task << " (whole)" << std::endl;
                    }
                    AMREX_ASSERT(forked[mff.owner_task].nComp() == orig.nComp());
                    orig.Redistribute(forked[mff.owner_task], 0, 0, orig.nComp(), mff.ngrow);
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
    AMREX_ASSERT(task_idx < NTasks());

    auto &dm_vec = dms[ba.getRefID()];
    if (dm_vec.size() == 0) {
        // new entry
        dm_vec.resize(NTasks());
    }
    AMREX_ASSERT(dm_vec.size() == NTasks());

    if (dm_vec[task_idx] == nullptr) {
        // create DM of current box array over current task's ranks
        int rank_lo = split_bounds[task_idx];  // note that these ranks are not necessarily global
        int nprocs_task = NProcsTask(task_idx);

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

// split top frame of stack
// TODO: write version that takes cached comm object as argument in case of repeated identical split calls
MPI_Comm
ForkJoin::split_tasks ()
{
    int myproc = ParallelContext::MyProcSub();
    for (task_me = 0; task_me < NTasks(); ++task_me) {
        int lo = split_bounds[task_me];
        int hi = split_bounds[task_me + 1];
        if (myproc >= lo && myproc < hi) {
            break;
        }
    }
    AMREX_ASSERT(task_me < NTasks());

#ifdef BL_USE_MPI
    MPI_Comm new_comm;
    MPI_Comm_split(ParallelContext::CommunicatorSub(), task_me, myproc, &new_comm);
#else
    MPI_Comm new_comm = ParallelContext::CommunicatorSub();
#endif

    return new_comm;
}

void ForkJoin::create_task_output_dir ()
{
    if (task_output_dir != "" && !amrex::FileExists(task_output_dir)) {
        if (flag_verbose) {
            Print() << "Creating task_output_dir: " << task_output_dir << std::endl;
        }
        if (ParallelContext::IOProcessorSub()) {
            amrex::UtilCreateDirectory(task_output_dir, 0755, flag_verbose);
        }
    }
}

std::string
ForkJoin::get_io_filename (bool flag_unique)
{
    std::string result = "";

    if (task_output_dir != "") {
        // build base filename
        std::string result_base = task_output_dir;
        result_base += "/T-" + str_join(get_frame_id_vec(), "-");
        result_base += ".R-" + std::to_string(ParallelContext::MyProcSub());

        if (flag_unique) {
            // concatenate an integer to the end to make unique
            int i = 0;
            do {
                result = result_base + ".I-" + std::to_string(i++) + ".out";
            } while (file_exists(result));
        } else {
            result = result_base + ".out";
        }
    }

    return result;
}

}
