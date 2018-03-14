#include <AMReX_ForkJoin.H>

namespace amrex {

ForkJoin::ForkJoin (Vector<int> trn)
    : task_rank_n(std::move(trn))
{
    auto rank_n = ParallelContext::NProcs(); // number of ranks in current frame
    AMREX_ASSERT(task_n() >= 2);
    AMREX_ASSERT(std::accumulate(task_rank_n.begin(),task_rank_n.end(),0) == rank_n);
}

ForkJoin::ForkJoin (const Vector<double> &task_rank_pct)
{
    auto rank_n = ParallelContext::NProcs(); // number of ranks in current frame
    auto ntasks = task_rank_pct.size();
    AMREX_ASSERT(ntasks >= 2);
    task_rank_n.resize(ntasks);
    int prev = 0;
    double accum = 0;
    for (int i = 0; i < ntasks; ++i) {
        accum += task_rank_pct[i];
        int cur = std::round(rank_n * accum);
        task_rank_n[i] = cur - prev;
            prev = cur;
    }
    AMREX_ASSERT(std::accumulate(task_rank_n.begin(),task_rank_n.end(),0) == rank_n);
}

void ForkJoin::reg_mf (MultiFab &mf, const std::string &name, int idx,
                     Strategy strategy, Intent intent, int owner)
{
    if (idx >= data[name].size()) {
        data[name].resize(idx + 1);
    }
    data[name][idx] = MFFork(&mf, strategy, intent, owner);
}

// this is called before ParallelContext::split
// the parent task is the top frame in ParallelContext's stack
void ForkJoin::copy_data_to_tasks ()
{
    if (flag_verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "Copying data into fork-join tasks ..." << std::endl;
    }
    for (auto &p : data) { // for each name
        const auto &mf_name = p.first;
        for (int idx = 0; idx < p.second.size(); ++idx) { // for each index
            auto &mff = p.second[idx];
            const MultiFab &orig = *mff.orig;
            const auto &ba = orig.boxArray();
            Vector<MultiFab> &forked = mff.forked;
            int comp_n = orig.nComp(); // number of components in original

            forked.reserve(task_n()); // does nothing if forked MFs already created
            for (int i = 0; i < task_n(); ++i) {
                // check if this task needs this MF
                if (mff.strategy != Strategy::single || i == mff.owner_task) {

                    // compute task's component lower and upper bound
                    int comp_lo, comp_hi;
                    if (mff.strategy == Strategy::split) {
                        // split components across tasks
                        comp_lo = comp_n *  i    / task_n();
                        comp_hi = comp_n * (i+1) / task_n();
                    } else {
                        // copy all components to task
                        comp_lo = 0;
                        comp_hi = comp_n;
                    }

                    // create task's MF if first time through
                    if (forked.size() <= i) {
                        if (flag_verbose && ParallelDescriptor::IOProcessor()) {
                            std::cout << "  Creating forked " << mf_name << "[" << idx << "] for task " << i
                                      << (mff.strategy == Strategy::split ? " (split)" : " (whole)") << std::endl;
                        }
                        // look up the distribution mapping for this (box array, task) pair
                        const DistributionMapping &dm = get_dm(ba, i);
                        forked.emplace_back(ba, dm, comp_hi - comp_lo, 0);
                    } else if (flag_verbose && ParallelDescriptor::IOProcessor()) {
                        std::cout << "  Forked " << mf_name << "[" << idx << "] for task " << i
                                  << " already created" << std::endl;
                    }
                    AMREX_ASSERT(i < forked.size());

                    // copy data if needed
                    if (mff.intent == Intent::in || mff.intent == Intent::inout) {
                        if (flag_verbose && ParallelDescriptor::IOProcessor()) {
                            std::cout << "    Copying " << mf_name << "[" << idx << "] into to task " << i << std::endl;
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
            AMREX_ASSERT(forked.size() == task_n());
        }
    }
}

// this is called after ParallelContext::unsplit
// the parent task is the top frame in ParallelContext's stack
void ForkJoin::copy_data_from_tasks () {
    if (flag_verbose && ParallelDescriptor::IOProcessor()) {
        std::cout << "Copying data out of fork-join tasks ..." << std::endl;
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
                    for (int i = 0; i < task_n(); ++i) {
                        if (flag_verbose && ParallelDescriptor::IOProcessor()) {
                            std::cout << "  Copying " << mf_name << "[" << idx << "] out from task " << i << "  (unsplit)" << std::endl;
                        }
                        int comp_lo = comp_n *  i    / task_n();
                        int comp_hi = comp_n * (i+1) / task_n();
                        orig.copy(forked[i], 0, comp_lo, comp_hi - comp_lo);
                    }
                } else { // mff.strategy == single or duplicate
                    // copy all components from owner_task
                    if (flag_verbose && ParallelDescriptor::IOProcessor()) {
                        std::cout << "Copying " << mf_name << " out from task " << mff.owner_task << "  (whole)" << std::endl;
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
const DistributionMapping & ForkJoin::get_dm (const BoxArray& ba, int task_idx)
{
    auto &dm_vec = dms[ba.getRefID()];

    if (dm_vec.size() == 0) {
        // new entry
        dm_vec.resize(task_n());
    }
    AMREX_ASSERT(task_idx < dm_vec.size());

    if (dm_vec[task_idx] == nullptr) {
        // create DM of current box array over current task's ranks
#if 0
        auto task_bounds = ParallelContext::get_split_bounds(task_rank_n);
        auto task_glo_rank_lo = ParallelContext::local_to_global_rank(task_bounds[task_idx].first);
        auto task_glo_rank_hi = ParallelContext::local_to_global_rank(task_bounds[task_idx].second);
        dm_vec[task_idx].reset(new DistributionMapping(ba, task_glo_rank_lo, task_glo_rank_hi));
#else
// xxxxx        // hard coded colors only right now
//        AMREX_ASSERT(task_rank_n.size() == ParallelDescriptor::NColors());
//        ParallelDescriptor::Color color = ParallelDescriptor::Color(task_idx);
        int nprocs = ParallelDescriptor::NProcs();
        dm_vec[task_idx].reset(new DistributionMapping(ba, nprocs));
#endif
        if (flag_verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "    Creating DM for (box array, task id) = ("
                      << ba.getRefID() << ", " << task_idx << ")" << std::endl;
        }
    } else {
        // DM has already been created
        if (flag_verbose && ParallelDescriptor::IOProcessor()) {
            std::cout << "    DM for (box array, task id) = (" << ba.getRefID() << ", " << task_idx
                      << ") already created" << std::endl;
        }
    }
    AMREX_ASSERT(dm_vec[task_idx] != nullptr);

    return *dm_vec[task_idx];
}

}
