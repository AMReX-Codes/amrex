#ifndef AMREX_FORKJOIN_H
#define AMREX_FORKJOIN_H

#include <AMReX_ParallelContext.H>
#include <AMReX_MultiFab.H>

namespace detail {

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

class ForkJoin
{
  public:

    enum class Strategy {
        SINGLE,     // one task gets a copy of whole MF
        DUPLICATE,  // all tasks get a copy of whole MF
        SPLIT,      // split MF components across tasks
    };
    enum class Access { RD, WR, RW };

    struct MFFork
    {
        MultiFab *orig;
        Strategy strategy;
        Access access;
        int owner_task; // only used if access == SINGLE or DUPLICATE
        Vector<MultiFab> forked; // holds new multifab for each task in fork

        MFFork() = default;
        MFFork(MultiFab *omf, Strategy s = Strategy::DUPLICATE,
               Access a = Access::RW, int own = -1)
            : orig(omf), strategy(s), access(a), owner_task(own) {}
    };

    ForkJoin(Vector<int> trn) : task_rank_n(std::move(trn)) {
        auto rank_n = ParallelContext::rank_n(); // number of ranks in current frame
        AMREX_ASSERT(task_n() >= 2);
        AMREX_ASSERT(detail::sum(task_rank_n) == rank_n);
    }

    ForkJoin(const Vector<double> &task_rank_pct) {
        auto rank_n = ParallelContext::rank_n(); // number of ranks in current frame
        auto task_n = task_rank_pct.size();
        AMREX_ASSERT(task_n >= 2);
        task_rank_n.resize(task_n);
        int prev = 0;
        double accum = 0;
        for (int i = 0; i < task_n; ++i) {
            accum += task_rank_pct[i];
            int cur = std::round(rank_n * accum);
            task_rank_n[i] = cur - prev;
            prev = cur;
        }
        AMREX_ASSERT(detail::sum(task_rank_n) == rank_n);
    }

    ForkJoin(int task_n = 2) : ForkJoin( Vector<double>(task_n, 1.0 / task_n) ) { }

    int task_n() const { return task_rank_n.size(); }

    void reg_mf(MultiFab &mf, std::string name, int idx,
                Strategy strategy, Access access, int owner = -1) {
        if (idx >= data[name].size()) {
            data[name].resize(idx + 1);
        }
        data[name][idx] = MFFork(&mf, strategy, access, owner);
    }

    void reg_mf(MultiFab &mf, std::string name,
                Strategy strategy, Access access, int owner = -1) {
        reg_mf(mf, name, 0, strategy, access, owner);
    }

    void reg_mf_vec(const Vector<MultiFab *> &mfs, std::string name,
                    Strategy strategy, Access access, int owner = -1) {
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, access, owner);
        }
    }

    // these overloads are for in case the MultiFab argument is const
    // access must be read-only
    void reg_mf(const MultiFab &mf, std::string name, int idx,
                Strategy strategy, Access access, int owner = -1) {
        AMREX_ASSERT_WITH_MESSAGE(access == Access::RD,
                                  "const MultiFab must be registered read-only");
        reg_mf(*const_cast<MultiFab *>(&mf), name, idx, strategy, access, owner);
    }
    void reg_mf(const MultiFab &mf, std::string name,
                Strategy strategy, Access access, int owner = -1) {
        reg_mf(mf, name, 0, strategy, access, owner);
    }

    // TODO: may want to de-register MFs if they change across invocations

    // this is called before ParallelContext::split
    // the parent task is the top frame in ParallelContext's stack
    void copy_data_to_tasks();

    // this is called after ParallelContext::unsplit
    // the parent task is the top frame in ParallelContext's stack
    void copy_data_from_tasks();

    template <class F>
    void fork_join(const F &fn)
    {
        copy_data_to_tasks(); // move data to local tasks
        task_me = ParallelContext::split(task_rank_n);
        fn(*this);
        ParallelContext::unsplit();
        copy_data_from_tasks(); // move local data back
    }

    int get_task_me() const { return task_me; }

    MultiFab &get_mf(std::string name, int idx = 0) {
        AMREX_ASSERT_WITH_MESSAGE(data.count(name) > 0 && idx < data[name].size(), "get_mf(): name or index not found");
        AMREX_ASSERT(task_me >= 0 && task_me < data[name][idx].forked.size());
        return data[name][idx].forked[task_me];
    }

    // vector of pointers to all MFs under a name
    Vector<MultiFab *> get_mf_vec(std::string name) {
        int dim = data.at(name).size();
        Vector<MultiFab *> result(dim);
        for (int idx = 0; idx < dim; ++idx) {
            result[idx] = &get_mf(name, idx);
        }
        return result;
    }

  private:
    Vector<int> task_rank_n; // number of ranks in each forked task
    int task_me = -1; // task the current rank belongs to
    std::map<BoxArray::RefID, Vector<std::unique_ptr<DistributionMapping>>> dms; // DM cache
    std::unordered_map<std::string, Vector<MFFork>> data;
    const static bool flag_verbose = false; // for debugging

    // multiple MultiFabs may share the same box array
    // only compute the DM once per unique (box array, task) pair and cache it
    // create map from box array RefID to vector of DistributionMapping indexed by task ID
    const DistributionMapping &get_dm(BoxArray ba, int task_idx);
};

}

#endif // AMREX_FORKJOIN_H
