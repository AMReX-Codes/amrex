#ifndef AMREX_FORKJOIN_H
#define AMREX_FORKJOIN_H

#include <AMReX_ParallelContext.H>
#include <AMReX_MultiFab.H>

namespace amrex {

class ForkJoin
{
  public:

    enum class Strategy {
        single,     // one task gets a copy of whole MF
        duplicate,  // all tasks get a copy of whole MF
        split,      // split MF components across tasks
    };
    enum class Intent { in, out, inout };

    ForkJoin (Vector<int> trn);

    ForkJoin (const Vector<double> &task_rank_pct);

    ForkJoin (int ntasks) : ForkJoin( Vector<double>(ntasks, 1.0 / ntasks) ) { }

    int task_n () const { return task_rank_n.size(); }

    void reg_mf (MultiFab &mf, const std::string &name, int idx,
                 Strategy strategy, Intent intent, int owner = -1);

    void reg_mf (MultiFab &mf, const std::string &name,
                 Strategy strategy, Intent intent, int owner = -1) {
        reg_mf(mf, name, 0, strategy, intent, owner);
    }

    // these overloads are for in case the MultiFab argument is const
    // intent must be in
    void reg_mf (const MultiFab &mf, const std::string &name, int idx,
                 Strategy strategy, Intent intent, int owner = -1) {
        AMREX_ASSERT_WITH_MESSAGE(intent == Intent::in,
                                  "const MultiFab must be registered read-only");
        reg_mf(const_cast<MultiFab&>(mf), name, idx, strategy, intent, owner);
    }
    void reg_mf (const MultiFab &mf, const std::string &name,
                 Strategy strategy, Intent intent, int owner = -1) {
        reg_mf(mf, name, 0, strategy, intent, owner);
    }

    void reg_mf_vec (const Vector<MultiFab *> &mfs, const std::string &name,
                     Strategy strategy, Intent intent, int owner = -1) {
        data[name].reserve(mfs.size());
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, intent, owner);
        }
    }

    // overload in case of vector of pointer to const MultiFab
    void reg_mf_vec (const Vector<MultiFab const *> &mfs, const std::string &name,
                     Strategy strategy, Intent intent, int owner = -1) {
        data[name].reserve(mfs.size());
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, intent, owner);
        }
    }

    // TODO: may want to de-register MFs if they change across invocations

    // this is called before ParallelContext::split
    // the parent task is the top frame in ParallelContext's stack
    void copy_data_to_tasks ();

    // this is called after ParallelContext::unsplit
    // the parent task is the top frame in ParallelContext's stack
    void copy_data_from_tasks ();

    template <class F>
    void fork_join (const F &fn)
    {
        copy_data_to_tasks(); // move data to local tasks
        task_me = ParallelContext::split(task_rank_n);
        fn(*this);
        ParallelContext::unsplit();
        copy_data_from_tasks(); // move local data back
    }

    int get_task_me () const { return task_me; }

    MultiFab &get_mf (const std::string &name, int idx = 0) {
        AMREX_ASSERT_WITH_MESSAGE(data.count(name) > 0 && idx < data[name].size(), "get_mf(): name or index not found");
        AMREX_ASSERT(task_me >= 0 && task_me < data[name][idx].forked.size());
        return data[name][idx].forked[task_me];
    }

    // vector of pointers to all MFs under a name
    Vector<MultiFab *> get_mf_vec (const std::string &name) {
        int dim = data.at(name).size();
        Vector<MultiFab *> result(dim);
        for (int idx = 0; idx < dim; ++idx) {
            result[idx] = &get_mf(name, idx);
        }
        return result;
    }

  private:

    struct MFFork
    {
        MultiFab *orig;
        Strategy strategy;
        Intent intent;
        int owner_task; // only used if strategy == single or duplicate
        Vector<MultiFab> forked; // holds new multifab for each task in fork

        MFFork () = default;
        MFFork (const MFFork&) = delete;
        MFFork& operator= (const MFFork&) = delete;
        MFFork (MFFork&&) = default;
        MFFork& operator= (MFFork&&) = default;
        MFFork (MultiFab *omf, Strategy s = Strategy::duplicate,
                Intent i = Intent::inout, int own = -1)
            : orig(omf), strategy(s), intent(i), owner_task(own) {}
    };

    Vector<int> task_rank_n; // number of ranks in each forked task
    int task_me = -1; // task the current rank belongs to
    std::map<BoxArray::RefID, Vector<std::unique_ptr<DistributionMapping>>> dms; // DM cache
    std::unordered_map<std::string, Vector<MFFork>> data;
    const static bool flag_verbose = false; // for debugging

    // multiple MultiFabs may share the same box array
    // only compute the DM once per unique (box array, task) pair and cache it
    // create map from box array RefID to vector of DistributionMapping indexed by task ID
    const DistributionMapping &get_dm (const BoxArray& ba, int task_idx);
};

}

#endif // AMREX_FORKJOIN_H
