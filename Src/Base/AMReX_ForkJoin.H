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

    // if strategy == split, use this to specify how to split components across tasks
    struct ComponentSet {
        // for now, just support contiguous range of components
        // only requires a single MultiFab::Copy() per task this way
        // in the future, consider changing this to a true set for arbitrary components
        ComponentSet() {}
        ComponentSet(int lo_, int hi_) : lo(lo_), hi(hi_) {}
        int lo; // inclusive bound
        int hi; // exclusive bound
    };

    ForkJoin (const Vector<int> &task_rank_n, const std::string &a_task_output_dir = "");

    ForkJoin (const Vector<double> &task_rank_pct, const std::string &a_task_output_dir = "");

    ForkJoin (int ntasks, const std::string& a_task_output_dir = "")
        : ForkJoin( Vector<double>(ntasks, 1.0 / ntasks), a_task_output_dir ) { }

    int NTasks () const { return split_bounds.size() - 1; }

    int MyTask () const { return task_me; }

    bool Verbose () const { return flag_verbose; }

    void SetVerbose (bool verbose_in) { flag_verbose = verbose_in; }

    ComponentSet ComponentBounds(const std::string& name, int idx=0) const;

    int NProcsTask (int task) const {
        AMREX_ASSERT(task + 1 < split_bounds.size());
        return split_bounds[task + 1] - split_bounds[task];
    }

    void reg_mf (MultiFab &mf, const std::string &name, int idx,
                 Strategy strategy, Intent intent, const IntVect ng = {AMREX_D_DECL(0,0,0)}, int owner = -1);

    void reg_mf (MultiFab &mf, const std::string &name,
                 Strategy strategy, Intent intent, const IntVect ng = {AMREX_D_DECL(0,0,0)}, int owner = -1) {
        reg_mf(mf, name, 0, strategy, intent, ng, owner);
    }

    // these overloads are for in case the MultiFab argument is const
    // intent must be in
    void reg_mf (const MultiFab &mf, const std::string &name, int idx,
                 Strategy strategy, Intent intent, const IntVect ng = {AMREX_D_DECL(0,0,0)}, int owner = -1) {
        AMREX_ASSERT_WITH_MESSAGE(intent == Intent::in,
                                  "const MultiFab must be registered read-only");
        reg_mf(const_cast<MultiFab&>(mf), name, idx, strategy, intent, ng, owner);
    }
    void reg_mf (const MultiFab &mf, const std::string &name,
                 Strategy strategy, Intent intent, const IntVect ng = {AMREX_D_DECL(0,0,0)}, int owner = -1) {
        reg_mf(mf, name, 0, strategy, intent, ng, owner);
    }

    void reg_mf_vec (const Vector<MultiFab *> &mfs, const std::string &name,
                     Strategy strategy, Intent intent, const IntVect ng = {AMREX_D_DECL(0,0,0)}, int owner = -1) {
        data[name].reserve(mfs.size());
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, intent, ng, owner);
        }
    }

    // overload in case of vector of pointer to const MultiFab
    void reg_mf_vec (const Vector<MultiFab const *> &mfs, const std::string &name,
                     Strategy strategy, Intent intent, const IntVect ng = {AMREX_D_DECL(0,0,0)}, int owner = -1) {
        data[name].reserve(mfs.size());
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, intent, ng, owner);
        }
    }

    // modify how the multifab is split along components across the tasks
    void modify_split (const std::string &name, int idx, Vector<ComponentSet> comp_split);

    // TODO: may want to de-register MFs if they change across invocations

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

    static void set_task_output_file (const std::string & filename)
    {
        ParallelContext::set_last_frame_ofs(filename);
    }

    template <class F>
    void fork_join (const F &fn)
    {
        const int io_rank = 0; // team's sub-rank 0 does IO
        MPI_Comm task_comm = split_tasks();
        copy_data_to_tasks(); // move data to local tasks
        ParallelContext::push(task_comm, task_me, io_rank);
        set_task_output_file(get_fresh_io_filename());
        fn(*this);
        ParallelContext::pop();
#ifdef BL_USE_MPI
        MPI_Comm_free(&task_comm);
#endif
        copy_data_from_tasks(); // move local data back
    }

  private:

    struct MFFork
    {
        MultiFab *orig = nullptr;
        Strategy strategy;
        Vector<ComponentSet> comp_split; // if strategy == split, how to split components to tasks
        Intent intent;
        IntVect nGrow;
        int owner_task; // only used if strategy == single or duplicate
        Vector<MultiFab> forked; // holds new multifab for each task in fork

        MFFork () = default;
        MFFork (const MFFork&) = delete;
        MFFork& operator= (const MFFork&) = delete;
        MFFork (MFFork&&) = default;
        MFFork& operator= (MFFork&&) = default;
        MFFork (MultiFab *omf, Strategy s = Strategy::duplicate,
                Intent i = Intent::inout, const IntVect ng = {AMREX_D_DECL(0,0,0)}, int own = -1)
            : orig(omf), strategy(s), intent(i), nGrow(ng), owner_task(own) {}

        bool empty() const { return orig == nullptr; }
    };

    Vector<int> split_bounds; // task i has ranks over the interval [result[i], result[i+1])
    int task_me = -1; // which forked task the rank belongs to
    std::map<BoxArray::RefID, Vector<std::unique_ptr<DistributionMapping>>> dms; // DM cache
    std::unordered_map<std::string, Vector<MFFork>> data;
    bool flag_verbose = false; // for debugging
    std::string task_output_dir;

    void init(const Vector<int> &task_rank_n, const std::string &task_output_dir);

    // multiple MultiFabs may share the same box array
    // only compute the DM once per unique (box array, task) pair and cache it
    // create map from box array RefID to vector of DistributionMapping indexed by task ID
    const DistributionMapping &get_dm (const BoxArray& ba, int task_idx,
                                       const DistributionMapping& dm_orig);

    // this is called before ParallelContext::split
    // the parent task is the top frame in ParallelContext's stack
    void copy_data_to_tasks ();

    // this is called after ParallelContext::unsplit
    // the parent task is the top frame in ParallelContext's stack
    void copy_data_from_tasks ();

    // split top frame of stack
    MPI_Comm split_tasks ();

    // unique output file for this sub-task
    std::string get_fresh_io_filename ();
};

}

#endif // AMREX_FORKJOIN_H
