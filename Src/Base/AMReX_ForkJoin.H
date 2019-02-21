#ifndef AMREX_FORKJOIN_H
#define AMREX_FORKJOIN_H

#include <AMReX_ParallelContext.H>
#include <AMReX_MultiFab.H>

namespace amrex {

class ForkJoin
{
  public:

    enum class Strategy {
        single,     //!< one task gets a copy of whole MF
        duplicate,  //!< all tasks get a copy of whole MF
        split,      //!< split MF components across tasks
    };
    enum class Intent { in, out, inout };

    //! if strategy == split, use this to specify how to split components across tasks
    struct ComponentSet {
        //! for now, just support contiguous range of components
        //! only requires a single MultiFab::Copy() per task this way
        //! in the future, consider changing this to a true set for arbitrary components
        ComponentSet() {}
        ComponentSet(int lo_, int hi_) : lo(lo_), hi(hi_) {}
        int lo; //!< inclusive bound
        int hi; //!< exclusive bound
    };

    ForkJoin (const Vector<int> &task_rank_n);

    ForkJoin (const Vector<double> &task_rank_pct);

    ForkJoin (int ntasks)
        : ForkJoin( Vector<double>(ntasks, 1.0 / ntasks) ) { }

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
                 Strategy strategy, Intent intent, int owner = -1);

    void reg_mf (MultiFab &mf, const std::string &name,
                 Strategy strategy, Intent intent, int owner = -1) {
        reg_mf(mf, name, 0, strategy, intent, owner);
    }

    //! these overloads are for in case the MultiFab argument is const
    //! intent must be in
    //!
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

    //! overload in case of vector of pointer to const MultiFab
    void reg_mf_vec (const Vector<MultiFab const *> &mfs, const std::string &name,
                     Strategy strategy, Intent intent, int owner = -1) {
        data[name].reserve(mfs.size());
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, intent, owner);
        }
    }

    //! modify the number of grow cells associated with the multifab
    void modify_ngrow (const std::string &name, int idx, IntVect ngrow);
    void modify_ngrow (const std::string &name, IntVect ngrow) {
        modify_ngrow(name, 0, ngrow);
    }

    //! modify how the multifab is split along components across the tasks
    void modify_split (const std::string &name, int idx, Vector<ComponentSet> comp_split);
    void modify_split (const std::string &name, Vector<ComponentSet> comp_split) {
        modify_split(name, 0, comp_split);
    }

    // TODO: may want to de-register MFs if they change across invocations

    MultiFab &get_mf (const std::string &name, int idx = 0) {
        AMREX_ASSERT_WITH_MESSAGE(data.count(name) > 0 && idx < data[name].size(), "get_mf(): name or index not found");
        AMREX_ASSERT(task_me >= 0 && task_me < data[name][idx].forked.size());
        return data[name][idx].forked[task_me];
    }

    //! vector of pointers to all MFs under a name
    Vector<MultiFab *> get_mf_vec (const std::string &name) {
        int dim = data.at(name).size();
        Vector<MultiFab *> result(dim);
        for (int idx = 0; idx < dim; ++idx) {
            result[idx] = &get_mf(name, idx);
        }
        return result;
    }

    void set_task_output_dir (const std::string & dir)
    {
        task_output_dir = dir;
    }

    static void set_task_output_file (const std::string & filename)
    {
        ParallelContext::set_last_frame_ofs(filename);
    }

    template <class F>
    void fork_join (const F &fn)
    {
        flag_invoked = true; // set invoked flag
        const int io_rank = 0; // team's sub-rank 0 does IO
        create_task_output_dir();
        MPI_Comm task_comm = split_tasks();
        copy_data_to_tasks(); // move data to local tasks
        ParallelContext::push(task_comm, task_me, io_rank);
        set_task_output_file(get_io_filename());
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
        Intent intent;
        int owner_task; //!< only used if strategy == single or duplicate
        IntVect ngrow;
        Vector<ComponentSet> comp_split; //!< if strategy == split, how to split components to tasks
        Vector<MultiFab> forked; //!< holds new multifab for each task in fork

        MFFork () = default;
        MFFork (const MFFork&) = delete;
        MFFork& operator= (const MFFork&) = delete;
        MFFork (MFFork&&) = default;
        MFFork& operator= (MFFork&&) = default;
        MFFork (MultiFab *omf, Strategy s, Intent i, int own,
                const IntVect &ng,
                Vector<ComponentSet> cs)
            : orig(omf), strategy(s), intent(i), owner_task(own), ngrow(ng), comp_split(std::move(cs)) {}

        bool empty() const { return orig == nullptr; }
    };

    bool flag_verbose = false; //!< for debugging
    bool flag_invoked = false; //!< track if object has been invoked yet
    Vector<int> split_bounds; //!< task i has ranks over the interval [result[i], result[i+1])
    int task_me = -1; //!< which forked task the rank belongs to
    std::map<BoxArray::RefID, Vector<std::unique_ptr<DistributionMapping>>> dms; //!< DM cache
    std::unordered_map<std::string, Vector<MFFork>> data;
    std::string task_output_dir = ""; //!< where to write task output

    void init(const Vector<int> &task_rank_n);

    //! multiple MultiFabs may share the same box array
    //! only compute the DM once per unique (box array, task) pair and cache it
    //! create map from box array RefID to vector of DistributionMapping indexed by task ID
    //!
    const DistributionMapping &get_dm (const BoxArray& ba, int task_idx,
                                       const DistributionMapping& dm_orig);

    //! this is called before ParallelContext::split
    //! the parent task is the top frame in ParallelContext's stack
    //!
    void copy_data_to_tasks ();

    //! this is called after ParallelContext::unsplit
    //! the parent task is the top frame in ParallelContext's stack
    //!
    void copy_data_from_tasks ();

    //! split top frame of stack
    MPI_Comm split_tasks ();

    //! create the task output directory
    void create_task_output_dir ();

    //! unique output file for this sub-task
    std::string get_io_filename (bool flag_unique = false);
};

}

#endif // AMREX_FORKJOIN_H
