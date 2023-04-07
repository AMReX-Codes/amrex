#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H
#include <AMReX_Config.H>

#include <AMReX_ccse-mpi.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_BLBackTrace.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLassert.H>
#include <AMReX_Extension.H>
#include <AMReX_INT.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_ValLocPair.H>

#ifndef BL_AMRPROF
#include <AMReX_Box.H>
#endif

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_PMI
#include <pmi.h>
#endif

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <algorithm>
#include <atomic>
#include <csignal>
#include <functional>
#include <limits>
#include <numeric>
#include <string>
#include <typeinfo>
#include <vector>

namespace amrex {

template <typename T> class LayoutData;

//! Parallel frontend that abstracts functionalities needed to spawn processes and handle communication
namespace ParallelDescriptor
{
    //
    //! Used as default argument to ParallelDescriptor::Barrier().
    const std::string Unnamed("Unnamed");

    //! Hold the description and status of communication data
    class Message
    {
    public:

        Message () :
            m_finished(true),
            m_type(MPI_DATATYPE_NULL),
            m_req(MPI_REQUEST_NULL) {}
        Message (MPI_Request req_, MPI_Datatype type_) :
            m_finished(false),
            m_type(type_),
            m_req(req_) {}
        Message (MPI_Status stat_, MPI_Datatype type_) :
            m_finished(true),
            m_type(type_),
            m_req(MPI_REQUEST_NULL), m_stat(stat_) {}
        void wait ();
        bool test ();
        size_t count () const;
        int tag () const;
        int pid () const;
        MPI_Datatype type () const { return m_type; }
        MPI_Request  req () const { return m_req; }
        MPI_Status   stat () const { return m_stat; }

    private:

        bool               m_finished; // NOLINT(modernize-use-default-member-init)
        MPI_Datatype       m_type;
        MPI_Request        m_req;
        mutable MPI_Status m_stat;
    };

#ifdef BL_USE_MPI
    void MPI_Error(const char* file, int line, const char* str, int rc);

#define BL_MPI_REQUIRE(x)                                                \
do                                                                        \
{                                                                        \
  if ( int l_status_ = (x) )                                                \
    {                                                                        \
        amrex::ParallelDescriptor::MPI_Error(__FILE__,__LINE__,#x, l_status_); \
    }                                                                        \
}                                                                        \
while ( false )
#endif

    /**
    * \brief Perform any needed parallel initialization.  This MUST be the
    * first routine in this class called from within a program.
    */
    void StartParallel (int*    argc = nullptr,
                        char*** argv = nullptr,
                        MPI_Comm mpi_comm = MPI_COMM_WORLD);

    void Initialize ();
    void Finalize ();

    extern AMREX_EXPORT int use_gpu_aware_mpi;
    inline bool UseGpuAwareMpi () { return use_gpu_aware_mpi; }

    //! Split the process pool into teams
    void StartTeams ();
    void EndTeams ();

    /**
    * \brief Perform any needed parallel finalization.  This MUST be the
    * last routine in this class called from within a program.
    */
    void EndParallel ();

    //! return the rank number local to the current Parallel Context
    inline int
    MyProc () noexcept
    {
        return ParallelContext::MyProcAll();
    }
    inline int
    MyProc (MPI_Comm comm) noexcept
    {
#ifdef BL_USE_MPI
        int r;
        MPI_Comm_rank(comm,&r);
        return r;
#else
        amrex::ignore_unused(comm);
        return 0;
#endif
    }

    //! Provide functionalities needed to construct a team of processes to perform a particular job
    struct ProcessTeam
    {
        typedef MPI_Comm team_t;

        //! synchronize processes within the team
        void Barrier () const {
            if (m_size > 1) {
#if defined(BL_USE_MPI3)
                MPI_Barrier(m_team_comm);
#endif
            }
        }

        //! memory fence
        void MemoryBarrier () const {
            if (m_size > 1) {
#ifdef AMREX_USE_OMP
                if (omp_in_parallel()) {
                    #pragma omp barrier
                }
                #pragma omp master
#endif
                {
#if defined(BL_USE_MPI3)
                    std::atomic_thread_fence(std::memory_order_release);
                    MPI_Barrier(m_team_comm);
                    std::atomic_thread_fence(std::memory_order_acquire);
#endif
                }
            }
        }

        //! free a communicator
        void clear () {
#if defined(BL_USE_MPI3)
            if (m_size > 1) {
                MPI_Comm_free(&m_team_comm);
                if (m_rankInTeam==0) MPI_Comm_free(&m_lead_comm);
            }
#endif
        }

        [[nodiscard]] const team_t& get() const {
            return m_team_comm;
        }
        //! return the communicator
        [[nodiscard]] const MPI_Comm& get_team_comm() const { return m_team_comm; }
        [[nodiscard]] const MPI_Comm& get_lead_comm() const { return m_lead_comm; }

        int    m_numTeams;
        int    m_size;
        int    m_color;
        int    m_lead;
        int    m_rankInTeam;
        int    m_do_team_reduce;

        MPI_Comm  m_team_comm;
        MPI_Comm  m_lead_comm;
    };

    extern AMREX_EXPORT ProcessTeam m_Team;

    extern AMREX_EXPORT int m_MinTag, m_MaxTag;
    inline int MinTag () noexcept { return m_MinTag; }
    inline int MaxTag () noexcept { return m_MaxTag; }

    extern AMREX_EXPORT MPI_Comm m_comm;
    inline MPI_Comm Communicator () noexcept { return m_comm; }

    extern AMREX_EXPORT int m_nprocs_per_node;
    //! Return the number of MPI ranks per node.  This may not be correct if
    //! MPI_COMM_TYPE_SHARED groups processes across nodes.
    inline int NProcsPerNode () noexcept { return m_nprocs_per_node; }

#ifdef AMREX_USE_MPI
    extern Vector<MPI_Datatype*> m_mpi_types;
    extern Vector<MPI_Op*> m_mpi_ops;
#endif

    //! return the number of MPI ranks local to the current Parallel Context
    inline int
    NProcs () noexcept
    {
        return ParallelContext::NProcsAll();
    }

    inline int
    NProcs (MPI_Comm comm) noexcept
    {
#ifdef BL_USE_MPI
        int s;
        BL_MPI_REQUIRE(MPI_Comm_size(comm, &s));
        return s;
#else
        amrex::ignore_unused(comm);
        return 1;
#endif
    }
    /**
    * \brief The MPI rank number of the I/O Processor (probably rank 0).
    * This rank is usually used to write to stdout.
    */
    extern AMREX_EXPORT const int ioProcessor;
    inline int
    IOProcessorNumber () noexcept
    {
        return ioProcessor;
    }
    /**
    * \brief Is this CPU the I/O Processor?
    * To get the rank number, call IOProcessorNumber()
    */
    inline bool
    IOProcessor () noexcept
    {
         return MyProc() == IOProcessorNumber();
    }

    inline int
    IOProcessorNumber (MPI_Comm comm) noexcept
    {
        return (comm == ParallelDescriptor::Communicator()) ? ioProcessor : 0;
    }

    inline bool
    IOProcessor (MPI_Comm comm) noexcept
    {
        return MyProc(comm) == IOProcessorNumber(comm);
    }

    //
    inline int
    TeamSize () noexcept
    {
        return m_Team.m_size;
    }
    inline int
    NTeams () noexcept
    {
        return m_Team.m_numTeams;
    }
    inline int
    MyTeamColor () noexcept
    {
        return m_Team.m_color;
    }
    inline int
    MyTeamLead () noexcept
    {
        return m_Team.m_lead;
    }
    inline int
    MyRankInTeam () noexcept
    {
        return m_Team.m_rankInTeam;
    }
    inline int
    TeamLead (int rank) noexcept
    {
        return (rank >= 0) ? (rank - rank % m_Team.m_size) : MPI_PROC_NULL;
    }
    inline bool
    isTeamLead () noexcept
    {
        return MyRankInTeam() == 0;
    }
    inline bool
    sameTeam (int rank) noexcept
    {
        return MyTeamLead() == TeamLead(rank);
    }
    inline bool
    sameTeam (int rankA, int rankB) noexcept
    {
        return TeamLead(rankA) == TeamLead(rankB);
    }
    inline int
    RankInLeadComm (int rank) noexcept
    {
        return (rank >= 0) ? (rank / m_Team.m_size) : MPI_PROC_NULL;
    }
    inline bool
    doTeamReduce () noexcept
    {
        return m_Team.m_do_team_reduce;
    }
    inline const ProcessTeam&
    MyTeam () noexcept
    {
        return m_Team;
    }
    inline std::pair<int,int>
    team_range (int begin, int end, int rit = -1, int nworkers = 0) noexcept
    {
        int rb, re;
        {
            if (rit < 0) rit = ParallelDescriptor::MyRankInTeam();
            if (nworkers == 0) nworkers = ParallelDescriptor::TeamSize();
            BL_ASSERT(rit<nworkers);
            if (nworkers == 1) {
                rb = begin;
                re = end;
            } else {
                int ntot = end - begin;
                int nr   = ntot / nworkers;
                int nlft = ntot - nr * nworkers;
                if (rit < nlft) {  // get nr+1 items
                    rb = begin + rit * (nr + 1);
                    re = rb + nr + 1;
                } else {           // get nr items
                    rb = begin + rit * nr + nlft;
                    re = rb + nr;
                }
            }
        }

#ifdef AMREX_USE_OMP
        int nthreads = omp_get_num_threads();
        if (nthreads > 1) {
            int tid = omp_get_thread_num();
            int ntot = re - rb;
            int nr   = ntot / nthreads;
            int nlft = ntot - nr * nthreads;
            if (tid < nlft) {  // get nr+1 items
                rb += tid * (nr + 1);
                re = rb + nr + 1;
            } else {           // get nr items
                rb += tid * nr + nlft;
                re = rb + nr;
            }
        }
#endif

        return std::make_pair(rb,re);
    }
    template <typename F>
    void team_for (int begin, int end, const F& f)
    {
        const auto& range = team_range(begin, end);
        for (int i = range.first; i < range.second; ++i) {
            f(i);
        }
    }
    template <typename F>               // rit: rank in team
    void team_for (int begin, int end, int rit, const F& f)
    {
        const auto& range = team_range(begin, end, rit);
        for (int i = range.first; i < range.second; ++i) {
            f(i);
        }
    }
    template <typename F>               // rit: rank in team
    void team_for (int begin, int end, int rit, int nworkers, const F& f)
    {
        const auto& range = team_range(begin, end, rit, nworkers);
        for (int i = range.first; i < range.second; ++i) {
            f(i);
        }
    }

    void Barrier (const std::string& message = Unnamed);
    void Barrier (const MPI_Comm &comm, const std::string& message = Unnamed);
    Message Abarrier ();
    Message Abarrier (const MPI_Comm &comm);

    void Test (MPI_Request& request, int& flag, MPI_Status& status);
    void Test (Vector<MPI_Request>& request, int& flag, Vector<MPI_Status>& status);

    void Comm_dup (MPI_Comm comm, MPI_Comm& newcomm);
    //! Abort with specified error code.
    void Abort (int errorcode = SIGABRT, bool backtrace = true);
    //! ErrorString return string associated with error internal error condition
    const char* ErrorString (int errorcode);
    //! Returns wall-clock seconds since start of execution.
    double second () noexcept;

    //! And-wise boolean reduction.
    void ReduceBoolAnd (bool& rvar);
    //! And-wise boolean reduction to specified cpu.
    void ReduceBoolAnd (bool& rvar, int cpu);

    //! Or-wise boolean reduction.
    void ReduceBoolOr  (bool& rvar);
    //! Or-wise boolean reduction to specified cpu.
    void ReduceBoolOr  (bool& rvar, int cpu);

    //! Real sum reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T& rvar);

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T* rvar, int cnt);

    // Having this for backward compatibility
    void ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar);
    //
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (Vector<std::reference_wrapper<T> >&& rvar);

    //! Real sum reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T& rvar, int cpu);

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T* rvar, int cnt, int cpu);

    // Having this for backward compatibility
    void ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar, int cpu);
    //
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (Vector<std::reference_wrapper<T> >&& rvar, int cpu);

    //! Real max reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T& rvar);

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T* rvar, int cnt);

    // Having this for backward compatibility
    void ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar);
    //
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (Vector<std::reference_wrapper<T> >&& rvar);

    //! Real max reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T& rvar, int cpu);

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T* rvar, int cnt, int cpu);

    // Having this for backward compatibility
    void ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar, int cpu);
    //
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (Vector<std::reference_wrapper<T> >&& rvar, int cpu);

    //! Real min reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T& rvar);

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T* rvar, int cnt);

    // Having this for backward compatibility
    void ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar);
    //
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (Vector<std::reference_wrapper<T> >&& rvar);

    //! Real min reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T& rvar, int cpu);

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T* rvar, int cnt, int cpu);

    // Having this for backward compatibility
    void ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar, int cpu);
    //
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (Vector<std::reference_wrapper<T> >&& rvar, int cpu);

    //! Integer sum reduction.
    void ReduceIntSum (int& rvar);
    void ReduceIntSum (int* rvar, int cnt);
    void ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar);
    //! Integer sum reduction to specified cpu.
    void ReduceIntSum (int& rvar, int cpu);
    void ReduceIntSum (int* rvar, int cnt, int cpu);
    void ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar, int cpu);

    //! Integer max reduction.
    void ReduceIntMax (int& rvar);
    void ReduceIntMax (int* rvar, int cnt);
    void ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar);
    //! Integer max reduction to specified cpu.
    void ReduceIntMax (int& rvar, int cpu);
    void ReduceIntMax (int* rvar, int cnt, int cpu);
    void ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar, int cpu);

    //! Integer min reduction.
    void ReduceIntMin (int& rvar);
    void ReduceIntMin (int* rvar, int cnt);
    void ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar);
    //! Integer min reduction to specified cpu.
    void ReduceIntMin (int& rvar, int cpu);
    void ReduceIntMin (int* rvar, int cnt, int cpu);
    void ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar, int cpu);

    //! Long sum reduction.
    void ReduceLongSum (Long& rvar);
    void ReduceLongSum (Long* rvar, int cnt);
    void ReduceLongSum (Vector<std::reference_wrapper<Long> >&& rvar);
    //! Long sum reduction to specified cpu.
    void ReduceLongSum (Long& rvar, int cpu);
    void ReduceLongSum (Long* rvar, int cnt, int cpu);
    void ReduceLongSum (Vector<std::reference_wrapper<Long> >&& rvar, int cpu);

    //! Long max reduction.
    void ReduceLongMax (Long& rvar);
    void ReduceLongMax (Long* rvar, int cnt);
    void ReduceLongMax (Vector<std::reference_wrapper<Long> >&& rvar);
    //! Long max reduction to specified cpu.
    void ReduceLongMax (Long& rvar, int cpu);
    void ReduceLongMax (Long* rvar, int cnt, int cpu);
    void ReduceLongMax (Vector<std::reference_wrapper<Long> >&& rvar, int cpu);

    //! Long min reduction.
    void ReduceLongMin (Long& rvar);
    void ReduceLongMin (Long* rvar, int cnt);
    void ReduceLongMin (Vector<std::reference_wrapper<Long> >&& rvar);
    //! Long min reduction to specified cpu.
    void ReduceLongMin (Long& rvar, int cpu);
    void ReduceLongMin (Long* rvar, int cnt, int cpu);
    void ReduceLongMin (Vector<std::reference_wrapper<Long> >&& rvar, int cpu);

    //! Long and-wise reduction.
    void ReduceLongAnd (Long& rvar);
    void ReduceLongAnd (Long* rvar, int cnt);
    void ReduceLongAnd (Vector<std::reference_wrapper<Long> >&& rvar);
    //! Long and-wise reduction to specified cpu.
    void ReduceLongAnd (Long& rvar, int cpu);
    void ReduceLongAnd (Long* rvar, int cnt, int cpu);
    void ReduceLongAnd (Vector<std::reference_wrapper<Long> >&& rvar, int cpu);

    //! Parallel gather.
    void Gather (Real const* sendbuf, int nsend, Real* recvbuf, int root);
    /**
    * \brief Returns sequential message sequence numbers, usually used as
    * tags for send/recv.
    */
    inline int SeqNum () noexcept { return ParallelContext::get_inc_mpi_tag(); }

    template <class T> Message Asend(const T*, size_t n, int pid, int tag);
    template <class T> Message Asend(const T*, size_t n, int pid, int tag, MPI_Comm comm);
    template <class T> Message Asend(const std::vector<T>& buf, int pid, int tag);

    template <class T> Message Arecv(T*, size_t n, int pid, int tag);
    template <class T> Message Arecv(T*, size_t n, int pid, int tag, MPI_Comm comm);
    template <class T> Message Arecv(std::vector<T>& buf, int pid, int tag);

    template <class T> Message Send(const T* buf, size_t n, int dst_pid, int tag);
    template <class T> Message Send(const T* buf, size_t n, int dst_pid, int tag, MPI_Comm comm);
    template <class T> Message Send(const std::vector<T>& buf, int dst_pid, int tag);

    template <class T> Message Recv(T*, size_t n, int pid, int tag);
    template <class T> Message Recv(T*, size_t n, int pid, int tag, MPI_Comm comm);
    template <class T> Message Recv(std::vector<T>& buf, int pid, int tag);

    template <class T> void Bcast(T*, size_t n, int root = 0);
    template <class T> void Bcast(T*, size_t n, int root, const MPI_Comm &comm);
    void Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

    template <class T, class T1> void Scatter(T*, size_t n, const T1*, size_t n1, int root);

    template <class T, class T1> void Gather(const T*, size_t n, T1*, size_t n1, int root);
    template <class T> std::vector<T> Gather(const T&, int root);

    template <class T> void Gatherv (const T* send, int sc,
                                     T* recv, const std::vector<int>& rc, const std::vector<int>& disp,
                                     int root);

    //! Gather LayoutData values to a vector on root
    template <class T> void GatherLayoutDataToVector (const LayoutData<T>& sendbuf,
                                                      Vector<T>& recvbuf,
                                                      int root);

    void Wait     (MPI_Request& req, MPI_Status& status);
    void Waitall  (Vector<MPI_Request>& reqs, Vector<MPI_Status>& status);
    void Waitany  (Vector<MPI_Request>& reqs, int &index, MPI_Status& status);
    void Waitsome (Vector<MPI_Request>&, int&, Vector<int>&, Vector<MPI_Status>&);

    void ReadAndBcastFile(const std::string &filename, Vector<char> &charBuf,
                          bool bExitOnError = true,
                          const MPI_Comm &comm = Communicator() );
    void IProbe(int src_pid, int tag, int &mflag, MPI_Status &status);
    void IProbe(int src_pid, int tag, MPI_Comm comm, int &mflag, MPI_Status &status);

    /** Convert an MPI_THREAD_<X> level to string
     *
     * @param mtlev MPI_THREAD_<X> level
     * @return string representation of the equivalent MPI macro name
     */
    std::string
    mpi_level_to_string (int mtlev);

    // PMI = Process Management Interface, available on Crays. Provides API to
    // query topology of the job.
#ifdef AMREX_PMI
    void PMI_Initialize();
    void PMI_PrintMeshcoords(const pmi_mesh_coord_t *pmi_mesh_coord);
#endif

#ifdef BL_USE_MPI
    int select_comm_data_type (std::size_t nbytes);
    std::size_t alignof_comm_data (std::size_t nbytes);
#endif
}
}

namespace amrex {

#ifdef BL_USE_MPI

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend (const T* buf, size_t n, int dst_pid, int tag)
{
    return Asend(buf, n, dst_pid, tag, Communicator());
}

namespace ParallelDescriptor {  // Have to use namespace here to work around a gcc bug
template <>
Message
Asend<char> (const char* buf, size_t n, int dst_pid, int tag, MPI_Comm comm);

template <class T>
Message
Asend (const T* buf, size_t n, int dst_pid, int tag, MPI_Comm comm)
{
    static_assert(!std::is_same<char,T>::value, "Asend: char version has been specialized");

    BL_PROFILE_T_S("ParallelDescriptor::Asend(TsiiM)", T);
    BL_COMM_PROFILE(BLProfiler::AsendTsiiM, n * sizeof(T), dst_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(buf),
                              n,
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              comm,
                              &req) );
    BL_COMM_PROFILE(BLProfiler::AsendTsiiM, BLProfiler::AfterCall(), dst_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend (const std::vector<T>& buf, int dst_pid, int tag)
{
    return Asend(buf.data(), buf.size(), dst_pid, tag, Communicator());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send (const T* buf, size_t n, int dst_pid, int tag)
{
    return Send(buf, n, dst_pid, tag, Communicator());
}

namespace ParallelDescriptor {  // Have to use namespace here to work around a gcc bug
template <>
Message
Send<char> (const char* buf, size_t n, int dst_pid, int tag, MPI_Comm comm);

template <class T>
Message
Send (const T* buf, size_t n, int dst_pid, int tag, MPI_Comm comm)
{
    static_assert(!std::is_same<char,T>::value, "Send: char version has been specialized");

    BL_PROFILE_T_S("ParallelDescriptor::Send(Tsii)", T);

#ifdef BL_COMM_PROFILING
    int dst_pid_world(-1);
    MPI_Group groupComm, groupWorld;
    BL_MPI_REQUIRE( MPI_Comm_group(comm, &groupComm) );
    BL_MPI_REQUIRE( MPI_Comm_group(Communicator(), &groupWorld) );
    BL_MPI_REQUIRE( MPI_Group_translate_ranks(groupComm, 1, &dst_pid, groupWorld, &dst_pid_world) );

    BL_COMM_PROFILE(BLProfiler::SendTsii, n * sizeof(T), dst_pid_world, tag);
#endif

    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(buf),
                             n,
                             Mpi_typemap<T>::type(),
                             dst_pid,
                             tag,
                             comm) );
    BL_COMM_PROFILE(BLProfiler::SendTsii, BLProfiler::AfterCall(), dst_pid, tag);
    return Message();
}
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send (const std::vector<T>& buf, int dst_pid, int tag)
{
    return Send(buf.data(), buf.size(), dst_pid, tag, Communicator());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (T* buf, size_t n, int src_pid, int tag)
{
    return Arecv(buf, n, src_pid, tag, Communicator());
}

namespace ParallelDescriptor {  // Have to use namespace here to work around a gcc bug
template <>
Message
Arecv<char> (char* buf, size_t n, int src_pid, int tag, MPI_Comm comm);

template <class T>
Message
Arecv (T* buf, size_t n, int src_pid, int tag, MPI_Comm comm)
{
    static_assert(!std::is_same<char,T>::value, "Arecv: char version has been specialized");

    BL_PROFILE_T_S("ParallelDescriptor::Arecv(TsiiM)", T);
    BL_COMM_PROFILE(BLProfiler::ArecvTsiiM, n * sizeof(T), src_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(buf,
                              n,
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              comm,
                              &req) );
    BL_COMM_PROFILE(BLProfiler::ArecvTsiiM, BLProfiler::AfterCall(), src_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (std::vector<T>& buf, int src_pid, int tag)
{
    return Arecv(buf.data(), buf.size(), src_pid, tag, Communicator());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv (T* buf, size_t n, int src_pid, int tag)
{
    return Recv(buf, n, src_pid, tag, Communicator());
}

namespace ParallelDescriptor {  // Have to use namespace here to work around a gcc bug
template <>
Message
Recv<char> (char* buf, size_t n, int src_pid, int tag, MPI_Comm comm);

template <class T>
Message
Recv (T* buf, size_t n, int src_pid, int tag, MPI_Comm comm)
{
    static_assert(!std::is_same<char,T>::value, "Recv: char version has been specialized");

    BL_PROFILE_T_S("ParallelDescriptor::Recv(Tsii)", T);
    BL_COMM_PROFILE(BLProfiler::RecvTsii, BLProfiler::BeforeCall(), src_pid, tag);

    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(buf,
                             n,
                             Mpi_typemap<T>::type(),
                             src_pid,
                             tag,
                             comm,
                             &stat) );
#ifdef BL_COMM_PROFILING
    int src_pid_comm(stat.MPI_SOURCE);
    int src_pid_world(stat.MPI_SOURCE);
    if(src_pid_comm != MPI_ANY_SOURCE) {
      MPI_Group groupComm, groupWorld;
      BL_MPI_REQUIRE( MPI_Comm_group(comm, &groupComm) );
      BL_MPI_REQUIRE( MPI_Comm_group(Communicator(), &groupWorld) );
      BL_MPI_REQUIRE( MPI_Group_translate_ranks(groupComm, 1, &src_pid_comm, groupWorld, &src_pid_world) );
    }

    BL_COMM_PROFILE(BLProfiler::RecvTsii, n * sizeof(T), src_pid_world, stat.MPI_TAG);
#endif
    return Message(stat, Mpi_typemap<T>::type());
}
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv (std::vector<T>& buf, int src_pid, int tag)
{
    return Recv(buf.data(), buf.size(), src_pid, tag, Communicator());
}

template <class T>
void
ParallelDescriptor::Bcast (T*     t,
                           size_t n,
                           int    root)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_ASSERT(n < static_cast<size_t>(std::numeric_limits<int>::max()));

    BL_PROFILE_T_S("ParallelDescriptor::Bcast(Tsi)", T);
    BL_COMM_PROFILE(BLProfiler::BCastTsi, BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    BL_MPI_REQUIRE( MPI_Bcast(t,
                              n,
                              Mpi_typemap<T>::type(),
                              root,
                              Communicator()) );
    BL_COMM_PROFILE(BLProfiler::BCastTsi, n * sizeof(T), root, BLProfiler::NoTag());
}

template <class T>
void
ParallelDescriptor::Bcast (T*     t,
                           size_t n,
                           int    root,
                           const MPI_Comm &comm)
{
#ifdef BL_LAZY
    int r;
    MPI_Comm_compare(comm, Communicator(), &r);
    if (r == MPI_IDENT)
        Lazy::EvalReduction();
#endif

    BL_ASSERT(n < static_cast<size_t>(std::numeric_limits<int>::max()));

    BL_PROFILE_T_S("ParallelDescriptor::Bcast(Tsi)", T);
    BL_COMM_PROFILE(BLProfiler::BCastTsi, BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    BL_MPI_REQUIRE( MPI_Bcast(t,
                              n,
                              Mpi_typemap<T>::type(),
                              root,
                              comm) );
    BL_COMM_PROFILE(BLProfiler::BCastTsi, n * sizeof(T), root, BLProfiler::NoTag());
}

template <class T, class T1>
void
ParallelDescriptor::Gather (const T* t,
                            size_t   n,
                            T1*      t1,
                            size_t   n1,
                            int      root)
{
    BL_PROFILE_T_S("ParallelDescriptor::Gather(TsT1si)", T);
    BL_COMM_PROFILE(BLProfiler::GatherTsT1Si, BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    BL_ASSERT(n  < static_cast<size_t>(std::numeric_limits<int>::max()));
    BL_ASSERT(n1 < static_cast<size_t>(std::numeric_limits<int>::max()));

    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(t),
                               n,
                               Mpi_typemap<T>::type(),
                               t1,
                               n1,
                               Mpi_typemap<T1>::type(),
                               root,
                               Communicator()) );
    BL_COMM_PROFILE(BLProfiler::GatherTsT1Si,  n * sizeof(T), root, BLProfiler::NoTag());
}

template <class T>
std::vector<T>
ParallelDescriptor::Gather (const T& t, int root)
{
    BL_PROFILE_T_S("ParallelDescriptor::Gather(Ti)", T);
    BL_COMM_PROFILE(BLProfiler::GatherTi, BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    std::vector<T> resl;
    if ( root == MyProc() ) resl.resize(NProcs());
    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(&t),
                               1,
                               Mpi_typemap<T>::type(),
                               resl.data(),
                               1,
                               Mpi_typemap<T>::type(),
                               root,
                               Communicator()) );
    BL_COMM_PROFILE(BLProfiler::GatherTi, sizeof(T), root, BLProfiler::NoTag());
    return resl;
}

template <class T>
void
ParallelDescriptor::Gatherv (const T* send, int sc,
                             T* recv, const std::vector<int>& rc, const std::vector<int>& disp,
                             int root)
{
    BL_PROFILE_T_S("ParallelDescriptor::Gatherv(Ti)", T);
    BL_COMM_PROFILE(BLProfiler::Gatherv,  BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    MPI_Gatherv(send, sc, ParallelDescriptor::Mpi_typemap<T>::type(),
                recv, &rc[0], &disp[0], ParallelDescriptor::Mpi_typemap<T>::type(),
                root, Communicator());

    BL_COMM_PROFILE(BLProfiler::Gatherv, std::accumulate(rc.begin(),rc.end(),0)*sizeof(T), root, BLProfiler::NoTag());
}

template <class T>
void
ParallelDescriptor::GatherLayoutDataToVector (const LayoutData<T>& sendbuf,
                                              Vector<T>& recvbuf, int root)
{
    BL_PROFILE_T_S("ParallelDescriptor::GatherLayoutData(Ti)", T);

    // Gather prelims
    Vector<T> T_to_send;
    T_to_send.reserve(sendbuf.local_size());

    for (int i : sendbuf.IndexArray())
    {
        T_to_send.push_back(sendbuf[i]);
    }

    int nprocs = ParallelContext::NProcsSub();
    Vector<int> recvcount(nprocs, 0);
    recvbuf.resize(sendbuf.size());
    const Vector<int>& old_pmap = sendbuf.DistributionMap().ProcessorMap();
    for (int i=0; i<old_pmap.size(); ++i)
    {
        ++recvcount[old_pmap[i]];
    }

    // Make a map from post-gather to pre-gather index
    Vector<Vector<int>> new_ind_to_old_ind(nprocs);
    for (int i=0; i<nprocs; ++i)
    {
        new_ind_to_old_ind[i].reserve(recvcount[i]);
    }
    for (int i=0; i<old_pmap.size(); ++i)
    {
        new_ind_to_old_ind[old_pmap[i]].push_back(i);
    }

    // Flatten
    Vector<int> new_index_to_old_index;
    new_index_to_old_index.reserve(old_pmap.size());
    for (const Vector<int>& v : new_ind_to_old_ind)
    {
        if (!v.empty())
        {
            for (int el : v)
            {
                new_index_to_old_index.push_back(el);
            }
        }
    }

    Vector<int> disp(nprocs);
    if (!disp.empty()) disp[0] = 0;
    std::partial_sum(recvcount.begin(), recvcount.end()-1, disp.begin()+1);
    Vector<T> new_index_to_T(sendbuf.size());

    MPI_Gatherv(T_to_send.data(), T_to_send.size(),
                ParallelDescriptor::Mpi_typemap<T>::type(),
                new_index_to_T.data(), recvcount.data(), disp.data(),
                ParallelDescriptor::Mpi_typemap<T>::type(),
                root, ParallelContext::CommunicatorSub());

    // Now work just on the root, which now has global information on the collected;
    // LayoutData information; with new_index_to_old_index and new_index_to_T,
    // sort the gathered vector in pre-gather index space
    if (ParallelContext::MyProcSub() == root)
    {
        // Invert the map (new_index) --> (old_index)
        Vector<int> old_index_to_new_index(sendbuf.size());

        for (int i=0; i<old_index_to_new_index.size(); ++i)
        {
            old_index_to_new_index[new_index_to_old_index[i]] = i;
        }

        for (int i=0; i<recvbuf.size(); ++i)
        {
            recvbuf[i] = new_index_to_T[old_index_to_new_index[i]];
        }
    }
}

template <class T, class T1>
void
ParallelDescriptor::Scatter (T*        t,
                             size_t    n,
                             const T1* t1,
                             size_t    n1,
                             int       root)
{
    BL_PROFILE_T_S("ParallelDescriptor::Scatter(TsT1si)", T);
    BL_COMM_PROFILE(BLProfiler::ScatterTsT1si,  BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    BL_MPI_REQUIRE( MPI_Scatter(const_cast<T1*>(t1),
                                n1,
                                Mpi_typemap<T1>::type(),
                                t,
                                n,
                                Mpi_typemap<T>::type(),
                                root,
                                Communicator()) );
    BL_COMM_PROFILE(BLProfiler::ScatterTsT1si, n * sizeof(T), root, BLProfiler::NoTag());
}

#else

namespace ParallelDescriptor
{
template <class T>
Message
Asend(const T* /*buf*/, size_t /*n*/, int /*dst_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
Message
Asend(const T* /*buf*/, size_t /*n*/, int /*dst_pid*/, int /*tag*/, MPI_Comm /*comm*/)
{
    return Message();
}

template <class T>
Message
Asend(const std::vector<T>& /*buf*/, int /*dst_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
Message
Send(const T* /*buf*/, size_t /*n*/, int /*dst_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
Message
Send(const T* /*buf*/, size_t /*n*/, int /*dst_pid*/, int /*tag*/, MPI_Comm /*comm*/)
{
    return Message();
}

template <class T>
Message
Send(const std::vector<T>& /*buf*/, int /*dst_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
Message
Arecv(T* /*buf*/, size_t /*n*/, int /*src_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
Message
Arecv(T* /*buf*/, size_t /*n*/, int /*src_pid*/, int /*tag*/, MPI_Comm /*comm*/)
{
    return Message();
}

template <class T>
Message
Arecv(std::vector<T>& /*buf*/, int /*src_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
Message
Recv(T* /*buf*/, size_t /*n*/, int /*src_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
Message
Recv(T* /*buf*/, size_t /*n*/, int /*src_pid*/, int /*tag*/, MPI_Comm /*comm*/)
{
    return Message();
}

template <class T>
Message
Recv(std::vector<T>& /*buf*/, int /*src_pid*/, int /*tag*/)
{
    return Message();
}

template <class T>
void
Bcast(T* /*t*/, size_t /*n*/, int /*root*/)
{}

template <class T>
void
Bcast(T* /*t*/, size_t /*n*/, int /*root*/, const MPI_Comm & /*comm*/)
{}

template <class T, class T1>
void
Gather (const T* t, size_t n, T1* t1, size_t n1, int /*root*/)
{
    BL_ASSERT(n == n1);
    amrex::ignore_unused(n1);

    int const sc = static_cast<int>(n);
    for (int j=0; j<sc; ++j) { t1[j] = t[j]; }
}

template <class T>
std::vector<T>
Gather(const T& t, int /*root*/)
{
    std::vector<T> resl(1);
    resl[0] = t;
    return resl;
}

template <class T>
void
Gatherv (const T* send, int sc,
         T* recv, const std::vector<int>& /*rc*/,
         const std::vector<int>& /*disp*/, int /*root*/)
{
    for (int j=0; j<sc; ++j) { recv[j] = send[j]; }
}

template <class T>
void
GatherLayoutDataToVector (const LayoutData<T>& sendbuf,
                          Vector<T>& recvbuf, int /*root*/)
{
    recvbuf.resize(sendbuf.size());

    for (int i=0; i<sendbuf.size(); ++i)
    {
        recvbuf[i] = sendbuf[i];
    }
}

template <class T, class T1>
void
Scatter(T* /*t*/, size_t /*n*/, const T1* /*t1*/, size_t /*n1*/, int /*root*/)
{}

}
#endif

namespace ParallelDescriptor {

#ifdef AMREX_USE_MPI

namespace detail {

template<typename T>
void DoAllReduce (T* r, MPI_Op op, int cnt)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_ASSERT(cnt > 0);

    BL_MPI_REQUIRE( MPI_Allreduce(MPI_IN_PLACE, r, cnt,
                                  Mpi_typemap<T>::type(), op,
                                  Communicator()) );
}

template<typename T>
void DoReduce (T* r, MPI_Op op, int cnt, int cpu)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_ASSERT(cnt > 0);

    if (MyProc() == cpu) {
        BL_MPI_REQUIRE( MPI_Reduce(MPI_IN_PLACE, r, cnt,
                                   Mpi_typemap<T>::type(), op,
                                   cpu, Communicator()) );
    } else {
        BL_MPI_REQUIRE( MPI_Reduce(r, r, cnt,
                                   Mpi_typemap<T>::type(), op,
                                   cpu, Communicator()) );
    }
}

}

    //! Real sum reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T& rvar) {
        detail::DoAllReduce<T>(&rvar,MPI_SUM,1);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T* rvar, int cnt) {
        detail::DoAllReduce<T>(rvar,MPI_SUM,cnt);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (Vector<std::reference_wrapper<T> >&& rvar)
    {
        int cnt = rvar.size();
        Vector<T> tmp{std::begin(rvar), std::end(rvar)};
        detail::DoAllReduce<T>(tmp.data(),MPI_SUM,cnt);
        for (int i = 0; i < cnt; ++i) {
            rvar[i].get() = tmp[i];
        }
    }

    //! Real sum reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T& rvar, int cpu) {
        detail::DoReduce<T>(&rvar,MPI_SUM,1,cpu);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T* rvar, int cnt, int cpu) {
        detail::DoReduce<T>(rvar,MPI_SUM,cnt,cpu);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (Vector<std::reference_wrapper<T> >&& rvar, int cpu)
    {
        int cnt = rvar.size();
        Vector<T> tmp{std::begin(rvar), std::end(rvar)};
        detail::DoReduce<T>(tmp.data(),MPI_SUM,cnt,cpu);
        for (int i = 0; i < cnt; ++i) {
            rvar[i].get() = tmp[i];
        }
    }

    //! Real max reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T& rvar) {
        detail::DoAllReduce<T>(&rvar,MPI_MAX,1);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T* rvar, int cnt) {
        detail::DoAllReduce<T>(rvar,MPI_MAX,cnt);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (Vector<std::reference_wrapper<T> >&& rvar)
    {
        int cnt = rvar.size();
        Vector<T> tmp{std::begin(rvar), std::end(rvar)};
        detail::DoAllReduce<T>(tmp.data(),MPI_MAX,cnt);
        for (int i = 0; i < cnt; ++i) {
            rvar[i].get() = tmp[i];
        }
    }

    //! Real max reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T& rvar, int cpu) {
        detail::DoReduce<T>(&rvar,MPI_MAX,1,cpu);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T* rvar, int cnt, int cpu) {
        detail::DoReduce<T>(rvar,MPI_MAX,cnt,cpu);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (Vector<std::reference_wrapper<T> >&& rvar, int cpu)
    {
        int cnt = rvar.size();
        Vector<T> tmp{std::begin(rvar), std::end(rvar)};
        detail::DoReduce<T>(tmp.data(),MPI_MAX,cnt,cpu);
        for (int i = 0; i < cnt; ++i) {
            rvar[i].get() = tmp[i];
        }
    }

    //! Real min reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T& rvar) {
        detail::DoAllReduce<T>(&rvar,MPI_MIN,1);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T* rvar, int cnt) {
        detail::DoAllReduce<T>(rvar,MPI_MIN,cnt);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (Vector<std::reference_wrapper<T> >&& rvar)
    {
        int cnt = rvar.size();
        Vector<T> tmp{std::begin(rvar), std::end(rvar)};
        detail::DoAllReduce<T>(tmp.data(),MPI_MIN,cnt);
        for (int i = 0; i < cnt; ++i) {
            rvar[i].get() = tmp[i];
        }
    }

    //! Real min reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T& rvar, int cpu) {
        detail::DoReduce<T>(&rvar,MPI_MIN,1,cpu);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T* rvar, int cnt, int cpu) {
        detail::DoReduce<T>(rvar,MPI_MIN,cnt,cpu);
    }

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (Vector<std::reference_wrapper<T> >&& rvar, int cpu)
    {
        int cnt = rvar.size();
        Vector<T> tmp{std::begin(rvar), std::end(rvar)};
        detail::DoReduce<T>(tmp.data(),MPI_MIN,cnt,cpu);
        for (int i = 0; i < cnt; ++i) {
            rvar[i].get() = tmp[i];
        }
    }

#else

    //! Real sum reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T& ) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T*, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (Vector<std::reference_wrapper<T> >&&) {}

    //! Real sum reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T&, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (T*, int, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealSum (Vector<std::reference_wrapper<T> >&&, int) {}

    //! Real max reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T&) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T*, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (Vector<std::reference_wrapper<T> >&&) {}

    //! Real max reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T&, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (T*, int, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMax (Vector<std::reference_wrapper<T> >&&, int) {}

    //! Real min reduction.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T&) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T*, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (Vector<std::reference_wrapper<T> >&&) {}

    //! Real min reduction to specified cpu.
    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T&, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (T*, int, int) {}

    template <typename T>
    typename std::enable_if<std::is_floating_point<T>::value>::type
    ReduceRealMin (Vector<std::reference_wrapper<T> >&&, int) {}

#endif
}

#ifdef AMREX_USE_MPI
namespace ParallelDescriptor {

template <class T>
struct Mpi_typemap<GpuComplex<T>>
{
    static MPI_Datatype type ()
    {
        static_assert(std::is_same<T,double>() ||
                      std::is_same<T,float >(),
                      "Unsupported type T for GpuComplex");
        if constexpr (std::is_same<T,double>()) {
            return MPI_C_DOUBLE_COMPLEX;
        } else {
            return MPI_C_FLOAT_COMPLEX;
        }
    }
};

template<typename TV, typename TI>
struct Mpi_typemap<ValLocPair<TV,TI>>
{
    static MPI_Datatype type ()
    {
        static MPI_Datatype mpi_type = MPI_DATATYPE_NULL;
        if (mpi_type == MPI_DATATYPE_NULL) {
            using T = ValLocPair<TV,TI>;
            static_assert(std::is_trivially_copyable<T>::value,
                          "To communicate with MPI, ValLocPair must be trivially copyable.");
            static_assert(std::is_standard_layout<T>::value,
                          "To communicate with MPI, ValLocPair must be standard layout");

            T vlp[2];
            MPI_Datatype types[] = {
                Mpi_typemap<TV>::type(),
                Mpi_typemap<TI>::type(),
            };
            int blocklens[] = { 1, 1 };
            MPI_Aint disp[2];
            BL_MPI_REQUIRE( MPI_Get_address(&vlp[0].value, &disp[0]) );
            BL_MPI_REQUIRE( MPI_Get_address(&vlp[0].index, &disp[1]) );
            disp[1] -= disp[0];
            disp[0] = 0;
            BL_MPI_REQUIRE( MPI_Type_create_struct(2, blocklens, disp, types,
                                               &mpi_type) );
            MPI_Aint lb, extent;
            BL_MPI_REQUIRE( MPI_Type_get_extent(mpi_type, &lb, &extent) );
            if (extent != sizeof(T)) {
                MPI_Datatype tmp = mpi_type;
                BL_MPI_REQUIRE( MPI_Type_create_resized(tmp, 0, sizeof(vlp[0]), &mpi_type) );
                BL_MPI_REQUIRE( MPI_Type_free(&tmp) );
            }
            BL_MPI_REQUIRE( MPI_Type_commit( &mpi_type ) );

            m_mpi_types.push_back(&mpi_type);
        }
        return mpi_type;
    }
};

template <typename T, typename F>
MPI_Op Mpi_op ()
{
    static MPI_Op mpi_op = MPI_OP_NULL;
    if (mpi_op == MPI_OP_NULL) {
        static auto user_fn = [] (void *invec, void *inoutvec, int* len, // NOLINT
                                  MPI_Datatype * /*datatype*/)
        {
            auto in = static_cast<T const*>(invec);
            auto out = static_cast<T*>(inoutvec);
            for (int i = 0; i < *len; ++i) {
                out[i] = F()(in[i],out[i]);
            }
        };
        BL_MPI_REQUIRE( MPI_Op_create(user_fn, 1, &mpi_op) );
        m_mpi_ops.push_back(&mpi_op);
    }
    return mpi_op;
}

}
#endif

}

#endif /*BL_PARALLELDESCRIPTOR_H*/
