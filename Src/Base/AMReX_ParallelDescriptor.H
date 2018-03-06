#ifndef BL_PARALLELDESCRIPTOR_H
#define BL_PARALLELDESCRIPTOR_H

#include <AMReX_BLBackTrace.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLassert.H>
#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#ifndef BL_AMRPROF
#include <AMReX_Box.H>
#endif
#include <AMReX_ccse-mpi.H>

#ifdef BL_USE_UPCXX
#include <upcxx.h>
#endif

#ifdef BL_USE_MPI3
#include <atomic>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef BL_LAZY
#include <AMReX_Lazy.H>
#endif

#ifdef AMREX_PMI
#include <pmi.h>
#endif

#include <iostream>
#include <vector>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>

namespace amrex {
//! Parallel frontend that abstracts functionalities needed to spawn processes and handle communication 
namespace ParallelDescriptor
{
#ifdef BL_USE_UPCXX
    struct UPCXX_MPI_Mode {
	unsigned short upcxx_cnt;
	bool mpi_ready;
	UPCXX_MPI_Mode () : upcxx_cnt(0), mpi_ready(true) { }
	void incr_upcxx () { ++upcxx_cnt; }
	void decr_upcxx () { --upcxx_cnt; }
	void set_upcxx_mode () { mpi_ready = false; }
	void set_mpi_mode (bool barrier_called=false) {
	    BL_ASSERT(upcxx_cnt == 0);
	    if (!mpi_ready) {
		if (!barrier_called) upcxx::barrier();
		mpi_ready = true;
	    }    
	}
    };
    extern UPCXX_MPI_Mode Mode;
#endif    

#ifdef BL_USE_MPI3
    extern MPI_Win cp_win;
    extern MPI_Win fb_win;
#endif

    //
    // Used as default argument to ParallelDescriptor::Barrier().
    //
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

	bool               m_finished;
	MPI_Datatype       m_type;
	MPI_Request        m_req;
	mutable MPI_Status m_stat;
    };
    /**
    * \brief Perform any needed parallel initialization.  This MUST be the
    * first routine in this class called from within a program.
    */
    void StartParallel (int*    argc = 0,
			char*** argv = 0,
                        MPI_Comm mpi_comm = MPI_COMM_WORLD);

    //! Split the process pool into teams
    void StartTeams ();
    void EndTeams ();

    //! Return true if MPI one sided is enabled
    bool MPIOneSided ();

    typedef int (*PTR_TO_SIGNAL_HANDLER)(int);
    void AddSignalHandler (PTR_TO_SIGNAL_HANDLER);
    /**
    * \brief Perform any needed parallel finalization.  This MUST be the
    * last routine in this class called from within a program.
    */
    void EndParallel ();

    extern const int myId_undefined;
    extern const int myId_notInGroup;
    extern int m_MyId_all;
    extern int m_MyId_comp;
    extern int m_MyId_sub;
    extern int m_MyId_sidecar;
    extern const int notInSidecar;

    extern const int nProcs_undefined;
    extern int m_nProcs_all;
    extern int m_nProcs_comp;
    extern int m_nProcs_sub;
    extern Vector<int> m_num_procs_clr;
    extern Vector<int> m_first_procs_clr;
    extern Vector<int> m_nProcs_sidecar;
    extern int nSidecars;
    extern int inWhichSidecar;

    //! return the rank number local to the comp or sidecar groups
    inline int
    MyProc ()
    {
        if (m_MyId_comp != myId_notInGroup)
        {
            BL_ASSERT(m_MyId_comp != myId_undefined);
            return m_MyId_comp;
        }
        else
        {
            BL_ASSERT(m_MyId_sidecar != myId_undefined);
            return m_MyId_sidecar;
        }
    }
    inline int
    MyProc (MPI_Comm comm)
    {
#ifdef BL_USE_MPI
	int r;
	MPI_Comm_rank(comm,&r);
	return r;
#else
	return 0;
#endif	
    }
    //! return the rank number for the global group
    inline int
    MyProcAll ()
    {
        BL_ASSERT(m_MyId_all != myId_undefined);
        return m_MyId_all;
    }
    inline int
    MyProcComp ()
    {
        BL_ASSERT(m_MyId_comp != myId_undefined);
        return m_MyId_comp;
    }
    inline int
    MyProcSidecar ()
    {
        BL_ASSERT(m_MyId_sidecar != myId_undefined);
        return m_MyId_sidecar;
    }
    inline bool
    InSidecarGroup ()
    {
        BL_ASSERT(m_MyId_sidecar != myId_undefined);
        BL_ASSERT( ( (inWhichSidecar == notInSidecar) &&
	             (m_MyId_sidecar == myId_notInGroup)
		   ) || (
		     (inWhichSidecar >= 0) &&
	             (m_MyId_sidecar != myId_notInGroup)
		   )
		);
        return (m_MyId_sidecar >= 0);
    }
    inline int
    InWhichSidecar ()
    {
        return inWhichSidecar;
    }
    inline bool
    InCompGroup ()
    {
        BL_ASSERT(m_MyId_comp != myId_undefined);
        return (m_MyId_comp >= 0);
    }

    /* Define variables `m_nProcs_sub' and `m_MyCommSubColor' */
    void
    init_clr_vars();

    /**
    * \brief This sets the number of MPI ranks in the sidecar and
    * initializes nProcs, MyId, groups, and communicators
    * for the comp and sidecar sets of ranks.
    * The user will need to handle any consequences of changing
    * the numbers of processors in each group (such as changing
    * distribution maps and moving data in FabArrays).
    *
    * The union of ranks in the comp and sidecar arrays must be
    * equal to the set of ranks in 'all' (the amrex global
    * communicator).  The sets must be disjoint.
    * Global rank zero is required to be in compRanksInAll[0],
    * this restriction may be removed in the future.
    */
    void
    SetNProcsSidecars (const Vector<int> &compRanksInAll,
                       const Vector<Vector<int> > &sidecarRanksInAll,
		       bool printRanks = false); 
    /**
    * \brief This version makes one sidecar (number zero in the sidecar array)
    * Sidecar ranks are set to the high end of the range [0, m_nProcs_all)
    */
    void
    SetNProcsSidecars (int nscp); 

    //! Provide functionalities needed to construct a team of processes to perform a particular job
    struct ProcessTeam
    {
#ifdef BL_USE_UPCXX
	typedef upcxx::team team_t;
#else
	typedef MPI_Comm team_t;
#endif

        //! synchronize processes within the team
	void Barrier () const { 
	    if (m_size > 1) {
#ifdef BL_USE_UPCXX
		m_upcxx_team.barrier();
#elif defined(BL_USE_MPI3)
		MPI_Barrier(m_team_comm);
#endif
	    }
	}

        //! memory fence
	void MemoryBarrier () const {
	    if (m_size > 1) {
#ifdef _OPENMP
		if (omp_in_parallel()) {
                    #pragma omp barrier
	        }
                #pragma omp master
#endif
		{
#ifdef BL_USE_UPCXX
                    m_upcxx_team.barrier();
#elif defined(BL_USE_MPI3)
		    std::atomic_thread_fence(std::memory_order_release);
		    MPI_Barrier(m_team_comm);
		    std::atomic_thread_fence(std::memory_order_acquire);
#endif
	        }
	    }
	}

        //! free a communicator
	void clear () {
#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
	    if (m_size > 1) {
		MPI_Comm_free(&m_team_comm);
		if (m_rankInTeam==0) MPI_Comm_free(&m_lead_comm);
	    }
#endif
	}
	
	const team_t& get() const {
#ifdef BL_USE_UPCXX
	    return m_upcxx_team; 
#else
	    return m_team_comm;
#endif
	}
        //! return the communicator
	const MPI_Comm& get_team_comm() const { return m_team_comm; }
	const MPI_Comm& get_lead_comm() const { return m_lead_comm; }

#ifdef BL_USE_UPCXX
	team_t m_upcxx_team;
#endif
	int    m_numTeams;
	int    m_size;
	int    m_color;
	int    m_lead;
	int    m_rankInTeam;
	int    m_do_team_reduce;

	MPI_Comm  m_team_comm;
	MPI_Comm  m_lead_comm;
    };

    extern ProcessTeam m_Team;

    extern int m_MinTag, m_MaxTag;
    inline int MinTag () { return m_MinTag; }
    inline int MaxTag () { return m_MaxTag; }

    //! return the number of MPI ranks local to the comp or sidecar groups
    inline int
    NProcs ()
    {
        if (m_MyId_comp != myId_notInGroup)
        {
            BL_ASSERT(m_nProcs_comp != nProcs_undefined);
            return m_nProcs_comp;
        }
        else if(inWhichSidecar != notInSidecar)
        {
            BL_ASSERT(m_nProcs_sidecar[inWhichSidecar] != nProcs_undefined);
            return m_nProcs_sidecar[inWhichSidecar];
        }
	else
	{
            return nProcs_undefined;
	}
    }
    inline int
    NProcsAll ()
    {
        BL_ASSERT(m_nProcs_all != nProcs_undefined);
        return m_nProcs_all;
    }
    inline int
    NProcsComp ()
    {
        BL_ASSERT(m_nProcs_comp != nProcs_undefined);
        return m_nProcs_comp;
    }
    inline int
    NProcsSidecar (int whichSidecar)
    {
	if(nSidecars == 0) {
	  return 0;
	}
        BL_ASSERT(whichSidecar >= 0 && whichSidecar < m_nProcs_sidecar.size());
        BL_ASSERT(m_nProcs_sidecar[whichSidecar] != nProcs_undefined);
        return m_nProcs_sidecar[whichSidecar];
    }
    /**
    * \brief The MPI rank number of the I/O Processor (probably rank 0).
    * This rank is usually used to write to stdout.
    */
    extern const int ioProcessor;
    inline int
    IOProcessorNumber ()
    {
        return ioProcessor;
    }
    inline int
    IOProcessorNumberAll ()
    {
        return ioProcessor;
    }
    /**
    * \brief Is this CPU the I/O Processor?
    * To get the rank number, call IOProcessorNumber()
    */
    inline bool
    IOProcessor ()
    {
         return MyProc() == IOProcessorNumber();
    }
    inline bool
    IOProcessorAll ()
    {
         return MyProcAll() == IOProcessorNumberAll();
    }
    //
    inline int
    TeamSize ()
    {
	return m_Team.m_size;
    }
    inline int
    NTeams ()
    {
	return m_Team.m_numTeams;
    }
    inline int
    MyTeamColor ()
    {
	return m_Team.m_color;
    }
    inline int
    MyTeamLead ()
    {
	return m_Team.m_lead;
    }
    inline int
    MyRankInTeam ()
    {
	return m_Team.m_rankInTeam;
    }
    inline int
    TeamLead (int rank)
    {
	return (rank >= 0) ? (rank - rank % m_Team.m_size) : MPI_PROC_NULL;
    }
    inline bool
    isTeamLead ()
    {
	return MyRankInTeam() == 0;
    }
    inline bool
    sameTeam (int rank)
    {
	return MyTeamLead() == TeamLead(rank);
    }
    inline bool
    sameTeam (int rankA, int rankB)
    {
	return TeamLead(rankA) == TeamLead(rankB);
    }
    inline int
    RankInLeadComm (int rank)
    {
	return (rank >= 0) ? (rank / m_Team.m_size) : MPI_PROC_NULL;
    }
    inline bool
    doTeamReduce () 
    { 
	return m_Team.m_do_team_reduce;
    }
    inline const ProcessTeam&
    MyTeam ()
    {
	return m_Team;
    }
    inline std::pair<int,int>
    team_range (int begin, int end, int rit = -1, int nworkers = 0)
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

#ifdef _OPENMP
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
    //
    // AMReX's Parallel Communicators
    //    m_comm_all       is for all ranks, probably duped from MPI_COMM_WORLD
    //    m_comm_comp      is for the ranks doing computations
    //    m_comm_sub       m_comm_comp is split into sub communicators
    //    m_comm_sidecar[] is for the ranks in the sidecar
    //    m_comm_inter[]   is for communicating between comp and sidecars
    //                     the user must create any sidecar to sidecar communicators
    //    m_comm_both[]    is for communicating within a group made from comp and
    //                     one of the sidecars
    //
    extern MPI_Comm m_comm_all;
    extern MPI_Comm m_comm_comp;
    extern MPI_Comm m_comm_sub;
    extern Vector<MPI_Comm> m_comm_sidecar;  // ---- [whichsidecar]
    extern Vector<MPI_Comm> m_comm_inter;    // ---- [whichsidecar]
    extern Vector<MPI_Comm> m_comm_both;     // ---- [whichsidecar]


    extern int m_nCommColors;
    extern int m_clr_map;

    //! A communicator can be split into multiple sub-communicators associated with different colors
    class Color
    {
    public:
	Color () : c(-100) {}
	explicit Color (int i) : c(i) {}
	int to_int () const { return c; }
	bool valid () const { return c >= 0 && c <= m_nCommColors; }
	friend bool operator== (const Color& lhs, const Color& rhs);
	friend bool operator!= (const Color& lhs, const Color& rhs);
	friend std::ostream& operator<< (std::ostream& os, const Color& color);
    private:
	int c;
    };
    inline bool operator== (const Color& lhs, const Color& rhs)
    {
	return lhs.c == rhs.c;
    }
    inline bool operator!= (const Color& lhs, const Color& rhs){
	return !(lhs == rhs);
    }
    inline std::ostream& operator<< (std::ostream& os, const Color& color)
    {
	os << color.c;
	if (os.fail())
	    amrex::Error("operator<<(ostream&,Const Color&) failed");
	return os;
    }

    extern Color m_MyCommSubColor;
    extern Color m_MyCommCompColor;

    void StartSubCommunicator ();
    void EndSubCommunicator ();

    //! Returns the "local" communicator
    inline MPI_Comm Communicator ()         
    {
        if (m_MyId_comp != myId_notInGroup)
        {
            return m_comm_comp;
        }
        else
        {
            BL_ASSERT(inWhichSidecar != notInSidecar && inWhichSidecar < m_comm_sidecar.size());
            return m_comm_sidecar[inWhichSidecar];
        }
    }
    inline MPI_Comm CommunicatorAll ()
    {
        return m_comm_all;
    }
    inline MPI_Comm CommunicatorComp ()
    {
        return m_comm_comp;
    }
    inline MPI_Comm CommunicatorSidecar ()
    {
	if(inWhichSidecar == notInSidecar) {
          return MPI_COMM_NULL;
	} else {
          return m_comm_sidecar[inWhichSidecar];
	}
    }
    inline MPI_Comm CommunicatorInter (int whichSidecar)
    {
        return m_comm_inter[whichSidecar];
    }

    inline MPI_Comm CommunicatorBoth (int whichSidecar)
    {
        return m_comm_both[whichSidecar];
    }

    inline int NColors () { return m_nCommColors; }
    inline Color DefaultColor () { return m_MyCommCompColor; }
    inline Color SubCommColor () { return m_MyCommSubColor; }
    inline bool isActive(Color color) 
    { 
	return color == DefaultColor() || color == SubCommColor();
    }
    inline int MyProc (Color color) 
    {
	if (color == DefaultColor()) {
	    return MyProc();
	} else if (color == SubCommColor()) {
	    return m_MyId_sub; 
	} else {
	    return MPI_PROC_NULL;
	}
    }
    inline int NProcs (Color color) 
    {
	if (color == DefaultColor()) {
	    return NProcs();
	} else if (color.valid()) {
	    return m_num_procs_clr[color.to_int()];
	} else {
	    return 0;
	}
    }
    inline MPI_Comm Communicator (Color color) 
    { 
	if (color == DefaultColor()) {
	    return Communicator();
	} else if (color == SubCommColor()) {
	    return m_comm_sub;
	} else {
	    return MPI_COMM_NULL;
	}
    }

    /* Given `rk_clrd', i.e. the rank in the colored `comm', return the
    rank in the global `CommComp' */
    int
    Translate(int rk_clrd,Color clr);

    inline bool
    IOProcessor (Color color)
    {
	if (color == DefaultColor()) {
	    return IOProcessor();
	} else if (color.valid()) {
	    return MyProc(color) == 0;
	} else {
	    return false;
	}
    }

    void Barrier (const std::string& message = Unnamed);
    void Barrier (const MPI_Comm &comm, const std::string& message = Unnamed);

    void Test (MPI_Request& request, int& flag, MPI_Status& status);

    void Comm_dup (MPI_Comm comm, MPI_Comm& newcomm);
    //! Abort with specified error code.
    void Abort (int errorcode = SIGABRT, bool backtrace = true);
    //! ErrorString return string associated with error internal error condition
    const char* ErrorString (int errcode);
    //! Returns wall-clock seconds since start of execution.
    double second ();

    //! And-wise boolean reduction.
    void ReduceBoolAnd (bool& rvar, Color color = DefaultColor());
    //! And-wise boolean reduction to specified cpu.
    void ReduceBoolAnd (bool& rvar, int cpu);

    //! Or-wise boolean reduction.
    void ReduceBoolOr  (bool& rvar, Color color = DefaultColor());
    //! Or-wise boolean reduction to specified cpu.
    void ReduceBoolOr  (bool& rvar, int cpu);

    //! Real sum reduction.
    void ReduceRealSum (Real& rvar, Color color = DefaultColor());
    void ReduceRealSum (Real* rvar, int cnt, Color color = DefaultColor());
    void ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar, Color color = DefaultColor());
    //! Real sum reduction to specified cpu.
    void ReduceRealSum (Real& rvar, int cpu);
    void ReduceRealSum (Real* rvar, int cnt, int cpu);
    void ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar, int cpu);

    //! Real max reduction.
    void ReduceRealMax (Real& rvar, Color color = DefaultColor());
    void ReduceRealMax (Real* rvar, int cnt, Color color = DefaultColor());
    void ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar, Color color = DefaultColor());
    //! Real max reduction to specified cpu.
    void ReduceRealMax (Real& rvar, int cpu);
    void ReduceRealMax (Real* rvar, int cnt, int cpu);
    void ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar, int cpu);

    //! Real min reduction.
    void ReduceRealMin (Real& rvar, Color color = DefaultColor());
    void ReduceRealMin (Real* rvar, int cnt, Color color = DefaultColor());
    void ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar, Color color = DefaultColor());
    //! Real min reduction to specified cpu.
    void ReduceRealMin (Real& rvar, int cpu);
    void ReduceRealMin (Real* rvar, int cnt, int cpu);
    void ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar, int cpu);

    //! Integer sum reduction.
    void ReduceIntSum (int& rvar, Color color = DefaultColor());
    void ReduceIntSum (int* rvar, int cnt, Color color = DefaultColor());
    void ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar, Color color = DefaultColor());
    //! Integer sum reduction to specified cpu.
    void ReduceIntSum (int& rvar, int cpu);
    void ReduceIntSum (int* rvar, int cnt, int cpu);
    void ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar, int cpu);

    //! Integer max reduction.
    void ReduceIntMax (int& rvar, Color color = DefaultColor());
    void ReduceIntMax (int* rvar, int cnt, Color color = DefaultColor());
    void ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar, Color color = DefaultColor());
    //! Integer max reduction to specified cpu.
    void ReduceIntMax (int& rvar, int cpu);
    void ReduceIntMax (int* rvar, int cnt, int cpu);
    void ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar, int cpu);

    //! Integer min reduction.
    void ReduceIntMin (int& rvar, Color color = DefaultColor());
    void ReduceIntMin (int* rvar, int cnt, Color color = DefaultColor());
    void ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar, Color color = DefaultColor());
    //! Integer min reduction to specified cpu.
    void ReduceIntMin (int& rvar, int cpu);
    void ReduceIntMin (int* rvar, int cnt, int cpu);
    void ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar, int cpu);

    //! Long sum reduction.
    void ReduceLongSum (long& rvar, Color color = DefaultColor());
    void ReduceLongSum (long* rvar, int cnt, Color color = DefaultColor());
    void ReduceLongSum (Vector<std::reference_wrapper<long> >&& rvar, Color color = DefaultColor());
    //! Long sum reduction to specified cpu.
    void ReduceLongSum (long& rvar, int cpu);
    void ReduceLongSum (long* rvar, int cnt, int cpu);
    void ReduceLongSum (Vector<std::reference_wrapper<long> >&& rvar, int cpu);

    //! Long max reduction.
    void ReduceLongMax (long& rvar, Color color = DefaultColor());
    void ReduceLongMax (long* rvar, int cnt, Color color = DefaultColor());
    void ReduceLongMax (Vector<std::reference_wrapper<long> >&& rvar, Color color = DefaultColor());
    //! Long max reduction to specified cpu.
    void ReduceLongMax (long& rvar, int cpu);
    void ReduceLongMax (long* rvar, int cnt, int cpu);
    void ReduceLongMax (Vector<std::reference_wrapper<long> >&& rvar, int cpu);

    //! Long min reduction.
    void ReduceLongMin (long& rvar, Color color = DefaultColor());
    void ReduceLongMin (long* rvar, int cnt, Color color = DefaultColor());
    void ReduceLongMin (Vector<std::reference_wrapper<long> >&& rvar, Color color = DefaultColor());
    //! Long min reduction to specified cpu.
    void ReduceLongMin (long& rvar, int cpu);
    void ReduceLongMin (long* rvar, int cnt, int cpu);
    void ReduceLongMin (Vector<std::reference_wrapper<long> >&& rvar, int cpu);

    //! Long and-wise reduction.
    void ReduceLongAnd (long& rvar, Color color = DefaultColor());
    void ReduceLongAnd (long* rvar, int cnt, Color color = DefaultColor());
    void ReduceLongAnd (Vector<std::reference_wrapper<long> >&& rvar, Color color = DefaultColor());
    //! Long and-wise reduction to specified cpu.
    void ReduceLongAnd (long& rvar, int cpu);
    void ReduceLongAnd (long* rvar, int cnt, int cpu);
    void ReduceLongAnd (Vector<std::reference_wrapper<long> >&& rvar, int cpu);

    // There are no color versions of reducion to specified cpu, because it could
    // be confusing what cpu means.  Is it in the global or colored communicator?

    //! Parallel gather.
    void Gather (Real* sendbuf,
                 int   sendcount,
                 Real* recvbuf,
                 int   root);
    /**
    * \brief Returns sequential message sequence numbers, usually used as
    * tags for send/recv.  getsetinc has these values:
    * 0:  return the current sequence number and increment
    * 1:  return the current sequence number
    * 2:  set the current sequence number to newvalue
    * 3:  jump the current sequence number and return the sequence number before the call.
    */
    int SeqNum (int getsetinc = 0, int newvalue = 0);
    int SubSeqNum (int getsetinc = 0, int newvalue = 0);

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
    template <class T> Message Recv(std::vector<T>& t, int pid, int tag);

    template <class T> void Bcast(T*, size_t n, int root = 0);
    template <class T> void Bcast(T*, size_t n, int root, const MPI_Comm &comm);
    void Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

    template <class T, class T1> void Scatter(T*, size_t n, const T1*, size_t n1, int root);

    template <class T, class T1> void Gather(const T*, size_t n, T1*, size_t n1, int root);
    template <class T> std::vector<T> Gather(const T&, int root);

    template <class T> void Gatherv (const T* send, int sc,
				     T* recv, const std::vector<int>& rc, const std::vector<int>& disp,
				     int root);
    template <class T> void Gatherv (const T* send, long sc,
				     T* recv, const std::vector<long>& rc, const std::vector<long>& disp,
				     int root);

    void Waitsome (Vector<MPI_Request>&, int&, Vector<int>&, Vector<MPI_Status>&);

    void MPI_Error(const char* file, int line, const char* msg, int rc);

    void ReadAndBcastFile(const std::string &filename, Vector<char> &charBuf,
                          bool bExitOnError = true,
			  const MPI_Comm &comm = Communicator() );
    void IProbe(int src_pid, int tag, int &mflag, MPI_Status &status);
    void IProbe(int src_pid, int tag, MPI_Comm comm, int &mflag, MPI_Status &status);

    // PMI = Process Management Interface, available on Crays. Provides API to
    // query topology of the job.
#ifdef AMREX_PMI
    void PMI_Initialize();
    void PMI_PrintMeshcoords(const pmi_mesh_coord_t *pmi_mesh_coord);
#endif
}
}

#define BL_MPI_REQUIRE(x)						\
do									\
{									\
  if ( int l_status_ = (x) )						\
    {									\
        amrex::ParallelDescriptor::MPI_Error(__FILE__,__LINE__,#x, l_status_); \
    }									\
}									\
while ( false )

namespace amrex {

#if BL_USE_MPI
template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend (const T* buf,
                           size_t   n,
                           int      dst_pid,
                           int      tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Asend(Tsii)", T);
    BL_COMM_PROFILE(BLProfiler::AsendTsii, n * sizeof(T), dst_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(buf),
                              n,
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(BLProfiler::AsendTsii, BLProfiler::AfterCall(), dst_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend (const T* buf,
                           size_t   n,
                           int      dst_pid,
                           int      tag,
                           MPI_Comm comm)
{
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

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Asend (const std::vector<T>& buf,
                           int                   dst_pid,
                           int                   tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Asend(vTii)", T);
    BL_COMM_PROFILE(BLProfiler::AsendvTii, buf.size() * sizeof(T), dst_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Isend(const_cast<T*>(&buf[0]),
                              buf.size(),
                              Mpi_typemap<T>::type(),
                              dst_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(BLProfiler::AsendvTii, BLProfiler::AfterCall(), dst_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send (const T* buf,
                          size_t   n,
                          int      dst_pid,
                          int      tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Send(Tsii)", T);
    BL_COMM_PROFILE(BLProfiler::SendTsii, n * sizeof(T), dst_pid, tag);

    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(buf),
                             n,
                             Mpi_typemap<T>::type(),
                             dst_pid,
                             tag,
                             Communicator()) );
    BL_COMM_PROFILE(BLProfiler::SendTsii, BLProfiler::AfterCall(), dst_pid, tag);
    return Message();
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send (const T* buf,
                          size_t   n,
                          int      dst_pid,
                          int      tag,
			  MPI_Comm comm)
{
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

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Send (const std::vector<T>& buf,
                          int                   dst_pid,
                          int                   tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Send(vTii)", T);
    BL_COMM_PROFILE(BLProfiler::SendvTii, buf.size() * sizeof(T), dst_pid, tag);

    BL_MPI_REQUIRE( MPI_Send(const_cast<T*>(&buf[0]),
                             buf.size(),
                             Mpi_typemap<T>::type(),
                             dst_pid,
                             tag,
                             Communicator()) );
    BL_COMM_PROFILE(BLProfiler::SendvTii, BLProfiler::AfterCall(), dst_pid, tag);
    return Message();
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (T*       buf,
                           size_t   n,
                           int      src_pid,
                           int      tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Arecv(Tsii)", T);
    BL_COMM_PROFILE(BLProfiler::ArecvTsii, n * sizeof(T), src_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(buf,
                              n,
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(BLProfiler::ArecvTsii, BLProfiler::AfterCall(), src_pid, tag);

    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (T*       buf,
                           size_t   n,
                           int      src_pid,
                           int      tag,
                           MPI_Comm comm)
{
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

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Arecv (std::vector<T>& buf,
                           int             src_pid,
                           int             tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Arecv(vTii)", T);
    BL_COMM_PROFILE(BLProfiler::ArecvvTii, buf.size() * sizeof(T), src_pid, tag);

    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Irecv(&buf[0],
                              buf.size(),
                              Mpi_typemap<T>::type(),
                              src_pid,
                              tag,
                              Communicator(),
                              &req) );
    BL_COMM_PROFILE(BLProfiler::ArecvvTii, BLProfiler::AfterCall(), src_pid, tag);
    return Message(req, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv (T*     buf,
                          size_t n,
                          int    src_pid,
                          int    tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Recv(Tsii)", T);
    BL_COMM_PROFILE(BLProfiler::RecvTsii, BLProfiler::BeforeCall(), src_pid, tag);

    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(buf,
                             n,
                             Mpi_typemap<T>::type(),
                             src_pid,
                             tag,
                             Communicator(),
                             &stat) );
    BL_COMM_PROFILE(BLProfiler::RecvTsii, n * sizeof(T), stat.MPI_SOURCE, stat.MPI_TAG);
    return Message(stat, Mpi_typemap<T>::type());
}

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv (T*     buf,
                          size_t n,
                          int    src_pid,
                          int    tag,
			  MPI_Comm comm)
{
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

template <class T>
ParallelDescriptor::Message
ParallelDescriptor::Recv (std::vector<T>& buf,
                          int             src_pid,
                          int             tag)
{
    BL_PROFILE_T_S("ParallelDescriptor::Recv(vTii)", T);
    BL_COMM_PROFILE(BLProfiler::RecvvTii, BLProfiler::BeforeCall(), src_pid, tag);

    MPI_Status stat;
    BL_MPI_REQUIRE( MPI_Recv(&buf[0],
                             buf.size(),
                             Mpi_typemap<T>::type(),
                             src_pid,
                             tag,
                             Communicator(),
                             &stat) );
    BL_COMM_PROFILE(BLProfiler::RecvvTii, buf.size() * sizeof(T), stat.MPI_SOURCE, stat.MPI_TAG);
    return Message(stat, Mpi_typemap<T>::type());
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

    BL_ASSERT(n < std::numeric_limits<int>::max());

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

    BL_ASSERT(n < std::numeric_limits<int>::max());

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

    BL_ASSERT(n  < std::numeric_limits<int>::max());
    BL_ASSERT(n1 < std::numeric_limits<int>::max());

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

    std::vector<T> resl(1);
    if ( root == MyProc() ) resl.resize(NProcs());
    BL_MPI_REQUIRE( MPI_Gather(const_cast<T*>(&t),
                               1,
                               Mpi_typemap<T>::type(),
                               &resl[0],
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
ParallelDescriptor::Gatherv (const T* send, long sc, 
			     T* recv, const std::vector<long>& rc, const std::vector<long>& disp,
			     int root)
{
    BL_PROFILE_T_S("ParallelDescriptor::Gatherv(Tlong)", T);
    BL_COMM_PROFILE(BLProfiler::Gatherv,  BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    int tag = ParallelDescriptor::SeqNum();

    if (MyProc() == root) {
	int nprocs = NProcs();
	BL_ASSERT(rc.size() == nprocs);
	BL_ASSERT(disp.size() == nprocs);

	std::vector<MPI_Request> request;
	for (int i = 0; i < nprocs; ++i) {
	    T* buf = recv + disp[i];
	    if (i == root) {
		for (int k = 0; k < rc[i]; ++k) {
		    *buf++ = send[k];
		}
	    } else {
		MPI_Request req;
		MPI_Irecv(buf, rc[i], ParallelDescriptor::Mpi_typemap<T>::type(), i,
			  tag, Communicator(), &req);
		request.push_back(req);
	    }
	}
	if (request.size() > 0) {
	    std::vector<MPI_Status> status(request.size());
	    MPI_Waitall(request.size(), &request[0], &status[0]);
	}
    } else {
	MPI_Send(const_cast<T*>(send), sc, ParallelDescriptor::Mpi_typemap<T>::type(), root, tag,
		 Communicator());
    }

    BL_COMM_PROFILE(BLProfiler::Gatherv, std::accumulate(rc.begin(),rc.end(),0L)*sizeof(T), root, BLProfiler::NoTag());
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
Asend(const T* buf, size_t n, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Asend(const T* buf, size_t n, int dst_pid, int tag, MPI_Comm comm)
{
    return Message();
}

template <class T>
Message
Asend(const std::vector<T>& buf, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Send(const T* buf, size_t n, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Send(const T* buf, size_t n, int dst_pid, int tag, MPI_Comm comm)
{
    return Message();
}

template <class T>
Message
Send(const std::vector<T>& buf, int dst_pid, int tag)
{
    return Message();
}

template <class T>
Message
Arecv(T* buf, size_t n, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Arecv(T* buf, size_t n, int src_pid, int tag, MPI_Comm comm)
{
    return Message();
}

template <class T>
Message
Arecv(std::vector<T>& buf, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Recv(T* buf, size_t n, int src_pid, int tag)
{
    return Message();
}

template <class T>
Message
Recv(T* buf, size_t n, int src_pid, int tag, MPI_Comm comm)
{
    return Message();
}

template <class T>
Message
Recv(std::vector<T>& buf, int src_pid, int tag)
{
    return Message();
}

template <class T>
void
Bcast(T* t, size_t n, int root)
{}

template <class T>
void
Bcast(T* t, size_t n, int root, const MPI_Comm &comm)
{}

template <class T, class T1>
void
Gather(const T* t, size_t n, T1* t1, size_t n1, int root)
{}

template <class T>
std::vector<T>
Gather(const T& t, int root)
{
    std::vector<T> resl(1);
    resl[0] = t;
    return resl;
}

template <class T, class T1>
void
Scatter(T* t, size_t n, const T1* t1, size_t n1, int root)
{}

}
#endif

}

#endif /*BL_PARALLELDESCRIPTOR_H*/
