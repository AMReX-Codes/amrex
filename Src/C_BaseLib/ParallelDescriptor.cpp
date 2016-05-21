
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <stack>
#include <list>

#include <Utility.H>
#include <BLProfiler.H>
#include <ParallelDescriptor.H>

#ifndef BL_AMRPROF
#include <ParmParse.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ParallelDescriptor
{
    //
    // My processor IDs.
    //
    const int myId_undefined  = -11;
    const int myId_notInGroup = -22;
    int m_MyId_all         = myId_undefined;
    int m_MyId_comp        = myId_undefined;
    int m_MyId_sub         = myId_undefined;
    int m_MyId_sidecar     = myId_undefined;
    //
    // The number of processors.
    //
    const int nProcs_undefined  = -33;
    int m_nProcs_all     = nProcs_undefined;
    int m_nProcs_comp    = nProcs_undefined;
    int m_nProcs_sub     = nProcs_undefined;
    int m_nProcs_sidecar = nProcs_undefined;
    int nSidecarProcs    = 0;
    //
    // Team
    //
    ProcessTeam m_Team;
    //
    // BoxLib's Communicators
    //
    MPI_Comm m_comm_all     = MPI_COMM_NULL;    // for all procs, probably MPI_COMM_WORLD
    MPI_Comm m_comm_comp    = MPI_COMM_NULL;    // for the computation procs
    MPI_Comm m_comm_sub     = MPI_COMM_NULL;    // for sub computation procs
    MPI_Comm m_comm_sidecar = MPI_COMM_NULL;    // for the in-situ performance monitor
    MPI_Comm m_comm_inter   = MPI_COMM_NULL;    // for communicating between comp and sidecar
    //
    // BoxLib's Groups
    //
    MPI_Group m_group_all     = MPI_GROUP_NULL;
    MPI_Group m_group_comp    = MPI_GROUP_NULL;
    MPI_Group m_group_sidecar = MPI_GROUP_NULL;

    int m_nCommColors = 1;
    Color m_MyCommSubColor;
    Color m_MyCommCompColor;

    int m_MinTag = 1000, m_MaxTag = -1;

    const int ioProcessor = 0;

#ifdef BL_USE_UPCXX
    UPCXX_MPI_Mode Mode;
#endif    

#ifdef BL_USE_MPI3
    MPI_Win cp_win;
    MPI_Win fb_win;
    MPI_Win fpb_win;
#endif
  
    namespace util
    {
	//
	// Reduce helper functons.
	//
	void DoAllReduceReal     (Real&      r, MPI_Op op, Color color = DefaultColor());
	void DoAllReduceLong     (long&      r, MPI_Op op, Color color = DefaultColor());
	void DoAllReduceInt      (int&       r, MPI_Op op, Color color = DefaultColor());

	void DoAllReduceReal     (Real*      r, MPI_Op op, int cnt, Color color = DefaultColor());
	void DoAllReduceLong     (long*      r, MPI_Op op, int cnt, Color color = DefaultColor());
	void DoAllReduceInt      (int*       r, MPI_Op op, int cnt, Color color = DefaultColor());

	void DoReduceReal     (Real&      r, MPI_Op op, int cpu);
	void DoReduceLong     (long&      r, MPI_Op op, int cpu);
	void DoReduceInt      (int&       r, MPI_Op op, int cpu);

	void DoReduceReal     (Real*      r, MPI_Op op, int cnt, int cpu);
	void DoReduceLong     (long*      r, MPI_Op op, int cnt, int cpu);
	void DoReduceInt      (int*       r, MPI_Op op, int cnt, int cpu);
    }

    typedef std::list<ParallelDescriptor::PTR_TO_SIGNAL_HANDLER> SH_LIST;
    SH_LIST The_Signal_Handler_List;
}

void
ParallelDescriptor::AddSignalHandler (PTR_TO_SIGNAL_HANDLER fp)
{
  The_Signal_Handler_List.push_back(fp);
}

#ifdef BL_USE_MPI

#include <ccse-mpi.H>

namespace
{
    const char*
    the_message_string (const char* file,
                        int         line,
                        const char* call,
                        int         status)
    {
	const int N = 512;
	static char buf[N];
	if ( status )
	{
	    snprintf(buf, N, "BoxLib MPI Error: File %s, line %d, %s: %s",
                     file, line, call, ParallelDescriptor::ErrorString(status));
	}
	else
	{
	    snprintf(buf, N, "BoxLib MPI Error: File %s, line %d, %s",
                     file, line, call);
	}
	return buf;
    }

}

namespace ParallelDescriptor
{
    void
    MPI_Error (const char* file, int line, const char* str, int rc)
    {
	BoxLib::Error(the_message_string(file, line, str, rc));
    }
}

void
ParallelDescriptor::Abort ()
{
#ifdef WIN32
    throw;
#endif
    MPI_Abort(CommunicatorAll(), -1);
}

void
ParallelDescriptor::Abort (int errorcode)
{
#ifdef BL_BGL
    MPI_Abort(CommunicatorAll(), errorcode);
#else
    BoxLib::Abort(ErrorString(errorcode));
#endif
}

const char*
ParallelDescriptor::ErrorString (int errorcode)
{
    BL_ASSERT(errorcode > 0 && errorcode <= MPI_ERR_LASTCODE);

    int len = 0;

    static char msg[MPI_MAX_ERROR_STRING+1];

    MPI_Error_string(errorcode, msg, &len);

    BL_ASSERT(len <= MPI_MAX_ERROR_STRING);

    return msg;
}

void
ParallelDescriptor::Message::wait ()
{
    BL_PROFILE_S("ParallelDescriptor::Message::wait()");

    BL_COMM_PROFILE(BLProfiler::Wait, sizeof(m_type), pid(), tag());
    BL_MPI_REQUIRE( MPI_Wait(&m_req, &m_stat) );
    BL_COMM_PROFILE(BLProfiler::Wait, sizeof(m_type), BLProfiler::AfterCall(), tag());
}

bool
ParallelDescriptor::Message::test ()
{
    int flag;
    BL_COMM_PROFILE(BLProfiler::Test, sizeof(m_type), pid(), tag());
    BL_MPI_REQUIRE( MPI_Test(&m_req, &flag, &m_stat) );
    BL_COMM_PROFILE(BLProfiler::Test, flag, BLProfiler::AfterCall(), tag());
    m_finished = flag != 0;
    return m_finished;
}

int
ParallelDescriptor::Message::tag () const
{
    if ( !m_finished ) BoxLib::Error("Message::tag: Not Finished!");
    return m_stat.MPI_TAG;
}

int
ParallelDescriptor::Message::pid () const
{
    if ( !m_finished ) BoxLib::Error("Message::pid: Not Finished!");
    return m_stat.MPI_SOURCE;
}

size_t
ParallelDescriptor::Message::count () const
{
    if ( m_type == MPI_DATATYPE_NULL ) BoxLib::Error("Message::count: Bad Type!");
    if ( !m_finished ) BoxLib::Error("Message::count: Not Finished!");
    int cnt;
    BL_MPI_REQUIRE( MPI_Get_count(&m_stat, m_type, &cnt) );
    return cnt;
}

void
ParallelDescriptor::StartParallel (int*    argc,
                                   char*** argv,
                                   MPI_Comm mpi_comm)
{
    m_comm_all = mpi_comm;

    int sflag(0);

    BL_MPI_REQUIRE( MPI_Initialized(&sflag) );

    if ( ! sflag) {
	BL_MPI_REQUIRE( MPI_Init(argc, argv) );
    }
    
    BL_MPI_REQUIRE( MPI_Comm_size(CommunicatorAll(), &m_nProcs_all) );
    BL_MPI_REQUIRE( MPI_Comm_rank(CommunicatorAll(), &m_MyId_all) );

#ifdef BL_USE_MPI3
    int mpi_version, mpi_subversion;
    BL_MPI_REQUIRE( MPI_Get_version(&mpi_version, &mpi_subversion) );
    if (mpi_version < 3) BoxLib::Abort("MPI 3 is needed because USE_MPI3=TRUE");
#endif

    if(m_MyId_all == 0 && nSidecarProcs > 0) {
      std::cout << "**** nSidecarProcs = " << nSidecarProcs << std::endl;
    }
    if(nSidecarProcs >= m_nProcs_all) {
      std::cerr << "**** nSidecarProcs >= m_nProcs_all:  " << nSidecarProcs
                << " >= " << m_nProcs_all << std::endl;
      BoxLib::Abort("Error in StartParallel:  bad nSidecarProcs.");
    }

    MPI_Comm_group(m_comm_all, &m_group_all);

    if(nSidecarProcs > 0) {
      Array<int> sidecarRanksInAll(nSidecarProcs, nProcs_undefined);
      for(int ip(0); ip < nSidecarProcs; ++ip) {
        sidecarRanksInAll[ip] = m_nProcs_all - nSidecarProcs + ip;
      }
      MPI_Group_excl(m_group_all, nSidecarProcs, sidecarRanksInAll.dataPtr(), &m_group_comp);
      MPI_Comm_create(m_comm_all, m_group_comp, &m_comm_comp);

      MPI_Group_incl(m_group_all, nSidecarProcs, sidecarRanksInAll.dataPtr(), &m_group_sidecar);
      MPI_Comm_create(m_comm_all, m_group_sidecar, &m_comm_sidecar);

      MPI_Group_size(m_group_sidecar, &m_nProcs_sidecar);
      MPI_Group_size(m_group_comp, &m_nProcs_comp);
    } else {
      m_comm_comp  = m_comm_all;
      m_group_comp = m_group_all;

      m_nProcs_comp = m_nProcs_all;
      m_nProcs_sidecar = 0;
    }

    // ---- find the maximum value for a tag
    int flag(0), *attrVal;
    BL_MPI_REQUIRE( MPI_Attr_get(m_comm_all, MPI_TAG_UB, &attrVal, &flag) );
    if(flag) {
      m_MaxTag = *attrVal;
      m_MaxTag -= 4;  // so we dont wrap if maxint
      m_MaxTag = std::max(m_MaxTag, 9000);
    } else {
      m_MaxTag = 9000;
    }
    BL_COMM_PROFILE_TAGRANGE(m_MinTag, m_MaxTag);

    if(nSidecarProcs > 0) {
      int tag(m_MaxTag + 1);


      if(m_MyId_all >= m_nProcs_all - nSidecarProcs) {  // ---- sidecar group
        MPI_Group_rank(m_group_sidecar, &m_MyId_sidecar);
        MPI_Intercomm_create(m_comm_sidecar, 0, m_comm_all, 0, tag, &m_comm_inter);
	m_MyId_comp = myId_notInGroup;

      } else {                        // ---- in computation group
        MPI_Group_rank(m_group_comp, &m_MyId_comp);
        MPI_Intercomm_create(m_comm_comp, 0, m_comm_all, m_nProcs_all - nSidecarProcs,
	                     tag, &m_comm_inter);
	m_MyId_sidecar = myId_notInGroup;
      }

    } else {
      m_MyId_comp   = m_MyId_all;
      m_MyId_sidecar   = myId_notInGroup;
    }

    //
    // Wait until all other processes are properly started.
    //
    BL_MPI_REQUIRE( MPI_Barrier(CommunicatorAll()) );

    if(m_MyId_all     == myId_undefined ||
       m_MyId_comp    == myId_undefined ||
       m_MyId_sidecar == myId_undefined)
    {
      std::cerr << "m_MyId_all m_MyId_comp m_MyId_sidecar = " << m_MyId_all << "  "
	          << m_MyId_comp << "  " << m_MyId_sidecar << std::endl;
      BoxLib::Abort("**** Error:  bad MyId in ParallelDescriptor::StartParallel()");
    }
    if(m_nProcs_all     == nProcs_undefined ||
       m_nProcs_comp    == nProcs_undefined ||
       m_nProcs_sidecar == nProcs_undefined ||
       (m_nProcs_comp + m_nProcs_sidecar != m_nProcs_all))
    {
      std::cerr << "m_nProcs_all m_nProcs_comp m_nProcs_sidecar = " << m_nProcs_all << "  "
	          << m_nProcs_comp << "  " << m_nProcs_sidecar << std::endl;
      BoxLib::Abort("**** Error:  bad nProcs in ParallelDescriptor::StartParallel()");
    }
}

void
ParallelDescriptor::EndParallel ()
{
    BL_ASSERT(m_MyId_all != -1);
    BL_ASSERT(m_nProcs_all != -1);

    MPI_Group_free(&m_group_all);

    BL_MPI_REQUIRE( MPI_Finalize() );
}

void
ParallelDescriptor::StartSubCommunicator ()
{
    ParmParse pp("boxlib");
    pp.query("ncolors", m_nCommColors);

    if (m_nCommColors > 1) {

#if defined(BL_USE_MPI3) || defined(BL_USE_UPCXX)
	//    m_nCommColors = 1;
	BoxLib::Abort("boxlib.ncolors > 1 not supported for MPI3 and UPCXX");
	if (doTeamReduce())
	    BoxLib::Abort("boxlib.ncolors > 1 not supported with team.reduce on");	    
#endif

#ifdef BL_LAZY
	BoxLib::Abort("boxlib.ncolors > 1 not supported for LAZY=TRUE");
#endif

	m_nProcs_sub = m_nProcs_comp / m_nCommColors;
	if (m_nProcs_sub * m_nCommColors != m_nProcs_comp) {
	    BoxLib::Abort("# of processors is not divisible by boxlib.ncolors");
	}
	m_MyCommSubColor  = Color(MyProc()/m_nProcs_sub);
	m_MyCommCompColor = Color(m_nCommColors);  // special color for CommComp color

	BL_MPI_REQUIRE( MPI_Comm_split(Communicator(), m_MyCommSubColor.to_int(), MyProc(), &m_comm_sub) );
	BL_MPI_REQUIRE( MPI_Comm_rank(m_comm_sub, &m_MyId_sub) );
    } else {
	m_nCommColors = 1;
	m_nProcs_sub  = NProcs();
	m_MyCommSubColor = Color(0);
	m_MyCommCompColor = Color(0);
	m_comm_sub    = Communicator();
	m_MyId_sub    = MyProc();
    }
}

void
ParallelDescriptor::EndSubCommunicator ()
{
    if (m_nCommColors > 1) {
	MPI_Comm_free(&m_comm_sub);
    }
}

double
ParallelDescriptor::second ()
{
    return MPI_Wtime();
}

void
ParallelDescriptor::Barrier (const std::string &message)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::Barrier()");
    BL_COMM_PROFILE_BARRIER(message, true);

    BL_MPI_REQUIRE( MPI_Barrier(ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE_BARRIER(message, false);
}

void
ParallelDescriptor::Barrier (MPI_Comm comm, const std::string &message)
{
#ifdef BL_LAZY
    int r;
    MPI_Comm_compare(comm, Communicator(), &r);
    if (r == MPI_IDENT)
	Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::Barrier(comm)");
    BL_COMM_PROFILE_BARRIER(message, true);

    BL_MPI_REQUIRE( MPI_Barrier(comm) );

    BL_COMM_PROFILE_BARRIER(message, false);
}

void
ParallelDescriptor::Test (MPI_Request& request, int& flag, MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Test()");
    BL_COMM_PROFILE(BLProfiler::Test, sizeof(char), status.MPI_SOURCE, status.MPI_TAG);

    BL_MPI_REQUIRE( MPI_Test(&request,&flag,&status) );

    BL_COMM_PROFILE(BLProfiler::Test, flag, BLProfiler::AfterCall(), status.MPI_TAG);
}

void
ParallelDescriptor::IProbe (int src_pid, int tag, int& flag, MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Iprobe()");
    BL_COMM_PROFILE(BLProfiler::Iprobe, sizeof(char), src_pid, tag);

    BL_MPI_REQUIRE( MPI_Iprobe(src_pid, tag, ParallelDescriptor::Communicator(),
                               &flag, &status) );

    BL_COMM_PROFILE(BLProfiler::Iprobe, flag, BLProfiler::AfterCall(), status.MPI_TAG);
}

void
ParallelDescriptor::IProbe (int src_pid, int tag, MPI_Comm comm, int& flag, MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Iprobe(comm)");
    BL_COMM_PROFILE(BLProfiler::Iprobe, sizeof(char), src_pid, tag);

    BL_MPI_REQUIRE( MPI_Iprobe(src_pid, tag, comm,
                               &flag, &status) );

    BL_COMM_PROFILE(BLProfiler::Iprobe, flag, BLProfiler::AfterCall(), status.MPI_TAG);
}

void
ParallelDescriptor::Comm_dup (MPI_Comm comm, MPI_Comm& newcomm)
{
    BL_PROFILE_S("ParallelDescriptor::Comm_dup()");
    BL_MPI_REQUIRE( MPI_Comm_dup(comm, &newcomm) );
}

void
ParallelDescriptor::util::DoAllReduceReal (Real&  r,
                                           MPI_Op op,
					   Color  color)
{
    if (!isActive(color)) return;

#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceReal()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, BLProfiler::BeforeCall(), true);

    Real recv;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Real recv_team;
	BL_MPI_REQUIRE( MPI_Reduce(&r, &recv_team, 1, Mpi_typemap<Real>::type(), op,
				   0, MyTeam().get_team_comm()) );
	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Allreduce(&recv_team, &recv, 1, Mpi_typemap<Real>::type(), op,
					  MyTeam().get_lead_comm()) );
	}
	BL_MPI_REQUIRE( MPI_Bcast(&recv, 1, Mpi_typemap<Real>::type(), 
				  0, MyTeam().get_team_comm()) );
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Allreduce(&r,
				      &recv,
				      1,
				      Mpi_typemap<Real>::type(),
				      op,
				      Communicator(color)) );
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, sizeof(Real), false);
    r = recv;
}

void
ParallelDescriptor::util::DoAllReduceReal (Real*  r,
                                           MPI_Op op,
                                           int    cnt,
					   Color  color)
{
    if (!isActive(color)) return;

#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceReal()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, BLProfiler::BeforeCall(), true);

    BL_ASSERT(cnt > 0);

    Array<Real> recv(cnt);

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Array<Real> recv_team(cnt);
	BL_MPI_REQUIRE( MPI_Reduce(r, recv_team.dataPtr(), cnt, Mpi_typemap<Real>::type(), op,
				   0, MyTeam().get_team_comm()) );
	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Allreduce(recv_team.dataPtr(), recv.dataPtr(), cnt, 
					  Mpi_typemap<Real>::type(), op,
					  MyTeam().get_lead_comm()) );
	}
	BL_MPI_REQUIRE( MPI_Bcast(recv.dataPtr(), cnt, Mpi_typemap<Real>::type(), 
				  0, MyTeam().get_team_comm()) );
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Allreduce(r,
				      recv.dataPtr(),
				      cnt,
				      Mpi_typemap<Real>::type(),
				      op,
				      Communicator(color)) );
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, cnt * sizeof(Real), false);
    for (int i = 0; i < cnt; i++)
        r[i] = recv[i];
}

void
ParallelDescriptor::util::DoReduceReal (Real&  r,
                                        MPI_Op op,
                                        int    cpu)
{
#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceReal()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceR, sizeof(Real), cpu);

    Real recv;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Real recv_team;
	BL_MPI_REQUIRE( MPI_Reduce(&r, &recv_team, 1, Mpi_typemap<Real>::type(), op,
				   0, MyTeam().get_team_comm()) );

	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Reduce(&recv_team, &recv, 1, Mpi_typemap<Real>::type(), op,
				       RankInLeadComm(cpu), MyTeam().get_lead_comm()) );
	}
	if (sameTeam(cpu)) {
	    BL_MPI_REQUIRE( MPI_Bcast(&recv, 1, Mpi_typemap<Real>::type(), 
				      0, MyTeam().get_team_comm()) );
	}
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Reduce(&r,
				   &recv,
				   1,
				   Mpi_typemap<Real>::type(),
				   op,
				   cpu,
				   Communicator()) );
    }
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceR, BLProfiler::AfterCall(), cpu);

    if (ParallelDescriptor::MyProc() == cpu)
	r = recv;
}

void
ParallelDescriptor::util::DoReduceReal (Real*  r,
                                        MPI_Op op,
                                        int    cnt,
                                        int    cpu)
{
#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceReal()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceR, cnt * sizeof(Real), cpu);

    BL_ASSERT(cnt > 0);

    Array<Real> recv(cnt);

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Array<Real> recv_team(cnt);
	BL_MPI_REQUIRE( MPI_Reduce(r, &recv_team[0], cnt, Mpi_typemap<Real>::type(), op,
				   0, MyTeam().get_team_comm()) );

	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Reduce(&recv_team[0], &recv[0], cnt, Mpi_typemap<Real>::type(), op,
				       RankInLeadComm(cpu), MyTeam().get_lead_comm()) );
	}
	if (sameTeam(cpu)) {
	    BL_MPI_REQUIRE( MPI_Bcast(&recv[0], cnt, Mpi_typemap<Real>::type(), 
				      0, MyTeam().get_team_comm()) );
	}
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Reduce(r,
				   recv.dataPtr(),
				   cnt,
				   Mpi_typemap<Real>::type(),
				   op,
				   cpu,
				   Communicator()) );
    }
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceR, BLProfiler::AfterCall(), cpu);

    if (ParallelDescriptor::MyProc() == cpu)
    {
        for (int i = 0; i < cnt; i++)
            r[i] = recv[i];
    }
}

void
ParallelDescriptor::ReduceRealMax (Real& r, Color color)
{
    BL_PROFILE("ReduceRealMax");
    util::DoAllReduceReal(r,MPI_MAX,color);
}

void
ParallelDescriptor::ReduceRealMin (Real& r, Color color)
{
    util::DoAllReduceReal(r,MPI_MIN,color);
}

void
ParallelDescriptor::ReduceRealSum (Real& r, Color color)
{
    util::DoAllReduceReal(r,MPI_SUM,color);
}

void
ParallelDescriptor::ReduceRealMax (Real* r, int cnt, Color color)
{
    BL_PROFILE("ReduceRealMax");
    util::DoAllReduceReal(r,MPI_MAX,cnt,color);
}

void
ParallelDescriptor::ReduceRealMin (Real* r, int cnt, Color color)
{
    util::DoAllReduceReal(r,MPI_MIN,cnt,color);
}

void
ParallelDescriptor::ReduceRealSum (Real* r, int cnt, Color color)
{
    util::DoAllReduceReal(r,MPI_SUM,cnt,color);
}

void
ParallelDescriptor::ReduceRealMax (Real& r, int cpu)
{
    BL_PROFILE("ReduceRealMax");
    util::DoReduceReal(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceRealMin (Real& r, int cpu)
{
    util::DoReduceReal(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceRealSum (Real& r, int cpu)
{
    util::DoReduceReal(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceRealMax (Real* r, int cnt, int cpu)
{
    BL_PROFILE("ReduceRealMax");
    util::DoReduceReal(r,MPI_MAX,cnt,cpu);
}


void
ParallelDescriptor::ReduceRealMin (Real* r, int cnt, int cpu)
{
    util::DoReduceReal(r,MPI_MIN,cnt,cpu);
}

void
ParallelDescriptor::ReduceRealSum (Real* r, int cnt, int cpu)
{
    util::DoReduceReal(r,MPI_SUM,cnt,cpu);
}

void
ParallelDescriptor::util::DoAllReduceLong (long&  r,
                                           MPI_Op op,
					   Color  color)
{
    if (!isActive(color)) return;

#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceLong()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceL, BLProfiler::BeforeCall(), true);

    long recv;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	long recv_team;
	BL_MPI_REQUIRE( MPI_Reduce(&r, &recv_team, 1, MPI_LONG, op,
				   0, MyTeam().get_team_comm()) );
	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Allreduce(&recv_team, &recv, 1, MPI_LONG, op,
					  MyTeam().get_lead_comm()) );
	}
	BL_MPI_REQUIRE( MPI_Bcast(&recv, 1, MPI_LONG, 
				  0, MyTeam().get_team_comm()) );
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Allreduce(&r,
				      &recv,
				      1,
				      MPI_LONG,
				      op,
				      Communicator(color)) );
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceL, sizeof(long), false);
    r = recv;
}

void
ParallelDescriptor::util::DoAllReduceLong (long*  r,
                                           MPI_Op op,
                                           int    cnt,
					   Color  color)
{
    if (!isActive(color)) return;

#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceLong()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceL, BLProfiler::BeforeCall(), true);

    BL_ASSERT(cnt > 0);

    Array<long> recv(cnt);

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Array<long> recv_team(cnt);
	BL_MPI_REQUIRE( MPI_Reduce(r, recv_team.dataPtr(), cnt, MPI_LONG, op,
				   0, MyTeam().get_team_comm()) );
	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Allreduce(recv_team.dataPtr(), recv.dataPtr(), cnt, 
					  MPI_LONG, op,
					  MyTeam().get_lead_comm()) );
	}
	BL_MPI_REQUIRE( MPI_Bcast(recv.dataPtr(), cnt, MPI_LONG,
				  0, MyTeam().get_team_comm()) );
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Allreduce(r,
				      recv.dataPtr(),
				      cnt,
				      MPI_LONG,
				      op,
				      Communicator(color)) );
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceL, cnt * sizeof(long), false);
    for (int i = 0; i < cnt; i++)
        r[i] = recv[i];
}

void
ParallelDescriptor::util::DoReduceLong (long&  r,
                                        MPI_Op op,
                                        int    cpu)
{
#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceLong()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceL, sizeof(long), cpu);

    long recv;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	long recv_team;
	BL_MPI_REQUIRE( MPI_Reduce(&r, &recv_team, 1, MPI_LONG, op,
				   0, MyTeam().get_team_comm()) );

	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Reduce(&recv_team, &recv, 1, MPI_LONG, op,
				       RankInLeadComm(cpu), MyTeam().get_lead_comm()) );
	}
	if (sameTeam(cpu)) {
	    BL_MPI_REQUIRE( MPI_Bcast(&recv, 1, MPI_LONG,
				      0, MyTeam().get_team_comm()) );
	}
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Reduce(&r,
				   &recv,
				   1,
				   MPI_LONG,
				   op,
				   cpu,
				   Communicator()));
    }
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceL, BLProfiler::AfterCall(), cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::util::DoReduceLong (long*  r,
                                        MPI_Op op,
                                        int    cnt,
                                        int    cpu)
{
#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceLong()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceL, cnt * sizeof(long), cpu);

    BL_ASSERT(cnt > 0);

    Array<long> recv(cnt);

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Array<long> recv_team(cnt);
	BL_MPI_REQUIRE( MPI_Reduce(r, &recv_team[0], cnt, MPI_LONG, op,
				   0, MyTeam().get_team_comm()) );

	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Reduce(&recv_team[0], &recv[0], cnt, MPI_LONG, op,
				       RankInLeadComm(cpu), MyTeam().get_lead_comm()) );
	}
	if (sameTeam(cpu)) {
	    BL_MPI_REQUIRE( MPI_Bcast(&recv[0], cnt, MPI_LONG, 
				      0, MyTeam().get_team_comm()) );
	}
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Reduce(r,
				   recv.dataPtr(),
				   cnt,
				   MPI_LONG,
				   op,
				   cpu,
				   Communicator()));
    }
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceL, BLProfiler::AfterCall(), cpu);

    if (ParallelDescriptor::MyProc() == cpu)
    {
        for (int i = 0; i < cnt; i++)
            r[i] = recv[i];
    }
}

void
ParallelDescriptor::ReduceLongAnd (long& r, Color color)
{
    util::DoAllReduceLong(r,MPI_LAND,color);
}

void
ParallelDescriptor::ReduceLongSum (long& r, Color color)
{
    util::DoAllReduceLong(r,MPI_SUM,color);
}

void
ParallelDescriptor::ReduceLongMax (long& r, Color color)
{
    util::DoAllReduceLong(r,MPI_MAX,color);
}

void
ParallelDescriptor::ReduceLongMin (long& r, Color color)
{
    util::DoAllReduceLong(r,MPI_MIN,color);
}

void
ParallelDescriptor::ReduceLongAnd (long* r, int cnt, Color color)
{
    util::DoAllReduceLong(r,MPI_LAND,cnt,color);
}

void
ParallelDescriptor::ReduceLongSum (long* r, int cnt, Color color)
{
    util::DoAllReduceLong(r,MPI_SUM,cnt,color);
}

void
ParallelDescriptor::ReduceLongMax (long* r, int cnt, Color color)
{
    util::DoAllReduceLong(r,MPI_MAX,cnt,color);
}

void
ParallelDescriptor::ReduceLongMin (long* r, int cnt, Color color)
{
    util::DoAllReduceLong(r,MPI_MIN,cnt,color);
}

void
ParallelDescriptor::ReduceLongAnd (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_LAND,cpu);
}

void
ParallelDescriptor::ReduceLongSum (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceLongMax (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceLongMin (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceLongAnd (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_LAND,cnt,cpu);
}

void
ParallelDescriptor::ReduceLongSum (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_SUM,cnt,cpu);
}

void
ParallelDescriptor::ReduceLongMax (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_MAX,cnt,cpu);
}

void
ParallelDescriptor::ReduceLongMin (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_MIN,cnt,cpu);
}

void
ParallelDescriptor::util::DoAllReduceInt (int&   r,
                                          MPI_Op op,
					  Color  color)
{
    if (!isActive(color)) return;

#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceInt()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceI, BLProfiler::BeforeCall(), true);

    int recv;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	int recv_team;
	BL_MPI_REQUIRE( MPI_Reduce(&r, &recv_team, 1, MPI_INT, op,
				   0, MyTeam().get_team_comm()) );
	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Allreduce(&recv_team, &recv, 1, MPI_INT, op,
					  MyTeam().get_lead_comm()) );
	}
	BL_MPI_REQUIRE( MPI_Bcast(&recv, 1, MPI_INT, 
				  0, MyTeam().get_team_comm()) );
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Allreduce(&r,
				      &recv,
				      1,
				      MPI_INT,
				      op,
				      Communicator(color)));
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceI, sizeof(int), false);
    r = recv;
}

void
ParallelDescriptor::util::DoAllReduceInt (int*   r,
                                          MPI_Op op,
                                          int    cnt,
					  Color  color)
{
    if (!isActive(color)) return;

#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceInt()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceI, BLProfiler::BeforeCall(), true);

    BL_ASSERT(cnt > 0);

    Array<int> recv(cnt);

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Array<int> recv_team(cnt);
	BL_MPI_REQUIRE( MPI_Reduce(r, recv_team.dataPtr(), cnt, MPI_INT, op,
				   0, MyTeam().get_team_comm()) );
	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Allreduce(recv_team.dataPtr(), recv.dataPtr(), cnt, 
					  MPI_INT, op,
					  MyTeam().get_lead_comm()) );
	}
	BL_MPI_REQUIRE( MPI_Bcast(recv.dataPtr(), cnt, MPI_INT,
				  0, MyTeam().get_team_comm()) );
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Allreduce(r,
				      recv.dataPtr(),
				      cnt,
				      MPI_INT,
				      op,
				      Communicator(color)));
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceI, cnt * sizeof(int), false);
    for (int i = 0; i < cnt; i++)
        r[i] = recv[i];
}

void
ParallelDescriptor::util::DoReduceInt (int&   r,
                                       MPI_Op op,
                                       int    cpu)
{
#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceInt()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceI, sizeof(int), cpu);

    int recv;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	int recv_team;
	BL_MPI_REQUIRE( MPI_Reduce(&r, &recv_team, 1, MPI_INT, op,
				   0, MyTeam().get_team_comm()) );

	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Reduce(&recv_team, &recv, 1, MPI_INT, op,
				       RankInLeadComm(cpu), MyTeam().get_lead_comm()) );
	}
	if (sameTeam(cpu)) {
	    BL_MPI_REQUIRE( MPI_Bcast(&recv, 1, MPI_INT,
				      0, MyTeam().get_team_comm()) );
	}
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Reduce(&r,
				   &recv,
				   1,
				   MPI_INT,
				   op,
				   cpu,
				   Communicator()));
    }
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceI, BLProfiler::AfterCall(), cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::util::DoReduceInt (int*   r,
                                       MPI_Op op,
                                       int    cnt,
                                       int    cpu)
{
#ifdef BL_USE_UPCXX
    Mode.set_mpi_mode();
#endif

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceInt()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceI, cnt * sizeof(int), cpu);

    BL_ASSERT(cnt > 0);

    Array<int> recv(cnt);

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Array<long> recv_team(cnt);
	BL_MPI_REQUIRE( MPI_Reduce(r, &recv_team[0], cnt, MPI_LONG, op,
				   0, MyTeam().get_team_comm()) );

	if (isTeamLead()) {
	    BL_MPI_REQUIRE( MPI_Reduce(&recv_team[0], &recv[0], cnt, MPI_LONG, op,
				       RankInLeadComm(cpu), MyTeam().get_lead_comm()) );
	}
	if (sameTeam(cpu)) {
	    BL_MPI_REQUIRE( MPI_Bcast(&recv[0], cnt, MPI_LONG, 
				      0, MyTeam().get_team_comm()) );
	}
    }
    else
#endif
    {
	BL_MPI_REQUIRE( MPI_Reduce(r,
				   recv.dataPtr(),
				   cnt,
				   MPI_INT,
				   op,
				   cpu,
				   Communicator()));
    }
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceI, BLProfiler::AfterCall(), cpu);

    if (ParallelDescriptor::MyProc() == cpu)
    {
        for (int i = 0; i < cnt; i++)
            r[i] = recv[i];
    }
}

void
ParallelDescriptor::ReduceIntSum (int& r, Color color)
{
    util::DoAllReduceInt(r,MPI_SUM,color);
}

void
ParallelDescriptor::ReduceIntMax (int& r, Color color)
{
    util::DoAllReduceInt(r,MPI_MAX,color);
}

void
ParallelDescriptor::ReduceIntMin (int& r, Color color)
{
    util::DoAllReduceInt(r,MPI_MIN,color);
}

void
ParallelDescriptor::ReduceIntSum (int* r, int cnt, Color color)
{
    util::DoAllReduceInt(r,MPI_SUM,cnt,color);
}

void
ParallelDescriptor::ReduceIntMax (int* r, int cnt, Color color)
{
    util::DoAllReduceInt(r,MPI_MAX,cnt,color);
}

void
ParallelDescriptor::ReduceIntMin (int* r, int cnt, Color color)
{
    util::DoAllReduceInt(r,MPI_MIN,cnt,color);
}

void
ParallelDescriptor::ReduceIntSum (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceIntMax (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceIntMin (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceIntSum (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_SUM,cnt,cpu);
}

void
ParallelDescriptor::ReduceIntMax (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_MAX,cnt,cpu);
}

void
ParallelDescriptor::ReduceIntMin (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_MIN,cnt,cpu);
}

void
ParallelDescriptor::ReduceBoolAnd (bool& r, Color color)
{
    if (!isActive(color)) return;

    int src = r; // src is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM,color);

    r = (src == ParallelDescriptor::NProcs(color)) ? true : false;
}

void
ParallelDescriptor::ReduceBoolOr (bool& r, Color color)
{
    if (!isActive(color)) return;

    int src = r; // src is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM,color);

    r = (src == 0) ? false : true;
}

void
ParallelDescriptor::ReduceBoolAnd (bool& r, int cpu)
{
    int src = r; // src is either 0 or 1.

    util::DoReduceInt(src,MPI_SUM,cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = (src == ParallelDescriptor::NProcs()) ? true : false;
}

void
ParallelDescriptor::ReduceBoolOr (bool& r, int cpu)
{
    int src = r; // src is either 0 or 1.

    util::DoReduceInt(src,MPI_SUM,cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = (src == 0) ? false : true;
}

void
ParallelDescriptor::Gather (Real* sendbuf,
                            int   nsend,
                            Real* recvbuf,
                            int   root)
{
    BL_PROFILE_S("ParallelDescriptor::Gather()");
    BL_COMM_PROFILE(BLProfiler::GatherRiRi, BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    BL_ASSERT(root >= 0);
    BL_ASSERT(nsend > 0);
    BL_ASSERT(!(sendbuf == 0));
    BL_ASSERT(!(recvbuf == 0));

    MPI_Datatype typ = Mpi_typemap<Real>::type();

    BL_MPI_REQUIRE( MPI_Gather(sendbuf,
                               nsend,
                               typ,
                               recvbuf,
                               nsend,
                               typ,
                               root,
                               Communicator()));
    BL_COMM_PROFILE(BLProfiler::GatherRiRi, nsend * sizeof(Real), root, BLProfiler::NoTag());
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<char>::type ()
{
    return  MPI_CHAR;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<short>::type ()
{
    return  MPI_SHORT;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<int>::type ()
{
    return  MPI_INT;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<long>::type ()
{
    return  MPI_LONG;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<unsigned char>::type ()
{
    return  MPI_UNSIGNED_CHAR;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<unsigned short>::type ()
{
    return  MPI_UNSIGNED_SHORT;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<unsigned int>::type ()
{
    return  MPI_UNSIGNED;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<unsigned long>::type ()
{
    return  MPI_UNSIGNED_LONG;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<float>::type ()
{
    return  MPI_FLOAT;
}

template <>
MPI_Datatype
ParallelDescriptor::Mpi_typemap<double>::type ()
{
    return  MPI_DOUBLE;
}

void
ParallelDescriptor::Waitsome (Array<MPI_Request>& reqs,
                              int&                completed,
                              Array<int>&         indx,
                              Array<MPI_Status>&  status)
{
    BL_PROFILE_S("ParallelDescriptor::Waitsome()");
#ifdef JEFF_TEST
    std::vector<MPI_Request> rq;
    for (int i = 0; i < reqs.size(); i++)
        if (reqs[i] != MPI_REQUEST_NULL)
            rq.push_back(reqs[i]);
    std::vector<MPI_Status> rst(rq.size());

    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitall, reqs, completed, indx, status, true);
    BL_MPI_REQUIRE( MPI_Waitall(rq.size(), &rq[0], &rst[0]) );
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitall, reqs, completed, indx, status, false);
    completed = rq.size();
    int c = 0;
    for ( int i = 0; i < reqs.size(); ++i )
        if (reqs[i] != MPI_REQUEST_NULL)
    {
	reqs[i] = rq[c];
	status[i] = rst[c];
	indx[c] = i;
	c++;
    }
#else
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitsome, reqs, completed, indx, status, true);
    BL_MPI_REQUIRE( MPI_Waitsome(reqs.size(),
                                 reqs.dataPtr(),
                                 &completed,
                                 indx.dataPtr(),
                                 status.dataPtr()));
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitsome, reqs, completed, indx, status, false);
#endif
}

void
ParallelDescriptor::Bcast(void *buf,
                          int count,
                          MPI_Datatype datatype,
                          int root,
                          MPI_Comm comm)
{
#ifdef BL_LAZY
    int r;
    MPI_Comm_compare(comm, Communicator(), &r);
    if (r == MPI_IDENT)
	Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::Bcast(viMiM)");
    BL_COMM_PROFILE(BLProfiler::BCastTsi, BLProfiler::BeforeCall(), root, BLProfiler::NoTag());

    BL_MPI_REQUIRE( MPI_Bcast(buf,
                              count,
                              datatype,
                              root,
                              comm) );
    int tsize(0);
    BL_MPI_REQUIRE( MPI_Type_size(datatype, &tsize) );
    BL_COMM_PROFILE(BLProfiler::BCastTsi, count * tsize, root, BLProfiler::NoTag());
}


#else /*!BL_USE_MPI*/

void
ParallelDescriptor::StartParallel (int*    argc,
                                   char*** argv,
                                   MPI_Comm)
{
    m_nProcs_all     = 1;
    m_nProcs_comp    = 1;
    m_nProcs_sidecar = 0;
    nSidecarProcs    = 0;

    m_MyId_all     = 0;
    m_MyId_comp    = 0;
    m_MyId_sidecar = myId_notInGroup;

    m_comm_all     = 0;
    m_comm_comp    = 0;
    m_comm_sidecar = 0;
    m_comm_inter   = 0;

    m_MaxTag    = 9000;
}

void
ParallelDescriptor::StartSubCommunicator ()
{
    m_nCommColors = 1;
    m_nProcs_sub  = 1;
    m_MyCommSubColor = Color(0);
    m_MyCommCompColor = Color(0);
    m_comm_sub    = 0;
    m_MyId_sub    = 0;
}

void
ParallelDescriptor::Gather (Real* sendbuf,
			    int   nsend,
			    Real* recvbuf,
			    int   root)
{
    BL_ASSERT(root == 0);
    BL_ASSERT(nsend > 0);
    BL_ASSERT(!(sendbuf == 0));
    BL_ASSERT(!(recvbuf == 0));

    for (int i = 0; i < nsend; ++i)
        recvbuf[i] = sendbuf[i];
}

void
ParallelDescriptor::Message::wait ()
{}

bool
ParallelDescriptor::Message::test ()
{
    return m_finished;
}

void ParallelDescriptor::EndParallel () {}

void ParallelDescriptor::EndSubCommunicator () {}

void ParallelDescriptor::Abort ()
{ 
#ifdef WIN32
    throw;
#else
    std::abort(); 
#endif
}
void ParallelDescriptor::Abort (int)
{ 
#ifdef WIN32
    throw;
#else
    std::abort(); 
#endif
}

const char* ParallelDescriptor::ErrorString (int) { return ""; }

void ParallelDescriptor::Barrier (const std::string &message) {}
void ParallelDescriptor::Barrier (MPI_Comm, const std::string &message) {}

void ParallelDescriptor::Test (MPI_Request&, int&, MPI_Status&) {}
void ParallelDescriptor::IProbe (int, int, int&, MPI_Status&) {}
void ParallelDescriptor::IProbe (int, int, MPI_Comm, int&, MPI_Status&) {}

void ParallelDescriptor::Comm_dup (MPI_Comm, MPI_Comm&) {}

void ParallelDescriptor::ReduceRealMax (Real&,Color) {}
void ParallelDescriptor::ReduceRealMin (Real&,Color) {}
void ParallelDescriptor::ReduceRealSum (Real&,Color) {}

void ParallelDescriptor::ReduceRealMax (Real&,int) {}
void ParallelDescriptor::ReduceRealMin (Real&,int) {}
void ParallelDescriptor::ReduceRealSum (Real&,int) {}

void ParallelDescriptor::ReduceRealMax (Real*,int,Color) {}
void ParallelDescriptor::ReduceRealMin (Real*,int,Color) {}
void ParallelDescriptor::ReduceRealSum (Real*,int,Color) {}

void ParallelDescriptor::ReduceRealMax (Real*,int,int) {}
void ParallelDescriptor::ReduceRealMin (Real*,int,int) {}
void ParallelDescriptor::ReduceRealSum (Real*,int,int) {}

void ParallelDescriptor::ReduceLongAnd (long&,Color) {}
void ParallelDescriptor::ReduceLongSum (long&,Color) {}
void ParallelDescriptor::ReduceLongMax (long&,Color) {}
void ParallelDescriptor::ReduceLongMin (long&,Color) {}

void ParallelDescriptor::ReduceLongAnd (long&,int) {}
void ParallelDescriptor::ReduceLongSum (long&,int) {}
void ParallelDescriptor::ReduceLongMax (long&,int) {}
void ParallelDescriptor::ReduceLongMin (long&,int) {}

void ParallelDescriptor::ReduceLongAnd (long*,int,Color) {}
void ParallelDescriptor::ReduceLongSum (long*,int,Color) {}
void ParallelDescriptor::ReduceLongMax (long*,int,Color) {}
void ParallelDescriptor::ReduceLongMin (long*,int,Color) {}

void ParallelDescriptor::ReduceLongAnd (long*,int,int) {}
void ParallelDescriptor::ReduceLongSum (long*,int,int) {}
void ParallelDescriptor::ReduceLongMax (long*,int,int) {}
void ParallelDescriptor::ReduceLongMin (long*,int,int) {}

void ParallelDescriptor::ReduceIntSum (int&,Color) {}
void ParallelDescriptor::ReduceIntMax (int&,Color) {}
void ParallelDescriptor::ReduceIntMin (int&,Color) {}

void ParallelDescriptor::ReduceIntSum (int&,int) {}
void ParallelDescriptor::ReduceIntMax (int&,int) {}
void ParallelDescriptor::ReduceIntMin (int&,int) {}

void ParallelDescriptor::ReduceIntSum (int*,int,Color) {}
void ParallelDescriptor::ReduceIntMax (int*,int,Color) {}
void ParallelDescriptor::ReduceIntMin (int*,int,Color) {}

void ParallelDescriptor::ReduceIntSum (int*,int,int) {}
void ParallelDescriptor::ReduceIntMax (int*,int,int) {}
void ParallelDescriptor::ReduceIntMin (int*,int,int) {}

void ParallelDescriptor::ReduceBoolAnd (bool&,Color) {}
void ParallelDescriptor::ReduceBoolOr  (bool&,Color) {}

void ParallelDescriptor::ReduceBoolAnd (bool&,int) {}
void ParallelDescriptor::ReduceBoolOr  (bool&,int) {}

void ParallelDescriptor::Bcast(void *, int, MPI_Datatype, int, MPI_Comm) {}

//
// Here so we don't need to include <Utility.H> in <ParallelDescriptor.H>.
//
double
ParallelDescriptor::second ()
{
    return BoxLib::wsecond();
}

void
ParallelDescriptor::Waitsome (Array<MPI_Request>& reqs,
                              int&                completed,
                              Array<int>&         indx,
                              Array<MPI_Status>&  status)
{}

#endif
//
// This function is the same whether or not we're using MPI.
//
int
ParallelDescriptor::SeqNum ()
{
    static int seqno = m_MinTag;
    int result = seqno;

    if (NColors() == 1) { 
	++seqno;
    } else {
	seqno += 2;
    } 

    if (seqno > m_MaxTag) {
	seqno = m_MinTag;
	BL_COMM_PROFILE_TAGWRAP();
    }

    return result;
}

// FIXME: Does COMM_PROFILE work with these two SeqNum functions?

int
ParallelDescriptor::SubSeqNum ()
{
    static int seqno = m_MinTag+1;    
    int result = seqno;

    seqno += 2;

    if (seqno > m_MaxTag) {
	seqno = m_MinTag+1;
	BL_COMM_PROFILE_TAGWRAP();
    }

    return result;
}


#include <BLFort.H>

BL_FORT_PROC_DECL(BL_PD_BARRIER,bl_pd_barrier)()
{
    ParallelDescriptor::Barrier();
}

BL_FORT_PROC_DECL(BL_PD_COMMUNICATOR,bl_pd_communicator)(void* vcomm)
{
    MPI_Comm* comm = reinterpret_cast<MPI_Comm*>(vcomm);

    *comm = ParallelDescriptor::Communicator();
}

BL_FORT_PROC_DECL(BL_PD_MYPROC,bl_pd_myproc)(int* myproc)
{
    *myproc = ParallelDescriptor::MyProc();
}

BL_FORT_PROC_DECL(BL_PD_NPROCS,bl_pd_nprocs)(int* nprocs)
{
    *nprocs = ParallelDescriptor::NProcs();
}

BL_FORT_PROC_DECL(BL_PD_IOPROC,bl_pd_ioproc)(int* ioproc)
{
    *ioproc = ParallelDescriptor::IOProcessorNumber();
}

BL_FORT_PROC_DECL(BL_PD_IS_IOPROC,bl_pd_is_ioproc)(int* ioproc)
{
    *ioproc = ParallelDescriptor::IOProcessor()?1:0;
}

BL_FORT_PROC_DECL(BL_PD_SECOND,bl_pd_second)(double* r)
{
    *r = ParallelDescriptor::second();
}

#ifdef BL_USE_FLOAT
BL_FORT_PROC_DECL(BL_PD_REDUCE_REAL_MAX_TO_IOPROC,bl_pd_reduce_real_max_to_ioproc)(float* r)
{
    ParallelDescriptor::ReduceRealMax(*r,ParallelDescriptor::IOProcessorNumber());
}

BL_FORT_PROC_DECL(BL_PD_REDUCE_REAL_SUM_TO_IOPROC,bl_pd_reduce_real_sum_to_ioproc)(float* r)
{
    ParallelDescriptor::ReduceRealSum(*r,ParallelDescriptor::IOProcessorNumber());
}
#else
BL_FORT_PROC_DECL(BL_PD_REDUCE_REAL_MAX_TO_IOPROC,bl_pd_reduce_real_max_to_ioproc)(double* r)
{
    ParallelDescriptor::ReduceRealMax(*r,ParallelDescriptor::IOProcessorNumber());
}

BL_FORT_PROC_DECL(BL_PD_REDUCE_REAL_SUM_TO_IOPROC,bl_pd_reduce_real_sum_to_ioproc)(double* r)
{
    ParallelDescriptor::ReduceRealSum(*r,ParallelDescriptor::IOProcessorNumber());
}
#endif

BL_FORT_PROC_DECL(BL_PD_ABORT,bl_pd_abort)()
{
    ParallelDescriptor::Abort();
}

#ifdef BL_USE_MPI
namespace ParallelDescriptor
{
#ifndef BL_AMRPROF
template <> MPI_Datatype Mpi_typemap<IntVect>::type()
{
    static MPI_Datatype mine(MPI_DATATYPE_NULL);
    if ( mine == MPI_DATATYPE_NULL )
    {
	IntVect iv[2];	// Used to construct the data types
	MPI_Datatype types[] = {
	    MPI_LB,
	    MPI_INT,
	    MPI_UB};
	int blocklens[] = { 1, BL_SPACEDIM, 1};
	MPI_Aint disp[3];
	int n = 0;
	BL_MPI_REQUIRE( MPI_Address(&iv[0],      &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[0].vect, &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[1],      &disp[n++]) );
	for ( int i = n-1; i >= 0; i-- )
	{
	    disp[i] -= disp[0];
	}
	BL_MPI_REQUIRE( MPI_Type_struct(n, blocklens, disp, types, &mine) );
	BL_MPI_REQUIRE( MPI_Type_commit( &mine ) );
    }
    return mine;
}

template <> MPI_Datatype Mpi_typemap<IndexType>::type()
{
    static MPI_Datatype mine(MPI_DATATYPE_NULL);
    if ( mine == MPI_DATATYPE_NULL )
    {
	IndexType iv[2];	// Used to construct the data types
	MPI_Datatype types[] = {
	    MPI_LB,
	    MPI_UNSIGNED,
	    MPI_UB};
	int blocklens[] = { 1, 1, 1};
	MPI_Aint disp[3];
	int n = 0;
	BL_MPI_REQUIRE( MPI_Address(&iv[0],       &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[0].itype, &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[1],       &disp[n++]) );
	for ( int i = n-1; i >= 0; i-- )
	{
	    disp[i] -= disp[0];
	}
	BL_MPI_REQUIRE( MPI_Type_struct(n, blocklens, disp, types, &mine) );
	BL_MPI_REQUIRE( MPI_Type_commit( &mine ) );
    }
    return mine;
}

template <> MPI_Datatype Mpi_typemap<Box>::type()
{
    static MPI_Datatype mine(MPI_DATATYPE_NULL);
    if ( mine == MPI_DATATYPE_NULL )
    {
	Box iv[2];	// Used to construct the data types
	MPI_Datatype types[] = {
	    MPI_LB,
	    Mpi_typemap<IntVect>::type(),
	    Mpi_typemap<IntVect>::type(),
	    Mpi_typemap<IndexType>::type(),
	    MPI_UB};
	int blocklens[] = { 1, 1, 1, 1, 1};
	MPI_Aint disp[5];
	int n = 0;
	BL_MPI_REQUIRE( MPI_Address(&iv[0],          &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[0].smallend, &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[0].bigend,   &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[0].btype,    &disp[n++]) );
	BL_MPI_REQUIRE( MPI_Address(&iv[1],          &disp[n++]) );
	for ( int i = n-1; i >= 0; i-- )
	{
	    disp[i] -= disp[0];
	}
	BL_MPI_REQUIRE( MPI_Type_struct(n, blocklens, disp, types, &mine) );
	BL_MPI_REQUIRE( MPI_Type_commit( &mine ) );
    }
    return mine;
}
#endif
}
#endif

void
ParallelDescriptor::ReadAndBcastFile (const std::string& filename,
                                      Array<char>&       charBuf,
				      bool               bExitOnError,
				      const MPI_Comm    &comm)
{
    enum { IO_Buffer_Size = 40960 * 32 };

#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif

    Array<Setbuf_Char_Type> io_buffer(IO_Buffer_Size);

    int fileLength = 0, fileLengthPadded;

    std::ifstream iss;

    if (ParallelDescriptor::IOProcessor())
    {
        iss.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        iss.open(filename.c_str(), std::ios::in);
        if ( ! iss.good())
        {
	  if(bExitOnError) {
            BoxLib::FileOpenFailed(filename);
	  } else {
            fileLength = -1;
	  }
        } else {
          iss.seekg(0, std::ios::end);
          fileLength = iss.tellg();
          iss.seekg(0, std::ios::beg);
	}
    }
    ParallelDescriptor::Bcast(&fileLength, 1,
                              ParallelDescriptor::IOProcessorNumber(), comm);

    if(fileLength == -1) {
      return;
    }

    fileLengthPadded = fileLength + 1;
    fileLengthPadded += fileLengthPadded % 8;
    charBuf.resize(fileLengthPadded);
    if (ParallelDescriptor::IOProcessor())
    {
        iss.read(charBuf.dataPtr(), fileLength);
        iss.close();
    }
    ParallelDescriptor::Bcast(charBuf.dataPtr(), fileLengthPadded,
                              ParallelDescriptor::IOProcessorNumber(), comm);
    charBuf[fileLength] = '\0';
}


void
ParallelDescriptor::SidecarProcess ()
{
#ifdef IN_TRANSIT
#ifdef BL_USE_MPI
    bool finished(false);
    int signal(-1);
    while (!finished)
    {
        // Receive the signal from the compute group.
        ParallelDescriptor::Bcast(&signal, 1, 0, ParallelDescriptor::CommunicatorInter());

        // Process the  signal with the set of user-provided signal handlers
        for (ParallelDescriptor::SH_LIST::iterator it=The_Signal_Handler_List.begin();
             it != The_Signal_Handler_List.end() && signal != ParallelDescriptor::SidecarQuitSignal;
             ++it)
        {
          signal = (* *it)(signal);
        }

        if (signal == ParallelDescriptor::SidecarQuitSignal)
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Sidecars received the quit signal." << std::endl;
            finished = true;
            break;
        }
    }
    if (ParallelDescriptor::IOProcessor())
      std::cout << "===== SIDECARS DONE. EXITING ... =====" << std::endl;
#endif
#endif
}


#ifndef BL_AMRPROF
void
ParallelDescriptor::StartTeams ()
{
    int team_size = 1;
    int do_team_reduce = 0;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    ParmParse pp("team");
    pp.query("size", team_size);
    pp.query("reduce", do_team_reduce);
#endif

    int nprocs = ParallelDescriptor::NProcs();
    int rank   = ParallelDescriptor::MyProc();

    if (nprocs % team_size != 0)
	BoxLib::Abort("Number of processes not divisible by team size");

    m_Team.m_numTeams    = nprocs / team_size;
    m_Team.m_size        = team_size;
    m_Team.m_color       = rank / team_size;
    m_Team.m_lead        = m_Team.m_color * team_size;
    m_Team.m_rankInTeam  = rank - m_Team.m_lead;

    m_Team.m_do_team_reduce = team_size > 0 && do_team_reduce;

#if defined(BL_USE_UPCXX) || defined(BL_USE_MPI3)
    {
	MPI_Group grp, team_grp, lead_grp;
	BL_MPI_REQUIRE( MPI_Comm_group(ParallelDescriptor::Communicator(), &grp) );
	int team_ranks[team_size];
	for (int i = 0; i < team_size; ++i) {
	    team_ranks[i] = MyTeamLead() + i;
	}
	BL_MPI_REQUIRE( MPI_Group_incl(grp, team_size, team_ranks, &team_grp) );
	BL_MPI_REQUIRE( MPI_Comm_create(ParallelDescriptor::Communicator(), 
					team_grp, &m_Team.m_team_comm) );

	std::vector<int>lead_ranks(m_Team.m_numTeams);
	for (int i = 0; i < lead_ranks.size(); ++i) {
	    lead_ranks[i] = i * team_size;
	}
	BL_MPI_REQUIRE( MPI_Group_incl(grp, lead_ranks.size(), &lead_ranks[0], &lead_grp) );
	BL_MPI_REQUIRE( MPI_Comm_create(ParallelDescriptor::Communicator(), 
					lead_grp, &m_Team.m_lead_comm) );

        BL_MPI_REQUIRE( MPI_Group_free(&grp) );
        BL_MPI_REQUIRE( MPI_Group_free(&team_grp) );
        BL_MPI_REQUIRE( MPI_Group_free(&lead_grp) );

#ifdef BL_USE_UPCXX
	upcxx::team* team;
	upcxx::team_all.split(MyTeamColor(), MyRankInTeam(), team);
        m_Team.m_upcxx_team = team;
#endif
    }
#endif
}
#endif

void
ParallelDescriptor::EndTeams ()
{
    m_Team.clear();
}


bool
ParallelDescriptor::MPIOneSided ()
{
    static bool do_onesided = false;

#ifndef BL_AMRPROF
#if defined(BL_USE_MPI3) && !defined(BL_USE_UPCXX)
    static bool first = true;
    if (first) {
	first = false;
	ParmParse pp("mpi");
	pp.query("onesided",do_onesided);
    }
#endif
#endif

    return do_onesided;
}
