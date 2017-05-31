
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <sstream>
#include <stack>
#include <list>
#include <chrono>

#include <AMReX_Utility.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLFort.H>
#include <AMReX_ParallelDescriptor.H>

#ifdef BL_USE_MPI
#include <AMReX_ccse-mpi.H>
#endif

#ifdef BL_USE_FORTRAN_MPI
extern "C" {
    void bl_fortran_mpi_comm_init (int fcomm);
    void bl_fortran_mpi_comm_free ();
    void bl_fortran_sidecar_mpi_comm_free (int fcomm);
    void bl_fortran_set_nprocs_sidecar(int nSidecarProcs, int inWhichSidecar,
                    int m_nProcs_all, int m_nProcs_comp, int *m_nProcs_sidecar,
                    int fcomma, int fcommc, int *fcomms,
                    int fgrpa, int fgrpc, int *fgrps,
                    int m_MyId_all, int m_MyId_comp, int m_MyId_sidecar);
    void bl_fortran_set_nprocs_sidecar_to_zero(int m_nProcs_all,
                    int fcomma, int fcommc,
                    int fgrpa, int fgrpc,
                    int m_MyId_all, int m_MyId_comp);
}
#endif

#ifndef BL_AMRPROF
#include <AMReX_ParmParse.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

namespace ParallelDescriptor
{
    const int myId_undefined   = -11;
    const int myId_notInGroup  = -22;
    const int nProcs_undefined = -33;
    const int notInSidecar  = -44;

#ifdef BL_USE_MPI
    //
    // My processor IDs.
    //
    int m_MyId_all         = myId_undefined;
    int m_MyId_comp        = myId_undefined;
    int m_MyId_sub         = myId_undefined;
    int m_MyId_sidecar     = myId_undefined;
    //
    // The number of processors.
    //
    int m_nProcs_all     = nProcs_undefined;
    int m_nProcs_comp    = nProcs_undefined;
    int m_nProcs_sub     = nProcs_undefined;
    Array<int> m_num_procs_clr;
    Array<int> m_first_procs_clr;
    Array<int> m_nProcs_sidecar;
    int nSidecars = 0;
    int inWhichSidecar = notInSidecar;
    //
    // Team
    //
    ProcessTeam m_Team;
    //
    // AMReX's Communicators
    //
    MPI_Comm m_comm_all     = MPI_COMM_NULL;    // for all ranks, probably MPI_COMM_WORLD
    MPI_Comm m_comm_comp    = MPI_COMM_NULL;    // for the ranks doing computations
    MPI_Comm m_comm_sub     = MPI_COMM_NULL;    // for sub computation procs
    Array<MPI_Comm> m_comm_sidecar;             // for the ranks in the sidecar
    Array<MPI_Comm> m_comm_inter;               // for communicating between comp and sidecar
    Array<MPI_Comm> m_comm_both;                // for communicating within comp and a sidecar
    //
    // AMReX's Groups
    //
    MPI_Group m_group_all     = MPI_GROUP_NULL;
    MPI_Group m_group_comp    = MPI_GROUP_NULL;
    Array<MPI_Group> m_group_sidecar;
#else
    //  Set these for non-mpi codes that do not call amrex::Initialize(...)
    int m_MyId_all         = 0;
    int m_MyId_comp        = 0;
    int m_MyId_sub         = 0;
    int m_MyId_sidecar     = myId_notInGroup;
    //
    int m_nProcs_all     = 1;
    int m_nProcs_comp    = 1;
    int m_nProcs_sub     = 1;
    Array<int> m_num_procs_clr;
    Array<int> m_first_procs_clr;
    Array<int> m_nProcs_sidecar;
    int nSidecars = 0;
    int inWhichSidecar = notInSidecar;
    //
    ProcessTeam m_Team;
    //
    MPI_Comm m_comm_all     = 0;
    MPI_Comm m_comm_comp    = 0;
    MPI_Comm m_comm_sub     = 0;
    Array<MPI_Comm> m_comm_sidecar;
    Array<MPI_Comm> m_comm_inter;
    //
    // AMReX's Groups
    //
    MPI_Group m_group_all     = 0;
    MPI_Group m_group_comp    = 0;
    MPI_Group m_group_sidecar = 0;
#endif

    int m_nCommColors = 1;
    Color m_MyCommSubColor;
    Color m_MyCommCompColor;

    int m_MinTag = 1000, m_MaxTag = -1, m_MaxTag_MPI = -1, tagBuffer = 32;

    const int ioProcessor = 0;

#ifdef BL_USE_UPCXX
    UPCXX_MPI_Mode Mode;
#endif    

#ifdef BL_USE_MPI3
    MPI_Win cp_win;
    MPI_Win fb_win;
#endif
  
    namespace util
    {
	//
	// Reduce helper functions.
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
	    snprintf(buf, N, "AMReX MPI Error: File %s, line %d, %s: %s",
                     file, line, call, ParallelDescriptor::ErrorString(status));
	}
	else
	{
	    snprintf(buf, N, "AMReX MPI Error: File %s, line %d, %s",
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
	amrex::Error(the_message_string(file, line, str, rc));
    }
}

void
ParallelDescriptor::Abort (int errorcode, bool backtrace)
{
    if (backtrace) {
	BLBackTrace::handler(errorcode);
    } else {
	MPI_Abort(CommunicatorAll(), errorcode);
    }
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
    if ( !m_finished ) amrex::Error("Message::tag: Not Finished!");
    return m_stat.MPI_TAG;
}

int
ParallelDescriptor::Message::pid () const
{
    if ( !m_finished ) amrex::Error("Message::pid: Not Finished!");
    return m_stat.MPI_SOURCE;
}

size_t
ParallelDescriptor::Message::count () const
{
    if ( m_type == MPI_DATATYPE_NULL ) amrex::Error("Message::count: Bad Type!");
    if ( !m_finished ) amrex::Error("Message::count: Not Finished!");
    int cnt;
    BL_MPI_REQUIRE( MPI_Get_count(&m_stat, m_type, &cnt) );
    return cnt;
}

void
ParallelDescriptor::StartParallel (int*    argc,
                                   char*** argv,
                                   MPI_Comm mpi_comm)
{
    int sflag(0);

    BL_MPI_REQUIRE( MPI_Initialized(&sflag) );

    if ( ! sflag) {
	BL_MPI_REQUIRE( MPI_Init(argc, argv) );
    }
    
    BL_MPI_REQUIRE( MPI_Comm_dup(mpi_comm, &m_comm_all) );

    // ---- find the maximum value for a tag
    int flag(0), *attrVal;
    BL_MPI_REQUIRE( MPI_Attr_get(m_comm_all, MPI_TAG_UB, &attrVal, &flag) );
    if(flag) {
      m_MaxTag_MPI = *attrVal;
      m_MaxTag = m_MaxTag_MPI - tagBuffer;  // ---- buffer for sidecar tags
      m_MaxTag = std::max(m_MaxTag, 9000);
    } else {
      m_MaxTag = 9000;
      m_MaxTag_MPI = m_MaxTag;
    }
    BL_COMM_PROFILE_TAGRANGE(m_MinTag, m_MaxTag);

    BL_MPI_REQUIRE( MPI_Comm_size(CommunicatorAll(), &m_nProcs_all) );
    BL_MPI_REQUIRE( MPI_Comm_rank(CommunicatorAll(), &m_MyId_all) );
    BL_MPI_REQUIRE( MPI_Comm_group(CommunicatorAll(), &m_group_all) );

#ifdef BL_USE_MPI3
    int mpi_version, mpi_subversion;
    BL_MPI_REQUIRE( MPI_Get_version(&mpi_version, &mpi_subversion) );
    if (mpi_version < 3) amrex::Abort("MPI 3 is needed because USE_MPI3=TRUE");
#endif

    SetNProcsSidecars(0);  // ---- users resize these later

    //
    // Wait until all other processes are properly started.
    //
    BL_MPI_REQUIRE( MPI_Barrier(CommunicatorAll()) );
}

void
ParallelDescriptor::EndParallel ()
{
    BL_ASSERT(m_MyId_all   != myId_undefined);
    BL_ASSERT(m_nProcs_all != nProcs_undefined);

    if(m_group_comp != MPI_GROUP_NULL && m_group_comp != m_group_all) {
      BL_MPI_REQUIRE( MPI_Group_free(&m_group_comp) );
    }
    if(m_comm_comp != MPI_COMM_NULL && m_comm_comp != m_comm_all) {
      BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_comp) );
    }
    if(m_group_all != MPI_GROUP_NULL) {
      BL_MPI_REQUIRE( MPI_Group_free(&m_group_all) );
    }
// bl_fortran_mpi_comm_free() has already freed the global communicator
#ifndef BL_USE_FORTRAN_MPI
    if(m_comm_all != MPI_COMM_NULL) {
      BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_all) );
    }
#endif

    BL_MPI_REQUIRE( MPI_Finalize() );
}

/* Given `rk_clrd', i.e. the rank in the colored `comm', return the
rank in the global `CommComp' */
int
ParallelDescriptor::Translate(int rk_clrd,Color clr)
{
  if (clr == DefaultColor())
  {
    return rk_clrd;
  }
  else if (clr.valid())
  {
    return m_first_procs_clr[clr.to_int()]+rk_clrd;
  }
  else
  {
    return MPI_PROC_NULL;
  }
}

/* Define variables `m_nProcs_sub', `m_MyCommSubColor',
`m_num_procs_clr', and `m_first_procs_clr' */
void
ParallelDescriptor::init_clr_vars()
{
  if (m_nProcs_comp < m_nCommColors)
  {
    amrex::Abort("Need `nProcs >= nColors'");
  }
  /* Number of proc.'s per color */
  m_nProcs_sub = m_nProcs_comp/m_nCommColors;
  /* Remaining proc.'s */
  int const numRemProcs = m_nProcs_comp-m_nCommColors*m_nProcs_sub;
  /* All numbers of proc.'s per color (clear) */
  m_num_procs_clr.clear();
  /* All numbers of proc.'s per color (init.) */
  m_num_procs_clr.resize(m_nCommColors,m_nProcs_sub);
  /* Distribute remaining proc.'s */
  for (int clr = 0; clr < numRemProcs; ++clr)
  {
    ++(m_num_procs_clr[clr]);
  }
  /* All first proc.'s (clear) */
  m_first_procs_clr.clear();
  /* All first proc.'s (init.) */
  m_first_procs_clr.resize(m_nCommColors,0);
  /* Add the previous `clr's `nProc's */
  for (int clr = 1; clr < m_nCommColors; ++clr)
  {
      m_first_procs_clr[clr] =
        m_first_procs_clr[clr-1]+m_num_procs_clr[clr-1];
  }
  /* My proc. must be larger than my first proc. */
  int clr = 0; int myClr = -1;
  while (clr < m_nCommColors && MyProc() > m_first_procs_clr[clr]-1)
  {
    ++clr; ++myClr;
  }
  /* Possibly adjust number of proc.'s per color */
  m_nProcs_sub = m_num_procs_clr[myClr];
  /* Define `Color' */
  m_MyCommSubColor = Color(myClr);

  return;
} /* `init_clr_vars( ...' */

void
ParallelDescriptor::SetNProcsSidecars (const Array<int> &compRanksInAll,
                                       const Array<Array<int> > &sidecarRanksInAll, bool printRanks)
{
    BL_ASSERT(compRanksInAll.size() > 0);
    BL_ASSERT(m_MyId_all != myId_undefined);

    bool inComp(false);
    nSidecars = sidecarRanksInAll.size();

    MPI_Comm commTemp = Communicator();  // ---- save to compare later

    // ---- check validity of the rank arrays and set inComp
    if(compRanksInAll[0] != 0) {  // ---- we require this for now
      amrex::Abort("**** Error in SetNProcsSidecars:  compRanksInAll[0] != 0");
    }
    std::set<int> rankSet;
    for(int i(0); i < compRanksInAll.size(); ++i) {
      rankSet.insert(compRanksInAll[i]);
      if(m_MyId_all == compRanksInAll[i]) {
        inComp = true;
      }
    }
    inWhichSidecar = notInSidecar;
    for(int i(0); i < sidecarRanksInAll.size(); ++i) {
      for(int j(0); j < sidecarRanksInAll[i].size(); ++j) {
        rankSet.insert(sidecarRanksInAll[i][j]);
	if(m_MyId_all == sidecarRanksInAll[i][j]) {
	  inWhichSidecar = i;
	}
      }
    }

    if(rankSet.size() != m_nProcs_all) {
      std::cerr << "**** rankSet.size() != m_nProcs_all:  " << rankSet.size()
                << " != " << m_nProcs_all << std::endl;
      amrex::Abort("**** Error in SetNProcsSidecars:  rankSet.size() != m_nProcs_all.");
    }
    int rankCheck(0);
    std::set<int>::const_iterator cit;
    for(cit = rankSet.begin(); cit != rankSet.end(); ++cit) {
      if(*cit != rankCheck) {
        amrex::Abort("**** Error in SetNProcsSidecars:  rankSet is not correct.");
      }
      ++rankCheck;
    }

    // ---- check to ensure ranks are only removed from comp or added to comp but not both
    // ---- during one call to this function.  this is currently disallowed because of
    // ---- data movement issues into and out of the computation, it is not a communicator problem
    // ---- this resctriction could be removed if there is a valid use case
    if(m_MyId_all == 0 && m_nProcs_comp > 0) {
      Array<int> oldCompRanks(m_nProcs_comp), oldCompRanksInAll(m_nProcs_comp);
      for(int i(0); i < oldCompRanksInAll.size(); ++i) {
        oldCompRanks[i] = i;
      }
      if(m_group_comp != MPI_GROUP_NULL && m_group_comp != m_group_all) {
        BL_MPI_REQUIRE( MPI_Group_translate_ranks(m_group_comp, oldCompRanks.size(), oldCompRanks.dataPtr(),
                                                  m_group_all,  oldCompRanksInAll.dataPtr()) );
      }
      std::set<int> ocrSet, criaSet;
      for(int i(0); i < oldCompRanksInAll.size(); ++i) {
	ocrSet.insert(oldCompRanksInAll[i]);
        if(printRanks) {
          std::cout << "oooo i oldCompRanks[i] oldCompRanksInAll[i] = " << i << "  "
                    << oldCompRanks[i] << "  " << oldCompRanksInAll[i] << std::endl;
	}
      }
      for(int i(0); i < compRanksInAll.size(); ++i) {
	criaSet.insert(compRanksInAll[i]);
      }

      for(int i(0); i < oldCompRanksInAll.size(); ++i) {
        if(criaSet.find(oldCompRanksInAll[i]) != criaSet.end()) {  // ---- erase from both sets
	  criaSet.erase(oldCompRanksInAll[i]);
	  ocrSet.erase(oldCompRanksInAll[i]);
	}
      }
      if(printRanks) {
        std::cout << "criaSet.size() ocrSet.size() = " << criaSet.size() << "  " << ocrSet.size() << std::endl;
        std::set<int>::iterator it;
        for(it = criaSet.begin(); it != criaSet.end(); ++it) {
          std::cout << "criaSet = " << *it << std::endl;
        }
        for(it = ocrSet.begin(); it != ocrSet.end(); ++it) {
          std::cout << "ocrSet = " << *it << std::endl;
        }
      }
      if(ocrSet.size() > 0 && criaSet.size() > 0) {  // ---- this is currently not allowed
        amrex::Abort("**** Error in SetNProcsSidecars:  adding and removing ranks from comp not supported.");
      }
    }

    // ---- print the ranks
    if(m_MyId_all == 0 && printRanks) {
      std::cout << "cccc nCompProcs = " << compRanksInAll.size() << std::endl;
      for(int i(0); i < compRanksInAll.size(); ++i) {
        std::cout << "cccc cccc compRanksInAll[" << i << "] = " << compRanksInAll[i] << std::endl;
      }
      std::cout << "ssss nSidecars = " << sidecarRanksInAll.size() << std::endl;
      for(int i(0); i < sidecarRanksInAll.size(); ++i) {
        std::cout << "ssss ssss sidecar[" << i << "].size() = " << sidecarRanksInAll[i].size() << std::endl;
        for(int j(0); j < sidecarRanksInAll[i].size(); ++j) {
          std::cout << "ssss ssss ssss sidecarRanksInAll[" << i << "][" << j << "] = "
	            << sidecarRanksInAll[i][j] << std::endl;
        }
      }
    }


    // ---- free existing groups and communicators and reinitialize values
    if(m_comm_comp != MPI_COMM_NULL && m_comm_comp != m_comm_all) {
      BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_comp) );
      m_comm_comp = MPI_COMM_NULL;
    }
    for(int i(0); i < m_comm_sidecar.size(); ++i) {
      if(m_comm_sidecar[i] != MPI_COMM_NULL) {
        BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_sidecar[i]) );
        m_comm_sidecar[i] = MPI_COMM_NULL;
      }
    }
    m_comm_sidecar.clear();

    for(int i(0); i < m_comm_inter.size(); ++i) {
      if(m_comm_inter[i] != MPI_COMM_NULL) {
        BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_inter[i]) );
      }
    }
    m_comm_inter.clear();

    for(int i(0); i < m_comm_both.size(); ++i) {
      if(m_comm_both[i] != MPI_COMM_NULL) {
        BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_both[i]) );
      }
    }
    m_comm_both.clear();

    if(m_comm_sub != MPI_COMM_NULL && m_comm_sub != commTemp) {
      BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_sub) );
      m_comm_sub = MPI_COMM_NULL;
    }

    if(m_group_comp != MPI_GROUP_NULL && m_group_comp != m_group_all) {
      BL_MPI_REQUIRE( MPI_Group_free(&m_group_comp) );
      m_group_comp = MPI_GROUP_NULL;
    }
    for(int i(0); i < m_group_sidecar.size(); ++i) {
      if(m_group_sidecar[i] != MPI_GROUP_NULL) {
        BL_MPI_REQUIRE( MPI_Group_free(&m_group_sidecar[i]) );
      }
    }
    m_group_sidecar.clear();

    m_nProcs_comp = nProcs_undefined;
    m_MyId_comp   = myId_undefined;
    m_nProcs_sidecar.clear();
    m_MyId_sidecar = myId_undefined;


    // ---- now create new groups and communicators
    if(nSidecars > 0) {
      m_nProcs_sidecar.resize(nSidecars, nProcs_undefined);
      m_group_sidecar.resize(nSidecars, MPI_GROUP_NULL);
      m_comm_sidecar.resize(nSidecars, MPI_COMM_NULL);
      int* castptr = (int*)compRanksInAll.dataPtr();
      BL_MPI_REQUIRE( MPI_Group_incl(m_group_all, compRanksInAll.size(), castptr, &m_group_comp) );
      BL_MPI_REQUIRE( MPI_Comm_create(m_comm_all, m_group_comp, &m_comm_comp) );

      for(int i(0); i < sidecarRanksInAll.size(); ++i) {
	if(sidecarRanksInAll[i].size() > 0) {
          int* castptr2 = (int*)sidecarRanksInAll[i].dataPtr();
          BL_MPI_REQUIRE( MPI_Group_incl(m_group_all, sidecarRanksInAll[i].size(), castptr2,
	                                 &m_group_sidecar[i]) );
          BL_MPI_REQUIRE( MPI_Comm_create(m_comm_all, m_group_sidecar[i], &m_comm_sidecar[i]) );
          BL_MPI_REQUIRE( MPI_Group_size(m_group_sidecar[i], &m_nProcs_sidecar[i]) );
	  BL_ASSERT(m_nProcs_sidecar[i] == sidecarRanksInAll[i].size());
	} else {
          m_nProcs_sidecar[i] = 0;
	}
      }
      BL_MPI_REQUIRE( MPI_Group_size(m_group_comp, &m_nProcs_comp) );
    } else {
      m_comm_comp   = m_comm_all;
      m_group_comp  = m_group_all;
      m_nProcs_comp = m_nProcs_all;
    }

    if(nSidecars > tagBuffer) {
      tagBuffer = nSidecars;
      m_MaxTag = m_MaxTag_MPI - tagBuffer;
      m_MaxTag = std::max(m_MaxTag, 9000);
      BL_COMM_PROFILE_TAGRANGE(m_MinTag, m_MaxTag);
    }

    // ---- create the inter communicators, but only between comp and each sidecar
    // ---- the user will have to create any sidecar to sidecar communicators
    if(nSidecars > 0) {
      m_comm_inter.resize(nSidecars, MPI_COMM_NULL);
      for(int i(0); i < sidecarRanksInAll.size(); ++i) {
	if(sidecarRanksInAll[i].size() > 0) {
          int tag(m_MaxTag + 1 + i);

          if(inComp) {                          // ---- in the computation group
	    if(inWhichSidecar >= 0) {
              amrex::Abort("**** Error 0:  bad inWhichSidecar in SetNProcsSidecars()");
	    }
            BL_MPI_REQUIRE( MPI_Group_rank(m_group_comp, &m_MyId_comp) );
            BL_MPI_REQUIRE( MPI_Intercomm_create(m_comm_comp, 0, m_comm_all, sidecarRanksInAll[i][0],
	                         tag, &m_comm_inter[i]) );
	    m_MyId_sidecar = myId_notInGroup;
          } else {                              // ---- in a sidecar group
	    if(inWhichSidecar < 0) {
              amrex::Abort("**** Error 1:  bad inWhichSidecar in SetNProcsSidecars()");
	    }
	    if(inWhichSidecar == i) {
              BL_MPI_REQUIRE( MPI_Group_rank(m_group_sidecar[i], &m_MyId_sidecar) );
              BL_MPI_REQUIRE( MPI_Intercomm_create(m_comm_sidecar[i], 0, m_comm_all, 0, tag, &m_comm_inter[i]) );
	    } else {
              if(m_MyId_sidecar == myId_undefined) {
	        m_MyId_sidecar = myId_notInGroup;
	      }
	    }
	    m_MyId_comp = myId_notInGroup;
          }
	}
      }
    } else {
      m_MyId_comp = m_MyId_all;
      m_MyId_sidecar = myId_notInGroup;
    }

    // ---- create communicator spanning comp and each sidecar
    if(nSidecars > 0) {
      m_comm_both.resize(nSidecars, MPI_COMM_NULL);
      for(int i(0); i < m_group_sidecar.size(); ++i) {
	MPI_Group groupBoth;
        BL_MPI_REQUIRE( MPI_Group_union(m_group_comp, m_group_sidecar[i], &groupBoth) );
        BL_MPI_REQUIRE( MPI_Comm_create(m_comm_all, groupBoth, &m_comm_both[i]) );
      }
    }

    // ---- recreate the color sub communicator
    if(m_nCommColors > 1) {
      /* Define variables `m_nProcs_sub', `m_MyCommSubColor',
      `m_num_procs_clr', and `m_first_procs_clr' */
      init_clr_vars();
      /* Special color for `CommComp's color */
      m_MyCommCompColor = Color(m_nCommColors); 

      BL_MPI_REQUIRE( MPI_Comm_split(Communicator(), m_MyCommSubColor.to_int(), MyProc(), &m_comm_sub) );
      BL_MPI_REQUIRE( MPI_Comm_rank(m_comm_sub, &m_MyId_sub) );
    } else {
      m_nProcs_sub  = NProcs();
      m_num_procs_clr.resize(1,0);
      m_first_procs_clr.resize(1,0);
      m_MyCommSubColor = Color(0);
      m_MyCommCompColor = Color(0);
      m_comm_sub    = Communicator();
      m_MyId_sub    = MyProc();
    }

    // ---- more error checking
    if(m_MyId_all     == myId_undefined ||
       m_MyId_comp    == myId_undefined ||
       m_MyId_sidecar == myId_undefined)
    {
      std::cerr << "m_MyId_all m_MyId_comp m_MyId_sidecar = " << m_MyId_all << "  "
	          << m_MyId_comp << "  " << m_MyId_sidecar << std::endl;
      amrex::Abort("**** Error 2:  bad MyId in ParallelDescriptor::SetNProcsSidecars()");
    }
    if(m_nProcs_all  == nProcs_undefined ||
       m_nProcs_comp == nProcs_undefined)
    {
      std::cerr << "m_nProcs_all m_nProcs_comp = " << m_nProcs_all << "  "
	          << m_nProcs_comp << std::endl;
      amrex::Abort("**** Error 3:  bad nProcs in ParallelDescriptor::SetNProcsSidecars()");
    }
    int nSCSum(0);
    for(int i(0); i < sidecarRanksInAll.size(); ++i) {
      nSCSum += sidecarRanksInAll[i].size();
      if(m_nProcs_sidecar[i] == nProcs_undefined) {
        std::cerr << "m_nProcs_sidecar[" << i << "] = " << m_nProcs_sidecar[i] << std::endl;
        amrex::Abort("**** Error 4:  bad m_nProcs_sidecar in ParallelDescriptor::SetNProcsSidecars()");
      }
    }
    if(m_nProcs_comp + nSCSum != m_nProcs_all) {
      std::cerr << "m_nProcs_all m_nProcs_comp + nSCSum = " << m_nProcs_all << "  "
                << m_nProcs_comp + nSCSum << std::endl;
      amrex::Abort("**** Error 5:  bad nProcs in ParallelDescriptor::SetNProcsSidecars()");
    }

#ifdef BL_USE_FORTRAN_MPI
    if(nSidecars > 0) {
      int fcomma = MPI_Comm_c2f(m_comm_all);
      int fcommc = MPI_Comm_c2f(m_comm_comp);
      Array<int> fcomms(nSidecars, -1);
      for(int i(0); i < fcomms.size(); ++i) {
        fcomms[i] = MPI_Comm_c2f(m_comm_sidecar[i]);
      }
      int fgrpa  = MPI_Group_c2f(m_group_all);
      int fgrpc  = MPI_Group_c2f(m_group_comp);
      Array<int> fgrps(nSidecars, -2);
      for(int i(0); i < fgrps.size(); ++i) {
        fgrps[i] = MPI_Group_c2f(m_group_sidecar[i]);
      }
      bl_fortran_set_nprocs_sidecar(nSidecars, inWhichSidecar,
                                    m_nProcs_all, m_nProcs_comp, m_nProcs_sidecar.dataPtr(),
                                    fcomma, fcommc, fcomms.dataPtr(),
                                    fgrpa, fgrpc, fgrps.dataPtr(),
                                    m_MyId_all, m_MyId_comp, m_MyId_sidecar);
    } else {
      int fcomma = MPI_Comm_c2f(m_comm_all);
      int fcommc = MPI_Comm_c2f(m_comm_comp);
      int fgrpa  = MPI_Group_c2f(m_group_all);
      int fgrpc  = MPI_Group_c2f(m_group_comp);
      bl_fortran_set_nprocs_sidecar_to_zero(m_nProcs_all, fcomma, fcommc, fgrpa, fgrpc,
                                            m_MyId_all, m_MyId_comp);
    }
#endif

    ParallelDescriptor::EndTeams();
    ParallelDescriptor::EndSubCommunicator();

    ParallelDescriptor::StartTeams();
    ParallelDescriptor::StartSubCommunicator();


}

void
ParallelDescriptor::SetNProcsSidecars (int nscp)
{
    BL_ASSERT(nscp >= 0);
    BL_ASSERT(nscp < m_nProcs_all);
    int nSidecarProcs(nscp);

    Array<int> compRanksInAll(m_nProcs_all - nSidecarProcs, nProcs_undefined);
    Array<Array<int> >sidecarRanksInAll;

    for(int ip(0); ip < compRanksInAll.size(); ++ip) {
      compRanksInAll[ip] = ip;
    }

    if(nSidecarProcs > 0) {
      int whichSidecar(0);
      sidecarRanksInAll.resize(1);
      sidecarRanksInAll[whichSidecar].resize(nSidecarProcs, nProcs_undefined);
      for(int ip(0); ip < nSidecarProcs; ++ip) {
        sidecarRanksInAll[whichSidecar][ip] = m_nProcs_all - nSidecarProcs + ip;
      }
    }

    SetNProcsSidecars(compRanksInAll, sidecarRanksInAll);
}

void
ParallelDescriptor::StartSubCommunicator ()
{
    ParmParse pp("amrex");
    pp.query("ncolors", m_nCommColors);

    if (m_nCommColors > 1)
    {

#if defined(BL_USE_MPI3) || defined(BL_USE_UPCXX)
	//    m_nCommColors = 1;
	amrex::Abort("amrex.ncolors > 1 not supported for MPI3 and UPCXX");
	if (doTeamReduce())
	    amrex::Abort("amrex.ncolors > 1 not supported with team.reduce on");	    
#endif

#ifdef BL_LAZY
	amrex::Abort("amrex.ncolors > 1 not supported for LAZY=TRUE");
#endif

      /* Define variables `m_nProcs_sub', `m_MyCommSubColor',
      `m_num_procs_clr', and `m_first_procs_clr' */
      init_clr_vars();
      /* Special color for `CommComp's color */
      m_MyCommCompColor = Color(m_nCommColors); 

      BL_MPI_REQUIRE( MPI_Comm_split(Communicator(), m_MyCommSubColor.to_int(), MyProc(), &m_comm_sub) );
      BL_MPI_REQUIRE( MPI_Comm_rank(m_comm_sub, &m_MyId_sub) );
    } else {
	m_nCommColors = 1;
	m_nProcs_sub  = NProcs();
  m_num_procs_clr.resize(1,0);
  m_first_procs_clr.resize(1,0);
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
      if(m_comm_sub != MPI_COMM_NULL && m_comm_sub != Communicator()) {
	BL_MPI_REQUIRE( MPI_Comm_free(&m_comm_sub) );
        m_comm_sub = MPI_COMM_NULL;
      }
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
ParallelDescriptor::Barrier (const MPI_Comm &comm, const std::string &message)
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
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitsome, reqs, completed, indx, status, true);
    BL_MPI_REQUIRE( MPI_Waitsome(reqs.size(),
                                 reqs.dataPtr(),
                                 &completed,
                                 indx.dataPtr(),
                                 status.dataPtr()));
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitsome, reqs, completed, indx, status, false);
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

    m_MyId_all     = 0;
    m_MyId_comp    = 0;
    m_MyId_sidecar = myId_notInGroup;

    m_comm_all     = 0;
    m_comm_comp    = 0;

    m_MaxTag    = 9000;
}

void
ParallelDescriptor::StartSubCommunicator ()
{
    m_nCommColors = 1;
    m_nProcs_sub  = 1;
    m_num_procs_clr.resize(1,0);
    m_first_procs_clr.resize(1,0);
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

void ParallelDescriptor::Abort (int s, bool backtrace)
{ 
    if (backtrace) {
	BLBackTrace::handler(s);
    } else {
	std::_Exit(EXIT_FAILURE);
    }
}

int
ParallelDescriptor::Translate(int rk_clrd,Color clr)
{
    return 0;
}

const char* ParallelDescriptor::ErrorString (int) { return ""; }

void ParallelDescriptor::Barrier (const std::string &message) {}
void ParallelDescriptor::Barrier (const MPI_Comm &comm, const std::string &message) {}

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

namespace {
    static auto clock_time_begin = std::chrono::steady_clock::now();
}

double
ParallelDescriptor::second ()
{
    auto t = std::chrono::steady_clock::now();
    using ds = std::chrono::duration<double>;
    return std::chrono::duration_cast<ds>(t-clock_time_begin).count();
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
ParallelDescriptor::SeqNum (int getsetinc, int newvalue)
{
    static int seqno = m_MinTag;
    int result = seqno;

    switch(getsetinc) {
      case 0:  // ---- increment and return result
      {
        if (NColors() == 1) { 
	    ++seqno;
        } else {
	    seqno += 2;
        } 
      }
      break;
      case 1:  // ---- get current seqno
      {
        result = seqno;
      }
      break;
      case 2:  // ---- set current seqno
      {
	seqno = newvalue;
        result = seqno;
      }
      break;
      default:  // ---- error
      {
	amrex::Abort("**** Error in ParallelDescriptor::SeqNum:  bad getsetinc.");
        result = -1;
      }
    }

    if (seqno > m_MaxTag) {
	seqno = m_MinTag;
	BL_COMM_PROFILE_TAGWRAP();
    }

    return result;
}


int
ParallelDescriptor::SubSeqNum (int getsetinc, int newvalue)
{
    static int seqno = m_MinTag + 1;    
    int result = seqno;

    switch(getsetinc) {
      case 0:  // ---- increment and return result
      {
	seqno += 2;
      }
      break;
      case 1:  // ---- get current seqno
      {
        result = seqno;
      }
      break;
      case 2:  // ---- set current seqno
      {
	seqno = newvalue;
        result = seqno;
      }
      break;
      default:  // ---- error
      {
	amrex::Abort("**** Error in ParallelDescriptor::SeqNum:  bad getsetinc.");
        result = -1;
      }
    }

    if (seqno > m_MaxTag) {
	seqno = m_MinTag + 1;
	BL_COMM_PROFILE_TAGWRAP();
    }

    return result;
}


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
    enum { IO_Buffer_Size = 262144 * 8 };

#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif

    Array<Setbuf_Char_Type> io_buffer(IO_Buffer_Size);

    int fileLength(0), fileLengthPadded(0);

    std::ifstream iss;

    if (ParallelDescriptor::IOProcessor()) {
        iss.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        iss.open(filename.c_str(), std::ios::in);
        if ( ! iss.good()) {
	  if(bExitOnError) {
            amrex::FileOpenFailed(filename);
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
    if (ParallelDescriptor::IOProcessor()) {
        iss.read(charBuf.dataPtr(), fileLength);
        iss.close();
    }
    ParallelDescriptor::Bcast(charBuf.dataPtr(), fileLengthPadded,
                              ParallelDescriptor::IOProcessorNumber(), comm);
    charBuf[fileLength] = '\0';
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
	amrex::Abort("Number of processes not divisible by team size");

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


}
