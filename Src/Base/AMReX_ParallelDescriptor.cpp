
#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <sstream>
#include <stack>
#include <list>
#include <chrono>

#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLFort.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_TypeTraits.H>

#ifdef BL_USE_MPI
#include <AMReX_ccse-mpi.H>
#endif

#ifdef AMREX_PMI
#include <pmi.h>
#include <unordered_set>
#endif

#ifndef BL_AMRPROF
#include <AMReX_ParmParse.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef BL_USE_MPI
namespace
{
    static int call_mpi_finalize = 0;
}
#endif

namespace amrex {

namespace ParallelDescriptor
{
#ifdef AMREX_USE_GPU
    int use_gpu_aware_mpi = false;
#else
    int use_gpu_aware_mpi = false;
#endif

    ProcessTeam m_Team;

    MPI_Comm m_comm = MPI_COMM_NULL;    // communicator for all ranks, probably MPI_COMM_WORLD

    int m_MinTag = 1000, m_MaxTag = -1;

    const int ioProcessor = 0;

    namespace util
    {
	//
	// Reduce helper functions.
	//
	void DoAllReduceReal     (Real&      r, MPI_Op op);
	void DoAllReduceLong     (long&      r, MPI_Op op);
	void DoAllReduceInt      (int&       r, MPI_Op op);

	void DoAllReduceReal     (Real*      r, MPI_Op op, int cnt);
	void DoAllReduceLong     (long*      r, MPI_Op op, int cnt);
	void DoAllReduceInt      (int*       r, MPI_Op op, int cnt);

	void DoReduceReal     (Real&      r, MPI_Op op, int cpu);
	void DoReduceLong     (long&      r, MPI_Op op, int cpu);
	void DoReduceInt      (int&       r, MPI_Op op, int cpu);

	void DoReduceReal     (Real*      r, MPI_Op op, int cnt, int cpu);
	void DoReduceLong     (long*      r, MPI_Op op, int cnt, int cpu);
	void DoReduceInt      (int*       r, MPI_Op op, int cnt, int cpu);
    }

#ifdef AMREX_PMI
    void PMI_Initialize()
    {
      int pmi_nid;
      pmi_mesh_coord_t pmi_mesh_coord;
      int PMI_stat;
      int spawned;
      int appnum;
      int pmi_size = ParallelDescriptor::NProcs();
      int pmi_rank = ParallelDescriptor::MyProc();

      PMI_stat = PMI2_Init(&spawned, &pmi_size, &pmi_rank, &appnum);
      if (PMI_stat != PMI_SUCCESS) {
        ParallelDescriptor::Abort();
      }
      PMI_stat = PMI_Get_nid(pmi_rank, &pmi_nid);
      if (PMI_stat != PMI_SUCCESS) {
        ParallelDescriptor::Abort();
      }
      PMI_stat = PMI_Get_meshcoord(pmi_nid, &pmi_mesh_coord);
      if (PMI_stat != PMI_SUCCESS) {
        ParallelDescriptor::Abort();
      }

      // Now each MPI Process knows where it lives in the network mesh.  On
      // Aries (the interconnect on the Cray XC40), the x-coord indicates the
      // electrical group (1 group = 1 pair of adjacent cabinets); the y-coord
      // indicates the chassis (3 chassis per cabinet; 6 chassis per group);
      // and the z-coord indicates the slot (blade) within each chassis (16 per
      // chassis). Each slot contains 4 nodes, so there are at most 64 nodes
      // per chassis, 192 per cabinet, and 384 per group. (Usually there are
      // fewer than this per cabinet, since cabinets are usually a mixture of
      // compute nodes, I/O nodes, service nodes, and other things.) The slots
      // within each chassis (same x and y) are connected all-to-all, as are
      // slots in the same row across chasses within each group (same x and z).

      // One can use this information apply any kind of optimization we like,
      // e.g., splitting MPI processes into separate communicators within which
      // all processes are connected all-to-all. The following is a placeholder
      // for such an optimization, in which we merely collect the mesh
      // coordinates onto the IOProcessor and print the unique number of
      // groups, chassis, and slots occupied by the job. (This is a crude
      // measure of how fragmented the job is across the network.)

      PMI_PrintMeshcoords(&pmi_mesh_coord);
    }

    void PMI_PrintMeshcoords(const pmi_mesh_coord_t *pmi_mesh_coord) {
      unsigned short all_x_meshcoords[ParallelDescriptor::NProcs()];
      unsigned short all_y_meshcoords[ParallelDescriptor::NProcs()];
      unsigned short all_z_meshcoords[ParallelDescriptor::NProcs()];
      MPI_Gather(&pmi_mesh_coord->mesh_x,
                 1,
                 MPI_UNSIGNED_SHORT,
                 all_x_meshcoords,
                 1,
                 MPI_UNSIGNED_SHORT,
                 ParallelDescriptor::IOProcessorNumber(),
                 ParallelDescriptor::Communicator());
      MPI_Gather(&pmi_mesh_coord->mesh_y,
                 1,
                 MPI_UNSIGNED_SHORT,
                 all_y_meshcoords,
                 1,
                 MPI_UNSIGNED_SHORT,
                 ParallelDescriptor::IOProcessorNumber(),
                 ParallelDescriptor::Communicator());
      MPI_Gather(&pmi_mesh_coord->mesh_z,
                 1,
                 MPI_UNSIGNED_SHORT,
                 all_z_meshcoords,
                 1,
                 MPI_UNSIGNED_SHORT,
                 ParallelDescriptor::IOProcessorNumber(),
                 ParallelDescriptor::Communicator());

      amrex::Print() << "PMI statistics:" << std::endl;

      std::vector<unsigned short> PMI_x_meshcoord(all_x_meshcoords, all_x_meshcoords + ParallelDescriptor::NProcs());
      std::vector<unsigned short> PMI_y_meshcoord(all_y_meshcoords, all_y_meshcoords + ParallelDescriptor::NProcs());
      std::vector<unsigned short> PMI_z_meshcoord(all_z_meshcoords, all_z_meshcoords + ParallelDescriptor::NProcs());

      std::sort(PMI_x_meshcoord.begin(), PMI_x_meshcoord.end());
      std::sort(PMI_y_meshcoord.begin(), PMI_y_meshcoord.end());
      std::sort(PMI_z_meshcoord.begin(), PMI_z_meshcoord.end());

      auto last = std::unique(PMI_x_meshcoord.begin(), PMI_x_meshcoord.end());
      amrex::Print() << "# of unique groups: " << std::distance(PMI_x_meshcoord.begin(), last) << std::endl;

      last = std::unique(PMI_y_meshcoord.begin(), PMI_y_meshcoord.end());
      amrex::Print() << "# of unique groups: " << std::distance(PMI_y_meshcoord.begin(), last) << std::endl;

      last = std::unique(PMI_z_meshcoord.begin(), PMI_z_meshcoord.end());
      amrex::Print() << "# of unique groups: " << std::distance(PMI_z_meshcoord.begin(), last) << std::endl;
    }
#endif
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
    if (backtrace && amrex::system::signal_handling) {
	BLBackTrace::handler(errorcode);
    } else {
	MPI_Abort(Communicator(), errorcode);
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
    BL_PROFILE_S("ParallelDescriptor::Message::test()");
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
                                   MPI_Comm a_mpi_comm)
{
    int sflag(0);
    MPI_Initialized(&sflag);

    if ( ! sflag) {
	MPI_Init(argc, argv);
        m_comm = MPI_COMM_WORLD;
        call_mpi_finalize = 1;
    } else {
        MPI_Comm_dup(a_mpi_comm, &m_comm);
        call_mpi_finalize = 0;
    }

    ParallelContext::push(m_comm);

    // ---- find the maximum value for a tag
    int flag(0), *p;
    // For Open MPI, calling this with subcommunicators will fail.
    // So we use MPI_COMM_WORLD here.
    BL_MPI_REQUIRE( MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &p, &flag) );
    m_MaxTag = *p;
    if(!flag) {
        amrex::Abort("MPI_Comm_get_attr() failed to get MPI_TAG_UB");
    }
    BL_COMM_PROFILE_TAGRANGE(m_MinTag, m_MaxTag);

#ifdef BL_USE_MPI3
    int mpi_version, mpi_subversion;
    BL_MPI_REQUIRE( MPI_Get_version(&mpi_version, &mpi_subversion) );
    if (mpi_version < 3) amrex::Abort("MPI 3 is needed because USE_MPI3=TRUE");
#endif

    // Wait until all other processes are properly started.
//    BL_MPI_REQUIRE( MPI_Barrier(Communicator()) );
}

void
ParallelDescriptor::EndParallel ()
{
    if (!call_mpi_finalize) {
        BL_MPI_REQUIRE( MPI_Comm_free(&m_comm) );
    }
    m_comm = MPI_COMM_NULL;

    ParallelContext::pop();

    if (call_mpi_finalize) {
        MPI_Finalize();
    }
}

double
ParallelDescriptor::second () noexcept
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
ParallelDescriptor::ReduceBoolAnd (bool& r)
{
    int src = r; // src is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM);

    r = (src == ParallelDescriptor::NProcs()) ? true : false;
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
ParallelDescriptor::ReduceBoolOr (bool& r)
{
    int src = r; // src is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM);

    r = (src == 0) ? false : true;
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
ParallelDescriptor::ReduceRealSum (Real& r)
{
    util::DoAllReduceReal(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceRealSum (Real* r, int cnt)
{
    util::DoAllReduceReal(r,MPI_SUM,cnt);
}

void
ParallelDescriptor::ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar)
{
    int cnt = rvar.size();
    Vector<Real> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceReal(tmp.data(),MPI_SUM,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceRealSum (Real& r, int cpu)
{
    util::DoReduceReal(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceRealSum (Real* r, int cnt, int cpu)
{
    util::DoReduceReal(r,MPI_SUM,cnt,cpu);
}

void
ParallelDescriptor::ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<Real> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceReal(tmp.data(),MPI_SUM,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceRealMax (Real& r)
{
    util::DoAllReduceReal(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceRealMax (Real* r, int cnt)
{
    util::DoAllReduceReal(r,MPI_MAX,cnt);
}

void
ParallelDescriptor::ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar)
{
    int cnt = rvar.size();
    Vector<Real> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceReal(tmp.data(),MPI_MAX,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceRealMax (Real& r, int cpu)
{
    util::DoReduceReal(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceRealMax (Real* r, int cnt, int cpu)
{
    util::DoReduceReal(r,MPI_MAX,cnt,cpu);
}

void
ParallelDescriptor::ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<Real> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceReal(tmp.data(),MPI_MAX,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceRealMin (Real& r)
{
    util::DoAllReduceReal(r,MPI_MIN);
}

void
ParallelDescriptor::ReduceRealMin (Real* r, int cnt)
{
    util::DoAllReduceReal(r,MPI_MIN,cnt);
}

void
ParallelDescriptor::ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar)
{
    int cnt = rvar.size();
    Vector<Real> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceReal(tmp.data(),MPI_MIN,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceRealMin (Real& r, int cpu)
{
    util::DoReduceReal(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceRealMin (Real* r, int cnt, int cpu)
{
    util::DoReduceReal(r,MPI_MIN,cnt,cpu);
}

void
ParallelDescriptor::ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<Real> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceReal(tmp.data(),MPI_MIN,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceIntSum (int& r)
{
    util::DoAllReduceInt(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceIntSum (int* r, int cnt)
{
    util::DoAllReduceInt(r,MPI_SUM,cnt);
}

void
ParallelDescriptor::ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceInt(tmp.data(),MPI_SUM,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceIntSum (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceIntSum (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_SUM,cnt,cpu);
}

void
ParallelDescriptor::ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceInt(tmp.data(),MPI_SUM,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceIntMax (int& r)
{
    util::DoAllReduceInt(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceIntMax (int* r, int cnt)
{
    util::DoAllReduceInt(r,MPI_MAX,cnt);
}

void
ParallelDescriptor::ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceInt(tmp.data(),MPI_MAX,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceIntMax (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceIntMax (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_MAX,cnt,cpu);
}

void
ParallelDescriptor::ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceInt(tmp.data(),MPI_MAX,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceIntMin (int& r)
{
    util::DoAllReduceInt(r,MPI_MIN);
}

void
ParallelDescriptor::ReduceIntMin (int* r, int cnt)
{
    util::DoAllReduceInt(r,MPI_MIN,cnt);
}

void
ParallelDescriptor::ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceInt(tmp.data(),MPI_MIN,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceIntMin (int& r, int cpu)
{
    util::DoReduceInt(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceIntMin (int* r, int cnt, int cpu)
{
    util::DoReduceInt(r,MPI_MIN,cnt,cpu);
}

void
ParallelDescriptor::ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceInt(tmp.data(),MPI_MIN,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongSum (long& r)
{
    util::DoAllReduceLong(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceLongSum (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_SUM,cnt);
}

void
ParallelDescriptor::ReduceLongSum (Vector<std::reference_wrapper<long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceLong(tmp.data(),MPI_SUM,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongSum (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceLongSum (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_SUM,cnt,cpu);
}

void
ParallelDescriptor::ReduceLongSum (Vector<std::reference_wrapper<long> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceLong(tmp.data(),MPI_SUM,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongMax (long& r)
{
    util::DoAllReduceLong(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceLongMax (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_MAX,cnt);
}

void
ParallelDescriptor::ReduceLongMax (Vector<std::reference_wrapper<long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceLong(tmp.data(),MPI_MAX,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongMax (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceLongMax (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_MAX,cnt,cpu);
}

void
ParallelDescriptor::ReduceLongMax (Vector<std::reference_wrapper<long> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceLong(tmp.data(),MPI_MAX,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongMin (long& r)
{
    util::DoAllReduceLong(r,MPI_MIN);
}

void
ParallelDescriptor::ReduceLongMin (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_MIN,cnt);
}

void
ParallelDescriptor::ReduceLongMin (Vector<std::reference_wrapper<long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceLong(tmp.data(),MPI_MIN,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongMin (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceLongMin (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_MIN,cnt,cpu);
}

void
ParallelDescriptor::ReduceLongMin (Vector<std::reference_wrapper<long> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceLong(tmp.data(),MPI_MIN,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongAnd (long& r)
{
    util::DoAllReduceLong(r,MPI_LAND);
}

void
ParallelDescriptor::ReduceLongAnd (long* r, int cnt)
{
    util::DoAllReduceLong(r,MPI_LAND,cnt);
}

void
ParallelDescriptor::ReduceLongAnd (Vector<std::reference_wrapper<long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoAllReduceLong(tmp.data(),MPI_LAND,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::ReduceLongAnd (long& r, int cpu)
{
    util::DoReduceLong(r,MPI_LAND,cpu);
}

void
ParallelDescriptor::ReduceLongAnd (long* r, int cnt, int cpu)
{
    util::DoReduceLong(r,MPI_LAND,cnt,cpu);
}

void
ParallelDescriptor::ReduceLongAnd (Vector<std::reference_wrapper<long> >&& rvar,int cpu)
{
    int cnt = rvar.size();
    Vector<long> tmp{std::begin(rvar), std::end(rvar)};
    util::DoReduceLong(tmp.data(),MPI_LAND,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ParallelDescriptor::util::DoAllReduceReal (Real&  r,
                                           MPI_Op op)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceReal()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, BLProfiler::BeforeCall(), true);

    Real recv;

#if defined(BL_USE_MPI3)
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
				      Communicator()) );
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, sizeof(Real), false);
    r = recv;
}

void
ParallelDescriptor::util::DoAllReduceReal (Real*  r,
                                           MPI_Op op,
                                           int    cnt)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceReal()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, BLProfiler::BeforeCall(), true);

    BL_ASSERT(cnt > 0);

    Vector<Real> recv(cnt);

#if defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Vector<Real> recv_team(cnt);
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
				      Communicator()) );
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
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceReal()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceR, sizeof(Real), cpu);

    Real recv;

#if defined(BL_USE_MPI3)
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
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceReal()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceR, cnt * sizeof(Real), cpu);

    BL_ASSERT(cnt > 0);

    Vector<Real> recv(cnt);

#if defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Vector<Real> recv_team(cnt);
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
ParallelDescriptor::util::DoAllReduceLong (long&  r,
                                           MPI_Op op)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceLong()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceL, BLProfiler::BeforeCall(), true);

    long recv;

#if defined(BL_USE_MPI3)
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
				      Communicator()) );
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceL, sizeof(long), false);
    r = recv;
}

void
ParallelDescriptor::util::DoAllReduceLong (long*  r,
                                           MPI_Op op,
                                           int    cnt)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceLong()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceL, BLProfiler::BeforeCall(), true);

    BL_ASSERT(cnt > 0);

    Vector<long> recv(cnt);

#if defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Vector<long> recv_team(cnt);
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
				      Communicator()) );
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
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceLong()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceL, sizeof(long), cpu);

    long recv;

#if defined(BL_USE_MPI3)
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
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceLong()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceL, cnt * sizeof(long), cpu);

    BL_ASSERT(cnt > 0);

    Vector<long> recv(cnt);

#if defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Vector<long> recv_team(cnt);
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
ParallelDescriptor::util::DoAllReduceInt (int&   r,
                                          MPI_Op op)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceInt()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceI, BLProfiler::BeforeCall(), true);

    int recv;

#if defined(BL_USE_MPI3)
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
				      Communicator()));
    }
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceI, sizeof(int), false);
    r = recv;
}

void
ParallelDescriptor::util::DoAllReduceInt (int*   r,
                                          MPI_Op op,
                                          int    cnt)
{
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoAllReduceInt()");
    BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceI, BLProfiler::BeforeCall(), true);

    BL_ASSERT(cnt > 0);

    Vector<int> recv(cnt);

#if defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Vector<int> recv_team(cnt);
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
				      Communicator()));
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
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceInt()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceI, sizeof(int), cpu);

    int recv;

#if defined(BL_USE_MPI3)
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
#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::util::DoReduceInt()");
    BL_COMM_PROFILE_REDUCE(BLProfiler::ReduceI, cnt * sizeof(int), cpu);

    BL_ASSERT(cnt > 0);

    Vector<int> recv(cnt);

#if defined(BL_USE_MPI3)
    if (doTeamReduce() > 1) {
	Vector<long> recv_team(cnt);
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
ParallelDescriptor::Wait (MPI_Request& req,
                          MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Wait()");
    BL_COMM_PROFILE_WAIT(BLProfiler::Wait, req, status, true);
    BL_MPI_REQUIRE( MPI_Wait(&req, &status) );
    BL_COMM_PROFILE_WAIT(BLProfiler::Wait, req, status, false);
}

void
ParallelDescriptor::Waitall (Vector<MPI_Request>& reqs,
                             Vector<MPI_Status>& status)
{
    BL_ASSERT(status.size() >= reqs.size());

    BL_PROFILE_S("ParallelDescriptor::Waitall()");
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitall, reqs, reqs.size(), status, true);
    BL_MPI_REQUIRE( MPI_Waitall(reqs.size(),
                                reqs.dataPtr(),
                                status.dataPtr()) );
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitall, reqs, status.size(), status, false);
}

void
ParallelDescriptor::Waitany (Vector<MPI_Request>& reqs,
                             int &index,
                             MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Waitany()");
    BL_COMM_PROFILE_WAIT(BLProfiler::Waitany, reqs[0], status, true);
    BL_MPI_REQUIRE( MPI_Waitany(reqs.size(),
                                reqs.dataPtr(),
                                &index,
                                &status) );
    BL_COMM_PROFILE_WAIT(BLProfiler::Waitany, reqs[index], status, false);
}

void
ParallelDescriptor::Waitsome (Vector<MPI_Request>& reqs,
                              int&                completed,
                              Vector<int>&         indx,
                              Vector<MPI_Status>&  status)
{
    BL_ASSERT(status.size() >= reqs.size());
    BL_ASSERT(indx.size() >= reqs.size());

    BL_PROFILE_S("ParallelDescriptor::Waitsome()");
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitsome, reqs, reqs.size(), status, true);
    BL_MPI_REQUIRE( MPI_Waitsome(reqs.size(),
                                 reqs.dataPtr(),
                                 &completed,
                                 indx.dataPtr(),
                                 status.dataPtr()));
    BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitsome, reqs, indx.size(), status, false);
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
    m_comm = 0;
    m_MaxTag = 9000;
    ParallelContext::push(m_comm);
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

void ParallelDescriptor::EndParallel () 
{
    ParallelContext::pop();
}

void ParallelDescriptor::Abort (int s, bool backtrace)
{ 
    if (backtrace && amrex::system::signal_handling) {
	BLBackTrace::handler(s);
    } else {
	std::_Exit(EXIT_FAILURE);
    }
}

const char* ParallelDescriptor::ErrorString (int) { return ""; }

void ParallelDescriptor::Barrier (const std::string &message) {}
void ParallelDescriptor::Barrier (const MPI_Comm &comm, const std::string &message) {}

void ParallelDescriptor::Test (MPI_Request&, int&, MPI_Status&) {}
void ParallelDescriptor::IProbe (int, int, int&, MPI_Status&) {}
void ParallelDescriptor::IProbe (int, int, MPI_Comm, int&, MPI_Status&) {}

void ParallelDescriptor::Comm_dup (MPI_Comm, MPI_Comm&) {}

void ParallelDescriptor::ReduceRealMax (Real&) {}
void ParallelDescriptor::ReduceRealMin (Real&) {}
void ParallelDescriptor::ReduceRealSum (Real&) {}

void ParallelDescriptor::ReduceRealMax (Real&,int) {}
void ParallelDescriptor::ReduceRealMin (Real&,int) {}
void ParallelDescriptor::ReduceRealSum (Real&,int) {}

void ParallelDescriptor::ReduceRealMax (Real*,int) {}
void ParallelDescriptor::ReduceRealMin (Real*,int) {}
void ParallelDescriptor::ReduceRealSum (Real*,int) {}

void ParallelDescriptor::ReduceRealMax (Real*,int,int) {}
void ParallelDescriptor::ReduceRealMin (Real*,int,int) {}
void ParallelDescriptor::ReduceRealSum (Real*,int,int) {}

void ParallelDescriptor::ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar) {}
void ParallelDescriptor::ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar) {}
void ParallelDescriptor::ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar) {}

void ParallelDescriptor::ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar, int cpu) {}
void ParallelDescriptor::ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar, int cpu) {}
void ParallelDescriptor::ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar, int cpu) {}

void ParallelDescriptor::ReduceLongAnd (long&) {}
void ParallelDescriptor::ReduceLongSum (long&) {}
void ParallelDescriptor::ReduceLongMax (long&) {}
void ParallelDescriptor::ReduceLongMin (long&) {}

void ParallelDescriptor::ReduceLongAnd (long&,int) {}
void ParallelDescriptor::ReduceLongSum (long&,int) {}
void ParallelDescriptor::ReduceLongMax (long&,int) {}
void ParallelDescriptor::ReduceLongMin (long&,int) {}

void ParallelDescriptor::ReduceLongAnd (long*,int) {}
void ParallelDescriptor::ReduceLongSum (long*,int) {}
void ParallelDescriptor::ReduceLongMax (long*,int) {}
void ParallelDescriptor::ReduceLongMin (long*,int) {}

void ParallelDescriptor::ReduceLongAnd (long*,int,int) {}
void ParallelDescriptor::ReduceLongSum (long*,int,int) {}
void ParallelDescriptor::ReduceLongMax (long*,int,int) {}
void ParallelDescriptor::ReduceLongMin (long*,int,int) {}

void ParallelDescriptor::ReduceLongAnd (Vector<std::reference_wrapper<long> >&& rvar) {}
void ParallelDescriptor::ReduceLongSum (Vector<std::reference_wrapper<long> >&& rvar) {}
void ParallelDescriptor::ReduceLongMax (Vector<std::reference_wrapper<long> >&& rvar) {}
void ParallelDescriptor::ReduceLongMin (Vector<std::reference_wrapper<long> >&& rvar) {}

void ParallelDescriptor::ReduceLongAnd (Vector<std::reference_wrapper<long> >&& rvar, int cpu) {}
void ParallelDescriptor::ReduceLongSum (Vector<std::reference_wrapper<long> >&& rvar, int cpu) {}
void ParallelDescriptor::ReduceLongMax (Vector<std::reference_wrapper<long> >&& rvar, int cpu) {}
void ParallelDescriptor::ReduceLongMin (Vector<std::reference_wrapper<long> >&& rvar, int cpu) {}

void ParallelDescriptor::ReduceIntSum (int&) {}
void ParallelDescriptor::ReduceIntMax (int&) {}
void ParallelDescriptor::ReduceIntMin (int&) {}

void ParallelDescriptor::ReduceIntSum (int&,int) {}
void ParallelDescriptor::ReduceIntMax (int&,int) {}
void ParallelDescriptor::ReduceIntMin (int&,int) {}

void ParallelDescriptor::ReduceIntSum (int*,int) {}
void ParallelDescriptor::ReduceIntMax (int*,int) {}
void ParallelDescriptor::ReduceIntMin (int*,int) {}

void ParallelDescriptor::ReduceIntSum (int*,int,int) {}
void ParallelDescriptor::ReduceIntMax (int*,int,int) {}
void ParallelDescriptor::ReduceIntMin (int*,int,int) {}

void ParallelDescriptor::ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar) {}
void ParallelDescriptor::ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar) {}
void ParallelDescriptor::ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar) {}

void ParallelDescriptor::ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar, int cpu) {}
void ParallelDescriptor::ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar, int cpu) {}
void ParallelDescriptor::ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar, int cpu) {}

void ParallelDescriptor::ReduceBoolAnd (bool&) {}
void ParallelDescriptor::ReduceBoolOr  (bool&) {}

void ParallelDescriptor::ReduceBoolAnd (bool&,int) {}
void ParallelDescriptor::ReduceBoolOr  (bool&,int) {}

void ParallelDescriptor::Bcast(void *, int, MPI_Datatype, int, MPI_Comm) {}

double
ParallelDescriptor::second () noexcept
{
    return amrex::second();
}

void
ParallelDescriptor::Wait (MPI_Request& req,
                          MPI_Status& status)
{}

void
ParallelDescriptor::Waitall (Vector<MPI_Request>& reqs,
                             Vector<MPI_Status>& status)
{}

void
ParallelDescriptor::Waitany (Vector<MPI_Request>& reqs,
                             int &index,
                             MPI_Status& status)
{}

void
ParallelDescriptor::Waitsome (Vector<MPI_Request>& reqs,
                              int&                completed,
                              Vector<int>&         indx,
                              Vector<MPI_Status>&  status)
{}

#endif

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
    static_assert(AMREX_IS_TRIVIALLY_COPYABLE(IntVect), "IntVect must be trivially copyable");
    static_assert(std::is_standard_layout<IntVect>::value, "IntVect must be standard layout");

    static MPI_Datatype mine(MPI_DATATYPE_NULL);
    if ( mine == MPI_DATATYPE_NULL )
    {
	MPI_Datatype types[] = { MPI_INT };
	int blocklens[] = { AMREX_SPACEDIM };
	MPI_Aint disp[] = { 0 };
	BL_MPI_REQUIRE( MPI_Type_create_struct(1, blocklens, disp, types, &mine) );
        MPI_Aint lb, extent;
        BL_MPI_REQUIRE( MPI_Type_get_extent(mine, &lb, &extent) );
        if (extent != sizeof(IntVect)) {
            MPI_Datatype tmp = mine;
            BL_MPI_REQUIRE( MPI_Type_create_resized(tmp, 0, sizeof(IntVect), &mine) );
            BL_MPI_REQUIRE( MPI_Type_free(&tmp) );
        }
	BL_MPI_REQUIRE( MPI_Type_commit( &mine ) );
    }
    return mine;
}

template <> MPI_Datatype Mpi_typemap<IndexType>::type()
{
    static_assert(AMREX_IS_TRIVIALLY_COPYABLE(IndexType), "IndexType must be trivially copyable");
    static_assert(std::is_standard_layout<IndexType>::value, "IndexType must be standard layout");

    static MPI_Datatype mine(MPI_DATATYPE_NULL);
    if ( mine == MPI_DATATYPE_NULL )
    {
	MPI_Datatype types[] = { MPI_UNSIGNED };
	int blocklens[] = { 1 };
	MPI_Aint disp[] = { 0 };
        BL_MPI_REQUIRE( MPI_Type_create_struct(1, blocklens, disp, types, &mine) );
        MPI_Aint lb, extent;
        BL_MPI_REQUIRE( MPI_Type_get_extent(mine, &lb, &extent) );
        if (extent != sizeof(IndexType)) {
            MPI_Datatype tmp = mine;
            BL_MPI_REQUIRE( MPI_Type_create_resized(tmp, 0, sizeof(IndexType), &mine) );
            BL_MPI_REQUIRE( MPI_Type_free(&tmp) );
        }
	BL_MPI_REQUIRE( MPI_Type_commit( &mine ) );
    }
    return mine;
}

template <> MPI_Datatype Mpi_typemap<Box>::type()
{
    static_assert(AMREX_IS_TRIVIALLY_COPYABLE(Box), "Box must be trivially copyable");
    static_assert(std::is_standard_layout<Box>::value, "Box must be standard layout");

    static MPI_Datatype mine(MPI_DATATYPE_NULL);
    if ( mine == MPI_DATATYPE_NULL )
    {
	Box bx[2];
	MPI_Datatype types[] = {
	    Mpi_typemap<IntVect>::type(),
	    Mpi_typemap<IntVect>::type(),
	    Mpi_typemap<IndexType>::type(),
        };
	int blocklens[] = { 1, 1, 1 };
	MPI_Aint disp[3];
	BL_MPI_REQUIRE( MPI_Get_address(&bx[0].smallend, &disp[0]) );
	BL_MPI_REQUIRE( MPI_Get_address(&bx[0].bigend,   &disp[1]) );
	BL_MPI_REQUIRE( MPI_Get_address(&bx[0].btype,    &disp[2]) );
        disp[2] -= disp[0];
        disp[1] -= disp[0];
        disp[0] = 0;
        BL_MPI_REQUIRE( MPI_Type_create_struct(3, blocklens, disp, types, &mine) );
        MPI_Aint lb, extent;
        BL_MPI_REQUIRE( MPI_Type_get_extent(mine, &lb, &extent) );
        if (extent != sizeof(bx[0])) {
            MPI_Datatype tmp = mine;
            BL_MPI_REQUIRE( MPI_Type_create_resized(tmp, 0, sizeof(bx[0]), &mine) );
            BL_MPI_REQUIRE( MPI_Type_free(&tmp) );
        }
	BL_MPI_REQUIRE( MPI_Type_commit( &mine ) );
    }
    return mine;
}
#endif
}
#endif

void
ParallelDescriptor::ReadAndBcastFile (const std::string& filename,
                                      Vector<char>&       charBuf,
				      bool               bExitOnError,
				      const MPI_Comm    &comm)
{
    enum { IO_Buffer_Size = 262144 * 8 };

#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif

    Vector<Setbuf_Char_Type> io_buffer(IO_Buffer_Size);

    long fileLength(0), fileLengthPadded(0);

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
          fileLength = static_cast<std::streamoff>(iss.tellg());
          iss.seekg(0, std::ios::beg);
	}
    }
    ParallelDescriptor::Bcast(&fileLength, 1,
                              ParallelDescriptor::IOProcessorNumber(), comm);

    if(fileLength == -1) {
      return;
    }

    fileLengthPadded = fileLength + 1;
//    fileLengthPadded += fileLengthPadded % 8;
    charBuf.resize(fileLengthPadded);
    if (ParallelDescriptor::IOProcessor()) {
        iss.read(charBuf.dataPtr(), fileLength);
        iss.close();
    }
    ParallelDescriptor::Bcast(charBuf.dataPtr(), fileLengthPadded,
                              ParallelDescriptor::IOProcessorNumber(), comm);
    charBuf[fileLength] = '\0';
}

void
ParallelDescriptor::Initialize ()
{
#ifndef BL_AMRPROF
    ParmParse pp("amrex");
    pp.query("use_gpu_aware_mpi", use_gpu_aware_mpi);

    StartTeams();
#endif
}

void
ParallelDescriptor::Finalize ()
{
#ifndef BL_AMRPROF
    EndTeams();
#endif
}

#ifndef BL_AMRPROF
void
ParallelDescriptor::StartTeams ()
{
    int team_size = 1;
    int do_team_reduce = 0;

#if defined(BL_USE_MPI3)
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

#if defined(BL_USE_MPI3)
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
    }
#endif
}
#endif

void
ParallelDescriptor::EndTeams ()
{
    m_Team.clear();
}

}


using namespace amrex;

extern "C" {
    int amrex_fi_pd_myproc () {
        return ParallelDescriptor::MyProc();
    }

    int amrex_fi_pd_nprocs () {
        return ParallelDescriptor::NProcs();
    }

    int amrex_fi_pd_ioprocessor () {
        return ParallelDescriptor::IOProcessor();
    }

    int amrex_fi_pd_ioprocessor_number () {
        return ParallelDescriptor::IOProcessorNumber();
    }

    void amrex_fi_pd_bcast_r (Real* x, int n, int root)
    {
        ParallelDescriptor::Bcast(x, n, root);
    }

    Real amrex_fi_pd_wtime ()
    {
        return ParallelDescriptor::second();
    }
}
