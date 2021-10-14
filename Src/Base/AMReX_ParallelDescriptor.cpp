
#include <AMReX.H>
#include <AMReX_Utility.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_BLFort.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Arena.H>

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

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <cstdio>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stack>
#include <list>
#include <chrono>

#ifdef BL_USE_MPI
namespace
{
    static int call_mpi_finalize = 0;
    static int num_startparallel_called = 0;
    static MPI_Datatype mpi_type_intvect   = MPI_DATATYPE_NULL;
    static MPI_Datatype mpi_type_indextype = MPI_DATATYPE_NULL;
    static MPI_Datatype mpi_type_box       = MPI_DATATYPE_NULL;
    static MPI_Datatype mpi_type_lull_t    = MPI_DATATYPE_NULL;
}
#endif

namespace amrex { namespace ParallelDescriptor {

#ifdef AMREX_USE_MPI
    template <> MPI_Datatype Mpi_typemap<IntVect>::type();
    template <> MPI_Datatype Mpi_typemap<IndexType>::type();
    template <> MPI_Datatype Mpi_typemap<Box>::type();
#endif

#ifdef AMREX_USE_GPU
    int use_gpu_aware_mpi = false;
#else
    int use_gpu_aware_mpi = false;
#endif

    ProcessTeam m_Team;

    MPI_Comm m_comm = MPI_COMM_NULL;    // communicator for all ranks, probably MPI_COMM_WORLD

    int m_MinTag = 1000, m_MaxTag = -1;

    const int ioProcessor = 0;

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

#ifdef BL_USE_MPI

namespace
{
    const char*
    the_message_string (const char* file,
                        int         line,
                        const char* call,
                        int         status)
    {
        constexpr int N = 1024;
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

void
MPI_Error (const char* file, int line, const char* str, int rc)
{
    amrex::Error(the_message_string(file, line, str, rc));
}

void
Abort (int errorcode, bool backtrace)
{
    if (backtrace && amrex::system::signal_handling) {
        BLBackTrace::handler(errorcode);
    } else {
        MPI_Abort(Communicator(), errorcode);
    }
}

const char*
ErrorString (int errorcode)
{
    BL_ASSERT(errorcode > 0 && errorcode <= MPI_ERR_LASTCODE);

    int len = 0;

    static char msg[MPI_MAX_ERROR_STRING+1];

    MPI_Error_string(errorcode, msg, &len);

    BL_ASSERT(len <= MPI_MAX_ERROR_STRING);

    return msg;
}

void
Message::wait ()
{
    BL_PROFILE_S("ParallelDescriptor::Message::wait()");

    BL_COMM_PROFILE(BLProfiler::Wait, sizeof(m_type), pid(), tag());
    BL_MPI_REQUIRE( MPI_Wait(&m_req, &m_stat) );
    BL_COMM_PROFILE(BLProfiler::Wait, sizeof(m_type), BLProfiler::AfterCall(), tag());
}

bool
Message::test ()
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
Message::tag () const
{
    if ( !m_finished ) amrex::Error("Message::tag: Not Finished!");
    return m_stat.MPI_TAG;
}

int
Message::pid () const
{
    if ( !m_finished ) amrex::Error("Message::pid: Not Finished!");
    return m_stat.MPI_SOURCE;
}

size_t
Message::count () const
{
    if ( m_type == MPI_DATATYPE_NULL ) amrex::Error("Message::count: Bad Type!");
    if ( !m_finished ) amrex::Error("Message::count: Not Finished!");
    int cnt;
    BL_MPI_REQUIRE( MPI_Get_count(&m_stat, m_type, &cnt) );
    return cnt;
}

void
StartParallel (int* argc, char*** argv, MPI_Comm a_mpi_comm)
{
    int sflag(0);
    MPI_Initialized(&sflag);

    if ( ! sflag) {

#ifdef AMREX_MPI_THREAD_MULTIPLE
        int requested = MPI_THREAD_MULTIPLE;
        int provided = -1;

        MPI_Init_thread(argc, argv, requested, &provided);
#else //
        MPI_Init(argc, argv);
#endif

        m_comm = MPI_COMM_WORLD;
        call_mpi_finalize = 1;
    } else {
        MPI_Comm_dup(a_mpi_comm, &m_comm);
        call_mpi_finalize = 0;
    }

    // It seems that for some MPI implementation, the first call to MPI_Wtime is always 0.  That
    // sometimes causes problems for amrex::UniqueString function.  So we call MPI_Wtime here.
    auto tfoo = MPI_Wtime();
    amrex::ignore_unused(tfoo);

#ifdef AMREX_MPI_THREAD_MULTIPLE
    if ( ! sflag) {  // we initialized
        int requested = MPI_THREAD_MULTIPLE;
        int provided = -1;
        MPI_Query_thread(&provided);

        if (provided < requested)
        {
            auto f = ParallelDescriptor::mpi_level_to_string;
            std::cout << "MPI provided < requested: " << f(provided) << " < "
                      << f(requested) << std::endl;;
            std::abort();
        }
    }
#endif

    ParallelContext::push(m_comm);

    // Create these types outside OMP parallel region
    auto t1 = Mpi_typemap<IntVect>::type();
    auto t2 = Mpi_typemap<IndexType>::type();
    auto t3 = Mpi_typemap<Box>::type();
    auto t4 = Mpi_typemap<ParallelDescriptor::lull_t>::type();
    amrex::ignore_unused(t1,t2,t3,t4);

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

    ++num_startparallel_called;
}

void
EndParallel ()
{
    --num_startparallel_called;
    if (num_startparallel_called == 0) {
        BL_MPI_REQUIRE( MPI_Type_free(&mpi_type_intvect) );
        BL_MPI_REQUIRE( MPI_Type_free(&mpi_type_indextype) );
        BL_MPI_REQUIRE( MPI_Type_free(&mpi_type_box) );
        BL_MPI_REQUIRE( MPI_Type_free(&mpi_type_lull_t) );
        mpi_type_intvect   = MPI_DATATYPE_NULL;
        mpi_type_indextype = MPI_DATATYPE_NULL;
        mpi_type_box       = MPI_DATATYPE_NULL;
        mpi_type_lull_t    = MPI_DATATYPE_NULL;
    }

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
second () noexcept
{
    return MPI_Wtime();
}

void
Barrier (const std::string &message)
{
    amrex::ignore_unused(message);

#ifdef BL_LAZY
    Lazy::EvalReduction();
#endif

    BL_PROFILE_S("ParallelDescriptor::Barrier()");
    BL_COMM_PROFILE_BARRIER(message, true);

    BL_MPI_REQUIRE( MPI_Barrier(ParallelDescriptor::Communicator()) );

    BL_COMM_PROFILE_BARRIER(message, false);
}

void
Barrier (const MPI_Comm &comm, const std::string &message)
{
    amrex::ignore_unused(message);

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


Message
Abarrier ()
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Ibarrier(ParallelDescriptor::Communicator(), &req) );

    return Message(req, MPI_DATATYPE_NULL);
}

Message
Abarrier (const MPI_Comm & comm)
{
    MPI_Request req;
    BL_MPI_REQUIRE( MPI_Ibarrier(comm, &req) );

    return Message(req, MPI_DATATYPE_NULL);
}


void
Test (MPI_Request& request, int& flag, MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Test()");
    BL_COMM_PROFILE(BLProfiler::Test, sizeof(char), status.MPI_SOURCE, status.MPI_TAG);

    BL_MPI_REQUIRE( MPI_Test(&request,&flag,&status) );

    BL_COMM_PROFILE(BLProfiler::Test, flag, BLProfiler::AfterCall(), status.MPI_TAG);
}

void
Test (Vector<MPI_Request>& request, int& flag, Vector<MPI_Status>& status)
{
    BL_MPI_REQUIRE( MPI_Testall(request.size(), request.data(), &flag, status.data()) );
}

void
IProbe (int src_pid, int tag, int& flag, MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Iprobe()");
    BL_COMM_PROFILE(BLProfiler::Iprobe, sizeof(char), src_pid, tag);

    BL_MPI_REQUIRE( MPI_Iprobe(src_pid, tag, ParallelDescriptor::Communicator(),
                               &flag, &status) );

    BL_COMM_PROFILE(BLProfiler::Iprobe, flag, BLProfiler::AfterCall(), status.MPI_TAG);
}

void
IProbe (int src_pid, int tag, MPI_Comm comm, int& flag, MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Iprobe(comm)");
    BL_COMM_PROFILE(BLProfiler::Iprobe, sizeof(char), src_pid, tag);

    BL_MPI_REQUIRE( MPI_Iprobe(src_pid, tag, comm,
                               &flag, &status) );

    BL_COMM_PROFILE(BLProfiler::Iprobe, flag, BLProfiler::AfterCall(), status.MPI_TAG);
}

void
Comm_dup (MPI_Comm comm, MPI_Comm& newcomm)
{
    BL_PROFILE_S("ParallelDescriptor::Comm_dup()");
    BL_MPI_REQUIRE( MPI_Comm_dup(comm, &newcomm) );
}

void
ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar)
{
    ReduceRealSum<Real>(std::move(rvar));
}

void
ReduceRealSum (Vector<std::reference_wrapper<Real> >&& rvar, int cpu)
{
    ReduceRealSum<Real>(std::move(rvar), cpu);
}

void
ReduceRealMax (Vector<std::reference_wrapper<Real> > && rvar)
{
    ReduceRealMax<Real>(std::move(rvar));
}

void
ReduceRealMax (Vector<std::reference_wrapper<Real> >&& rvar, int cpu)
{
    ReduceRealMax<Real>(std::move(rvar), cpu);
}

void
ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar)
{
    ReduceRealMin<Real>(std::move(rvar));
}

void
ReduceRealMin (Vector<std::reference_wrapper<Real> >&& rvar, int cpu)
{
    ReduceRealMin<Real>(std::move(rvar), cpu);
}

void
ReduceBoolAnd (bool& r)
{
    int src = r; // src is either 0 or 1.

    detail::DoAllReduce<int>(&src,MPI_SUM,1);

    r = (src == ParallelDescriptor::NProcs()) ? true : false;
}

void
ReduceBoolAnd (bool& r, int cpu)
{
    int src = r; // src is either 0 or 1.

    detail::DoReduce<int>(&src,MPI_SUM,1,cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = (src == ParallelDescriptor::NProcs()) ? true : false;
}

void
ReduceBoolOr (bool& r)
{
    int src = r; // src is either 0 or 1.

    detail::DoAllReduce<int>(&src,MPI_SUM,1);

    r = (src == 0) ? false : true;
}

void
ReduceBoolOr (bool& r, int cpu)
{
    int src = r; // src is either 0 or 1.

    detail::DoReduce<int>(&src,MPI_SUM,1,cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = (src == 0) ? false : true;
}

void
ReduceIntSum (int& r)
{
    detail::DoAllReduce<int>(&r,MPI_SUM,1);
}

void
ReduceIntSum (int* r, int cnt)
{
    detail::DoAllReduce<int>(r,MPI_SUM,cnt);
}

void
ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoAllReduce<int>(tmp.data(),MPI_SUM,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceIntSum (int& r, int cpu)
{
    detail::DoReduce<int>(&r,MPI_SUM,1,cpu);
}

void
ReduceIntSum (int* r, int cnt, int cpu)
{
    detail::DoReduce<int>(r,MPI_SUM,cnt,cpu);
}

void
ReduceIntSum (Vector<std::reference_wrapper<int> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoReduce<int>(tmp.data(),MPI_SUM,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceIntMax (int& r)
{
    detail::DoAllReduce<int>(&r,MPI_MAX,1);
}

void
ReduceIntMax (int* r, int cnt)
{
    detail::DoAllReduce<int>(r,MPI_MAX,cnt);
}

void
ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoAllReduce<int>(tmp.data(),MPI_MAX,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceIntMax (int& r, int cpu)
{
    detail::DoReduce<int>(&r,MPI_MAX,1,cpu);
}

void
ReduceIntMax (int* r, int cnt, int cpu)
{
    detail::DoReduce<int>(r,MPI_MAX,cnt,cpu);
}

void
ReduceIntMax (Vector<std::reference_wrapper<int> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoReduce<int>(tmp.data(),MPI_MAX,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceIntMin (int& r)
{
    detail::DoAllReduce<int>(&r,MPI_MIN,1);
}

void
ReduceIntMin (int* r, int cnt)
{
    detail::DoAllReduce<int>(r,MPI_MIN,cnt);
}

void
ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoAllReduce<int>(tmp.data(),MPI_MIN,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceIntMin (int& r, int cpu)
{
    detail::DoReduce<int>(&r,MPI_MIN,1,cpu);
}

void
ReduceIntMin (int* r, int cnt, int cpu)
{
    detail::DoReduce<int>(r,MPI_MIN,cnt,cpu);
}

void
ReduceIntMin (Vector<std::reference_wrapper<int> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<int> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoReduce<int>(tmp.data(),MPI_MIN,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongSum (Long& r)
{
    detail::DoAllReduce<Long>(&r,MPI_SUM,1);
}

void
ReduceLongSum (Long* r, int cnt)
{
    detail::DoAllReduce<Long>(r,MPI_SUM,cnt);
}

void
ReduceLongSum (Vector<std::reference_wrapper<Long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoAllReduce<Long>(tmp.data(),MPI_SUM,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongSum (Long& r, int cpu)
{
    detail::DoReduce<Long>(&r,MPI_SUM,1,cpu);
}

void
ReduceLongSum (Long* r, int cnt, int cpu)
{
    detail::DoReduce<Long>(r,MPI_SUM,cnt,cpu);
}

void
ReduceLongSum (Vector<std::reference_wrapper<Long> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoReduce<Long>(tmp.data(),MPI_SUM,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongMax (Long& r)
{
    detail::DoAllReduce<Long>(&r,MPI_MAX,1);
}

void
ReduceLongMax (Long* r, int cnt)
{
    detail::DoAllReduce<Long>(r,MPI_MAX,cnt);
}

void
ReduceLongMax (Vector<std::reference_wrapper<Long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoAllReduce<Long>(tmp.data(),MPI_MAX,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongMax (Long& r, int cpu)
{
    detail::DoReduce<Long>(&r,MPI_MAX,1,cpu);
}

void
ReduceLongMax (Long* r, int cnt, int cpu)
{
    detail::DoReduce<Long>(r,MPI_MAX,cnt,cpu);
}

void
ReduceLongMax (Vector<std::reference_wrapper<Long> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoReduce<Long>(tmp.data(),MPI_MAX,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongMin (Long& r)
{
    detail::DoAllReduce<Long>(&r,MPI_MIN,1);
}

void
ReduceLongMin (Long* r, int cnt)
{
    detail::DoAllReduce<Long>(r,MPI_MIN,cnt);
}

void
ReduceLongMin (Vector<std::reference_wrapper<Long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoAllReduce<Long>(tmp.data(),MPI_MIN,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongMin (Long& r, int cpu)
{
    detail::DoReduce<Long>(&r,MPI_MIN,1,cpu);
}

void
ReduceLongMin (Long* r, int cnt, int cpu)
{
    detail::DoReduce<Long>(r,MPI_MIN,cnt,cpu);
}

void
ReduceLongMin (Vector<std::reference_wrapper<Long> >&& rvar, int cpu)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoReduce<Long>(tmp.data(),MPI_MIN,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongAnd (Long& r)
{
    detail::DoAllReduce<Long>(&r,MPI_LAND,1);
}

void
ReduceLongAnd (Long* r, int cnt)
{
    detail::DoAllReduce<Long>(r,MPI_LAND,cnt);
}

void
ReduceLongAnd (Vector<std::reference_wrapper<Long> >&& rvar)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoAllReduce<Long>(tmp.data(),MPI_LAND,cnt);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
ReduceLongAnd (Long& r, int cpu)
{
    detail::DoReduce<Long>(&r,MPI_LAND,1,cpu);
}

void
ReduceLongAnd (Long* r, int cnt, int cpu)
{
    detail::DoReduce<Long>(r,MPI_LAND,cnt,cpu);
}

void
ReduceLongAnd (Vector<std::reference_wrapper<Long> >&& rvar,int cpu)
{
    int cnt = rvar.size();
    Vector<Long> tmp{std::begin(rvar), std::end(rvar)};
    detail::DoReduce<Long>(tmp.data(),MPI_LAND,cnt,cpu);
    for (int i = 0; i < cnt; ++i) {
        rvar[i].get() = tmp[i];
    }
}

void
Gather (Real* sendbuf, int nsend, Real* recvbuf, int root)
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
Mpi_typemap<char>::type ()
{
    return  MPI_CHAR;
}

template <>
MPI_Datatype
Mpi_typemap<short>::type ()
{
    return  MPI_SHORT;
}

template <>
MPI_Datatype
Mpi_typemap<int>::type ()
{
    return  MPI_INT;
}

template <>
MPI_Datatype
Mpi_typemap<long>::type ()
{
    return  MPI_LONG;
}

template <>
MPI_Datatype
Mpi_typemap<long long>::type ()
{
    return  MPI_LONG_LONG;
}

template <>
MPI_Datatype
Mpi_typemap<unsigned char>::type ()
{
    return  MPI_UNSIGNED_CHAR;
}

template <>
MPI_Datatype
Mpi_typemap<unsigned short>::type ()
{
    return  MPI_UNSIGNED_SHORT;
}

template <>
MPI_Datatype
Mpi_typemap<unsigned int>::type ()
{
    return  MPI_UNSIGNED;
}

template <>
MPI_Datatype
Mpi_typemap<unsigned long>::type ()
{
    return  MPI_UNSIGNED_LONG;
}

template <>
MPI_Datatype
Mpi_typemap<unsigned long long>::type ()
{
    return  MPI_UNSIGNED_LONG_LONG;
}

template <>
MPI_Datatype
Mpi_typemap<float>::type ()
{
    return  MPI_FLOAT;
}

template <>
MPI_Datatype
Mpi_typemap<double>::type ()
{
    return  MPI_DOUBLE;
}

template <>
MPI_Datatype
Mpi_typemap<ParallelDescriptor::lull_t>::type ()
{
    if (mpi_type_lull_t == MPI_DATATYPE_NULL)
    {
        BL_MPI_REQUIRE( MPI_Type_contiguous(sizeof(lull_t), MPI_CHAR, &mpi_type_lull_t) );
        BL_MPI_REQUIRE( MPI_Type_commit(&mpi_type_lull_t) );
    }
    return mpi_type_lull_t;
}

void
Wait (MPI_Request& req, MPI_Status& status)
{
    BL_PROFILE_S("ParallelDescriptor::Wait()");
    BL_COMM_PROFILE_WAIT(BLProfiler::Wait, req, status, true);
    BL_MPI_REQUIRE( MPI_Wait(&req, &status) );
    BL_COMM_PROFILE_WAIT(BLProfiler::Wait, req, status, false);
}

void
Waitall (Vector<MPI_Request>& reqs, Vector<MPI_Status>& status)
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
Waitany (Vector<MPI_Request>& reqs, int &index, MPI_Status& status)
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
Waitsome (Vector<MPI_Request>& reqs, int& completed,
          Vector<int>& indx, Vector<MPI_Status>& status)
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
Bcast(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
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
StartParallel (int* /*argc*/, char*** /*argv*/, MPI_Comm)
{
    m_comm = 0;
    m_MaxTag = 9000;
    ParallelContext::push(m_comm);
}

void
Gather (Real* sendbuf, int nsend, Real* recvbuf, int root)
{
    amrex::ignore_unused(root);
    BL_ASSERT(root == 0);
    BL_ASSERT(nsend > 0);
    BL_ASSERT(!(sendbuf == 0));
    BL_ASSERT(!(recvbuf == 0));

    for (int i = 0; i < nsend; ++i)
        recvbuf[i] = sendbuf[i];
}

void
Message::wait ()
{}

bool
Message::test ()
{
    return m_finished;
}

void
EndParallel ()
{
    ParallelContext::pop();
}

void
Abort (int s, bool backtrace)
{
    if (backtrace && amrex::system::signal_handling) {
        BLBackTrace::handler(s);
    } else {
        std::_Exit(EXIT_FAILURE);
    }
}

const char* ErrorString (int) { return ""; }

void Barrier (const std::string &/*message*/) {}
void Barrier (const MPI_Comm &/*comm*/, const std::string &/*message*/) {}
Message Abarrier () { return Message(); }
Message Abarrier (const MPI_Comm &/*comm*/) { return Message(); }

void Test (MPI_Request&, int&, MPI_Status&) {}
void Test (Vector<MPI_Request>&, int&, Vector<MPI_Status>&) {}
void IProbe (int, int, int&, MPI_Status&) {}
void IProbe (int, int, MPI_Comm, int&, MPI_Status&) {}

void Comm_dup (MPI_Comm, MPI_Comm&) {}

void ReduceRealSum (Vector<std::reference_wrapper<Real> >&& /*rvar*/) {}
void ReduceRealMax (Vector<std::reference_wrapper<Real> >&& /*rvar*/) {}
void ReduceRealMin (Vector<std::reference_wrapper<Real> >&& /*rvar*/) {}

void ReduceRealSum (Vector<std::reference_wrapper<Real> >&& /*rvar*/, int /*cpu*/) {}
void ReduceRealMax (Vector<std::reference_wrapper<Real> >&& /*rvar*/, int /*cpu*/) {}
void ReduceRealMin (Vector<std::reference_wrapper<Real> >&& /*rvar*/, int /*cpu*/) {}

void ReduceLongAnd (Long&) {}
void ReduceLongSum (Long&) {}
void ReduceLongMax (Long&) {}
void ReduceLongMin (Long&) {}

void ReduceLongAnd (Long&,int) {}
void ReduceLongSum (Long&,int) {}
void ReduceLongMax (Long&,int) {}
void ReduceLongMin (Long&,int) {}

void ReduceLongAnd (Long*,int) {}
void ReduceLongSum (Long*,int) {}
void ReduceLongMax (Long*,int) {}
void ReduceLongMin (Long*,int) {}

void ReduceLongAnd (Long*,int,int) {}
void ReduceLongSum (Long*,int,int) {}
void ReduceLongMax (Long*,int,int) {}
void ReduceLongMin (Long*,int,int) {}

void ReduceLongAnd (Vector<std::reference_wrapper<Long> >&& /*rvar*/) {}
void ReduceLongSum (Vector<std::reference_wrapper<Long> >&& /*rvar*/) {}
void ReduceLongMax (Vector<std::reference_wrapper<Long> >&& /*rvar*/) {}
void ReduceLongMin (Vector<std::reference_wrapper<Long> >&& /*rvar*/) {}

void ReduceLongAnd (Vector<std::reference_wrapper<Long> >&& /*rvar*/, int /*cpu*/) {}
void ReduceLongSum (Vector<std::reference_wrapper<Long> >&& /*rvar*/, int /*cpu*/) {}
void ReduceLongMax (Vector<std::reference_wrapper<Long> >&& /*rvar*/, int /*cpu*/) {}
void ReduceLongMin (Vector<std::reference_wrapper<Long> >&& /*rvar*/, int /*cpu*/) {}

void ReduceIntSum (int&) {}
void ReduceIntMax (int&) {}
void ReduceIntMin (int&) {}

void ReduceIntSum (int&,int) {}
void ReduceIntMax (int&,int) {}
void ReduceIntMin (int&,int) {}

void ReduceIntSum (int*,int) {}
void ReduceIntMax (int*,int) {}
void ReduceIntMin (int*,int) {}

void ReduceIntSum (int*,int,int) {}
void ReduceIntMax (int*,int,int) {}
void ReduceIntMin (int*,int,int) {}

void ReduceIntSum (Vector<std::reference_wrapper<int> >&& /*rvar*/) {}
void ReduceIntMax (Vector<std::reference_wrapper<int> >&& /*rvar*/) {}
void ReduceIntMin (Vector<std::reference_wrapper<int> >&& /*rvar*/) {}

void ReduceIntSum (Vector<std::reference_wrapper<int> >&& /*rvar*/, int /*cpu*/) {}
void ReduceIntMax (Vector<std::reference_wrapper<int> >&& /*rvar*/, int /*cpu*/) {}
void ReduceIntMin (Vector<std::reference_wrapper<int> >&& /*rvar*/, int /*cpu*/) {}

void ReduceBoolAnd (bool&) {}
void ReduceBoolOr  (bool&) {}

void ReduceBoolAnd (bool&,int) {}
void ReduceBoolOr  (bool&,int) {}

void Bcast(void *, int, MPI_Datatype, int, MPI_Comm) {}

double
second () noexcept
{
    return amrex::second();
}

void
Wait (MPI_Request& /*req*/, MPI_Status& /*status*/)
{}

void
Waitall (Vector<MPI_Request>& /*reqs*/, Vector<MPI_Status>& /*status*/)
{}

void
Waitany (Vector<MPI_Request>& /*reqs*/, int &/*index*/, MPI_Status& /*status*/)
{}

void
Waitsome (Vector<MPI_Request>& /*reqs*/, int& /*completed*/,
          Vector<int>& /*indx*/, Vector<MPI_Status>& /*status*/)
{}

#endif

#ifndef BL_NO_FORT

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

#endif

#if defined(BL_USE_MPI) && !defined(BL_AMRPROF)
template <> MPI_Datatype Mpi_typemap<IntVect>::type()
{
    static_assert(std::is_trivially_copyable<IntVect>::value, "IntVect must be trivially copyable");
    static_assert(std::is_standard_layout<IntVect>::value, "IntVect must be standard layout");

    if ( mpi_type_intvect == MPI_DATATYPE_NULL )
    {
        MPI_Datatype types[] = { MPI_INT };
        int blocklens[] = { AMREX_SPACEDIM };
        MPI_Aint disp[] = { 0 };
        BL_MPI_REQUIRE( MPI_Type_create_struct(1, blocklens, disp, types, &mpi_type_intvect) );
        MPI_Aint lb, extent;
        BL_MPI_REQUIRE( MPI_Type_get_extent(mpi_type_intvect, &lb, &extent) );
        if (extent != sizeof(IntVect)) {
            MPI_Datatype tmp = mpi_type_intvect;
            BL_MPI_REQUIRE( MPI_Type_create_resized(tmp, 0, sizeof(IntVect), &mpi_type_intvect) );
            BL_MPI_REQUIRE( MPI_Type_free(&tmp) );
        }
        BL_MPI_REQUIRE( MPI_Type_commit( &mpi_type_intvect ) );
    }
    return mpi_type_intvect;
}

template <> MPI_Datatype Mpi_typemap<IndexType>::type()
{
    static_assert(std::is_trivially_copyable<IndexType>::value, "IndexType must be trivially copyable");
    static_assert(std::is_standard_layout<IndexType>::value, "IndexType must be standard layout");

    if ( mpi_type_indextype == MPI_DATATYPE_NULL )
    {
        MPI_Datatype types[] = { MPI_UNSIGNED };
        int blocklens[] = { 1 };
        MPI_Aint disp[] = { 0 };
        BL_MPI_REQUIRE( MPI_Type_create_struct(1, blocklens, disp, types, &mpi_type_indextype) );
        MPI_Aint lb, extent;
        BL_MPI_REQUIRE( MPI_Type_get_extent(mpi_type_indextype, &lb, &extent) );
        if (extent != sizeof(IndexType)) {
            MPI_Datatype tmp = mpi_type_indextype;
            BL_MPI_REQUIRE( MPI_Type_create_resized(tmp, 0, sizeof(IndexType), &mpi_type_indextype) );
            BL_MPI_REQUIRE( MPI_Type_free(&tmp) );
        }
        BL_MPI_REQUIRE( MPI_Type_commit( &mpi_type_indextype ) );
    }
    return mpi_type_indextype;
}

template <> MPI_Datatype Mpi_typemap<Box>::type()
{
    static_assert(std::is_trivially_copyable<Box>::value, "Box must be trivially copyable");
    static_assert(std::is_standard_layout<Box>::value, "Box must be standard layout");

    if ( mpi_type_box == MPI_DATATYPE_NULL )
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
        BL_MPI_REQUIRE( MPI_Type_create_struct(3, blocklens, disp, types, &mpi_type_box) );
        MPI_Aint lb, extent;
        BL_MPI_REQUIRE( MPI_Type_get_extent(mpi_type_box, &lb, &extent) );
        if (extent != sizeof(bx[0])) {
            MPI_Datatype tmp = mpi_type_box;
            BL_MPI_REQUIRE( MPI_Type_create_resized(tmp, 0, sizeof(bx[0]), &mpi_type_box) );
            BL_MPI_REQUIRE( MPI_Type_free(&tmp) );
        }
        BL_MPI_REQUIRE( MPI_Type_commit( &mpi_type_box ) );
    }
    return mpi_type_box;
}
#endif

void
ReadAndBcastFile (const std::string& filename, Vector<char>& charBuf,
                  bool bExitOnError, const MPI_Comm&comm)
{
    enum { IO_Buffer_Size = 262144 * 8 };

#ifdef BL_SETBUF_SIGNED_CHAR
    typedef signed char Setbuf_Char_Type;
#else
    typedef char Setbuf_Char_Type;
#endif

    Vector<Setbuf_Char_Type> io_buffer(IO_Buffer_Size);

    Long fileLength(0), fileLengthPadded(0);

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
Initialize ()
{
#ifndef BL_AMRPROF
    ParmParse pp("amrex");
    pp.query("use_gpu_aware_mpi", use_gpu_aware_mpi);

    StartTeams();
#endif
}

void
Finalize ()
{
#ifndef BL_AMRPROF
    EndTeams();
#endif
}

#ifndef BL_AMRPROF
void
StartTeams ()
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
EndTeams ()
{
    m_Team.clear();
}

std::string
mpi_level_to_string (int mtlev)
{
    switch (mtlev) {
#ifdef AMREX_USE_MPI
        case MPI_THREAD_SINGLE:
            return std::string("MPI_THREAD_SINGLE");
        case MPI_THREAD_FUNNELED:
            return std::string("MPI_THREAD_FUNNELED");
        case MPI_THREAD_SERIALIZED:
            return std::string("MPI_THREAD_SERIALIZED");
        case MPI_THREAD_MULTIPLE:
            return std::string("MPI_THREAD_MULTIPLE");
#endif
        default:
            return std::string("UNKNOWN");
    }
}

#ifdef BL_USE_MPI

int
select_comm_data_type (std::size_t nbytes)
{
    if (nbytes <= std::size_t(std::numeric_limits<int>::max()))
    {
        return 1;
    }
    else if (amrex::aligned_size(sizeof(unsigned long long), nbytes) <=
             sizeof(unsigned long long)*std::size_t(std::numeric_limits<int>::max()))
    {
        return 2;
    } else if (amrex::aligned_size(sizeof(ParallelDescriptor::lull_t), nbytes) <=
               sizeof(ParallelDescriptor::lull_t)*std::size_t(std::numeric_limits<int>::max()))
    {
        return 3;
    }
    else
    {
        return 0;
    }
}

std::size_t
alignof_comm_data (std::size_t nbytes)
{
    const int t = select_comm_data_type(nbytes);
    if (t == 1) {
        return 1;
    } else if (t == 2) {
        return sizeof(unsigned long long);
    } else if (t == 3) {
        return sizeof(ParallelDescriptor::lull_t);
    } else {
        amrex::Abort("TODO: message size is too big");
        return 0;
    }
}

template <>
Message
Asend<char> (const char* buf, size_t n, int pid, int tag, MPI_Comm comm)
{
    BL_PROFILE_T_S("ParallelDescriptor::Asend(TsiiM)", char);
    BL_COMM_PROFILE(BLProfiler::AsendTsiiM, n * sizeof(char), pid, tag);

    MPI_Request req;
    Message msg;
    const int comm_data_type = ParallelDescriptor::select_comm_data_type(n);
    if (comm_data_type == 1) {
        BL_MPI_REQUIRE( MPI_Isend(const_cast<char*>(buf),
                                  n,
                                  Mpi_typemap<char>::type(),
                                  pid, tag, comm, &req) );
        msg = Message(req, Mpi_typemap<char>::type());
    } else if (comm_data_type == 2) {
        if (!amrex::is_aligned(buf, alignof(unsigned long long))
            || (n % sizeof(unsigned long long)) != 0) {
            amrex::Abort("Message size is too big as char, and it cannot be sent as unsigned long long.");
        }
        BL_MPI_REQUIRE( MPI_Isend(const_cast<unsigned long long*>
                                     (reinterpret_cast<unsigned long long const*>(buf)),
                                  n/sizeof(unsigned long long),
                                  Mpi_typemap<unsigned long long>::type(),
                                  pid, tag, comm, &req) );
        msg = Message(req, Mpi_typemap<unsigned long long>::type());
    } else if (comm_data_type == 3) {
        if (!amrex::is_aligned(buf, alignof(ParallelDescriptor::lull_t))
            || (n % sizeof(ParallelDescriptor::lull_t)) != 0) {
            amrex::Abort("Message size is too big as char or unsigned long long, and it cannot be sent as ParallelDescriptor::lull_t");
        }
        BL_MPI_REQUIRE( MPI_Isend(const_cast<ParallelDescriptor::lull_t*>
                                     (reinterpret_cast<ParallelDescriptor::lull_t const*>(buf)),
                                  n/sizeof(ParallelDescriptor::lull_t),
                                  Mpi_typemap<ParallelDescriptor::lull_t>::type(),
                                  pid, tag, comm, &req) );
        msg = Message(req, Mpi_typemap<ParallelDescriptor::lull_t>::type());
    } else {
        amrex::Abort("TODO: message size is too big");
    }

    BL_COMM_PROFILE(BLProfiler::AsendTsiiM, BLProfiler::AfterCall(), pid, tag);
    return msg;
}

template <>
Message
Send<char> (const char* buf, size_t n, int pid, int tag, MPI_Comm comm)
{
    BL_PROFILE_T_S("ParallelDescriptor::Send(Tsii)", char);
    BL_COMM_PROFILE(BLProfiler::SendTsii, n * sizeof(char), pid, tag);

    const int comm_data_type = ParallelDescriptor::select_comm_data_type(n);
    if (comm_data_type == 1) {
        BL_MPI_REQUIRE( MPI_Send(const_cast<char*>(buf),
                                 n,
                                 Mpi_typemap<char>::type(),
                                 pid, tag, comm) );
    } else if (comm_data_type == 2) {
        if (!amrex::is_aligned(buf, alignof(unsigned long long))
            || (n % sizeof(unsigned long long)) != 0) {
            amrex::Abort("Message size is too big as char, and it cannot be sent as unsigned long long.");
        }
        BL_MPI_REQUIRE( MPI_Send(const_cast<unsigned long long*>
                                     (reinterpret_cast<unsigned long long const*>(buf)),
                                 n/sizeof(unsigned long long),
                                 Mpi_typemap<unsigned long long>::type(),
                                 pid, tag, comm) );
    } else if (comm_data_type == 3) {
        if (!amrex::is_aligned(buf, alignof(ParallelDescriptor::lull_t))
            || (n % sizeof(ParallelDescriptor::lull_t)) != 0) {
            amrex::Abort("Message size is too big as char or unsigned long long, and it cannot be sent as ParallelDescriptor::lull_t");
        }
        BL_MPI_REQUIRE( MPI_Send(const_cast<ParallelDescriptor::lull_t*>
                                     (reinterpret_cast<ParallelDescriptor::lull_t const*>(buf)),
                                 n/sizeof(ParallelDescriptor::lull_t),
                                 Mpi_typemap<ParallelDescriptor::lull_t>::type(),
                                 pid, tag, comm) );
    } else {
        amrex::Abort("TODO: message size is too big");
    }

    BL_COMM_PROFILE(BLProfiler::SendTsii, BLProfiler::AfterCall(), pid, tag);
    return Message();
}

template <>
Message
Arecv<char> (char* buf, size_t n, int pid, int tag, MPI_Comm comm)
{
    BL_PROFILE_T_S("ParallelDescriptor::Arecv(TsiiM)", char);
    BL_COMM_PROFILE(BLProfiler::ArecvTsiiM, n * sizeof(char), pid, tag);

    MPI_Request req;
    Message msg;
    const int comm_data_type = ParallelDescriptor::select_comm_data_type(n);
    if (comm_data_type == 1) {
        BL_MPI_REQUIRE( MPI_Irecv(buf,
                                  n,
                                  Mpi_typemap<char>::type(),
                                  pid, tag, comm, &req) );
        msg = Message(req, Mpi_typemap<char>::type());
    } else if (comm_data_type == 2) {
        if (!amrex::is_aligned(buf, alignof(unsigned long long))
            || (n % sizeof(unsigned long long)) != 0) {
            amrex::Abort("Message size is too big as char, and it cannot be received as unsigned long long.");
        }
        BL_MPI_REQUIRE( MPI_Irecv((unsigned long long *)buf,
                                  n/sizeof(unsigned long long),
                                  Mpi_typemap<unsigned long long>::type(),
                                  pid, tag, comm, &req) );
        msg = Message(req, Mpi_typemap<unsigned long long>::type());
    } else if (comm_data_type == 3) {
        if (!amrex::is_aligned(buf, alignof(ParallelDescriptor::lull_t))
            || (n % sizeof(ParallelDescriptor::lull_t)) != 0) {
            amrex::Abort("Message size is too big as char or unsigned long long, and it cannot be received as ParallelDescriptor::lull_t");
        }
        BL_MPI_REQUIRE( MPI_Irecv((ParallelDescriptor::lull_t *)buf,
                                  n/sizeof(ParallelDescriptor::lull_t),
                                  Mpi_typemap<ParallelDescriptor::lull_t>::type(),
                                  pid, tag, comm, &req) );
        msg = Message(req, Mpi_typemap<ParallelDescriptor::lull_t>::type());
    } else {
        amrex::Abort("Message size is too big");
    }

    BL_COMM_PROFILE(BLProfiler::ArecvTsiiM, BLProfiler::AfterCall(), pid, tag);
    return msg;
}

template <>
Message
Recv<char> (char* buf, size_t n, int pid, int tag, MPI_Comm comm)
{
    BL_PROFILE_T_S("ParallelDescriptor::Recv(Tsii)", char);
    BL_COMM_PROFILE(BLProfiler::RecvTsii, BLProfiler::BeforeCall(), pid, tag);

    MPI_Status stat;
    Message msg;
    const int comm_data_type = ParallelDescriptor::select_comm_data_type(n);
    if (comm_data_type == 1) {
        BL_MPI_REQUIRE( MPI_Recv(buf,
                                 n,
                                 Mpi_typemap<char>::type(),
                                 pid, tag, comm, &stat) );
        msg = Message(stat, Mpi_typemap<char>::type());
    } else if (comm_data_type == 2) {
        if (!amrex::is_aligned(buf, alignof(unsigned long long))
            || (n % sizeof(unsigned long long)) != 0) {
            amrex::Abort("Message size is too big as char, and it cannot be received as unsigned long long.");
        }
        BL_MPI_REQUIRE( MPI_Recv((unsigned long long *)buf,
                                 n/sizeof(unsigned long long),
                                 Mpi_typemap<unsigned long long>::type(),
                                 pid, tag, comm, &stat) );
        msg = Message(stat, Mpi_typemap<unsigned long long>::type());
    } else if (comm_data_type == 3) {
        if (!amrex::is_aligned(buf, alignof(ParallelDescriptor::lull_t))
            || (n % sizeof(ParallelDescriptor::lull_t)) != 0) {
            amrex::Abort("Message size is too big as char or unsigned long long, and it cannot be received as ParallelDescriptor::lull_t");
        }
        BL_MPI_REQUIRE( MPI_Recv((ParallelDescriptor::lull_t *)buf,
                                 n/sizeof(ParallelDescriptor::lull_t),
                                 Mpi_typemap<ParallelDescriptor::lull_t>::type(),
                                 pid, tag, comm, &stat) );
        msg = Message(stat, Mpi_typemap<ParallelDescriptor::lull_t>::type());
    } else {
        amrex::Abort("Message size is too big");
    }

    BL_COMM_PROFILE(BLProfiler::RecvTsii, n * sizeof(char), pid, stat.MPI_TAG);
    return msg;
}

#endif

}}


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
        return static_cast<Real>(ParallelDescriptor::second());
    }
}
