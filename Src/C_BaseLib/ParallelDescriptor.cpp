//
// $Id: ParallelDescriptor.cpp,v 1.96 2002-10-31 18:09:00 car Exp $
//
#include <cstdio>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <ParmParse.H>

FabComTag::FabComTag ()
{
    fromProc          = 0;
    toProc            = 0;
    fabIndex          = 0;
    fineIndex         = 0;
    srcComp           = 0;
    destComp          = 0;
    nComp             = 0;
    face              = 0;
    fabArrayId        = 0;
    fillBoxId         = 0;
    procThatNeedsData = 0;
    procThatHasData   = 0;
}
//
// Definition of non-inline members of CommData.
//
namespace ParallelDescriptor
{
    //
    // My processor ID.
    //
    int m_MyId = -1;
    //
    // The number of processors.
    //
    int m_nProcs = -1;
    //
    // The number of processors in CFD part of computation.
    //
    int m_nProcsCFD = -1;
    //
    // BoxLib's Communicator
    //
    MPI_Comm m_comm;

    const int ioProcessor = 0;

    namespace util
    {
	//
	// Reduce helper functons.
	//
	void DoAllReduceReal (Real& r, MPI_Op op);
	void DoAllReduceLong (long& r, MPI_Op op);
	void DoAllReduceInt (int& r, MPI_Op op);
	void DoReduceReal (Real& r, MPI_Op op, int cpu);
	void DoReduceLong (long& r, MPI_Op op, int cpu);
	void DoReduceInt (int& r, MPI_Op op, int cpu);
	//
	// Sets number of CPUs to use in CFD portion of computation via ParmParse.
	//
	void SetNProcsCFD ();
    }
}

int
ParallelDescriptor::NProcsCFD ()
{
    if (m_nProcsCFD == -1)
        util::SetNProcsCFD();

    BL_ASSERT(m_nProcsCFD != -1);

    return m_nProcsCFD;
}

CommData::CommData ()
{
    for (int i = 0; i < length(); i++)
        m_data[i] = 0;
}

CommData::CommData (int        face,
                    int        fabindex,
                    int        fromproc,
                    int        id,
                    int        ncomp,
                    int        srccomp,
                    int        fabarrayid,
                    const Box& box)
{
    m_data[0] = face;
    m_data[1] = fabindex;
    m_data[2] = fromproc;
    m_data[3] = id;
    m_data[4] = ncomp;
    m_data[5] = srccomp;
    m_data[6] = fabarrayid;

    const IntVect& sm = box.smallEnd();

    for (int i = 0; i < BL_SPACEDIM; i++)
        m_data[7+i] = sm[i];

    const IntVect& bg = box.bigEnd();

    for (int i = 0; i < BL_SPACEDIM; i++)
        m_data[7+BL_SPACEDIM+i] = bg[i];

    IntVect typ = box.type();

    for (int i = 0; i < BL_SPACEDIM; i++)
        m_data[7+2*BL_SPACEDIM+i] = typ[i];
}

CommData::CommData (const CommData& rhs)
{
    for (int i = 0; i < length(); i++)
        m_data[i] = rhs.m_data[i];
}

CommData&
CommData::operator= (const CommData& rhs)
{
    if (!(this == &rhs))
    {
        for (int i = 0; i < length(); i++)
            m_data[i] = rhs.m_data[i];
    }
    return *this;
}

bool
CommData::operator== (const CommData& rhs) const
{
    for (int i = 0; i < length(); i++)
        if (!(m_data[i] == rhs.m_data[i]))
            return false;

    return true;
}

std::ostream&
operator<< (std::ostream&   os,
            const CommData& cd)
{
    os << cd.face()       << ' '
       << cd.fabindex()   << ' '
       << cd.fromproc()   << ' '
       << cd.id()         << ' '
       << cd.nComp()      << ' '
       << cd.srcComp()    << ' '
       << cd.fabarrayid() << ' '
       << cd.box();

    return os;
}

std::ostream&
operator<< (std::ostream&          os,
            const Array<CommData>& cd)
{
    for (int i = 0; i < cd.size(); i++)
        os << cd[i] << '\n';
    return os;
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
	//
	// Should be large enough.
	//
	const int DIM = 1024;
	static char buf[DIM];
	if ( status )
	{
	    std::sprintf(buf, "BoxLib MPI Error: File %s, line %d, %s: %s",
			 file, line, call, ParallelDescriptor::ErrorString(status));
	}
	else
	{
	    std::sprintf(buf, "BoxLib MPI Error: File %s, line %d, %s",
			 file, line, call);
	}
	buf[DIM-1] = '\0';		// Just to be safe.
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
    MPI_Abort(Communicator(), -1);
}

void
ParallelDescriptor::Abort (int errorcode)
{
    BoxLib::Abort(ErrorString(errorcode));
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

ParallelDescriptor::Message::Message ()
    :
    m_finished(true),
    m_type(MPI_DATATYPE_NULL),
    m_req(MPI_REQUEST_NULL)
{}

ParallelDescriptor::Message::Message (MPI_Request req_, MPI_Datatype type_)
    :
    m_finished(false),
    m_type(type_),
    m_req(req_)
{}

ParallelDescriptor::Message::Message (MPI_Status stat_, MPI_Datatype type_)
    :
    m_finished(true),
    m_type(type_),
    m_req(MPI_REQUEST_NULL), m_stat(stat_)
{}

void
ParallelDescriptor::Message::wait ()
{
    BL_MPI_REQUIRE( MPI_Wait(&m_req, &m_stat) );
}

bool
ParallelDescriptor::Message::test ()
{
    int flag;
    BL_MPI_REQUIRE( MPI_Test(&m_req, &flag, &m_stat) );
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

MPI_Datatype
ParallelDescriptor::Message::type () const
{
    return m_type;
}

MPI_Request
ParallelDescriptor::Message::req () const
{
    return m_req;
}

void
ParallelDescriptor::util::SetNProcsCFD ()
{
    BL_ASSERT(m_nProcs != -1);
    BL_ASSERT(m_nProcsCFD == -1);

    m_nProcsCFD = m_nProcs;

    ParmParse pp("ParallelDescriptor");

    if (pp.query("nProcsCFD",m_nProcsCFD))
    {
        if (!(m_nProcsCFD > 0 && m_nProcsCFD <= m_nProcs))
            m_nProcsCFD = m_nProcs;

        if (ParallelDescriptor::IOProcessor())
            std::cout << "--> Running job with NProcsCFD = "
                      << m_nProcsCFD
                      << std::endl;
    }
}

void
ParallelDescriptor::StartParallel (int*    argc,
                                   char*** argv)
{
    BL_ASSERT(m_MyId == -1);
    BL_ASSERT(m_nProcs == -1);
    BL_ASSERT(m_nProcsCFD == -1);

    m_comm = MPI_COMM_WORLD;

    int sflag;

    BL_MPI_REQUIRE( MPI_Initialized(&sflag) );

    if (!sflag)
	BL_MPI_REQUIRE( MPI_Init(argc, argv) );
    
    BL_MPI_REQUIRE( MPI_Comm_size(Communicator(), &m_nProcs) );

    BL_MPI_REQUIRE( MPI_Comm_rank(Communicator(), &m_MyId) );

    //
    // Wait till all other processes are properly started.
    //
    BL_MPI_REQUIRE( MPI_Barrier(Communicator()) );
}

void
ParallelDescriptor::EndParallel ()
{
    BL_ASSERT(m_MyId != -1);
    BL_ASSERT(m_nProcs != -1);

    BL_MPI_REQUIRE( MPI_Finalize() );
}

double
ParallelDescriptor::second ()
{
    return MPI_Wtime();
}

void
ParallelDescriptor::Barrier ()
{
    BL_MPI_REQUIRE( MPI_Barrier(ParallelDescriptor::Communicator()) );
}

void
ParallelDescriptor::Barrier (MPI_Comm comm)
{
    BL_MPI_REQUIRE( MPI_Barrier(comm) );
}

void
ParallelDescriptor::Test (MPI_Request& request, int& flag, MPI_Status& status)
{
    BL_MPI_REQUIRE( MPI_Test(&request,&flag,&status) );
}

void
ParallelDescriptor::Comm_dup (MPI_Comm comm, MPI_Comm& newcomm)
{
    BL_MPI_REQUIRE( MPI_Comm_dup(comm, &newcomm) );
}

void
ParallelDescriptor::util::DoAllReduceReal (Real&  r,
                                           MPI_Op op)
{
    Real recv;

    BL_MPI_REQUIRE( MPI_Allreduce(&r,
                                  &recv,
                                  1,
                                  Mpi_typemap<Real>::type(),
                                  op,
                                  Communicator()) );
    r = recv;
}

void
ParallelDescriptor::util::DoReduceReal (Real&  r,
                                        MPI_Op op,
                                        int    cpu)
{
    Real recv;

    BL_MPI_REQUIRE( MPI_Reduce(&r,
                               &recv,
                               1,
                               Mpi_typemap<Real>::type(),
                               op,
                               cpu,
                               Communicator()) );

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::ReduceRealMax (Real& r)
{
    util::DoAllReduceReal(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceRealMin (Real& r)
{
    util::DoAllReduceReal(r,MPI_MIN);
}

void
ParallelDescriptor::ReduceRealSum (Real& r)
{
    util::DoAllReduceReal(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceRealMax (Real& r, int cpu)
{
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
ParallelDescriptor::util::DoAllReduceLong (long&  r,
                                           MPI_Op op)
{
    long recv;

    BL_MPI_REQUIRE( MPI_Allreduce(&r,
                                  &recv,
                                  1,
                                  MPI_LONG,
                                  op,
                                  Communicator()) );
    r = recv;
}

void
ParallelDescriptor::util::DoReduceLong (long&  r,
                                        MPI_Op op,
                                        int    cpu)
{
    long recv;

    BL_MPI_REQUIRE( MPI_Reduce(&r,
                               &recv,
                               1,
                               MPI_LONG,
                               op,
                               cpu,
                               Communicator()));

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::ReduceLongAnd (long& r)
{
    util::DoAllReduceLong(r,MPI_LAND);
}

void
ParallelDescriptor::ReduceLongSum (long& r)
{
    util::DoAllReduceLong(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceLongMax (long& r)
{
    util::DoAllReduceLong(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceLongMin (long& r)
{
    util::DoAllReduceLong(r,MPI_MIN);
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
ParallelDescriptor::util::DoAllReduceInt (int&   r,
                                          MPI_Op op)
{
    int recv;

    BL_MPI_REQUIRE( MPI_Allreduce(&r,
                                  &recv,
                                  1,
                                  MPI_INT,
                                  op,
                                  Communicator()));
    r = recv;
}

void
ParallelDescriptor::util::DoReduceInt (int&   r,
                                       MPI_Op op,
                                       int    cpu)
{
    int recv;

    BL_MPI_REQUIRE( MPI_Reduce(&r,
                               &recv,
                               1,
                               MPI_INT,
                               op,
                               cpu,
                               Communicator()));

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::ReduceIntSum (int& r)
{
    util::DoAllReduceInt(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceIntMax (int& r)
{
    util::DoAllReduceInt(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceIntMin (int& r)
{
    util::DoAllReduceInt(r,MPI_MIN);
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
ParallelDescriptor::ReduceBoolAnd (bool& r)
{
    int src = r; // `src' is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM);

    r = (src == ParallelDescriptor::NProcs()) ? true : false;
}

void
ParallelDescriptor::ReduceBoolOr (bool& r)
{
    int src = r; // `src' is either 0 or 1.

    util::DoAllReduceInt(src,MPI_SUM);

    r = (src == 0) ? false : true;
}

void
ParallelDescriptor::ReduceBoolAnd (bool& r, int cpu)
{
    int src = r; // `src' is either 0 or 1.

    util::DoReduceInt(src,MPI_SUM,cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = (src == ParallelDescriptor::NProcs()) ? true : false;
}

void
ParallelDescriptor::ReduceBoolOr (bool& r, int cpu)
{
    int src = r; // `src' is either 0 or 1.

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
}

MPI_Op
ParallelDescriptor::Max::op ()
{
    return  MPI_MAX;
}

MPI_Op
ParallelDescriptor::Min::op ()
{
    return  MPI_MIN;
}

MPI_Op
ParallelDescriptor::Sum::op ()
{
    return  MPI_SUM;
}

MPI_Op
ParallelDescriptor::Prod::op ()
{
    return  MPI_PROD;
}

MPI_Op
ParallelDescriptor::Logical_And::op ()
{
    return  MPI_LAND;
}

MPI_Op
ParallelDescriptor::Boolean_And::op ()
{
    return  MPI_BAND;
}

MPI_Op
ParallelDescriptor::Logical_Or::op ()
{
    return  MPI_LOR;
}

MPI_Op
ParallelDescriptor::Boolean_Or::op ()
{
    return  MPI_BOR;
}

MPI_Op
ParallelDescriptor::Logical_XOr::op ()
{
    return  MPI_LXOR;
}

MPI_Op
ParallelDescriptor::Boolean_XOr::op ()
{
    return  MPI_BXOR;
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
    BL_MPI_REQUIRE( MPI_Waitsome(reqs.size(),
                                 reqs.dataPtr(),
                                 &completed,
                                 indx.dataPtr(),
                                 status.dataPtr()));
}

#else /*!BL_USE_MPI*/

void ParallelDescriptor::StartParallel (int*, char***)
{
    m_nProcs    = 1;
    m_nProcsCFD = 1;
    m_MyId      = 0;
    m_comm = 0;
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

ParallelDescriptor::Message::Message ()
    :
    m_finished(true)
{}

void
ParallelDescriptor::Message::wait ()
{}

bool
ParallelDescriptor::Message::test ()
{
    return m_finished;
}

MPI_Request
ParallelDescriptor::Message::req () const
{
    return m_req;
}

void ParallelDescriptor::util::SetNProcsCFD () {}

void ParallelDescriptor::EndParallel () {}

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

void ParallelDescriptor::Barrier () {}
void ParallelDescriptor::Barrier (MPI_Comm) {}

void ParallelDescriptor::Test (MPI_Request&, int&, MPI_Status&) {}

void ParallelDescriptor::Comm_dup (MPI_Comm, MPI_Comm&) {}

void ParallelDescriptor::ReduceRealMax (Real&) {}
void ParallelDescriptor::ReduceRealMin (Real&) {}
void ParallelDescriptor::ReduceRealSum (Real&) {}

void ParallelDescriptor::ReduceRealMax (Real&,int) {}
void ParallelDescriptor::ReduceRealMin (Real&,int) {}
void ParallelDescriptor::ReduceRealSum (Real&,int) {}

void ParallelDescriptor::ReduceLongAnd (long&) {}
void ParallelDescriptor::ReduceLongSum (long&) {}
void ParallelDescriptor::ReduceLongMax (long&) {}
void ParallelDescriptor::ReduceLongMin (long&) {}

void ParallelDescriptor::ReduceLongAnd (long&,int) {}
void ParallelDescriptor::ReduceLongSum (long&,int) {}
void ParallelDescriptor::ReduceLongMax (long&,int) {}
void ParallelDescriptor::ReduceLongMin (long&,int) {}

void ParallelDescriptor::ReduceIntSum (int&) {}
void ParallelDescriptor::ReduceIntMax (int&) {}
void ParallelDescriptor::ReduceIntMin (int&) {}

void ParallelDescriptor::ReduceIntSum (int&,int) {}
void ParallelDescriptor::ReduceIntMax (int&,int) {}
void ParallelDescriptor::ReduceIntMin (int&,int) {}

void ParallelDescriptor::ReduceBoolAnd (bool&) {}
void ParallelDescriptor::ReduceBoolOr  (bool&) {}

void ParallelDescriptor::ReduceBoolAnd (bool&,int) {}
void ParallelDescriptor::ReduceBoolOr  (bool&,int) {}
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
    const int BEG = 1000;
    const int END = 9000;

    static int seqno = BEG;

    int result = seqno++;

    if (seqno > END) seqno = BEG;

    return result;
}

#if  defined(BL_FORT_USE_UPPERCASE)
#define FORT_BL_PD_BARRIER 	BL_PD_BARRIER
#define FORT_BL_PD_COMMUNICATOR BL_PD_COMMUNICATOR
#define FORT_BL_PD_MYPROC 	BL_PD_MYPROC
#define FORT_BL_PD_NPROCS 	BL_PD_NPROCS
#define FORT_BL_PD_IOPROC 	BL_PD_IOPROC
#define FORT_BL_PD_ABORT  	BL_PD_ABORT
#elif defined(BL_FORT_USE_LOWERCASE)
#define FORT_BL_PD_BARRIER 	bl_pd_barrier
#define FORT_BL_PD_COMMUNICATOR bl_pd_communicator
#define FORT_BL_PD_MYPROC 	bl_pd_myproc
#define FORT_BL_PD_NPROCS 	bl_pd_nprocs
#define FORT_BL_PD_IOPROC 	bl_pd_ioproc
#define FORT_BL_PD_ABORT  	bl_pd_abort
#elif defined(BL_FORT_USE_UNDERSCORE)
#define FORT_BL_PD_BARRIER 	bl_pd_barrier_
#define FORT_BL_PD_COMMUNICATOR bl_pd_communicator_
#define FORT_BL_PD_MYPROC 	bl_pd_myproc_
#define FORT_BL_PD_NPROCS 	bl_pd_nprocs_
#define FORT_BL_PD_IOPROC 	bl_pd_ioproc_
#define FORT_BL_PD_ABORT  	bl_pd_abort_
#endif

extern  "C"
void
FORT_BL_PD_BARRIER ()
{
    ParallelDescriptor::Barrier();
}

extern  "C"
void
FORT_BL_PD_COMMUNICATOR (void* vcomm)
{
    MPI_Comm* comm = reinterpret_cast<MPI_Comm*>(vcomm);

    *comm = ParallelDescriptor::Communicator();
}

extern  "C"
void
FORT_BL_PD_MYPROC (int* myproc)
{
    *myproc = ParallelDescriptor::MyProc();
}

extern  "C"
void
FORT_BL_PD_NPROCS (int* nprocs)
{
    *nprocs = ParallelDescriptor::NProcs();
}

extern  "C"
void
FORT_BL_PD_IOPROC (int* ioproc)
{
    *ioproc = ParallelDescriptor::IOProcessorNumber();
}

extern  "C"
void
FORT_BL_PD_ABORT ()
{
    ParallelDescriptor::Abort();
}

#ifdef BL_USE_MPI
namespace ParallelDescriptor
{
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
}
#endif
