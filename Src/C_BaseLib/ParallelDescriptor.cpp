//BL_COPYRIGHT_NOTICE

//
// $Id: ParallelDescriptor.cpp,v 1.50 1998-11-24 01:00:54 lijewski Exp $
//

#include <Utility.H>
#include <ParallelDescriptor.H>

//
// Definition of non-inline members of CommData.
//

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

#ifdef BL_USE_MPI

#include <ccse-mpi.H>
#include <RunStats.H>

static const aString REDUCE("mpi_reduce");

int ParallelDescriptor::m_nProcs = -1;
int ParallelDescriptor::m_MyId   = -1;

void
ParallelDescriptor::Abort ()
{
#ifdef WIN32
    // if ( m_nProcs <= 1 ) throw;
    throw;
#endif
    MPI_Abort(MPI_COMM_WORLD, -1);
}

void
ParallelDescriptor::Abort (int errorcode)
{
    BoxLib::Abort(ErrorString(errorcode));
}

const char* ParallelDescriptor::ErrorString (int errorcode)
{
    assert(errorcode > 0 && errorcode <= MPI_ERR_LASTCODE);

    int len = 0;

    static char msg[MPI_MAX_ERROR_STRING+1];

    MPI_Error_string(errorcode, msg, &len);

    assert(len <= MPI_MAX_ERROR_STRING);

    return msg;
}

void
ParallelDescriptor::StartParallel (int*    argc,
                                   char*** argv)
{
    assert(m_MyId == -1);
    assert(m_nProcs == -1);

    int rc;

    if ((rc = MPI_Init(argc, argv)) != MPI_SUCCESS)
        ParallelDescriptor::Abort(rc);

    if ((rc = MPI_Comm_size(MPI_COMM_WORLD, &m_nProcs)) != MPI_SUCCESS)
        ParallelDescriptor::Abort(rc);

    if ((rc = MPI_Comm_rank(MPI_COMM_WORLD, &m_MyId)) != MPI_SUCCESS)
        ParallelDescriptor::Abort(rc);
    //
    // Now wait till all other processes are properly started.
    //
    if ((rc = MPI_Barrier(MPI_COMM_WORLD)) != MPI_SUCCESS)
        ParallelDescriptor::Abort(rc);
}

void
ParallelDescriptor::EndParallel ()
{
    assert(m_MyId != -1);
    assert(m_nProcs != -1);

    int rc = MPI_Finalize();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
}

double
ParallelDescriptor::second ()
{
    return MPI_Wtime();
}

void
ParallelDescriptor::Barrier ()
{
    static RunStats mpi_stats("mpi_barrier");

    mpi_stats.start();

    int rc = MPI_Barrier(MPI_COMM_WORLD);

    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
}

void
ParallelDescriptor::DoAllReduceReal (Real& r,
                                     int   op)
{
    Real recv;

    static RunStats mpi_stats(REDUCE);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           mpi_data_type(&recv),
                           op,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
    
    r = recv;
}

void
ParallelDescriptor::DoReduceReal (Real& r,
                                  int   op,
                                  int   cpu)
{
    Real recv;

    static RunStats mpi_stats(REDUCE);

    mpi_stats.start();

    int rc = MPI_Reduce(&r,
                        &recv,
                        1,
                        mpi_data_type(&recv),
                        op,
                        cpu,
                        MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::ReduceRealMax (Real& r)
{
    DoAllReduceReal(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceRealMin (Real& r)
{
    DoAllReduceReal(r,MPI_MIN);
}

void
ParallelDescriptor::ReduceRealSum (Real& r)
{
    DoAllReduceReal(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceRealMax (Real& r, int cpu)
{
    DoReduceReal(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceRealMin (Real& r, int cpu)
{
    DoReduceReal(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceRealSum (Real& r, int cpu)
{
    DoReduceReal(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::DoAllReduceLong (long& r,
                                     int   op)
{
    long recv;

    static RunStats mpi_stats(REDUCE);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           MPI_LONG,
                           op,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
    
    r = recv;
}

void
ParallelDescriptor::DoReduceLong (long& r,
                                  int   op,
                                  int   cpu)
{
    long recv;

    static RunStats mpi_stats(REDUCE);

    mpi_stats.start();

    int rc = MPI_Reduce(&r,
                        &recv,
                        1,
                        MPI_LONG,
                        op,
                        cpu,
                        MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::ReduceLongAnd (long& r)
{
    DoAllReduceLong(r,MPI_LAND);
}

void
ParallelDescriptor::ReduceLongSum (long& r)
{
    DoAllReduceLong(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceLongMax (long& r)
{
    DoAllReduceLong(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceLongMin (long& r)
{
    DoAllReduceLong(r,MPI_MIN);
}

void
ParallelDescriptor::ReduceLongAnd (long& r, int cpu)
{
    DoReduceLong(r,MPI_LAND,cpu);
}

void
ParallelDescriptor::ReduceLongSum (long& r, int cpu)
{
    DoReduceLong(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceLongMax (long& r, int cpu)
{
    DoReduceLong(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceLongMin (long& r, int cpu)
{
    DoReduceLong(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::DoAllReduceInt (int& r,
                                    int  op)
{
    int recv;

    static RunStats mpi_stats(REDUCE);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           MPI_INT,
                           op,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
    
    r = recv;
}

void
ParallelDescriptor::DoReduceInt (int& r,
                                 int  op,
                                 int  cpu)
{
    int recv;

    static RunStats mpi_stats(REDUCE);

    mpi_stats.start();

    int rc = MPI_Reduce(&r,
                        &recv,
                        1,
                        MPI_INT,
                        op,
                        cpu,
                        MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    if (ParallelDescriptor::MyProc() == cpu)
        r = recv;
}

void
ParallelDescriptor::ReduceIntSum (int& r)
{
    DoAllReduceInt(r,MPI_SUM);
}

void
ParallelDescriptor::ReduceIntMax (int& r)
{
    DoAllReduceInt(r,MPI_MAX);
}

void
ParallelDescriptor::ReduceIntMin (int& r)
{
    DoAllReduceInt(r,MPI_MIN);
}

void
ParallelDescriptor::ReduceIntSum (int& r, int cpu)
{
    DoReduceInt(r,MPI_SUM,cpu);
}

void
ParallelDescriptor::ReduceIntMax (int& r, int cpu)
{
    DoReduceInt(r,MPI_MAX,cpu);
}

void
ParallelDescriptor::ReduceIntMin (int& r, int cpu)
{
    DoReduceInt(r,MPI_MIN,cpu);
}

void
ParallelDescriptor::ReduceBoolAnd (bool& r)
{
    int src = r; // `src' is either 0 or 1.

    DoAllReduceInt(src,MPI_SUM);

    r = (src == ParallelDescriptor::NProcs()) ? true : false;
}

void
ParallelDescriptor::ReduceBoolOr (bool& r)
{
    int src = r; // `src' is either 0 or 1.

    DoAllReduceInt(src,MPI_SUM);

    r = (src == 0) ? false : true;
}

void
ParallelDescriptor::ReduceBoolAnd (bool& r, int cpu)
{
    int src = r; // `src' is either 0 or 1.

    DoReduceInt(src,MPI_SUM,cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = (src == ParallelDescriptor::NProcs()) ? true : false;
}

void
ParallelDescriptor::ReduceBoolOr (bool& r, int cpu)
{
    int src = r; // `src' is either 0 or 1.

    DoReduceInt(src,MPI_SUM,cpu);

    if (ParallelDescriptor::MyProc() == cpu)
        r = (src == 0) ? false : true;
}

void
ParallelDescriptor::Broadcast (int   fromproc,
                               void* src,
                               void* dest,
                               int   nbytes)
{
    assert(src == dest);

    static RunStats mpi_stats("mpi_broadcast");

    mpi_stats.start();

    MPI_Bcast(src, nbytes, MPI_BYTE, fromproc, MPI_COMM_WORLD);

    mpi_stats.end();
}

void
ParallelDescriptor::Gather (Real* sendbuf,
                            int   nsend,
                            Real* recvbuf,
                            int   root)
{
    assert(root >= 0);
    assert(nsend > 0);
    assert(!(sendbuf == 0));
    assert(!(recvbuf == 0));

    MPI_Datatype typ = mpi_data_type(sendbuf);

    static RunStats mpi_stats("mpi_gather");

    mpi_stats.start();

    int rc = MPI_Gather(sendbuf,
                        nsend,
                        typ,
                        recvbuf,
                        nsend,
                        typ,
                        root,
                        MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
}

#else

int ParallelDescriptor::m_nProcs = 1;
int ParallelDescriptor::m_MyId   = 0;

void
ParallelDescriptor::Gather (Real* sendbuf,
			    int   nsend,
			    Real* recvbuf,
			    int   root)
{
    assert(root == 0);
    assert(nsend > 0);
    assert(!(sendbuf == 0));
    assert(!(recvbuf == 0));

    for (int i = 0; i < nsend; ++i)
        recvbuf[i] = sendbuf[i];
}

void ParallelDescriptor::StartParallel (int*, char***) {}
void ParallelDescriptor::EndParallel () {}

void ParallelDescriptor::Abort () { ::abort(); }
void ParallelDescriptor::Abort (int) { ::abort(); }
const char* ParallelDescriptor::ErrorString (int) { return ""; }
void ParallelDescriptor::Barrier () {}

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

void ParallelDescriptor::Broadcast (int,void*,void*,int) {}
//
// Here so we don't need to include <Utility.H> in <ParallelDescriptor.H>.
//
double
ParallelDescriptor::second ()
{
    return Utility::second();
}

#endif
