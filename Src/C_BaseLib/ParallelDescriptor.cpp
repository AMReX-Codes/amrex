//BL_COPYRIGHT_NOTICE

//
// $Id: ParallelDescriptor.cpp,v 1.39 1998-08-11 20:33:40 lijewski Exp $
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

#ifdef BL_USE_MPI

#include <ccse-mpi.H>
#include <RunStats.H>

static const aString MPI_Stats("mpi");

int ParallelDescriptor::m_nProcs = -1;
int ParallelDescriptor::m_MyId   = -1;

void
ParallelDescriptor::Abort (const char* msg)
{
    BoxLib::Warning(msg);

    MPI_Abort(MPI_COMM_WORLD, -1);
}

void
ParallelDescriptor::Abort (int errorcode)
{
    assert(errorcode > 0 && errorcode <= MPI_ERR_LASTCODE);

    int len = 0;

    char msg[MPI_MAX_ERROR_STRING+1];

    MPI_Error_string(errorcode, msg, &len);

    assert(len <= MPI_MAX_ERROR_STRING);

    ParallelDescriptor::Abort(msg);
}

void
ParallelDescriptor::StartParallel (int,
                                   int*    argc,
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
    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Barrier(MPI_COMM_WORLD);

    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
}

void
ParallelDescriptor::ReduceBoolAnd (bool& r)
{
    int src = r, recv; // `src' is either 0 or 1.

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&src,
                           &recv,
                           1,
                           MPI_INT,
                           MPI_SUM,
                           MPI_COMM_WORLD);

    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    r = (recv == ParallelDescriptor::NProcs()) ? true : false;
}

void
ParallelDescriptor::ReduceBoolOr  (bool& r)
{
    int src = r, recv; // `src' is either 0 or 1.

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&src,
                           &recv,
                           1,
                           MPI_INT,
                           MPI_SUM,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    r = (recv == 0) ? false : true;
}

void
ParallelDescriptor::ReduceRealMax (Real& r)
{
    Real recv;

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           mpi_data_type(&recv),
                           MPI_MAX,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);
    
    r = recv;
}

void
ParallelDescriptor::ReduceRealMin (Real& r)
{
    Real recv;

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           mpi_data_type(&recv),
                           MPI_MIN,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    r = recv;
}

void
ParallelDescriptor::ReduceRealSum (Real& r)
{
    Real recv;

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           mpi_data_type(&recv),
                           MPI_SUM,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    r = recv;
}

void
ParallelDescriptor::ReduceLongAnd (long& r)
{
    long recv;

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           MPI_LONG,
                           MPI_LAND,
                           MPI_COMM_WORLD);
    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    r = recv;
}

void
ParallelDescriptor::ReduceLongSum (long& r)
{
    long recv;

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           MPI_LONG,
                           MPI_SUM,
                           MPI_COMM_WORLD);

    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    r = recv;
}

void
ParallelDescriptor::ReduceIntSum (int& r)
{
    int recv;

    static RunStats mpi_stats(MPI_Stats);

    mpi_stats.start();

    int rc = MPI_Allreduce(&r,
                           &recv,
                           1,
                           MPI_INT,
                           MPI_SUM,
                           MPI_COMM_WORLD);

    mpi_stats.end();

    if (!(rc == MPI_SUCCESS))
        ParallelDescriptor::Abort(rc);

    r = recv;
}

void
ParallelDescriptor::Broadcast (int   fromproc,
                               void* src,
                               void* dest,
                               int   nbytes)
{
    assert(src == dest);

    static RunStats mpi_stats(MPI_Stats);

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

    static RunStats mpi_stats(MPI_Stats);

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

void ParallelDescriptor::StartParallel(int, int*, char***) {}
void ParallelDescriptor::EndParallel() {}

void ParallelDescriptor::Abort (const char* str) { BoxLib::Abort(str); }
void ParallelDescriptor::Barrier () {}

void ParallelDescriptor::ReduceBoolAnd (bool& rvar) {}
void ParallelDescriptor::ReduceBoolOr  (bool& rvar) {}

void ParallelDescriptor::ReduceRealSum (Real& rvar) {}
void ParallelDescriptor::ReduceRealMax (Real& rvar) {}
void ParallelDescriptor::ReduceRealMin (Real& rvar) {}

void ParallelDescriptor::ReduceIntSum (int& rvar) {}
void ParallelDescriptor::ReduceIntMax (int& rvar) {}
void ParallelDescriptor::ReduceIntMin (int& rvar) {}

void ParallelDescriptor::ReduceLongSum (long& rvar) {}
void ParallelDescriptor::ReduceLongMax (long& rvar) {}
void ParallelDescriptor::ReduceLongMin (long& rvar) {}
void ParallelDescriptor::ReduceLongAnd (long& rvar) {}

void ParallelDescriptor::Broadcast (int    fromproc,
                                    void*  src,
                                    void*  dest,
                                    int    nbytes) {}
//
// Here so we don't need to include <Utility.H> in <ParallelDescriptor.H>.
//
double
ParallelDescriptor::second ()
{
    return Utility::second();
}

#endif
