//BL_COPYRIGHT_NOTICE

//
// $Id: ParallelDescriptor.cpp,v 1.23 1998-04-09 22:37:27 car Exp $
//

#include <Utility.H>
#include <ParallelDescriptor.H>

#if defined(BL_USE_BSP)

#include "bsp.h"

//
// Type of function pointer required by bsp_fold().
//
typedef void (*VFVVVI)(void*,void*,void*,int*);

#ifdef FIXBSPLIBLEVEL1HEADER
extern "C"
{
  extern void bsp_fold (VFVVVI,void*,void*,int);
  extern void bsp_fold_cpp (VFVVVI,void*,void*,int,int,char*);
}
#endif /*FIXBSPLIBLEVEL1HEADER*/

#include "bsp_level1.h"

void
ParallelDescriptor::StartParallel (int nprocs, int*, char***)
{
    bsp_begin(nprocs);
}

void
ParallelDescriptor::EndParallel ()
{
    bsp_end();
}

int
ParallelDescriptor::MyProc ()
{ 
    return bsp_pid();                
}

int
ParallelDescriptor::NProcs ()
{ 
    return bsp_nprocs();             
}

bool
ParallelDescriptor::IOProcessor ()
{ 
    return bsp_pid() == ioProcessor; 
}

int
ParallelDescriptor::IOProcessorNumber ()
{ 
    return ioProcessor;              
}

void
ParallelDescriptor::Abort (const char* msg)  
{ 
    bsp_abort((char*)msg);           
}

double
ParallelDescriptor::second ()              
{ 
    return bsp_time();               
}

void
ParallelDescriptor::ShareVar (const void* var,
                              int         bytes)
{
    bsp_pushregister(var, bytes);
}

void
ParallelDescriptor::UnshareVar (const void* var)
{
    bsp_popregister(var);
}

void
ParallelDescriptor::WriteData (int         procnum,
                               const void* src,
                               void*       dest,
                               int         offset,
                               int         bytes)
{
    bsp_hpput(procnum, src, dest, offset, bytes);
}

void
ParallelDescriptor::ReadData (int         procnum,
                              const void* src,
                              int         offset,
                              void*       dest,
                              int         bytes)
{
    bsp_hpget(procnum, src, offset, dest, bytes);
}

void
ParallelDescriptor::SetMessageHeaderSize (int messageHeaderSize)
{
    bsp_set_tag_size(&messageHeaderSize);
} 

bool
ParallelDescriptor::GetMessageHeader (int&  dataSize,
                                      void* messageHeader)
{
    bsp_get_tag(&dataSize, messageHeader);
    
    return dataSize == -1
	? false  // By bsp definition
	: true;  // A message is waiting
} 

void
ParallelDescriptor::SendData (int         toproc,
                              const void* messageHeader,
                              const void* data,
                              int         datasizeinbytes)
{
    bsp_send(toproc, messageHeader, data, datasizeinbytes);
}

void
ParallelDescriptor::ReceiveData (void* data,
                                 int   datasizeinbytes)
{
    bsp_move(data, datasizeinbytes);
}

void
ParallelDescriptor::Synchronize ()
{
    bsp_sync();
}

void
ParallelDescriptor::Barrier ()
{
    ParallelDescriptor::Synchronize();
}

void
ParallelDescriptor::ReduceBoolAnd (bool& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(bool));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpBoolAnd, &rvar, &rvar, sizeof(bool));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceBoolOr  (bool& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(bool));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpBoolOr , &rvar, &rvar, sizeof(bool));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealSum (Real& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpRealSum, &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealMax (Real& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpRealMax, &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealMin (Real& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpRealMin, &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntSum (int& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpIntSum, &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntMax (int& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpIntMax, &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntMin (int& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpIntMin, &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongSum (long& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpLongSum, &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongMax (long& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpLongMax, &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongMin (long& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpLongMin, &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongAnd (long& rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI) Utility::OpLongAnd, &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

bool
ParallelDescriptor::MessageQueueEmpty ()
{
    int dataWaitingSize;
    FabComTag fabComTag;
    ParallelDescriptor::SetMessageHeaderSize(sizeof(FabComTag));
    return ParallelDescriptor::GetMessageHeader(dataWaitingSize,&fabComTag);
}

void
ParallelDescriptor::Broadcast (int   fromproc,
                               void* src,
                               void* dest,
                               int   nbytes)
{
    bsp_bcast(fromproc, src, dest, nbytes);
}

void
ParallelDescriptor::Gather (Real* sendbuf,
                            int   sendcount,
                            Real* recvbuf,
                            int   root)
{
    assert(root >= 0);
    assert(sendcount > 0);
    assert(!(sendbuf == 0));
    assert(!(recvbuf == 0));

    int myproc = ParallelDescriptor::MyProc();

    ParallelDescriptor::SetMessageHeaderSize(sizeof(int));

    if (!(myproc == root))
    {
        ParallelDescriptor::SendData(root,
                                     &myproc,
                                     sendbuf,
                                     sizeof(Real)*sendcount);
    }
    ParallelDescriptor::Synchronize();

    if (myproc == root)
    {
        int len, nproc;

        for ( ; ParallelDescriptor::GetMessageHeader(len, &nproc); )
        {
            assert(len == sizeof(Real)*sendcount);

            ParallelDescriptor::ReceiveData(recvbuf+nproc, len);
        }
        memcpy(recvbuf+root, sendbuf+root, sizeof(Real)*sendcount);
    }
}

#elif defined(BL_USE_MPI)

#include "mpi.h"

static int numprocs = -1;
static int myid     = -1;

inline
MPI_Datatype
TheRealType ()
{
    return sizeof(Real) == sizeof(float) ? MPI_FLOAT : MPI_DOUBLE;
}

void
ParallelDescriptor::StartParallel (int,
                                   int*    argc,
                                   char*** argv)
{
    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

void
ParallelDescriptor::EndParallel ()
{
    MPI_Finalize();
}

int
ParallelDescriptor::MyProc ()
{
    assert(myid != -1);
    return myid;
}

int
ParallelDescriptor::NProcs ()
{
    assert(numprocs != -1);
    return numprocs;
}

bool
ParallelDescriptor::IOProcessor ()
{
    assert(myid == -1);
    return myid == ioProcessor;
}

int
ParallelDescriptor::IOProcessorNumber ()
{
    return ioProcessor;
}

double
ParallelDescriptor::second ()
{
    return MPI_Wtime();
}

void
ParallelDescriptor::Barrier ()
{
    MPI_Barrier(MPI_COMM_WORLD);
}

void ParallelDescriptor::Synchronize () {}

void
ParallelDescriptor::ReduceRealMax (Real& r)
{
    Real recv;
    MPI_Allreduce(&r, &recv, 1, TheRealType(), MPI_MAX, MPI_COMM_WORLD);
    r = recv;
}

void
ParallelDescriptor::ReduceRealMin (Real& r)
{
    Real recv;
    MPI_Allreduce(&r, &recv, 1, TheRealType(), MPI_MIN, MPI_COMM_WORLD);
    r = recv;
}

void
ParallelDescriptor::ReduceRealSum (Real& r)
{
    Real recv;
    MPI_Allreduce(&r, &recv, 1, TheRealType(), MPI_SUM, MPI_COMM_WORLD);
    r = recv;
}

void
ParallelDescriptor::ReduceLongAnd (long& r)
{
    long recv;
    MPI_Allreduce(&r, &recv, 1, MPI_LONG, MPI_LAND, MPI_COMM_WORLD);
    r = recv;
}

void
ParallelDescriptor::ReduceLongSum (long& r)
{
    long recv;
    MPI_Allreduce(&r, &recv, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
    r = recv;
}

void
ParallelDescriptor::ReduceIntSum (int& r)
{
    int recv;
    MPI_Allreduce(&r, &recv, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    r = recv;
}

void
ParallelDescriptor::Abort (const char* msg)
{
    BoxLib::Warning(msg);
    MPI_Abort(MPI_COMM_WORLD, -1);
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

    MPI_Datatype typ = TheRealType();

    MPI_Gather(sendbuf, nsend, typ, recvbuf, nsend, typ, root, MPI_COMM_WORLD);
}

#else

void
ParallelDescriptor::Gather (Real* sendbuf,
			    int nsend,
			    Real* recvbuf,
			    int root)
{
    assert(root == 0);
    assert(nsend > 0);
    assert(!(sendbuf == 0));
    assert(!(recvbuf == 0));

    for (int i = 0; i < nsend; ++i) recvbuf[i] = sendbuf[i];
}

void ParallelDescriptor::StartParallel(int, int*, char***) {}
void ParallelDescriptor::EndParallel() {}

void ParallelDescriptor::Abort (const char* str) { BoxLib::Abort(str); }
int ParallelDescriptor::MyProc () { return 0; }
int ParallelDescriptor::NProcs () { return 1; }
void ParallelDescriptor::Barrier () {}
void ParallelDescriptor::Synchronize () {}
bool ParallelDescriptor::IOProcessor () { return true; }
int  ParallelDescriptor::IOProcessorNumber () { return 0; }

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


void ParallelDescriptor::ShareVar (const void* var, int bytes) {}
void ParallelDescriptor::UnshareVar (const void* var) {}
void ParallelDescriptor::WriteData (int         procnum,
                                    const void* src,
                                    void*       dest,
                                    int         offset,
                                    int         bytes) {}

void ParallelDescriptor::ReadData (int         procnum,
                                   const void* src,
                                   int         offset,
                                   void*       dest,
                                   int         bytes) {}

void ParallelDescriptor::SetMessageHeaderSize (int messageHeaderSize) {} 

bool ParallelDescriptor::GetMessageHeader (int&  dataSize,
                                           void* messageHeader)
{
    return false;  // No messages waiting.
}
bool ParallelDescriptor::MessageQueueEmpty ()
{
    return true;  // No messages waiting.
} 
void ParallelDescriptor::SendData (int         toproc,
                                   const void* messageHeader,
                                   const void* data,
                                   int         datasizeinbytes) {}

void ParallelDescriptor::ReceiveData (void* data,
                                      int   datasizeinbytes) {}

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
