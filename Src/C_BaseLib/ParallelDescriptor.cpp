//BL_COPYRIGHT_NOTICE

//
// $Id: ParallelDescriptor.cpp,v 1.12 1998-03-25 18:23:08 car Exp $
//
#include <Utility.H>
#include <ParallelDescriptor.H>

#ifdef BL_USE_BSP
#include "bsp.h"

#ifdef FIXBSPLIBLEVEL1HEADER
extern "C"
{
  extern void bsp_fold(void (*)(void *, void *, void *, int *),
                       void*,void*,int);
  extern void bsp_fold_cpp(void (*)(void *, void *, void *, int *),
                           void*,void* , int, int, char *);
}
#endif /*FIXBSPLIBLEVEL1HEADER*/

#include "bsp_level1.h"
#endif /*BL_USE_BSP*/


#if defined(BL_USE_BSP)


void StartParallel(int nprocs)
{
    bsp_begin(nprocs);
}

void StartParallelAllProcs() 
{
    bsp_begin(bsp_nprocs());
}

void EndParallel()            
{
    bsp_end();
}

//
// Type of function pointer required by bsp_fold().
//
typedef void (*VFVVVI)(void*,void*,void*,int*);

int ParallelDescriptor::MyProc ()                 
{ 
    return bsp_pid();                
}
int ParallelDescriptor::NProcs ()                 
{ 
    return bsp_nprocs();             
}
bool ParallelDescriptor::IOProcessor ()           
{ 
    return bsp_pid() == ioProcessor; 
}
int  ParallelDescriptor::IOProcessorNumber ()     
{ 
    return ioProcessor;              
}
void ParallelDescriptor::Abort (const char* msg)  
{ 
    bsp_abort((char*)msg);           
}
double ParallelDescriptor::second ()              
{ 
    return bsp_time();               
}

void ParallelDescriptor::ShareVar (const void* var,
				   int         bytes)
{
    bsp_pushregister(var, bytes);
}

void ParallelDescriptor::UnshareVar (const void* var)
{
    bsp_popregister(var);
}
void ParallelDescriptor::WriteData (int         procnum,
				    const void* src,
				    void*       dest,
				    int         offset,
				    int         bytes)
{
    bsp_hpput(procnum, src, dest, offset, bytes);
}
void ParallelDescriptor::ReadData (int         procnum,
				   const void* src,
				   int         offset,
				   void*       dest,
				   int         bytes)
{
    bsp_hpget(procnum, src, offset, dest, bytes);
}
void ParallelDescriptor::SetMessageHeaderSize (int messageHeaderSize)
{
    bsp_set_tag_size(&messageHeaderSize);
} 
bool ParallelDescriptor::GetMessageHeader (int& dataSize,
					   void* messageHeader)
{
    bsp_get_tag(&dataSize, messageHeader);
    
    return dataSize == -1
	? false  // By bsp definition
	: true;  // A message is waiting
} 

void ParallelDescriptor::SendData (int         toproc,
				   const void* messageHeader,
				   const void* data,
				   int         datasizeinbytes)
{
    bsp_send(toproc, messageHeader, data, datasizeinbytes);
}
void ParallelDescriptor::ReceiveData (void* data,
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
ParallelDescriptor::Synchronize (const char* msg)
{
    cout << "----- " << bsp_pid() << " :  about to sync:  " << msg << endl;
    bsp_sync();
}

void
ParallelDescriptor::ReduceBoolAnd (bool &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(bool));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpBoolAnd, &rvar, &rvar, sizeof(bool));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceBoolOr  (bool &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(bool));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpBoolOr , &rvar, &rvar, sizeof(bool));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealSum (Real &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpRealSum, &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealMax (Real &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpRealMax, &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealMin (Real &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpRealMin, &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntSum (int &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpIntSum, &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntMax (int &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpIntMax, &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntMin (int &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpIntMin, &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongSum (long &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpLongSum, &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongMax (long &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpLongMax, &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongMin (long &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((VFVVVI)Utility::OpLongMin, &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongAnd (long &rvar)
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

void ParallelDescriptor::Broadcast (int fromproc, void*  src, void*  dest, int nbytes)
{
    bsp_bcast(fromproc, (void *) src, (void *) dest, nbytes);
}

#elif defined(BL_USE_MPI)


#else


void StartParallel(int ) {}
void StartParallelAllProcs()  {}
void EndParallel() {}

void ParallelDescriptor::Abort (const char* str)
{
    BoxLib::Abort(str);
}
int ParallelDescriptor::MyProc () { return 0; }
int ParallelDescriptor::NProcs () { return 1; }
void ParallelDescriptor::Synchronize () {}
void ParallelDescriptor::Synchronize (const char* msg)
{
    cout << "----- " << 0 << " :  about to sync:  " << msg << endl;
}
bool ParallelDescriptor::IOProcessor () { return true; }
int  ParallelDescriptor::IOProcessorNumber () { return 0; }

// reduction operations
//template<class T> static void ReduceMin (T &rvar) {}

// bool
void ParallelDescriptor::ReduceBoolAnd (bool& rvar) {}
void ParallelDescriptor::ReduceBoolOr  (bool& rvar) {}
// Real
void ParallelDescriptor::ReduceRealSum (Real& rvar) {}
void ParallelDescriptor::ReduceRealMax (Real& rvar) {}
void ParallelDescriptor::ReduceRealMin (Real& rvar) {}
// int
void ParallelDescriptor::ReduceIntSum (int& rvar) {}
void ParallelDescriptor::ReduceIntMax (int& rvar) {}
void ParallelDescriptor::ReduceIntMin (int& rvar) {}
// long
void ParallelDescriptor::ReduceLongSum (long& rvar) {}
void ParallelDescriptor::ReduceLongMax (long& rvar) {}
void ParallelDescriptor::ReduceLongMin (long& rvar) {}
void ParallelDescriptor::ReduceLongAnd (long& rvar) {}

// data transfer functions
void ParallelDescriptor::ShareVar (const void* var, int         bytes) {}
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
    return false;  // no messages waiting
} 
bool ParallelDescriptor::MessageQueueEmpty ()
{
    return true;  // no messages waiting
} 
void ParallelDescriptor::SendData (int         toproc,
	       const void* messageHeader,
	       const void* data,
	       int         datasizeinbytes) {}

void ParallelDescriptor::ReceiveData (void* data,
		  int   datasizeinbytes) {}

void ParallelDescriptor::Broadcast (int fromproc, void*  src, void*  dest, int nbytes) {}

//
// Here so we don't need to include <Utility.H> in <ParallelDescriptor.H>.
//

double
ParallelDescriptor::second ()
{
    return Utility::second();
}

#endif /*BL_USE_BSP*/
