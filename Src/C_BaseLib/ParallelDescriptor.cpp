//BL_COPYRIGHT_NOTICE

//
// $Id: ParallelDescriptor.cpp,v 1.9 1998-02-17 22:42:05 lijewski Exp $
//

#include <Utility.H>
#include <ParallelDescriptor.H>

#ifdef BL_USE_BSP

#include <RunStats.H>

//
// Type of function pointer required by bsp_fold().
//
typedef void (*VFVVVI)(void*,void*,void*,int*);

void
ParallelDescriptor::Synchronize ()
{
    RunStats stats("parallel_sync");
    stats.start();
    bsp_sync();
    stats.end();
}

void
ParallelDescriptor::Synchronize (const char* msg)
{
    cout << "----- " << bsp_pid() << " :  about to sync:  " << msg << endl;
    RunStats stats("parallel_sync");
    stats.start();
    bsp_sync();
    stats.end();
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

#else

//
// Here so we don't need to include <Utility.H> in <ParallelDescriptor.H>.
//

double
ParallelDescriptor::second ()
{
    return Utility::second();
}

#endif /*BL_USE_BSP*/
