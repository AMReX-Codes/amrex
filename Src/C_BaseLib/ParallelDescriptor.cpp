//BL_COPYRIGHT_NOTICE

//
// $Id: ParallelDescriptor.cpp,v 1.5 1997-11-20 00:49:44 lijewski Exp $
//

#ifdef BL_USE_BSP

#include <ParallelDescriptor.H>
#include <Utility.H>

void
ParallelDescriptor::ReduceBoolAnd (bool &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(bool));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpBoolAnd,
              &rvar, &rvar, sizeof(bool));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealSum (Real &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpRealSum,
              &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealMax (Real &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpRealMax,
              &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceRealMin (Real &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(Real));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpRealMin,
              &rvar, &rvar, sizeof(Real));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntSum (int &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpIntSum,
              &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntMax (int &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpIntMax,
              &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceIntMin (int &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(int));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpIntMin,
              &rvar, &rvar, sizeof(int));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongSum (long &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpLongSum,
              &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongMax (long &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpLongMax,
              &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongMin (long &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpLongMin,
              &rvar, &rvar, sizeof(long));
    ParallelDescriptor::UnshareVar(&rvar);
}

void
ParallelDescriptor::ReduceLongAnd (long &rvar)
{
    ParallelDescriptor::ShareVar(&rvar, sizeof(long));
    ParallelDescriptor::Synchronize();
    bsp_fold((void (*)(void *, void *, void *, int *)) Utility::OpLongAnd,
              &rvar, &rvar, sizeof(long));
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

#endif
