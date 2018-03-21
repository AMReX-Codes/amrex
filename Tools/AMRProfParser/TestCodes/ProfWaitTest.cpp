// --------------------------------------------------------------
// AMRProfTest.cpp
// --------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <math.h>
using std::cout;
using std::endl;

#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>

using namespace amrex;

//---------------------------------------------------------------
void WaitRegion() {
  BL_PROFILE_REGION_START("R::Wait");
  BL_PROFILE("Wait");
  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());

  Vector<long>        snds(nProcs,0);
  Vector<MPI_Status>  stats(nProcs-1);
  Vector<MPI_Request> rreqs(nProcs-1);
  const int SeqNum = ParallelDescriptor::SeqNum();
  Vector<int> data(nProcs, myProc);
  int nSends = nProcs-1;
  int j;

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      rreqs[j] = ParallelDescriptor::Arecv(&data[i], 1, i, SeqNum).req();
      j++;
    }
  }

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      ParallelDescriptor::Send(&data[j], 1, i, SeqNum);
      j++;
    } 
  }

  for (int i=0; i<nSends; ++i)
  {
    ParallelDescriptor::Wait(rreqs[i], stats[i]);
  }

  if (myProc == 0)
  {
    cout << "Wait: " << endl;
    for (int i=0; i<nProcs; i++)
    { cout << "data[" << i << "] = " << data[i] << endl; } 
  }

  BL_PROFILE_REGION_STOP("R::Wait");
}

//---------------------------------------------------------------
void WaitAllRegion() {
  BL_PROFILE_REGION_START("R::WaitAll");
  BL_PROFILE("WaitAll");
  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());

  Vector<long>        snds(nProcs,0);
  Vector<MPI_Status>  stats(nProcs-1);
  Vector<MPI_Request> rreqs(nProcs-1);
  const int SeqNum = ParallelDescriptor::SeqNum();
  Vector<int> data(nProcs, myProc);
  int nSends = nProcs-1;
  int j;

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      rreqs[j] = ParallelDescriptor::Arecv(&data[i], 1, i, SeqNum).req();
      j++;
    }
  }

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      ParallelDescriptor::Send(&data[j], 1, i, SeqNum);
      j++;
    } 
  }

  ParallelDescriptor::Waitall(rreqs, stats);

  if (myProc == 0)
  {
    cout << "Waitall: " << endl;
    for (int i=0; i<nProcs; i++)
    { cout << "data[" << i << "] = " << data[i] << endl; } 
  }

  BL_PROFILE_REGION_STOP("R::WaitAll");
}

//---------------------------------------------------------------
void WaitSomeRegion() {
  BL_PROFILE_REGION_START("R::WaitSome");
  BL_PROFILE("WaitSome");
  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());

  Vector<long>        snds(nProcs,0);
  Vector<MPI_Status>  stats(nProcs-1);
  Vector<MPI_Request> rreqs(nProcs-1);
  const int SeqNum = ParallelDescriptor::SeqNum();
  Vector<int> data(nProcs, myProc);
  int nSends = nProcs-1;
  int j;

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      rreqs[j] = ParallelDescriptor::Arecv(&data[i], 1, i, SeqNum).req();
      j++;
    }
  }

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      ParallelDescriptor::Send(&data[j], 1, i, SeqNum);
      j++;
    } 
  }

  int completed = -1; 
  Vector<int> indx;
  ParallelDescriptor::Waitsome(rreqs, completed, indx, stats);

  if (myProc == 0)
  {
    cout << "Waitsome: " << endl;
    for (int i=0; i<nProcs; i++)
    { cout << "data[" << i << "] = " << data[i] << endl; } 
  }
  BL_PROFILE_REGION_STOP("R::WaitSome");
}

//---------------------------------------------------------------
void WaitAnyRegion() {
  BL_PROFILE_REGION_START("R::WaitAny");
  BL_PROFILE("WaitAny");
  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());

  Vector<long>        snds(nProcs,0);
  Vector<MPI_Status>  stats(nProcs-1);
  Vector<MPI_Request> rreqs(nProcs-1);
  const int SeqNum = ParallelDescriptor::SeqNum();
  Vector<int> data(nProcs, myProc);
  int nSends = nProcs-1;
  int j;

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      rreqs[j] = ParallelDescriptor::Arecv(&data[i], 1, i, SeqNum).req();
      j++;
    }
  }

  j=0;
  for (int i=0; i<nProcs; ++i)
  {
    if (i != myProc)
    {
      ParallelDescriptor::Send(&data[j], 1, i, SeqNum);
      j++;
    } 
  }

  int indx = -1;
  ParallelDescriptor::Waitany(rreqs, indx, stats);

  if (myProc == 0)
  {
    cout << "Waitall: " << endl;
    for (int i=0; i<nProcs; i++)
    { cout << "data[" << i << "] = " << data[i] << endl; } 
  }
  BL_PROFILE_REGION_STOP("R::WaitAny");
}


// --------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  BL_PROFILE_REGION_START("main()");
  BL_PROFILE_VAR("main()", pmain);

  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());

/*
  amrex::ParallelDescriptor::Barrier();
  {
    BL_PROFILE_REGION_START("USleep1");
    BL_PROFILE("USleep1()");
    amrex::Print() << "USleep1." << std::endl;
    amrex::USleep(1);
    BL_PROFILE_REGION_STOP("USleep1");
  }
  amrex::ParallelDescriptor::Barrier();
  {
    BL_PROFILE_REGION_START("USleep2");
    BL_PROFILE("USleep2()");
    amrex::Print() << "USleep2." << std::endl;
    amrex::USleep(2);
    BL_PROFILE_REGION_STOP("USleep2");
  }
  amrex::ParallelDescriptor::Barrier();
  {
    BL_PROFILE_REGION_START("USleep3");
    BL_PROFILE("USleep3()");
    amrex::Print() << "USleep3." << std::endl;
    amrex::USleep(3);
    BL_PROFILE_REGION_STOP("USleep3");
  }
  amrex::ParallelDescriptor::Barrier();
*/

  WaitRegion();
  amrex::ParallelDescriptor::Barrier("Wait");

  WaitAllRegion();
  amrex::ParallelDescriptor::Barrier("WaitAll");

  WaitAnyRegion();
  amrex::ParallelDescriptor::Barrier("WaitAny");

  WaitSomeRegion();
  amrex::ParallelDescriptor::Barrier("WaitSome");

  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_REGION_STOP("main()");

  amrex::Finalize();
}

// --------------------------------------------------------------
// --------------------------------------------------------------

