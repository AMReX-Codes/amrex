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

    amrex::Print() << "Wait: " << endl;
    for (int i=0; i<nProcs; i++)
    { amrex::Print() << "data[" << i << "] = " << data[i] << endl; } 

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

    amrex::Print() << "Waitall: " << endl;
    for (int i=0; i<nProcs; i++)
    { amrex::Print() << "data[" << i << "] = " << data[i] << endl; } 

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

  int finished = 0;
  while (finished < nSends)
  {
    int completed = -1;
    Vector<int> indx(rreqs.size());
    ParallelDescriptor::Waitsome(rreqs, completed, indx, stats);
    finished += completed;

    amrex::Print() << "Waitsome: " << endl;
    for (int i=0; i<indx.size(); i++)
    { amrex::Print() << myProc << " :data[" << i+1 << "] = " << data[i+1] << endl; } 
  }
  amrex::Print() << myProc << " :data[" << myProc+1 << "] = " << data[myProc] << endl; 

  BL_PROFILE_REGION_STOP("R::WaitSome");
}

//---------------------------------------------------------------
void WaitAnyRegion() {
  BL_PROFILE_REGION_START("R::WaitAny");
  BL_PROFILE("WaitAny");
  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());

  Vector<long>        snds(nProcs,0);
  MPI_Status stats;
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

  amrex::Print() << "Waitany: " << endl;
  int indx = -1;
  for (int i=0; i<nSends; ++i)
  {
    ParallelDescriptor::Waitany(rreqs, indx, stats);
    amrex::Print() << myProc << " :data[" << indx+1 << "] = " << data[indx+1] << endl;
  }
  amrex::Print() << myProc << " :data[" << myProc << "] = " << data[myProc] << endl;

  BL_PROFILE_REGION_STOP("R::WaitAny");
}


// --------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  BL_PROFILE_REGION_START("main()");
  BL_PROFILE_VAR("main()", pmain);

/*
  amrex::ParallelDescriptor::Barrier();
  {
    BL_PROFILE_REGION_START("Sleep1");
    BL_PROFILE("Sleep1()");
    amrex::Print() << "Sleep1." << std::endl;
    amrex::Sleep(1);
    BL_PROFILE_REGION_STOP("Sleep1");
  }
  amrex::ParallelDescriptor::Barrier();
  {
    BL_PROFILE_REGION_START("Sleep2");
    BL_PROFILE("Sleep2()");
    amrex::Print() << "Sleep2." << std::endl;
    amrex::Sleep(2);
    BL_PROFILE_REGION_STOP("Sleep2");
  }
  amrex::ParallelDescriptor::Barrier();
  {
    BL_PROFILE_REGION_START("Sleep3");
    BL_PROFILE("Sleep3()");
    amrex::Print() << "Sleep3." << std::endl;
    amrex::Sleep(3);
    BL_PROFILE_REGION_STOP("Sleep3");
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

