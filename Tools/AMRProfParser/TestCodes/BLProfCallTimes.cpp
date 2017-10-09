// --------------------------------------------------------------
// BLProfCallTimes.cpp
// --------------------------------------------------------------
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <math.h>
using std::cout;
using std::endl;

#include <AMReX_Utility.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;


//---------------------------------------------------------------
void TestCall(const std::string &fname) {
  BL_PROFILE(fname);
}


//---------------------------------------------------------------
void TestFunctionCalls() {
  BL_PROFILE("TestFunctionCalls()");
  BL_PROFILE_REGION_START("TestFunctionCalls()");

  //cout << "Testing function calls." << endl;

  //long nChunks(100), chunkSize(1000000), nFuncs(12);
  long nChunks(8), chunkSize(100000), nFuncs(12);
  Real totalTime(0.0);
  std::string fNameBase("fName_");

  for(long t(0); t < nChunks; ++t) {
    std::stringstream fNameStr;
    fNameStr << fNameBase << amrex::Random_int(nFuncs);
    std::string fName(fNameStr.str());
    Real startTime(ParallelDescriptor::second());
    for(long n(0); n < chunkSize; ++n) {
      TestCall(fName);
    }
    totalTime += ParallelDescriptor::second() - startTime;
  }

  long nCalls(nChunks * chunkSize);
  long millionCalls(nCalls / 1000000);
  cout << ParallelDescriptor::MyProc()
       << "::  TestFunctionCalls:  NCalls/proc = " << millionCalls
       << " million   totalTime = " << totalTime
       << "  mCPS = " << millionCalls / totalTime << endl;
  BL_PROFILE_REGION_STOP("TestFunctionCalls()");
}


//---------------------------------------------------------------
void TestReductionCalls() {
  BL_PROFILE("TestReductionCalls()");

  if(ParallelDescriptor::IOProcessor()) {
    cout << "Testing reduction calls." << endl;
  }

  //long nChunks(100), chunkSize(1000000);
  long nChunks(10), chunkSize(1000000);
  Real totalTime(0.0);

  for(long t(0); t < nChunks; ++t) {
    Real startTime(ParallelDescriptor::second());
    for(long n(0); n < chunkSize; ++n) {
      BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, BLProfiler::BeforeCall(), true);
      BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, sizeof(Real), false);
    }
    totalTime += ParallelDescriptor::second() - startTime;
  }

  long nCalls(nChunks * chunkSize * 2);  // ---- two calls in loop above
  long millionCalls(nCalls / 1000000);
  cout << ParallelDescriptor::MyProc()
       << "::  TestFunctionCalls:  NCalls/proc = " << millionCalls
       << " million   totalTime = " << totalTime
       << "  mCPS = " << millionCalls / totalTime << endl;
}


//---------------------------------------------------------------
void TestCommCalls() {
  BL_PROFILE("TestCommCalls()");
  BL_PROFILE_REGION_START("TestCommCalls()");

  if(ParallelDescriptor::IOProcessor()) {
    cout << "Testing mixed communication calls." << endl;
  }

//long nChunks(100), chunkSize(1000000);
  long nChunks(4), chunkSize(1000);
  Real totalTime(0.0);
  long reqsSize(37);
  Vector<MPI_Request> reqs(reqsSize);
  Vector<MPI_Status> status(reqsSize);
  Vector<int> index(reqsSize);
  int completed(0), sendSize(32), destPid(42), tag(127);

  for(long t(0); t < nChunks; ++t) {
    BL_PROFILE_REGION_START("Chunks()");
    Real startTime(ParallelDescriptor::second());
    for(long n(0); n < chunkSize; ++n) {
      BL_COMM_PROFILE_BARRIER("Barrier0", true);
      BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, BLProfiler::BeforeCall(), true);
      BL_COMM_PROFILE_ALLREDUCE(BLProfiler::AllReduceR, sizeof(Real), false);
      BL_COMM_PROFILE_BARRIER("Barrier0", false);
      BL_COMM_PROFILE_BARRIER("Barrier1", true);
      BL_COMM_PROFILE_REDUCE(BLProfiler::AllReduceR, BLProfiler::BeforeCall(), true);
      BL_COMM_PROFILE(BLProfiler::AsendTsii, sendSize * sizeof(Real), destPid, tag);
      BL_COMM_PROFILE_WAITSOME(BLProfiler::Waitsome, reqs, completed, index, status, true);
      BL_COMM_PROFILE(BLProfiler::ArecvTsiiM, sendSize * sizeof(Real), destPid, tag);
      BL_COMM_PROFILE_BARRIER("Barrier1", false);
    }
    totalTime += ParallelDescriptor::second() - startTime;
    BL_PROFILE_REGION_STOP("Chunks()");
  }

  long nCalls(nChunks * chunkSize * 10);  // ---- ten calls in loop above
  long millionCalls(nCalls / 1000000);
  cout << ParallelDescriptor::MyProc()
       << "::  TestFunctionCalls:  NCalls/proc = " << millionCalls
       << " million   totalTime = " << totalTime
       << "  mCPS = " << millionCalls / totalTime << endl;

  BL_PROFILE_REGION_STOP("TestCommCalls()");
}


//---------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv, false);

  BL_PROFILE_REGION_START("main()");
  BL_PROFILE_VAR("main()", pmain);

  bool testFilter(false);
  if(testFilter) {
    BL_COMM_PROFILE_FILTER(BLProfiler::AllReduceR);
  }
  //TestReductionCalls();

  ParallelDescriptor::Barrier();

  TestFunctionCalls();

  ParallelDescriptor::Barrier();

  if(argc == 2) {
    int oldcsfs(BLProfiler::GetCSFlushSize());
    int csfs(atoi(argv[1]));
    BLProfiler::SetCSFlushSize(csfs);
    int newcsfs(BLProfiler::GetCSFlushSize());
    if(ParallelDescriptor::IOProcessor()) {
      std::cout << "---- Setting CSFlushSize to:  " << csfs
                << "  :: old new = " << oldcsfs << "  " << newcsfs << std::endl;
    }
  }

  TestCommCalls();

  ParallelDescriptor::Barrier();

  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_REGION_STOP("main()");

  amrex::Finalize();
}

// --------------------------------------------------------------
// --------------------------------------------------------------

