// --------------------------------------------------------------
// AMRProfTest.cpp
// --------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <unistd.h>

#include <AMReX_BLProfiler.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMRPROFTEST_F.H>

using std::cout;
using std::endl;
using namespace amrex;

//===============================================================
namespace {
  const unsigned int msps(1000000);
}


// --------------------------------------------------------------
void ReportSleep(double s, bool us = false) {
  if(ParallelDescriptor::IOProcessor()) {
    if(us) {
      s /= msps;
    }
    std::cout << "Sleeping " << s << " s." << std::endl;
  }
}


// --------------------------------------------------------------
void Nap(double sleepsec) {
  BL_PROFILE_REGION_START("R::Nap");
  BL_PROFILE("Nap()");
  ReportSleep(sleepsec);
  usleep(sleepsec * msps);
  BL_PROFILE_REGION_STOP("R::Nap");
}


// --------------------------------------------------------------
void Nap1234() {
  BL_PROFILE_REGION_START("R::Nap1234");
  BL_PROFILE("Nap1234()");
  double sleepsec(0.1234);
  usleep(sleepsec * msps);
  BL_PROFILE_REGION_STOP("R::Nap1234");
}


// --------------------------------------------------------------
void Sleep(double sleepsec) {
  BL_PROFILE_REGION_START("R::Sleep");
  BL_PROFILE("Sleep()");
  ReportSleep(sleepsec);
  usleep(sleepsec * msps);
  BL_PROFILE_REGION_STOP("R::Sleep");
}


// --------------------------------------------------------------
void SleepNoProf(double sleepsec) {
  ReportSleep(sleepsec);
  usleep(sleepsec * msps);
}


// --------------------------------------------------------------
void SleepOneProc(double sleepsec) {
  BL_PROFILE("SleepOneProc()");
  ReportSleep(sleepsec);
  usleep(sleepsec * msps);
}


// --------------------------------------------------------------
void SleepProcTimes(double sleepsec) {
  BL_PROFILE("SleepProcTimes()");
  ReportSleep(sleepsec);
  usleep(sleepsec * msps);
}


// --------------------------------------------------------------
void NestSleep(double s) {
  BL_PROFILE_REGION_START("R::NestSleep");
  BL_PROFILE("NestSleep()");
  Sleep(0.1 * s);
  SleepNoProf(0.15 * s);
  BL_PROFILE_REGION_STOP("R::NestSleep");
}


// --------------------------------------------------------------
void NestSleep2(double s) {
  BL_PROFILE_REGION_START("R::NestSleep2");
  BL_PROFILE("NestSleep2()");
  SleepNoProf(0.2 * s);
  NestSleep(0.3 * s);
  BL_PROFILE_REGION_STOP("R::NestSleep2");
}


// --------------------------------------------------------------
void NoSleep() {
  BL_PROFILE("NoSleep()");
}


// --------------------------------------------------------------
void RecursiveSleep(double s, int nr) {
  BL_PROFILE_REGION_START("R::RecursiveSleep");
  BL_PROFILE("RecursiveSleep()");
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "_in RecursiveSleep:  nr = " << nr << std::endl;
  }
  Sleep(s);
  if(nr > 0) {
    RecursiveSleep(s, nr - 1);
  }
  BL_PROFILE_REGION_STOP("R::RecursiveSleep");
}


// --------------------------------------------------------------
void SameProfName(double sleepsec) {
  BL_PROFILE_VAR("SameProfName()", SameProfName);
  SleepProcTimes(sleepsec);
  BL_PROFILE_VAR_STOP(SameProfName);
}


// --------------------------------------------------------------
void CProfInt() {
  int largeNCalls(10000);
  BL_PROFILE_VAR_NS("CProfInt()", cprofint);
  for(int i(0); i < largeNCalls; ++i) {
    BL_PROFILE_VAR_START(cprofint);
    //usleep(0.2);
    BL_PROFILE_VAR_STOP(cprofint);
  }
}


// --------------------------------------------------------------
void FlushTest() {
  BL_PROFILE_VAR("FlushTest20()", FlushTest20);
    BL_PROFILE_VAR("FlushTest3()", FlushTest3);
    usleep(0.2);
    BL_PROFILE_VAR_STOP(FlushTest3);

    BL_TRACE_PROFILE_FLUSH();

  BL_PROFILE_VAR_STOP(FlushTest20);

  BL_PROFILE_VAR("FlushTest21()", FlushTest21);
  BL_PROFILE_VAR_STOP(FlushTest21);
}


// --------------------------------------------------------------
int main(int argc, char *argv[]) {

#ifdef BL_USE_MPI
  MPI_Init(&argc, &argv);
#endif

  amrex::Initialize(argc, argv);
  BL_PROFILE_INIT_PARAMS(3.0, true, true);
  BL_PROFILE_REGION_START("R::main");
  BL_PROFILE_VAR("main()", pmain);

  // ---- test the profiling timer
  double tpStart(ParallelDescriptor::second());
  Nap1234();
  if(ParallelDescriptor::IOProcessor()) {
    cout << "Test profiling time = " << ParallelDescriptor::second() - tpStart << endl;
  }

  // ---- test simple functions
  Nap(0.321);  // ---- this one contains profiling
  BL_PROFILE_VAR("SimpleSleepTest()", SimpleSleepTest);
  SleepNoProf(0.42);  // ---- this one contains no profiling
  BL_PROFILE_VAR_STOP(SimpleSleepTest);

  // ---- test nested functions
  NestSleep(0.333);
  NestSleep2(0.555);

  // ---- test recursive function
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "RecursiveNapTest." << std::endl;
  }
  RecursiveSleep(0.2, 4);

  // ---- test function called on only one proc
  if(ParallelDescriptor::IOProcessor()) {
    SleepOneProc(0.4567);
  }

  // ---- test different function times
  int myProc(ParallelDescriptor::MyProc());
  SleepProcTimes(myProc * 0.1111);


/*
  // ---- test fortran functions
  BL_PROFILE_VAR("TESTFORTPROF", testfortprof);
  FORT_AMRPROFTEST();
  BL_PROFILE_VAR_STOP(testfortprof);

  // ---- test fortran functions with a large number of calls
  BL_PROFILE_CHANGE_FORT_INT_NAME("fort_amrproftestint", 8);
  BL_PROFILE_VAR("TESTPROFINT", testfortprofint);
  FORT_AMRPROFTESTINT();
  BL_PROFILE_VAR_STOP(testfortprofint);
*/


/*
  // ---- test c++ functions with a large number of calls
  BL_PROFILE_VAR("CProfIntCall", cprofintcall);
  SleepNoProf(0.42);  // ---- this one contains no profiling
  CProfInt();
  BL_PROFILE_VAR_STOP(cprofintcall);
*/


  // ---- test regions
  BL_PROFILE_REGION_START("R::Part 1");
  SleepNoProf(0.77);
  NoSleep();
  BL_PROFILE_REGION_START("R::P12 overlap");
  SleepNoProf(0.88);
  BL_PROFILE_REGION_STOP("R::Part 1");
  BL_PROFILE_REGION_START("R::Part 2");
  SleepNoProf(0.99);
  BL_PROFILE_REGION_STOP("R::P12 overlap");
  NoSleep();
  BL_PROFILE_REGION_STOP("R::Part 2");


  // ---- test using the same profile name in a called function
  BL_PROFILE_VAR("SameProfName()", SameProfName);
  SameProfName(0.42);
  BL_PROFILE_VAR_STOP(SameProfName);
  Nap(0.2);

  // ---- test unmatched start (an error)
  //BL_PROFILE_VAR("unmatchedstart()", unmatchedstart);

/*
  // ---- test flushing for html output
  BL_TRACE_PROFILE_SETFLUSHSIZE(0);
  BL_PROFILE_VAR("CallingFlushTest", CallingFlushTest);
  usleep(0.3);
  FlushTest();
  BL_PROFILE_VAR_STOP(CallingFlushTest);
  BL_PROFILE_VAR("AfterCallingFlushTest", AfterCallingFlushTest);
  BL_PROFILE_VAR_STOP(AfterCallingFlushTest);
*/

  ParallelDescriptor::Barrier("EndOfMain");


  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_REGION_STOP("R::main");
  usleep(0.1 * msps);

  bool finalizeMPI(false);
  amrex::Finalize(finalizeMPI);

#ifdef BL_USE_MPI
  MPI_Finalize();
#endif

}

// --------------------------------------------------------------
// --------------------------------------------------------------
