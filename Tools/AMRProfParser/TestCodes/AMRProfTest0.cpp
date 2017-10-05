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


// --------------------------------------------------------------
void ReportSleep(double s) {
  amrex::Print() << "Sleeping " << s << " s." << std::endl;
}

// --------------------------------------------------------------
void SleepProcTimes(double sleepsec) {
  ReportSleep(sleepsec);
  amrex::USleep(sleepsec);
}

//---------------------------------------------------------------
void InitRegion() {
  BL_PROFILE_REGION_START("R::InitRegion");
  BL_PROFILE("InitRegion()");
  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());
  double sleepTime(1 + myProc * 0.1);

  amrex::Print(Print::AllProcs) << myProc << "::InitRegion = "
    << sleepTime << " s." << endl;
  
  SleepProcTimes(sleepTime);

  BL_PROFILE_REGION_STOP("R::InitRegion");

}

//---------------------------------------------------------------
void ComputeRegion() {
  BL_PROFILE_REGION_START("R::ComputeRegion");
  BL_PROFILE("ComputeRegion()");

  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());

  double sleeptime(2 + (nProcs - 1) - myProc);

  amrex::Print() << "Compute Region." << std::endl;

  SleepProcTimes(sleeptime);
  
  amrex::Print(Print::AllProcs) << myProc << "::ComputeRegion myProc = "
    <<  sleeptime << " s." << endl;

 BL_PROFILE_REGION_STOP("R::ComputeRegion");
}

//---------------------------------------------------------------
void ConcludeRegion() {
  BL_PROFILE_REGION_START("R::ConcludeRegion");
  BL_PROFILE_VAR("ConcludeRegion()", concluderegion);
  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());
  double sleepTime(1 + (nProcs - myProc) * 0.1);

  amrex::Print() << "Conclude Region" << std::endl;

  SleepProcTimes(sleepTime);

  amrex::Print(Print::AllProcs) << myProc << "::ConcludeRegion = "
    << sleepTime << " s." << endl;
  
  BL_PROFILE_VAR_STOP(concluderegion);
  BL_PROFILE_REGION_STOP("R::ConcludeRegion");
}


// --------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  BL_PROFILE_REGION_START("main()");
  BL_PROFILE_VAR("main()", pmain);

  int myProc(ParallelDescriptor::MyProc());

  InitRegion();

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
  /*
  {
    BL_PROFILE_REGION_START("USleepRplus1");
    BL_PROFILE("USleepRplus1()");
    amrex::Print() << "USleepRplus1." << std::endl;
    amrex::USleep(myProc+1);
    BL_PROFILE_REGION_STOP("USleepRplus1");
  }
  amrex::ParallelDescriptor::Barrier();
  */

  ComputeRegion();

  ConcludeRegion();

  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_REGION_STOP("main()");

  amrex::Finalize();
}

// --------------------------------------------------------------
// --------------------------------------------------------------

