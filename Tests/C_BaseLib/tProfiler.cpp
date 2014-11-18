// --------------------------------------------------------------
// tProfiler.cpp
// --------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <string>
#include <sstream>
#include <unistd.h>
using std::cout;
using std::endl;

#include <Profiler.H>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <TPROFILER_F.H>


// --------------------------------------------------------------
void Sleep(unsigned int sleeptime) {
  BL_PROFILE_REGION_START("R::Sleep");
  if(ParallelDescriptor::IOProcessor()) {
    cout << "Sleeping " << sleeptime << endl;
  }
  sleep(sleeptime);
  ParallelDescriptor::Barrier();
  BL_PROFILE_REGION_STOP("R::Sleep");
}


// --------------------------------------------------------------
void nap(unsigned int sleeptime) {
  BL_PROFILE("nap()");
  BL_PROFILE_REGION_START("R::nap");
  sleep(sleeptime);
  BL_PROFILE_REGION_STOP("R::nap");
}


// --------------------------------------------------------------
void napabort(unsigned int sleeptime) {
  BL_PROFILE("napabort()");
  Sleep(sleeptime);
  BoxLib::Finalize();
  BoxLib::Abort("From napabort");
}


// --------------------------------------------------------------
void nestnap(int s) {
  BL_PROFILE("nestnap()");
  BL_PROFILE_REGION_START("R::nestnap");
  Sleep(1 * s);
  nap(2 * s);
  nap(1 * s);
  BL_PROFILE_REGION_STOP("R::nestnap");
}


// --------------------------------------------------------------
void nestnap2(int s) {
  BL_PROFILE("nestnap2()");
  BL_PROFILE_REGION_START("R::nestnap2");
  Sleep(1 * s);
  nestnap(2 * s);
  //Sleep(2 * s);
  nap(3 * s);
  BL_PROFILE_REGION_STOP("R::nestnap2");
}


// --------------------------------------------------------------
void nestnapabort(int s) {
  BL_PROFILE("nestnapabort()");
  Sleep(1 * s);
  napabort(2 * s);
}


// --------------------------------------------------------------
void nonap() {
  BL_PROFILE("nonap()");
}


// --------------------------------------------------------------
int main(int argc, char *argv[]) {
  BoxLib::Initialize(argc, argv);

  sleep(1);
  BL_PROFILE_INIT_PARAMS(3.0, true, true);
  BL_PROFILE_REGION_START("R::main");
  BL_PROFILE_VAR("main()", pmain);

  int myProc(ParallelDescriptor::MyProc());
  int nProcs(ParallelDescriptor::NProcs());

/*
{
  Array<std::string> localStrings, syncedStrings;
  bool alreadySynced;

  localStrings.push_back("allString 0");
  localStrings.push_back("allString 1");
  if(ParallelDescriptor::IOProcessor()) {
    localStrings.push_back("proc_0");
  }
  std::stringstream sstr;
  sstr << "proc_";
  for(int i(0); i < myProc; ++i) {
    sstr << myProc;
  }
  for(int i(0); i < myProc; ++i) {
    localStrings.push_back(sstr.str());
  }
  localStrings.push_back("allString zzz");
  ParallelDescriptor::Barrier();


  BoxLib::SyncStrings(localStrings, syncedStrings, alreadySynced);

  if( ! alreadySynced) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << "not already synced." << endl;
    }
    std::ofstream ofs("syncedstrings.txt");
    ParallelDescriptor::Barrier();
    if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < syncedStrings.size(); ++i) {
        cout << myProc << "::ss[" << i << "] = " << syncedStrings[i] << endl;
        ofs << myProc << "::ss[" << i << "] = " << syncedStrings[i] << endl;
      }
    }
    ofs.close();
  } else {
    if(ParallelDescriptor::IOProcessor()) {
      cout << "already synced." << endl;
    }
  }

}

{
  Array<std::string> localStrings, syncedStrings;
  bool alreadySynced;
  
  localStrings.push_back("samestrings 0");
  localStrings.push_back("samestrings 1");
  std::stringstream sstr;
  localStrings.push_back("allString zzz");
  ParallelDescriptor::Barrier();
  
  BoxLib::SyncStrings(localStrings, syncedStrings, alreadySynced);
  
  if( ! alreadySynced) {
    if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < syncedStrings.size(); ++i) {
        cout << myProc << "::ss[" << i << "] = " << syncedStrings[i] << endl;
      }
    }
  }

  if(ParallelDescriptor::IOProcessor()) {
    if(alreadySynced) {
      cout << "already synced." << endl;
    } else {
      cout << "not already synced." << endl;
    }
  }

}
*/



/*
  if(myProc == nProcs - 1) {
    BL_PROFILE_VAR("NProcs-1 func(2)", pnpm1);
    sleep(2);
    BL_PROFILE_VAR_STOP(pnpm1);
  }
  if(myProc == 0) {
    BL_PROFILE_VAR("Proc 0 func(1)", p0fun);
    sleep(1);
    BL_PROFILE_VAR_STOP(p0fun);
  }
*/

  //sleep(2);
  //nap(3);

/*
  BL_PROFILE_REGION_START("R::Part 1");
  //nap(1);
  //nonap();
  BL_PROFILE_REGION_START("R::P12 overlap");
  nestnap(1);
  BL_PROFILE_REGION_STOP("R::Part 1");
  BL_PROFILE_REGION_START("R::Part 2");
  //nap(1);
  //nestnap2(1);
  //nap(3);
  BL_PROFILE_REGION_STOP("R::P12 overlap");
  //sleep(2);
  //nestnapabort(1);
  nonap();
  BL_PROFILE_REGION_STOP("R::Part 2");
  //Sleep(2);
*/

  Real tpStart(ParallelDescriptor::second());
  BL_PROFILE_VAR("TESTPROF", ptp);
  FORT_TESTPROFILER();
  BL_PROFILE_VAR_STOP(ptp);
  cout << "Test fort time = " << ParallelDescriptor::second() - tpStart << endl;

  Real tpiStart(ParallelDescriptor::second());
  BL_PROFILE_CHANGE_FORT_INT_NAME("eos3", 3);
  BL_PROFILE_VAR("TESTPROFINT", ptpi);
  FORT_TESTPROFILERINT();
  BL_PROFILE_VAR_STOP(ptpi);
  cout << "Test fort int time = " << ParallelDescriptor::second() - tpiStart << endl;

  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_REGION_STOP("R::main");
  //sleep(1);

  BL_PROFILE_FINALIZE();
  BoxLib::Finalize();
}
// --------------------------------------------------------------
// --------------------------------------------------------------
