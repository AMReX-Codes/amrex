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

#include <AMReX_BLProfiler.H>
#include <AMReX_Utility.H>
#include <AMReX_ParallelDescriptor.H>
#include <TPROFILER_F.H>

#include <AMReX_FArrayBox.H>
#include <AMReX_FabConv.H>
#include <AMReX_FPC.H>

using namespace amrex;

// --------------------------------------------------------------
void Sleep(unsigned int sleeptime) {
  BL_PROFILE("Sleep()");
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
  amrex::Finalize();
  amrex::Abort("From napabort");
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
  amrex::Initialize(argc, argv);

//  sleep(1);
  BL_PROFILE_INIT_PARAMS(3.0, true, true);
  BL_PROFILE_REGION_START("R::main");
  BL_PROFILE_VAR("main()", pmain);

  int myProc(ParallelDescriptor::MyProc());
  int nProcs(ParallelDescriptor::NProcs());

{ // ---- test sync strings
  Vector<std::string> localStrings, syncedStrings;
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


  amrex::SyncStrings(localStrings, syncedStrings, alreadySynced);

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

{ // ---- test already synced strings
  Vector<std::string> localStrings, syncedStrings;
  bool alreadySynced;
  
  localStrings.push_back("samestrings 0");
  localStrings.push_back("samestrings 1");
  std::stringstream sstr;
  localStrings.push_back("allString zzz");
  ParallelDescriptor::Barrier();
  
  amrex::SyncStrings(localStrings, syncedStrings, alreadySynced);
  
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



    BL_PROFILE_VAR("NProcs-1 func(2)", pnpm1);
    sleep(2);
    BL_PROFILE_VAR_STOP(pnpm1);

  if(myProc == 0) {
    BL_PROFILE_VAR("Proc 0 func(1)", p0fun);
    sleep(1);
    BL_PROFILE_VAR_STOP(p0fun);
  }

  BL_PROFILE_REGION_START("R::Sleepnap");
  Sleep(2);
  nap(3);
  BL_PROFILE_REGION_STOP("R::Sleepnap");

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


  //sleep(1);

/*
  Box b(IntVect(0,0), IntVect(80,80));
  Box b0(IntVect(0,0), IntVect(39,80));
  Box b1(IntVect(40,0), IntVect(80,80));
  FArrayBox dfab(b);
  dfab.setVal(42.0e-128, b0, 0, 1);
  dfab.setVal(42.0e+128, b1, 0, 1);
  std::ofstream dfabos("fab.fab");
  dfab.writeOn(dfabos);
  dfabos.close();

  FArrayBox::setFormat(FABio::FAB_IEEE_32);
  FArrayBox sfab(b);
  sfab.copy(dfab);
  std::ofstream sfabos("fab_ieee_32.fab");
  sfab.writeOn(sfabos);
  sfabos.close();

  FArrayBox::setFormat(FABio::FAB_NATIVE);
  FArrayBox nfab(b);
  nfab.copy(dfab);
  std::ofstream nfabos("fab_native.fab");
  nfab.writeOn(nfabos);
  nfabos.close();

  FArrayBox::setFormat(FABio::FAB_NATIVE_32);
  FArrayBox n32fab(b);
  n32fab.copy(dfab);
  std::ofstream n32fabos("fab_native32.fab");
  n32fab.writeOn(n32fabos);
  n32fabos.close();

  cout << "minmax double = " << std::numeric_limits<double>::min()
       << "  " << std::numeric_limits<double>::max() << endl;
  cout << "minmax float  = " << std::numeric_limits<float>::min()
       << "  " << std::numeric_limits<float>::max() << endl;

  int nTimes(100);
  long nPts(128 * 128 * 128);
  Vector<Real> nativeVals(nPts);
  Vector<float> floatVals(nPts);
  Real nMin(std::numeric_limits<Real>::max());
  Real nMax(-std::numeric_limits<Real>::max());
  float fMin(std::numeric_limits<float>::max());
  float fMax(-std::numeric_limits<float>::max());
  unsigned int rint(std::numeric_limits<unsigned int>::max());
  cout << "rint = " << rint << endl;

  for(int i(0); i < nPts; ++i) {
    int ri(amrex::Random_int(rint));
    nativeVals[i] = ri;
    if(i == 0) {
      nativeVals[i] = 1.234e+123;
    }
    if(i == 1) {
      nativeVals[i] = 2.341e-231;
    }
    if(i < 256) {
      cout << "nativeVals[" << i << "] = " << nativeVals[i] << endl;
    }
    nMin = std::min(nMin, nativeVals[i]);
    nMax = std::max(nMax, nativeVals[i]);
  }
  cout << "nMin nMax = " << nMin << "  " << nMax << endl;

  RealDescriptor rdNative(FPC::NativeRealDescriptor());
  cout << "rdNative = " << rdNative << endl;
  RealDescriptor rdIEEE32(FPC::Ieee32NormalRealDescriptor());
  cout << "rdIEEE32 = " << rdIEEE32 << endl;
  RealDescriptor rdNative32(FPC::ieee_float, FPC::reverse_float_order, 4);
  cout << "rdNative32 = " << rdNative32 << endl;
  RealDescriptor n32RD(FPC::Native32RealDescriptor());
  cout << "n32RD = " << n32RD << endl;
*/

/*
// ----
  BL_PROFILE_VAR("TestPD_convert_native", tpdcnative);
  for(int nt(0); nt < nTimes; ++nt) {
    RealDescriptor::convertFromNativeFormat(floatVals.dataPtr(), nPts,
                                            nativeVals.dataPtr(), rdNative32);
  }
  BL_PROFILE_VAR_STOP(tpdcnative);

  fMin = std::numeric_limits<float>::max();
  fMax = -std::numeric_limits<float>::max();
  for(int i(0); i < nPts; ++i) {
    fMin = std::min(fMin, floatVals[i]);
    fMax = std::max(fMax, floatVals[i]);
  }
  cout << "after TestPD_convert_native:  fMin fMax = " << fMin << "  " << fMax << endl;

  // ----
  BL_PROFILE_VAR("TestPD_convert", tpdc);
  for(int nt(0); nt < nTimes; ++nt) {
    RealDescriptor::convertFromNativeFormat(floatVals.dataPtr(), nPts,
                                            nativeVals.dataPtr(), rdIEEE32);
  }
  BL_PROFILE_VAR_STOP(tpdc);

  fMin = std::numeric_limits<float>::max();
  fMax = -std::numeric_limits<float>::max();
  for(int i(0); i < nPts; ++i) {
    fMin = std::min(fMin, floatVals[i]);
    fMax = std::max(fMax, floatVals[i]);
  }
  cout << "after TestPD_convert:  fMin fMax = " << fMin << "  " << fMax << endl;
*/

/*
  // ----
  BL_PROFILE_VAR("TestPD_convert_native32", tpdcn32);
  for(int nt(0); nt < nTimes; ++nt) {
    RealDescriptor::convertFromNativeFormat(floatVals.dataPtr(), nPts,
                                            nativeVals.dataPtr(), n32RD);
  }
  BL_PROFILE_VAR_STOP(tpdcn32);

  fMin = std::numeric_limits<float>::max();
  fMax = -std::numeric_limits<float>::max();
  for(int i(0); i < nPts; ++i) {
    fMin = std::min(fMin, floatVals[i]);
    fMax = std::max(fMax, floatVals[i]);
  }
  cout << "after TestPD_convert_native32:  fMin fMax = " << fMin << "  " << fMax << endl;

  // ----
  BL_PROFILE_VAR("TestCast_convert", tcc);
  for(int nt(0); nt < nTimes; ++nt) {
    for(int i(0); i < nPts; ++i) {
      floatVals[i] = nativeVals[i];
    }
  }
  BL_PROFILE_VAR_STOP(tcc);

  fMin = std::numeric_limits<float>::max();
  fMax = -std::numeric_limits<float>::max();
  for(int i(0); i < nPts; ++i) {
    fMin = std::min(fMin, floatVals[i]);
    fMax = std::max(fMax, floatVals[i]);
  }
  cout << "after TestCast_convert:  fMin fMax = " << fMin << "  " << fMax << endl;
*/



  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_REGION_STOP("R::main");

  BL_PROFILE_FINALIZE();
  amrex::Finalize();
}

// --------------------------------------------------------------
// --------------------------------------------------------------



