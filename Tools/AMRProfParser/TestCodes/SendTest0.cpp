// -------------------------------------------------------------
// SendTest0.cpp
// -------------------------------------------------------------
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <iostream>
#include <strstream>
#include <fstream>
#include <iomanip>

#include <unistd.h>

using std::cout;
using std::endl;
using std::ends;
using std::ostrstream;
using std::ofstream;
using std::streamoff;

using namespace amrex;

const int XDIR(0);
const int YDIR(1);
const int ZDIR(2);

// -------------------------------------------------------------
BoxArray MakeBoxArray(int maxgrid,  int nboxes) {
#if (BL_SPACEDIM == 2)
  IntVect ivlo(0, 0);
  IntVect ivhi(maxgrid - 1, maxgrid - 1);
#else
  IntVect ivlo(0, 0, 0);
  IntVect ivhi(maxgrid - 1, maxgrid - 1, maxgrid - 1);
#endif
  int iSide(pow(static_cast<Real>(nboxes), 1.0/3.0));
  Box tempBox(ivlo, ivhi);
  BoxArray bArray(nboxes);
  int ix(0), iy(0), iz(0);
  for(int ibox(0); ibox < nboxes; ++ibox) {
    Box sBox(tempBox);
    sBox.shift(XDIR, ix * maxgrid);
    sBox.shift(YDIR, iy * maxgrid);
#if (BL_SPACEDIM == 3)
    sBox.shift(ZDIR, iz * maxgrid);
#endif
    bArray.set(ibox, sBox);
    ++ix;
    if(ix > iSide) {
      ix = 0;
      ++iy;
    }
    if(iy > iSide) {
      iy = 0;
      ++iz;
    }
  }
  return bArray;
}


// -------------------------------------------------------------
int main(int argc, char *argv[]) {

  amrex::Initialize(argc,argv);

  ParallelDescriptor::Barrier("Starting main.");

  BL_PROFILE_VAR("main()", pmain);

  int myproc(ParallelDescriptor::MyProc());
  int nprocs(ParallelDescriptor::NProcs());
  int maxgrid(64), ncomps(16), nboxes(nprocs), ntimes(6);
  bool raninit(false);
  
  for(int itimes(0); itimes < ntimes; ++itimes) {
    double wallTimeStart(ParallelDescriptor::second());
    if(ParallelDescriptor::IOProcessor()) {
      cout << endl << "--------------------------------------------------" << endl;
    }

    int bSize(20000);
    if(myproc == ParallelDescriptor::IOProcessorNumber()) {
      Vector<int> iBuff(bSize, 0);
      int wakeUpPID = nprocs / 2;
      int tag       = itimes;
      ParallelDescriptor::Send(iBuff.dataPtr(), bSize, wakeUpPID, tag);
								      
    }
    if(myproc == (nprocs / 2)) {
      sleep(1);
      Vector<int> iBuff(bSize, 0);
      int waitForPID = ParallelDescriptor::IOProcessorNumber();
      int tag        = itimes;
      ParallelDescriptor::Recv(iBuff.dataPtr(), bSize, waitForPID, tag);
    }

    ParallelDescriptor::Barrier();

    double wallTime(ParallelDescriptor::second() - wallTimeStart);
    double wallTimeMax(wallTime);
    double wallTimeMin(wallTime);

    ParallelDescriptor::ReduceRealMin(wallTimeMin);
    ParallelDescriptor::ReduceRealMax(wallTimeMax);


    if(ParallelDescriptor::IOProcessor()) {
      cout << "  Wall clock time = " << wallTimeMax << endl;
      cout << "  Min wall clock time = " << wallTimeMin << endl;
      cout << "  Max wall clock time = " << wallTimeMax << endl;
      cout << "==================================================" << endl;
    }


  }

  BL_PROFILE_VAR_STOP(pmain);

  amrex::Finalize();
  return 0;
}
// -------------------------------------------------------------
// -------------------------------------------------------------
