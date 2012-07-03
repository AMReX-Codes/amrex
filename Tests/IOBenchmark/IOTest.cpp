// -------------------------------------------------------------
// IOTest.cpp
// -------------------------------------------------------------
#include <Array.H>
#include <IntVect.H>
#include <Box.H>
#include <BoxArray.H>
#include <FArrayBox.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <VisMF.H>
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

const int XDIR(0);

// -------------------------------------------------------------
BoxArray MakeBoxArray(int maxgrid,  int nboxes) {
#if (BL_SPACEDIM == 2)
  IntVect ivlo(0, 0);
  IntVect ivhi(maxgrid - 1, maxgrid - 1);
#else
  IntVect ivlo(0, 0, 0);
  IntVect ivhi(maxgrid - 1, maxgrid - 1, maxgrid - 1);
#endif
  Box tempBox(ivlo, ivhi);
  BoxArray bArray(nboxes);
  for(int ix(0); ix < nboxes; ++ix) {
    Box sBox(tempBox);
    sBox.shift(XDIR, ix * maxgrid);
    bArray.set(ix, sBox);
  }

  return bArray;
}  // end MakeBoxArray()


// -------------------------------------------------------------
void TestWriteNFiles(int nfiles, int maxgrid, int ncomps, int nboxes) {
  int myProc(ParallelDescriptor::MyProc());

  VisMF::SetNOutFiles(nfiles);

  BoxArray bArray(MakeBoxArray(maxgrid, nboxes));
  if(ParallelDescriptor::IOProcessor()) {
    cout << "  Timings for writing to " << nfiles << " files:" << endl;
  }

  // make a MultiFab
  MultiFab mfout(bArray, ncomps, 0);
  for(MFIter mfiset(mfout); mfiset.isValid(); ++mfiset) {
    for(int invar(0); invar < ncomps; ++invar) {
      mfout[mfiset].setVal((100.0 * mfiset.index()) + invar, invar);
    }
  }

  long npts(bArray[0].numPts());
  long totalNBytes(npts * ncomps * nboxes *sizeof(Real));
  std::string mfName("TestMF");

  ParallelDescriptor::Barrier();
  double wallTimeStart(ParallelDescriptor::second());

  VisMF::Write(mfout, mfName); 

  double wallTime(ParallelDescriptor::second() - wallTimeStart);

  double wallTimeMax(wallTime);
  double wallTimeMin(wallTime);

  ParallelDescriptor::ReduceRealMin(wallTimeMin);
  ParallelDescriptor::ReduceRealMax(wallTimeMax);

  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(5);
    cout << "  Total megabytes = " << ((Real) totalNBytes/1000000.0) << endl;
    cout << "  Megabytes/sec   = "
	 << ((Real) totalNBytes/wallTimeMax)/1000000.0 << endl;
    cout << "  Wall clock time = " << wallTimeMax << endl;
    cout << "  Min wall clock time = " << wallTimeMin << endl;
    cout << "  Max wall clock time = " << wallTimeMax << endl;
  }
}


// -------------------------------------------------------------
void TestReadMF() {
  int myProc(ParallelDescriptor::MyProc());

  MultiFab mfin;
  std::string mfName("TestMF");

  ParallelDescriptor::Barrier();
  double wallTimeStart(ParallelDescriptor::second());

  VisMF::Read(mfin, mfName); 

  double wallTime(ParallelDescriptor::second() - wallTimeStart);

  double wallTimeMax(wallTime);
  double wallTimeMin(wallTime);

  ParallelDescriptor::ReduceRealMin(wallTimeMin);
  ParallelDescriptor::ReduceRealMax(wallTimeMax);

  long npts(mfin.boxArray()[0].numPts());
  int  ncomps(mfin.nComp());
  int  nboxes(mfin.boxArray().size());
  long totalNBytes(npts * ncomps * nboxes *sizeof(Real));

  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(5);
    cout << "  ncomps = " << ncomps << endl;
    cout << "  nboxes = " << nboxes << endl;
    cout << "  Total megabytes = " << ((Real) totalNBytes/1000000.0) << endl;
    cout << "  Megabytes/sec   = "
	 << ((Real) totalNBytes/wallTimeMax)/1000000.0 << endl;
    cout << "  Wall clock time = " << wallTimeMax << endl;
    cout << "  Min wall clock time = " << wallTimeMin << endl;
    cout << "  Max wall clock time = " << wallTimeMax << endl;
  }
}
// -------------------------------------------------------------
// -------------------------------------------------------------




