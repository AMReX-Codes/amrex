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
#include <Utility.H>
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
const int YDIR(1);
const int ZDIR(2);
Real bytesPerMB(1.0e+06);

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
void TestWriteNFiles(int nfiles, int maxgrid, int ncomps, int nboxes,
                     bool raninit, bool mb2)
{
  int myProc(ParallelDescriptor::MyProc());

  VisMF::SetNOutFiles(nfiles);
  if(mb2) {
    bytesPerMB = pow(2.0, 20);
  }

  BoxArray bArray(MakeBoxArray(maxgrid, nboxes));
  if(ParallelDescriptor::IOProcessor()) {
    cout << "  Timings for writing to " << nfiles << " files:" << endl;
  }

  // make a MultiFab
  MultiFab mfout(bArray, ncomps, 0);
  for(MFIter mfiset(mfout); mfiset.isValid(); ++mfiset) {
    for(int invar(0); invar < ncomps; ++invar) {
      if(raninit) {
        Real *dp = mfout[mfiset].dataPtr(invar);
	for(int i(0); i < mfout[mfiset].box().numPts(); ++i) {
	  dp[i] = BoxLib::Random() + (1.0 + static_cast<Real> (invar));
	}
      } else {
        mfout[mfiset].setVal((100.0 * mfiset.index()) + invar, invar);
      }
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
  Real megabytes((static_cast<Real> (totalNBytes)) / bytesPerMB);

  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(5);
    cout << "  Total megabytes = " << megabytes << endl;
    cout << "  Write:  Megabytes/sec   = " << megabytes/wallTimeMax << endl;
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
  Real megabytes((static_cast<Real> (totalNBytes)) / bytesPerMB);

  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(5);
    cout << "  ncomps = " << ncomps << endl;
    cout << "  nboxes = " << nboxes << endl;
    cout << "  Total megabytes = " << megabytes << endl;
    cout << "  Read:  Megabytes/sec   = " << megabytes/wallTimeMax << endl;
    cout << "  Wall clock time = " << wallTimeMax << endl;
    cout << "  Min wall clock time = " << wallTimeMin << endl;
    cout << "  Max wall clock time = " << wallTimeMax << endl;
  }
}
// -------------------------------------------------------------
// -------------------------------------------------------------


