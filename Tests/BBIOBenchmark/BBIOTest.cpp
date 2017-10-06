// -------------------------------------------------------------
// BBIOTest.cpp
// -------------------------------------------------------------
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_BoxArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_VisMF.H>
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
using std::ifstream;
using std::streamoff;

using namespace amrex;

const int XDIR(0);
const int YDIR(1);
const int ZDIR(2);
Real bytesPerMB(1.0e+06);

std::string dirName(".");

// -------------------------------------------------------------
void SetDirName(const std::string &dirname) {
  dirName = dirname;
}


// -------------------------------------------------------------
void TestWriteNFiles(int nfiles, int nMB, bool raninit, bool mb2)
{
  int myProc(ParallelDescriptor::MyProc());
  int nProcs(ParallelDescriptor::NProcs());
  bool initData(true);

  if(mb2) {
    bytesPerMB = pow(2.0, 20);
  }

  if(ParallelDescriptor::IOProcessor()) {
    cout << "  Timings for writing to " << nfiles << " files:" << endl;
  }

  // make the data array
  Vector<long> dataArray(nMB * bytesPerMB / sizeof(long), 0);
  if(initData) {
    long *dp = dataArray.dataPtr();
    for(long i(0); i < dataArray.size(); ++i) {
      dp[i] = i;
    }
  }

  Vector<Real> writeSize, writeRate;
  long npts(dataArray.size()), nItemsToWrite(0);
  long totalNBytes(npts * sizeof(long) * nProcs);
  std::string fileName(dirName + "/TestArray_");
  fileName = amrex::Concatenate(fileName, myProc, 4);
  cout << myProc << "::fileName = " << fileName << endl << endl;

  ParallelDescriptor::Barrier();

  for(int buffItems(1024); buffItems < 160000000; buffItems *= 2) {
    if(ParallelDescriptor::IOProcessor()) {
      writeSize.push_back(buffItems * sizeof(long) / bytesPerMB);
    }
    double wallTimeStart(ParallelDescriptor::second());
    long *dp = dataArray.dataPtr();

    ofstream os(fileName.c_str());
    int nItems(dataArray.size());
    while(nItems > 0) {
      nItemsToWrite = nItems > buffItems ? buffItems : nItems;
      os.write((char *) dp, nItemsToWrite * sizeof(long));
      nItems -= nItemsToWrite;
      dp += nItemsToWrite;
    }
    os.close();

    double wallTime(ParallelDescriptor::second() - wallTimeStart);

    double wallTimeMax(wallTime);
    double wallTimeMin(wallTime);

    ParallelDescriptor::ReduceRealMin(wallTimeMin);
    ParallelDescriptor::ReduceRealMax(wallTimeMax);
    Real megabytes((static_cast<Real> (totalNBytes)) / bytesPerMB);

    dp = dataArray.dataPtr();
    for(long i(0); i < dataArray.size(); ++i) {
      dp[i] = -1;
    }
    ifstream is(fileName.c_str());
    is.read((char *) dataArray.dataPtr(), dataArray.size() * sizeof(long));
    is.close();
    int badCheck(0);
    for(long i(0); i < dataArray.size(); ++i) {
      if(dp[i] != i) {
        ++badCheck;
      }
    }
    if(badCheck > 0) {
      cout << myProc << "::**** ****  Error:  BadCheck = " << badCheck << endl;
    }

    if(ParallelDescriptor::IOProcessor()) {
      cout << std::setprecision(5);
      cout << "  Buffer Items    = " << buffItems << endl;
      cout << "  Write Size MB   = " << buffItems * sizeof(long) / bytesPerMB << endl;
      cout << "  BadCheck        = " << badCheck << endl;
      cout << "  Total megabytes = " << megabytes << endl;
      cout << "  Write:  MB/sec  = " << megabytes/wallTimeMax << endl;
      cout << "  Wall clock time = " << wallTimeMax << endl;
      cout << "  Min wall time   = " << wallTimeMin << endl;
      cout << "  Max wall time   = " << wallTimeMax << endl;
      cout << endl << endl;

      writeRate.push_back(megabytes/wallTimeMax);
    }

  }
  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(4) << std::fixed << std::setw(12);
    cout << "Write Size (MB) Write Rate (MB/s)" << endl;
    if(writeSize.size() != writeRate.size()) {
      cout << "**** Error:  writeSize.size() != writeRate.size()" << endl;
    } else {
      for(int i(0); i < writeSize.size(); ++i) {
        cout << std::setprecision(4) << writeSize[i] << "    \t"
	     << std::setprecision(2) << writeRate[i] << endl;
      }
    }
    cout << endl << endl;
  }
}


// -------------------------------------------------------------
void TestReadNFiles(int nfiles, int nMB, bool raninit, bool mb2)
{
  int myProc(ParallelDescriptor::MyProc());
  int nProcs(ParallelDescriptor::NProcs());

  if(mb2) {
    bytesPerMB = pow(2.0, 20);
  }

  if(ParallelDescriptor::IOProcessor()) {
    cout << "  Timings for reading from " << nfiles << " files:" << endl;
  }

  // make the data array
  Vector<long> dataArray(nMB * bytesPerMB / sizeof(long), 0);

  Vector<Real> readSize, readRate;
  long npts(dataArray.size()), nItemsToRead(0);
  long totalNBytes(npts * sizeof(long) * nProcs);
  std::string fileName(dirName + "/TestArray_");
  fileName = amrex::Concatenate(fileName, myProc, 4);
  cout << myProc << "::fileName = " << fileName << endl << endl;

  ParallelDescriptor::Barrier();

  for(int buffItems(1024); buffItems < 160000000; buffItems *= 2) {
    if(ParallelDescriptor::IOProcessor()) {
      readSize.push_back(buffItems * sizeof(long) / bytesPerMB);
    }
    double wallTimeStart(ParallelDescriptor::second());
    long *dp = dataArray.dataPtr();

    ifstream is(fileName.c_str());
    int nItems(dataArray.size());
    while(nItems > 0) {
      nItemsToRead = nItems > buffItems ? buffItems : nItems;
      is.read((char *) dp, nItemsToRead * sizeof(long));
      nItems -= nItemsToRead;
      dp += nItemsToRead;
    }
    is.close();

    double wallTime(ParallelDescriptor::second() - wallTimeStart);

    double wallTimeMax(wallTime);
    double wallTimeMin(wallTime);

    ParallelDescriptor::ReduceRealMin(wallTimeMin);
    ParallelDescriptor::ReduceRealMax(wallTimeMax);
    Real megabytes((static_cast<Real> (totalNBytes)) / bytesPerMB);

    dp = dataArray.dataPtr();
    int badCheck(0);
    for(long i(0); i < dataArray.size(); ++i) {
      if(dp[i] != i) {
        ++badCheck;
      }
    }
    if(badCheck > 0) {
      cout << myProc << "::**** ****  Error:  BadCheck = " << badCheck << endl;
    }

    if(ParallelDescriptor::IOProcessor()) {
      cout << std::setprecision(5);
      cout << "  Buffer Items    = " << buffItems << endl;
      cout << "  Read Size MB    = " << buffItems * sizeof(long) / bytesPerMB << endl;
      cout << "  BadCheck        = " << badCheck << endl;
      cout << "  Total megabytes = " << megabytes << endl;
      cout << "  Read:  MB/sec   = " << megabytes/wallTimeMax << endl;
      cout << "  Wall clock time = " << wallTimeMax << endl;
      cout << "  Min wall time   = " << wallTimeMin << endl;
      cout << "  Max wall time   = " << wallTimeMax << endl;
      cout << endl << endl;

      readRate.push_back(megabytes/wallTimeMax);
    }

  }
  if(ParallelDescriptor::IOProcessor()) {
    cout << std::setprecision(4) << std::fixed << std::setw(12);
    cout << "Read Size (MB) Read Rate (MB/s)" << endl;
    if(readSize.size() != readRate.size()) {
      cout << "**** Error:  readSize.size() != readRate.size()" << endl;
    } else {
      for(int i(0); i < readSize.size(); ++i) {
        cout << std::setprecision(4) << readSize[i] << "    \t"
	     << std::setprecision(2) << readRate[i] << endl;
      }
    }
    cout << endl << endl;
  }
}
// -------------------------------------------------------------
// -------------------------------------------------------------


