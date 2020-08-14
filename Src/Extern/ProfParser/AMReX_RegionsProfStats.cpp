// ----------------------------------------------------------------------
//  RegionsProfStats.cpp
// ----------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <map>
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <limits>
#include <iomanip>
#include <set>
#include <sys/time.h>
#include<unistd.h>

using std::cout;
using std::endl;
using std::flush;
using std::string;
using std::ifstream;
using std::ofstream;
using std::map;
using std::unordered_map;
using std::unordered_multimap;
using std::vector;
using std::pair;

#include <AMReX_RegionsProfStats.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_BLProfUtilities.H>

using namespace amrex;


#ifdef _OPENMP
#include <omp.h>
#endif

#define SHOWVAL(val) { cout << #val << " = " << val << endl; }

const int XDIR(0);
const int YDIR(1);
const int ZDIR(2);

bool RegionsProfStats::bInitDataBlocks(true);
bool RegionsProfStats::persistentStreams(true);

std::map<Real, std::string, std::greater<Real> > RegionsProfStats::mTimersTotalsSorted;
Vector<std::string> RegionsProfStats::regHeaderFileNames;
std::map<std::string, int> RegionsProfStats::regDataFileNames;
Vector<std::ifstream *> RegionsProfStats::regDataStreams;

extern std::string SanitizeName(const std::string &s);
extern void amrex::MakeFuncPctTimesMF(const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                                      const Vector<std::string> &blpFNames,
			              const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
			              Real runTime, int dataNProcs);
extern void amrex::CollectMProfStats(std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                                     const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                                     const Vector<std::string> &fNames,
                                     Real runTime, int whichProc);
extern void amrex::GraphTopPct(const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                               const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                               const Vector<std::string> &fNames,
                               Real runTime, int dataNProcs, Real gPercent);


#define PRINTCS(CS) CS.csFNameNumber << " :: " << fNumberNames[CS.csFNameNumber] << " :: " \
                 << CS.totalTime << " :: " << CS.stackTime << " :: " << \
                 ((CS.totalTime > 0.0) ? (( 1.0 - (CS.stackTime / CS.totalTime)) * 100.0) :  ( 0.0 )) \
		 << " % :: " \
		 << CS.nCSCalls  << " :: " << CS.callStackDepth << " :: " \
		 << CS.callTime

#define PRINTCSNC(CS) CS.csFNameNumber << " :: " << fNumberNames[CS.csFNameNumber] << " :: " \
                   << CS.totalTime << " :: " << CS.stackTime << " :: " << \
                   ((CS.totalTime > 0.0) ? (( 1.0 - (CS.stackTime / CS.totalTime)) * 100.0) :  ( 0.0 )) \
		   << " % :: " \
		   << CS.nCSCalls  << " :: " << CS.callStackDepth


// ----------------------------------------------------------------------
RegionsProfStats::RegionsProfStats()
  :  currentDataBlock(0)
{
  maxFNumber = -1;
}


// ----------------------------------------------------------------------
RegionsProfStats::~RegionsProfStats() {
}


// ----------------------------------------------------------------------
void RegionsProfStats::AddCStatsHeaderFileName(const string &hfn) {
  if(std::find(regHeaderFileNames.begin(),
     regHeaderFileNames.end(), hfn) == regHeaderFileNames.end())
  {
    regHeaderFileNames.push_back(hfn);
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::SyncFNamesAndNumbers() {
  if(mFNameNumbersPerProc.size() == 0) {
    return;
  }
  std::set<std::string> localNames;
  for(int p(0); p < dataNProcs; ++p) {
    std::map<std::string, int>::iterator mfnnit;
    for(mfnnit = mFNameNumbersPerProc[p].begin();
        mfnnit != mFNameNumbersPerProc[p].end(); ++mfnnit)
    {
      localNames.insert(mfnnit->first);
    }
  }
  Vector<std::string> localStrings, syncedStrings;
  bool alreadySynced;
  for(std::set<std::string>::iterator lsit = localNames.begin();
      lsit != localNames.end(); ++lsit)
  {
    localStrings.push_back(*lsit);
  }

  amrex::SyncStrings(localStrings, syncedStrings, alreadySynced);
  fnameRemap.resize(dataNProcs);

  for(int p(0); p < dataNProcs; ++p) {
    fnameRemap[p].resize(syncedStrings.size(), -1);  // -1 are names not on proc p
    Vector<int> foundStrings(syncedStrings.size(), -1);
    std::map<std::string, int>::iterator mfnnit;
    for(mfnnit = mFNameNumbersPerProc[p].begin();
        mfnnit != mFNameNumbersPerProc[p].end(); ++mfnnit)
    {
      const std::string &findName = mfnnit->first;
      int localIndex(mfnnit->second);
      for(int n(0); n < syncedStrings.size(); ++n) {
        if(findName == syncedStrings[n]) {
	  fnameRemap[p][localIndex] = n;
	  foundStrings[n] = n;
	  //cout << "      p fname localIndex n = " << p << "  " << findName << "  "
	       //<< localIndex << "  " << n << endl;
	}
      }
    }
    for(int n(0); n < foundStrings.size(); ++n) {  // fill in unfound strings
      if(foundStrings[n] < 0) {
        for(int ii(0); ii < fnameRemap[p].size(); ++ii) {
	  if(fnameRemap[p][ii] < 0) {
	    fnameRemap[p][ii] = n;
	    break;
	  }
	}
      }
    }
  }

  // ---- the index here is the remapped one, must access numbersToFName with
  // ---- the data processors local index in fnameRemap
  numbersToFName = syncedStrings;

}


// ----------------------------------------------------------------------
void RegionsProfStats::AddFunctionName(const std::string &fname, int fnumber) {
  if(mFNameNumbersPerProc.size() != dataNProcs) {
    mFNameNumbersPerProc.resize(dataNProcs);
  }
  std::size_t found;
  std::string fnameNQ(fname.substr(1, fname.length() - 2));  // ---- remove quotes
  while((found = fnameNQ.find(" ")) != std::string::npos) {  // ---- replace spaces
    fnameNQ.replace(found, 1, "_");
  }

  std::map<std::string, int>::iterator mfnnit = mFNameNumbersPerProc[currentProc].find(fnameNQ);
  if(mfnnit == mFNameNumbersPerProc[currentProc].end()) {
    mFNameNumbersPerProc[currentProc].insert(std::pair<std::string, int>(fnameNQ, fnumber));
  } else {
    if(mfnnit->first != fnameNQ || mfnnit->second != fnumber) {
      cout << "************ conflict:  fname fnum = " << fnameNQ << "  " << fnumber << endl;
    }
  }

  maxFNumber = std::max(maxFNumber, fnumber);
}


// ----------------------------------------------------------------------
void RegionsProfStats::AddTimeMinMax(double tmin, double tmax) {
  dataBlocks[currentDataBlock].timeMin = tmin;
  dataBlocks[currentDataBlock].timeMax = tmax;
}


// ----------------------------------------------------------------------
BLProfStats::TimeRange RegionsProfStats::MakeRegionPlt(FArrayBox &rFab, int noregionnumber,
                                     int width, int height,
				     Vector<Vector<Box>> &regionBoxes)
{
  amrex::ignore_unused(noregionnumber);

#if (BL_SPACEDIM != 2)
  cout << "**** Error:  RegionsProfStats::MakeRegionPlt only supported for 2D" << endl;
  return TimeRange(0, 0);
#else
  BL_PROFILE("RegionsProfStats::MakeRegionPlt()");
  int xLength(width), yHeight(height);
  int nRegions(maxRNumber + 1);
  int whichProc(ParallelDescriptor::IOProcessorNumber());
  Real notInRegionValue(-1.0);
  Box b(IntVect(0, 0), IntVect(xLength - 1, (yHeight * nRegions) - 1));
  rFab.resize(b, 1);
  rFab.setVal<RunOn::Host>(notInRegionValue);

  Vector<Real> rStartTime(nRegions, -1.0);
  regionBoxes.clear();
  regionBoxes.resize(nRegions);

  // need a better way to get the real minmax time
  Real timeMax(-std::numeric_limits<Real>::max());
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    timeMax = std::max(timeMax, dBlock.timeMax);
    if(dBlock.proc == whichProc) {
      ReadBlock(dBlock, true, false);  // dont need to read the trace data
      for(int iss(0); iss < dBlock.rStartStop.size(); ++iss) {
        BLProfiler::RStartStop &rss = dBlock.rStartStop[iss];
        timeMax = std::max(timeMax, rss.rssTime);
      }
      ClearBlock(dBlock);
    }
  }

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(dBlock.proc != whichProc) {
      continue;
    }
    ReadBlock(dBlock, true, false);  // dont need to read the trace data

    for(int i(0); i < dBlock.rStartStop.size(); ++i) {
      BLProfiler::RStartStop &rss = dBlock.rStartStop[i];
      if(rss.rssStart) {     // start region
        if(rStartTime[rss.rssRNumber] < 0.0) {  // not started yet
          rStartTime[rss.rssRNumber] = rss.rssTime;
        } else {                             // already started, mismatched start/stop
        }
      } else {            // stop region
        if(rStartTime[rss.rssRNumber] < 0.0) {  // not started yet, mismatched start/stop
        } else {                             // stopping
          Real rtStart(rStartTime[rss.rssRNumber]), rtStop(rss.rssTime);
          rStartTime[rss.rssRNumber] = -1.0;
          int xStart = int(xLength * rtStart / timeMax);
          int xStop = int(xLength * rtStop / timeMax);
	  xStop = std::min(xStop, xLength - 1);
          int yLo(rss.rssRNumber * yHeight), yHi(((rss.rssRNumber + 1) *  yHeight) - 1);
          Box rBox(IntVect(xStart, yLo), IntVect(xStop, yHi));
	  regionBoxes[rss.rssRNumber].push_back(rBox);
          rFab.setVal<RunOn::Host>(rss.rssRNumber, rBox, 0);
        }
      }
    }
    ClearBlock(dBlock);
  }

  Real timeMin(0.0);
  return TimeRange(timeMin, timeMax);
#endif
}


// ----------------------------------------------------------------------
void RegionsProfStats::FillRegionTimeRanges(Vector<Vector<TimeRange>> &rtr,
                                            int whichProc)
{
  int nRegions(maxRNumber + 1);
  Vector<Real> rStartTime(nRegions, -1.0);

  rtr.resize(nRegions);

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(dBlock.proc != whichProc) {
      continue;
    }
    ReadBlock(dBlock, true, false);        // dont need to read the trace data

    for(int i(0); i < dBlock.rStartStop.size(); ++i) {
      BLProfiler::RStartStop &rss = dBlock.rStartStop[i];
      if(rss.rssStart) {     // start region
        if(rStartTime[rss.rssRNumber] < 0.0) {  // not started yet
          rStartTime[rss.rssRNumber] = rss.rssTime;
        } else {                             // already started, mismatched start/stop
        }
      } else {            // stop region
        if(rStartTime[rss.rssRNumber] < 0.0) {  // not started yet, mismatched start/stop
        } else {                             // stopping
	  rtr[rss.rssRNumber].push_back(TimeRange(rStartTime[rss.rssRNumber], rss.rssTime));
          rStartTime[rss.rssRNumber] = -1.0;
        }
      }
    }
    ClearBlock(dBlock);
  }
}


// ----------------------------------------------------------------------
bool RegionsProfStats::InitRegionTimeRanges(const Box &procBox) {
  BL_PROFILE("RegionsProfStats::InitRegionTimeRanges()");

  int nProcs(ParallelDescriptor::NProcs());
  int myProc(ParallelDescriptor::MyProc());
  int smallY(procBox.smallEnd(YDIR)), bigY(procBox.bigEnd(YDIR));

//#define BL_SERIAL_INIT
#ifdef BL_SERIAL_INIT
  bool bIOP(ParallelDescriptor::IOProcessor());
  if(bIOP) cout << "Starting serial InitRegionTimeRanges." << endl;
  BL_PROFILE_VAR("RegionsProfStats::InitRegionTimeRanges_Serial()", RPSIRTR_S);
  Vector<Vector<Vector<TimeRange> > > checkRegionTimeRanges;  // [proc][rnum][range]
  checkRegionTimeRanges.resize(dataNProcs);
  for(int p(0); p < checkRegionTimeRanges.size(); ++p) {
    checkRegionTimeRanges[p].resize(maxRNumber + 1);
  }

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
      if (persistentStreams){
         ReadBlockNoOpen(dBlock, true, false);  // dont need to read the trace data
      }
      else{
         ReadBlock(dBlock, true, false); // dont need to read the trace data
      }

      for(int i(0); i < dBlock.rStartStop.size(); ++i) {
        BLProfiler::RStartStop &rss = dBlock.rStartStop[i];
        if(rss.rssStart) {     // start region
          checkRegionTimeRanges[dBlock.proc][rss.rssRNumber].push_back(TimeRange(rss.rssTime, -1.0));
        } else {            // stop region
          checkRegionTimeRanges[dBlock.proc][rss.rssRNumber].back().stopTime = rss.rssTime;
        }
      }
      ClearBlock(dBlock);
  }
  BL_PROFILE_VAR_STOP(RPSIRTR_S);
  if(bIOP) cout << "Finished serial InitRegionTimeRanges." << endl;
#endif

  cout << myProc << " InitRegionTimeRanges Box = " << procBox << endl;

  BL_PROFILE_VAR("RegionsProfStats::InitRegionTimeRanges_Parallel()", RPSIRTR_P);
  regionTimeRanges.resize(dataNProcs);
  for(int p(0); p < regionTimeRanges.size(); ++p) {
    regionTimeRanges[p].resize(maxRNumber + 1);
  }

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(dBlock.proc >= smallY && dBlock.proc <= bigY) {    // ---- within myproc range
      if (persistentStreams){
         ReadBlockNoOpen(dBlock, true, false);  // dont need to read the trace data
      }
      else{
         ReadBlock(dBlock, true, false); // dont need to read the trace data
      }


      for(int i(0); i < dBlock.rStartStop.size(); ++i) {
        BLProfiler::RStartStop &rss = dBlock.rStartStop[i];
        if(rss.rssStart) {     // start region
	  if(rss.rssRNumber >= 0) {
            regionTimeRanges[dBlock.proc][rss.rssRNumber].push_back(TimeRange(rss.rssTime, -1.0));
	  }
        } else {            // stop region
	  if(rss.rssRNumber >= 0) {
            regionTimeRanges[dBlock.proc][rss.rssRNumber].back().stopTime = rss.rssTime;
	  }
        }
      }
      ClearBlock(dBlock);
    }
  }

  Vector<int> nRanges(maxRNumber + 1, 0);

  for(int p(0); p < regionTimeRanges.size(); ++p) {
    for(int r(0); r < regionTimeRanges[p].size(); ++r) {
      nRanges[r] = std::max(nRanges[r], static_cast<int> (regionTimeRanges[p][r].size()));
    }
  }
  ParallelDescriptor::ReduceIntMax(nRanges.dataPtr(), nRanges.size());
  long totalRanges(0);
  for(int r(0); r < nRanges.size(); ++r) {
    totalRanges += nRanges[r];
  }

  Vector<int> gSmallY(nProcs, -1);
  Vector<int> gBigY(nProcs, -1);
  gSmallY[myProc] = smallY;
  gBigY[myProc]   = bigY;
  ParallelDescriptor::ReduceIntMax(gSmallY.dataPtr(), gSmallY.size());
  ParallelDescriptor::ReduceIntMax(gBigY.dataPtr(),   gBigY.size());

  // ---- now resize ranges not on this processor and collect them
  for(int p(0); p < regionTimeRanges.size(); ++p) {
    for(int r(0); r < regionTimeRanges[p].size(); ++r) {
      if( ! (p >= smallY && p <= bigY)) {    // ---- not within myproc range
//        if(regionTimeRanges[p][r].size() > 0) {
//	  amrex::Abort("regionTimeRanges size error 0");
//	}
	regionTimeRanges[p][r].resize(nRanges[r]);
      }
    }
  }

  long allRanges = totalRanges * dataNProcs * (maxRNumber + 1) * 2;

  Vector<Real> gAllRanges(allRanges);
  for(int p(0); p < regionTimeRanges.size(); ++p) {
    for(int r(0); r < regionTimeRanges[p].size(); ++r) {
      if(p >= smallY && p <= bigY) {    // ---- within myproc range
        for(int t(0); t < regionTimeRanges[p][r].size(); ++t) {
/*	  int index((p * (maxRNumber + 1) * totalRanges * 2) +
	            (r * totalRanges * 2) + (t * 2)); */
	  long index((r * totalRanges * 2) + (t * 2) +
                     (p * (maxRNumber + 1) * totalRanges * 2));

          gAllRanges[index]     = regionTimeRanges[p][r][t].startTime;
          gAllRanges[index + 1] = regionTimeRanges[p][r][t].stopTime;
        }
      }
    }
  }

  Vector<int> recvDispl(nProcs, 0), recvCounts(nProcs, 0);
  for(int p(0); p < nProcs; ++p) {
    recvCounts[p] = (long(gBigY[p]) - gSmallY[p] + 1) * (maxRNumber + 1) * totalRanges * 2;
    recvDispl[p]  = long(gSmallY[p]) * (maxRNumber + 1) * totalRanges * 2;
  }


#ifdef BL_USE_MPI
  int myStartIndex(gSmallY[myProc] * (maxRNumber + 1) * totalRanges * 2);
  int sendCount((gBigY[myProc] - gSmallY[myProc] + 1) * (maxRNumber + 1) * totalRanges * 2);
  Vector<Real> localGAllRanges(sendCount, 0.0);
  for(int i(0); i < sendCount; ++i) {
    localGAllRanges[i] = gAllRanges[myStartIndex + i];
  }
  MPI_Allgatherv(localGAllRanges.dataPtr(), sendCount, ParallelDescriptor::Mpi_typemap<Real>::type(),
                 gAllRanges.dataPtr(), recvCounts.dataPtr(), recvDispl.dataPtr(),
		 ParallelDescriptor::Mpi_typemap<Real>::type(),
                 ParallelDescriptor::Communicator());
#endif

  for(int p(0); p < regionTimeRanges.size(); ++p) {
    for(int r(0); r < regionTimeRanges[p].size(); ++r) {
      if( ! (p >= smallY && p <= bigY)) {    // ---- not within myproc range
        for(int t(0); t < regionTimeRanges[p][r].size(); ++t) {
/*        int index((p * (maxRNumber + 1) * totalRanges * 2) +
                    (r * totalRanges * 2) + (t * 2));       */
	  long index((r * totalRanges * 2) + (t * 2) +
                     (p * (maxRNumber + 1) * totalRanges * 2));
          regionTimeRanges[p][r][t].startTime = gAllRanges[index];
          regionTimeRanges[p][r][t].stopTime  = gAllRanges[index + 1];
        }
      }
    }
  }

  BL_PROFILE_VAR_STOP(RPSIRTR_P);

#ifdef BL_SERIAL_INIT
  // ---- check vs serial
  for(int p(0); p < regionTimeRanges.size(); ++p) {
    for(int r(0); r < regionTimeRanges[p].size(); ++r) {
      for(int t(0); t < regionTimeRanges[p][r].size(); ++t) {
        if(regionTimeRanges[p][r][t].startTime != checkRegionTimeRanges[p][r][t].startTime) {
	  amrex::Abort("Bad checkRegionTimeRanges startTime");
	}
        if(regionTimeRanges[p][r][t].stopTime != checkRegionTimeRanges[p][r][t].stopTime) {
	  amrex::Abort("Bad checkRegionTimeRanges stopTime");
	}
      }
    }
  }
#endif

  // have to remove the last noRegion
  for(int p(0); p < regionTimeRanges.size(); ++p) {
    if(regionTimeRanges[p][0].size() > 0) {
      if(regionTimeRanges[p][0].back().stopTime < 0.0) {
        regionTimeRanges[p][0].pop_back();
      }
    }
  }
  cout << "----:: regionTimeRanges.size() = " << regionTimeRanges.size() << endl;
  return(maxRNumber > 0);
}


// ----------------------------------------------------------------------
bool RegionsProfStats::Include(const FuncStat &/*fs*/) {
  std::set<int>::iterator it;
  bool binclude(bDefaultInclude);
  return binclude;
}


// ----------------------------------------------------------------------
bool RegionsProfStats::AllCallTimesFAB(FArrayBox &actFab,
                                       const std::string &whichFuncName)
{
#if (BL_SPACEDIM != 2)
  cout << "**** Error:  RegionsProfStats::AllCallTimesFAB only supported in 2D." << endl;
  return false;
#else

  int whichFuncNameInt(-1);
  for(int i(0); i < numbersToFName.size(); ++i) {
    if(numbersToFName[i] == whichFuncName) {
      whichFuncNameInt = i;
    }
  }
  cout << "**** whichFuncName whichFuncNameInt = " << whichFuncName << "  " <<  whichFuncNameInt << endl;
  Vector<Vector<Real> > whichFuncAllTimes(dataNProcs);  // [proc][functime]
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if (persistentStreams){
       ReadBlockNoOpen(dBlock);
    }
    else{
       ReadBlock(dBlock);
    }
    for(int i(0); i < dBlock.vCallStats.size(); ++i) {
      BLProfiler::CallStats &cs = dBlock.vCallStats[i];
      if(cs.csFNameNumber < 0) {  // ---- the unused cs
        continue;
      }
      if(InTimeRange(dBlock.proc, cs.callTime)) {
        int remappedIndex(fnameRemap[dBlock.proc][cs.csFNameNumber]);
        if(remappedIndex == whichFuncNameInt) {
          if(numbersToFName[remappedIndex] != whichFuncName) {
            SHOWVAL(numbersToFName[remappedIndex]);
            amrex::Abort("**** Error 0:  fab whichFuncName");
          }
          whichFuncAllTimes[dBlock.proc].push_back(cs.stackTime);
        }
      }
    }
    ClearBlock(dBlock);
  }

  // make a fab with xdir == processor  ydir == each call time
  // needs to be parallelized and plotfiled
  int whichFuncNCalls(whichFuncAllTimes[0].size());
  bool bSameNCalls(true);
  for(int p(0); p < whichFuncAllTimes.size(); ++p) {
    if(whichFuncAllTimes[p].size() != whichFuncNCalls) {
      cout << "==== bSameNCalls = false" << endl;
      bSameNCalls = false;
    }
  }
  if(bSameNCalls) {
    Box actBox(IntVect(0,0), IntVect(dataNProcs - 1, whichFuncNCalls - 1));
    actFab.resize(actBox, 1);
    actFab.setVal<RunOn::Host>(0.0);
    Real *dptr = actFab.dataPtr(0);
    int nX(actBox.length(XDIR)), nY(actBox.length(YDIR));
    for(int p(0); p < nX; ++p) {
      for(int cnum(0); cnum < nY; ++cnum) {
        int index((cnum * nX) + p);
        dptr[index] = whichFuncAllTimes[p][cnum];
      }
    }
    std::ofstream actfout("whichFuncAllTimes.fab");
    actFab.writeOn(actfout);
    actfout.close();
  }

  return bSameNCalls;
#endif
}


// ----------------------------------------------------------------------
void RegionsProfStats::FillAllCallTimes(Vector<Vector<Real> > &allCallTimes,
                                        const std::string whichFuncName,
					int whichFuncNumber, const Box &procBox)
{
  BL_PROFILE("RegionsProfStats::FillAllCallTimes");

  int smallY(procBox.smallEnd(YDIR)), bigY(procBox.bigEnd(YDIR));
  int proc;

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    proc = dBlock.proc;

    if(proc >= smallY && proc <= bigY) {    // ---- within from proc range
      if (persistentStreams){
        ReadBlockNoOpen(dBlock, false, true);  // ---- only read trace data
      }
      else{
        ReadBlock(dBlock, false, true);  // ------ only read trace data
      }
      for(int i(0); i < dBlock.vCallStats.size(); ++i) {
        BLProfiler::CallStats &cs = dBlock.vCallStats[i];
        if(cs.csFNameNumber < 0) {  // ---- the unused cs
          continue;
        }
        if(InTimeRange(dBlock.proc, cs.callTime)) {
          int remappedIndex(fnameRemap[dBlock.proc][cs.csFNameNumber]);
          if(remappedIndex == whichFuncNumber) {
            if(numbersToFName[remappedIndex] != whichFuncName) {
              SHOWVAL(numbersToFName[remappedIndex]);
              amrex::Abort("**** Error 0:  fab whichFuncName");
            }
            allCallTimes[dBlock.proc].push_back(cs.stackTime);
          }
        }
      }
      ClearBlock(dBlock);
    }
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::CollectFuncStats(Vector<Vector<FuncStat> > &funcStats)
{
  funcStats.resize(numbersToFName.size());  // [fnum][proc]
  for(int n(0); n < funcStats.size(); ++n) {
    funcStats[n].resize(dataNProcs);
  }

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    ReadBlock(dBlock);

    for(int i(0); i < dBlock.vCallStats.size(); ++i) {
      BLProfiler::CallStats &cs = dBlock.vCallStats[i];
      if(cs.csFNameNumber < 0) {  // ---- the unused cs
        continue;
      }
      if(InTimeRange(dBlock.proc, cs.callTime)) {
        int remappedIndex(fnameRemap[dBlock.proc][cs.csFNameNumber]);
        funcStats[remappedIndex][dBlock.proc].totalTime += cs.stackTime;
        funcStats[remappedIndex][dBlock.proc].nCalls  += 1;
      }
    }
    ClearBlock(dBlock);
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::WriteSummary(std::ostream &ios, bool /*bwriteavg*/,
                                    int whichProc, bool graphTopPct)
{
  if( ! ParallelDescriptor::IOProcessor()) {
    return;
  }

  Real timeMin(std::numeric_limits<Real>::max());
  Real timeMax(-std::numeric_limits<Real>::max());
	 
  Vector<std::string> fNames(numbersToFName.size());
  for(int i(0); i < fNames.size(); ++i) {
    if(i >= 0) {
      fNames[i] = numbersToFName[i];
    }
  }

  Vector<BLProfiler::CallStats> vCallStatsAllOneProc;
  Vector<Vector<FuncStat> > funcStats(fNames.size());  // [fnum][proc]
  for(int n(0); n < funcStats.size(); ++n) {
    funcStats[n].resize(dataNProcs);
  }

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    ReadBlock(dBlock);
    if(dBlock.proc == whichProc) {
      for(int i(0); i < dBlock.vCallStats.size(); ++i) {
	// ---- here we have to add only the part of this
	// ---- callstat that intersects the region time range
        BLProfiler::CallStats &cs = dBlock.vCallStats[i];
	TimeRange tRangeFull(cs.callTime, cs.callTime + cs.totalTime);
	std::list<TimeRange> intersectList =
	    RegionsProfStats::RangeIntersection(filterTimeRanges[whichProc], tRangeFull);
	std::list<TimeRange>::iterator tri;
	for(tri = intersectList.begin(); tri != intersectList.end(); ++tri) {
          BLProfiler::CallStats csis(dBlock.vCallStats[i]);
	  csis.callTime  = tri->startTime;
	  csis.totalTime = tri->stopTime - tri->startTime;
          if(InTimeRange(dBlock.proc, cs.callTime)) {
            vCallStatsAllOneProc.push_back(csis);
          }
	}
      }
    }

    for(int i(0); i < dBlock.vCallStats.size(); ++i) {
      BLProfiler::CallStats &cs = dBlock.vCallStats[i];
      if(cs.csFNameNumber < 0) {  // ---- the unused cs
        continue;
      }
      if(InTimeRange(dBlock.proc, cs.callTime)) {
        int remappedIndex(fnameRemap[dBlock.proc][cs.csFNameNumber]);
        funcStats[remappedIndex][dBlock.proc].totalTime += cs.stackTime;
        funcStats[remappedIndex][dBlock.proc].nCalls  += 1;
      }
    }
    timeMin = std::min(timeMin, dBlock.timeMin);
    timeMax = std::max(timeMax, dBlock.timeMax);
    ClearBlock(dBlock);
  }

  Real calcRunTime(timeMax - timeMin);
  Real rangeRunTime(0.0);
  BLProfiler::SetRunTime(calcRunTime);
  for(std::list<TimeRange>::iterator trisum = filterTimeRanges[whichProc].begin();
      trisum != filterTimeRanges[whichProc].end(); ++trisum)
  {
    rangeRunTime += trisum->stopTime - trisum->startTime;
  }
  cout << "+++++++++++++++ calcRunTime rangeRunTime = " << calcRunTime << "  "
       << rangeRunTime << endl;

  std::map<std::string, BLProfiler::ProfStats> mProfStats;  // [fname, pstats]
  CollectMProfStats(mProfStats, funcStats, fNames, calcRunTime, whichProc);

  if(graphTopPct) {
    GraphTopPct(mProfStats, funcStats, fNames, rangeRunTime, dataNProcs, gPercent);
  }

#ifdef BL_TRACE_PROFILING
  MakeFuncPctTimesMF(funcStats, fNames, mProfStats, rangeRunTime, dataNProcs);
#endif

  if(ParallelDescriptor::IOProcessor()) {
    bool writeAvg(true), writeInclusive(true);
    BLProfilerUtils::WriteStats(ios, mProfStats, mFNameNumbersPerProc[whichProc],
                              vCallStatsAllOneProc, writeAvg, writeInclusive);
  }
}
// ----------------------------------------------------------------------
void RegionsProfStats::CheckRegionsData()
{
  bool bIOP(ParallelDescriptor::IOProcessor());
  int myProc(ParallelDescriptor::MyProc());
  cout << myProc << ":  " << "---------------------- checking regions data." << endl;
  SHOWVAL(dataNProcs);
  SHOWVAL(dataBlocks.size());
  SHOWVAL(maxRNumber);
  cout << myProc << ":  " << "----" << endl;

  Vector<Vector<Vector<TimeRange> > > checkRegionTimeRanges;  // [proc][rnum][range]
  Vector<Vector<int>> regionTimeRangesCount;                  // [region][proc] = count
  checkRegionTimeRanges.resize(dataNProcs);
  regionTimeRangesCount.resize(maxRNumber+1);
  for(int p(0); p < checkRegionTimeRanges.size(); ++p) {
    checkRegionTimeRanges[p].resize(maxRNumber + 1);
  }
  for(int p(0); p < regionTimeRangesCount.size(); ++p) {
    regionTimeRangesCount[p].resize(dataNProcs, 0);
  }
 
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(verbose) {
      cout << myProc << ":  " << "RegionsProfProc  " << dBlock.proc << "  nTraceStats  "
           << dBlock.nTraceStats << " nRSS " << dBlock.nRSS << " fileName "
           << dBlock.fileName << "  seekpos  "
	   << dBlock.seekpos << endl;
    }
    ReadBlock(dBlock, true, false); // dont need to read the trace data

    for(int i(0); i < dBlock.rStartStop.size(); ++i) {
      BLProfiler::RStartStop &rss = dBlock.rStartStop[i];
      if(rss.rssRNumber > (maxRNumber + 1))
        if(bIOP)
        {
          cerr << "***RegionsProfStats::CheckRegionsData: region number is greater than max number: "
               << rss.rssRNumber << " > " << maxRNumber + 1 << endl; 
        }
      if(rss.rssStart) {     // start region
        regionTimeRangesCount[rss.rssRNumber][dBlock.proc]++;
        checkRegionTimeRanges[dBlock.proc][rss.rssRNumber].push_back(TimeRange(rss.rssTime, -1.0));
      } else {            // stop region
        regionTimeRangesCount[rss.rssRNumber][dBlock.proc]++;
        checkRegionTimeRanges[dBlock.proc][rss.rssRNumber].back().stopTime = rss.rssTime;
      }
    }
    ClearBlock(dBlock);
  }

  cout << myProc << ":  " << "---------------------- checking regions consistency." << endl;
 
  for (int r(0); r<regionTimeRangesCount.size(); ++r) {
    if (verbose)
    {
      cout << "Region # " << r << " has " << regionTimeRangesCount[r][0] << " piece(s) of time data. " << endl;
    }
    for (int n(0); n<regionTimeRangesCount[r].size(); ++n) {
      if (regionTimeRangesCount[r][0] != regionTimeRangesCount[r][n])
      {
        cerr << "***Region " << r << " was called a different number of times on processor " << n << " : "
             << regionTimeRangesCount[r][0] << " != " << regionTimeRangesCount[r][n] << endl;
      }
    } 
  }

  cout << myProc << ":  " << "---------------------- checking time range consistency." << endl;

  for (int n(0); n<checkRegionTimeRanges.size(); ++n) {
    for (int r(0); r<checkRegionTimeRanges[n].size(); ++r) {
      for (int t(0); t<checkRegionTimeRanges[n][r].size()-1; ++t) {
        if ((verbose) && (n == 0))
        {
          cout << "RTR[" << n << "][" << r << "][" << t << "] = " << checkRegionTimeRanges[n][r][t] << endl;
        }
        if (checkRegionTimeRanges[n][r][t].startTime > checkRegionTimeRanges[n][r][t].stopTime)
        {
          cerr << "***Start time for RTR[" << n << "][" << r << "][" << t << "] is greater than stop time "
               << checkRegionTimeRanges[n][r][t] << endl;
        }
      }
    }
  }
}

// ----------------------------------------------------------------------
void RegionsProfStats::WriteHTML(std::ostream &csHTMLFile,
                                 bool simpleCombine, int whichProc)
{
  BLProfiler::CallStats *combCallStats = 0;

  Vector<std::string> fNumberNames(mFNameNumbersPerProc[whichProc].size());
  for(std::map<std::string, int>::const_iterator it = mFNameNumbersPerProc[whichProc].begin();
      it != mFNameNumbersPerProc[whichProc].end(); ++it)
  {
    fNumberNames[it->second] = it->first;
  }

  Real rangeRunTime(0.0);
  for(std::list<TimeRange>::iterator trisum = filterTimeRanges[whichProc].begin();
      trisum != filterTimeRanges[whichProc].end(); ++trisum)
  {
    rangeRunTime += trisum->stopTime - trisum->startTime;
  }
  cout << "rangeRunTime = " << rangeRunTime << endl;
  Real colorLinkPct(0.05), colorLinkTime(rangeRunTime);
  std::stack<Real> colorLinkTimeStack;
  colorLinkTimeStack.push(rangeRunTime);

  // write to html file
  std::stack<std::string> listEnds;

  csHTMLFile << "<!DOCTYPE html>" << '\n';
  csHTMLFile << "<html>" << '\n';
  csHTMLFile << "<head>" << '\n';
  csHTMLFile << "<title>Call Tree</title>" << '\n';
  csHTMLFile << "</head>" << '\n';
  csHTMLFile << '\n';

  csHTMLFile << "<body>" << '\n';
  csHTMLFile << '\n';
  csHTMLFile << "<script type=\"text/javascript\">" << '\n';
  csHTMLFile << "function collapse(id) {" << '\n';
  csHTMLFile << "  var elem = document.getElementById(id);" << '\n';
  csHTMLFile << "  if(elem.style.display == '') {" << '\n';
  csHTMLFile << "    elem.style.display = 'none';" << '\n';
  csHTMLFile << "  } else {" << '\n';
  csHTMLFile << "    elem.style.display = '';" << '\n';
  csHTMLFile << "  }" << '\n';
  csHTMLFile << "}" << '\n';
  csHTMLFile << "</script>" << '\n';
  csHTMLFile << '\n';

  csHTMLFile << "<h3>Function call times  "
             << "(function number :: function name :: inclusive time :: exclusive time :: 1-e/i % :: ncalls :: callstackdepth :: call time)</h3>"
	     << '\n';

  csHTMLFile << "<ul>" << '\n';
  listEnds.push("</ul>");

// the next two lines will indent the html
//#define IcsHTMLFile for(int id(0); id <= listEnds.size(); ++id) csHTMLFile << "  "; csHTMLFile
//#define IIcsHTMLFile for(int id(0); id < listEnds.size(); ++id) csHTMLFile << "  "; csHTMLFile
#define IcsHTMLFile csHTMLFile
#define IIcsHTMLFile csHTMLFile

  int nBlocks(0);
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(dBlock.proc == whichProc) {
      ++nBlocks;
    }
  }
  std::cout << "************  nBlocks = " << nBlocks << std::endl;
  bool bFirstBlock(true), bLastBlock(false);
  int nodeNumber(-1);
  BLProfiler::CallStats lastFlushedCS;

  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(dBlock.proc == whichProc) {
      if(--nBlocks == 0) {
        bLastBlock = true;
      }
      ReadBlock(dBlock, false, true);  // read only the traces
      Vector<BLProfiler::CallStats> &vCallTrace = dBlock.vCallStats;

      std::cout << "vCallTrace.size() = " << vCallTrace.size() << std::endl;
      if( ! bFirstBlock) {
        vCallTrace[0] = lastFlushedCS;  // copy to the unused cs
      }
      for(int iCT(0); iCT < vCallTrace.size(); ++iCT) {
	++nodeNumber;
        BLProfiler::CallStats &cs = vCallTrace[iCT];
        if(cs.callStackDepth < 0 || cs.csFNameNumber < 0) {  // ---- the unused cs
	  continue;
        }
        if(cs.nCSCalls > 1) {
	  static int count(0);
	  if(count++ < 8) {
            std::cout << "DDDDDDDDDD cs.nCalls = " << cs.nCSCalls << "  "
	              << fNumberNames[cs.csFNameNumber] << std::endl;
	  }
        }

        if(iCT == vCallTrace.size() - 1) {
	  if(combCallStats != 0) {
            IcsHTMLFile << "<li>" << PRINTCS((*combCallStats)) << "</li>" << '\n';
	    delete combCallStats;
	    combCallStats = 0;
	  }
	  if(bLastBlock) {
            IcsHTMLFile << "<li>" << PRINTCS(cs) << "</li>" << '\n';
            for(int n(0); n < cs.callStackDepth; ++n) {
	      if( ! listEnds.empty()) {
                IIcsHTMLFile << listEnds.top() << '\n';
                listEnds.pop();
	      } else {
	        std::cout << "WriteHTML::0:  listEnds.empty():  csd n = "
		          << cs.callStackDepth << "  " << n << std::endl;
	      }
	      if( ! listEnds.empty()) {
                IIcsHTMLFile << listEnds.top() << '\n';
                listEnds.pop();
	      } else {
	        std::cout << "WriteHTML::1:  listEnds.empty():  csd n = "
		          << cs.callStackDepth << "  " << n << std::endl;
	      }
            }
	  } else {  // ---- save the last calltrace
	    lastFlushedCS = cs;
	  }
        } else {
          BLProfiler::CallStats &csNext = vCallTrace[iCT + 1];
          if(csNext.callStackDepth > cs.callStackDepth) {
	    if(combCallStats != 0) {
              IcsHTMLFile << "<li>" << PRINTCS((*combCallStats)) << "</li>" << '\n';
	      delete combCallStats;
	      combCallStats = 0;
	    }
            IcsHTMLFile << "<li>" << '\n';
            listEnds.push("</li>");
	    colorLinkTime = colorLinkTimeStack.top();
	    colorLinkTimeStack.push(cs.totalTime);
	    if(cs.totalTime > colorLinkTime * colorLinkPct) {
              IcsHTMLFile << "<a style=\"color:#800000\" href=\"javascript:void(0)\" onclick=\"collapse('node"
	                  << nodeNumber << "')\">" << PRINTCS(cs) << "</a>" << '\n';
	    } else {
              IcsHTMLFile << "<a href=\"javascript:void(0)\" onclick=\"collapse('node"
	                  << nodeNumber << "')\">" << PRINTCS(cs) << "</a>" << '\n';
	    }
	    if(cs.totalTime > colorLinkTime * colorLinkPct) {  // ---- expand link
              IcsHTMLFile << "<ul id=\"node" << nodeNumber << "\" style=\"display:\">" << '\n';
            } else {
              IcsHTMLFile << "<ul id=\"node" << nodeNumber << "\" style=\"display:none\">" << '\n';
            }
            listEnds.push("</ul>");
          } else  if(csNext.callStackDepth == cs.callStackDepth) {
	    if(simpleCombine) {
	      if(iCT <  vCallTrace.size() - 2 && cs.csFNameNumber == csNext.csFNameNumber) {
	        if(combCallStats == 0) {
		  combCallStats = new BLProfiler::CallStats(cs);
		} else {
		  combCallStats->nCSCalls  += cs.nCSCalls;
		  combCallStats->totalTime += cs.totalTime;
		  combCallStats->stackTime += cs.stackTime;
		}
	      } else {
	        if(combCallStats != 0) {
		  if(cs.csFNameNumber == combCallStats->csFNameNumber) {
		    combCallStats->nCSCalls  += cs.nCSCalls;
		    combCallStats->totalTime += cs.totalTime;
		    combCallStats->stackTime += cs.stackTime;
		  }
                  IcsHTMLFile << "<li>" << PRINTCS((*combCallStats)) << "</li>" << '\n';
	          delete combCallStats;
	          combCallStats = 0;
		} else {
                  IcsHTMLFile << "<li>" << PRINTCS(cs) << "</li>" << '\n';
	        }
	      }
	    } else {
              IcsHTMLFile << "<li>" << PRINTCS(cs) << "</li>" << '\n';
	    }
          } else {
	    if(combCallStats != 0) {
	      if(cs.csFNameNumber == combCallStats->csFNameNumber) {
	        combCallStats->nCSCalls  += cs.nCSCalls;
	        combCallStats->totalTime += cs.totalTime;
	        combCallStats->stackTime += cs.stackTime;
	      }
              IcsHTMLFile << "<li>" << PRINTCS((*combCallStats)) << "</li>" << '\n';
	      delete combCallStats;
	      combCallStats = 0;
	    } else {
              IcsHTMLFile << "<li>" << PRINTCS(cs) << "</li>" << '\n';
	    }
            for(int n(0); n < cs.callStackDepth - csNext.callStackDepth; ++n) {
              IIcsHTMLFile << listEnds.top() << '\n';
              listEnds.pop();
              IIcsHTMLFile << listEnds.top() << '\n';
              listEnds.pop();
	      colorLinkTimeStack.pop();
            }
	    colorLinkTime = colorLinkTimeStack.top();
          }
        }
      }
      bFirstBlock = false;
    }
  }


  if(listEnds.size() != 1) {
    std::cout << "**** Error:  listEnds.size() = " << listEnds.size() << std::endl;
  } else {
    csHTMLFile << listEnds.top() << '\n';
    listEnds.pop();
  }

  csHTMLFile << "</body>" << '\n';
  csHTMLFile << "</html>" << '\n';
}


// ----------------------------------------------------------------------
void RegionsProfStats::CreateVCallStats(CallTreeNode &callTree,
                                        Vector<BLProfiler::CallStats> &vCallStatsNC)
{
  std::map<int, CallTreeNode>::iterator miter;

  if(callTree.fnameNumber >= 0) {
    vCallStatsNC.push_back(BLProfiler::CallStats(callTree.stackDepth,
                                                 callTree.fnameNumber, callTree.nCalls,
                                                 callTree.totalTime, callTree.stackTime, 0.0));
  }
  for(miter = callTree.calledFunctions.begin(); miter != callTree.calledFunctions.end(); ++miter) {
    CreateVCallStats(miter->second, vCallStatsNC);
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::PrintCallTreeNode(CallTreeNode &callTree,
                                         Vector<std::string> &fNumberNames)
{
  std::map<int, CallTreeNode>::iterator miter;

  if(callTree.fnameNumber >= 0) {
    std::cout << "PCTN:  " << fNumberNames[callTree.fnameNumber]
              << "  stackDepth = " << callTree.stackDepth
              << "  nCalls = " << callTree.nCalls << "  stackTime = "
	      << callTree.stackTime << std::endl;
  }
  for(miter = callTree.calledFunctions.begin(); miter != callTree.calledFunctions.end(); ++miter) {
    PrintCallTreeNode(miter->second, fNumberNames);
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::WriteHTMLNC(std::ostream &csHTMLFile, int whichProc)
{
  Vector<std::string> fNumberNames(mFNameNumbersPerProc[whichProc].size());
  for(std::map<std::string, int>::const_iterator it = mFNameNumbersPerProc[whichProc].begin();
      it != mFNameNumbersPerProc[whichProc].end(); ++it)
  {
    fNumberNames[it->second] = it->first;
  }

  // write to html file
  std::stack<std::string> listEnds;

  csHTMLFile << "<!DOCTYPE html>" << '\n';
  csHTMLFile << "<html>" << '\n';
  csHTMLFile << "<head>" << '\n';
  csHTMLFile << "<title>Call Tree</title>" << '\n';
  csHTMLFile << "</head>" << '\n';
  csHTMLFile << '\n';

  csHTMLFile << "<body>" << '\n';
  csHTMLFile << '\n';
  csHTMLFile << "<script type=\"text/javascript\">" << '\n';
  csHTMLFile << "function collapse(id) {" << '\n';
  csHTMLFile << "  var elem = document.getElementById(id);" << '\n';
  csHTMLFile << "  if(elem.style.display == '') {" << '\n';
  csHTMLFile << "    elem.style.display = 'none';" << '\n';
  csHTMLFile << "  } else {" << '\n';
  csHTMLFile << "    elem.style.display = '';" << '\n';
  csHTMLFile << "  }" << '\n';
  csHTMLFile << "}" << '\n';
  csHTMLFile << "</script>" << '\n';
  csHTMLFile << '\n';

  csHTMLFile << "<h3>Function calls "
             << "(function number :: function name :: inclusive time :: exclusive time :: 1-e/i % :: ncalls :: callstackdepth)</h3>"
	     << '\n';

  csHTMLFile << "<ul>" << '\n';
  listEnds.push("</ul>");

// the next two lines will indent the html
//#define IcsHTMLFile for(int id(0); id <= listEnds.size(); ++id) csHTMLFile << "  "; csHTMLFile
//#define IIcsHTMLFile for(int id(0); id < listEnds.size(); ++id) csHTMLFile << "  "; csHTMLFile
#define IcsHTMLFile csHTMLFile
#define IIcsHTMLFile csHTMLFile

  CallTreeNode callTree;
  callTree.stackDepth = -1;
  CallTreeNode *currentCTN(&callTree);
  std::map<int, CallTreeNode>::iterator miter;
  std::stack<CallTreeNode *> ctStack;
  ctStack.push(currentCTN);

  int nBlocks(0);
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(dBlock.proc == whichProc) {
      ++nBlocks;
    }
  }
  std::cout << "************  nBlocks = " << nBlocks << std::endl;
  bool bFirstBlock(true), bLastBlock(false), bAddFirstNode(true);
  BLProfiler::CallStats lastFlushedCS;

  int totalCalls(0);
  
  for(int idb(0); idb < dataBlocks.size(); ++idb) {
    DataBlock &dBlock = dataBlocks[idb];
    if(dBlock.proc == whichProc) {
      if(--nBlocks == 0) {
        bLastBlock = true;
      }
      ReadBlock(dBlock, false, true);  // read only the traces
      Vector<BLProfiler::CallStats> &vCallTrace = dBlock.vCallStats;

      std::cout << "vCallTrace.size() = " << vCallTrace.size() << std::endl;
      if( ! bFirstBlock) {
	vCallTrace[0] = lastFlushedCS;  // copy to the unused cs
      }
      for(int iCT(0); iCT < vCallTrace.size(); ++iCT) {
        BLProfiler::CallStats &cs = vCallTrace[iCT];
        if(cs.callStackDepth < 0 || cs.csFNameNumber < 0) {  // ---- the unused cs
	  continue;
        }

	if(bAddFirstNode) {
	// ---- add the node
	std::pair<std::map<int, CallTreeNode>::iterator, bool> retval;
	retval = currentCTN->calledFunctions.insert(
	               std::pair<int, CallTreeNode>(cs.csFNameNumber, CallTreeNode()));
	miter = retval.first;
	if(retval.second) {  // ---- value inserted
	  miter->second.fnameNumber = cs.csFNameNumber;
	  miter->second.stackDepth  = cs.callStackDepth;
	  miter->second.nCalls = 1;
	  miter->second.totalTime = cs.totalTime;
	  miter->second.stackTime = cs.stackTime;
	} else {      // ---- value existed
	  miter->second.nCalls += 1;
	  miter->second.totalTime += cs.totalTime;
	  miter->second.stackTime += cs.stackTime;
	}
	++totalCalls;
	}
	bAddFirstNode = true;
	    
        if(iCT == vCallTrace.size() - 1) {
	  if(bLastBlock) {
	  } else {
	    lastFlushedCS = cs;
	    bAddFirstNode = false;
	  }
	} else {

          BLProfiler::CallStats &csNext = vCallTrace[iCT + 1];

          if(csNext.callStackDepth > cs.callStackDepth) {
	    currentCTN = &(miter->second);
	    if(currentCTN != 0) {
	      ctStack.push(currentCTN);
	    }
	  } else  if(csNext.callStackDepth == cs.callStackDepth) {
	    // ---- node added above
	  } else {
	    for(int s(0); s < cs.callStackDepth  - csNext.callStackDepth; ++s) {
	      ctStack.pop();
	    }
	    currentCTN = ctStack.top();
	  }
	}
      }
      bFirstBlock = false;
    }
  }

  cout << "++++++++++ totalCalls = " << totalCalls << endl;

  Vector<BLProfiler::CallStats> vCallStatsNC;
  CreateVCallStats(callTree, vCallStatsNC);

      Vector<BLProfiler::CallStats> &vCallTrace = vCallStatsNC;
      int nodeNumber(-1);

      for(int iCT(0); iCT < vCallTrace.size(); ++iCT) {
	++nodeNumber;
        BLProfiler::CallStats &cs = vCallTrace[iCT];
        if(cs.callStackDepth < 0 || cs.csFNameNumber < 0) {  // ---- the unused cs
	  continue;
        }
        if(iCT == vCallTrace.size() - 1) {
          IcsHTMLFile << "<li>" << PRINTCSNC(cs) << "</li>" << '\n';
          for(int n(0); n < cs.callStackDepth; ++n) {
            if( ! listEnds.empty()) {
              IIcsHTMLFile << listEnds.top() << '\n';
              listEnds.pop();
            } else {
              std::cout << "WriteHTMLNC::0:  listEnds.empty():  csd n = "
	                << cs.callStackDepth << "  " << n << std::endl;
            }
            if( ! listEnds.empty()) {
              IIcsHTMLFile << listEnds.top() << '\n';
              listEnds.pop();
            } else {
              std::cout << "WriteHTMLNC::1:  listEnds.empty():  csd n = "
	                << cs.callStackDepth << "  " << n << std::endl;
            }
          }
        } else {
          BLProfiler::CallStats &csNext = vCallTrace[iCT + 1];
          if(csNext.callStackDepth > cs.callStackDepth) {
            IcsHTMLFile << "<li>" << '\n';
            listEnds.push("</li>");
            IcsHTMLFile << "<a href=\"javascript:void(0)\" onclick=\"collapse('node"
	                << nodeNumber << "')\">" << PRINTCSNC(cs) << "</a>" << '\n';
            if(cs.callStackDepth < 100) {
              IcsHTMLFile << "<ul id=\"node" << nodeNumber << "\" style=\"display:\">" << '\n';
            } else {
              IcsHTMLFile << "<ul id=\"node" << nodeNumber << "\" style=\"display:none\">" << '\n';
            }
            listEnds.push("</ul>");
          } else  if(csNext.callStackDepth == cs.callStackDepth) {
            IcsHTMLFile << "<li>" << PRINTCSNC(cs) << "</li>" << '\n';
          } else {
            IcsHTMLFile << "<li>" << PRINTCSNC(cs) << "</li>" << '\n';
            for(int n(0); n < cs.callStackDepth - csNext.callStackDepth; ++n) {
              if( ! listEnds.empty()) {
                IIcsHTMLFile << listEnds.top() << '\n';
                listEnds.pop();
              } else {
                std::cout << "WriteHTMLNC::2:  listEnds.empty():  csd n = "
		          << cs.callStackDepth << "  " << n << std::endl;
              }
              if( ! listEnds.empty()) {
                IIcsHTMLFile << listEnds.top() << '\n';
                listEnds.pop();
              } else {
                std::cout << "WriteHTMLNC::3:  listEnds.empty():  csd n = "
		          << cs.callStackDepth << "  " << n << std::endl;
              }
            }
          }
        }
      }


  if(listEnds.size() != 1) {
    std::cout << "**** Error:  listEnds.size() = " << listEnds.size() << std::endl;
  } else {
    csHTMLFile << listEnds.top() << '\n';
    listEnds.pop();
  }
  

  csHTMLFile << "</body>" << '\n';
  csHTMLFile << "</html>" << '\n';
}


#define PRINTCSTT(IOS, CS, DEL) for(int indent(0); indent < CS.callStackDepth; ++indent) \
    { IOS << "---|"; } IOS << DEL << CS.csFNameNumber \
    << DEL << fNumberNames[CS.csFNameNumber] << DEL << CS.totalTime << DEL << CS.stackTime \
    << DEL << CS.nCSCalls  << DEL << CS.callStackDepth << DEL << CS.callTime << '\n';

// ----------------------------------------------------------------------
void RegionsProfStats::WriteTextTrace(std::ostream &ios, bool simpleCombine,
                                      int whichProc, std::string delimString)
{
  Vector<std::string> fNumberNames(mFNameNumbersPerProc[whichProc].size());
  for(std::map<std::string, int>::const_iterator it = mFNameNumbersPerProc[whichProc].begin();
      it != mFNameNumbersPerProc[whichProc].end(); ++it)
  {
    fNumberNames[it->second] = it->first;
  }

  ios << "Text Trace  (function number :: function name :: inclusive time :: exclusive time :: ncalls :: callstackdepth :: call time)\n\n";
  ios << std::setprecision(16);

  if(simpleCombine) {
    BLProfiler::CallStats *combCallStats = 0;

    for(int idb(0); idb < dataBlocks.size(); ++idb) {
      DataBlock &dBlock = dataBlocks[idb];
      if(dBlock.proc == whichProc) {
        ReadBlock(dBlock, false, true);  // read only the traces
        Vector<BLProfiler::CallStats> &vCallTrace = dBlock.vCallStats;

        for(int i(0); i < vCallTrace.size(); ++i) {
          BLProfiler::CallStats &cs = vCallTrace[i];
	  if(cs.callStackDepth < 0 || cs.csFNameNumber < 0) {  // ---- the unused cs
	    continue;
	  }

	  if(i < vCallTrace.size() - 2) {
            BLProfiler::CallStats &csNext = vCallTrace[i + 1];
	    if(csNext.callStackDepth == cs.callStackDepth &&
	       cs.csFNameNumber == csNext.csFNameNumber)
	    {
              if(combCallStats == 0) {
                combCallStats = new BLProfiler::CallStats(cs);
              } else {
                combCallStats->nCSCalls  += cs.nCSCalls;
                combCallStats->totalTime += cs.totalTime;
                combCallStats->stackTime += cs.stackTime;
              }
	    } else {
              if(combCallStats == 0) {
	        PRINTCSTT(ios, cs, delimString);
	      } else {
                combCallStats->nCSCalls  += cs.nCSCalls;
                combCallStats->totalTime += cs.totalTime;
                combCallStats->stackTime += cs.stackTime;
	        PRINTCSTT(ios, (*combCallStats), delimString);
                delete combCallStats;
                combCallStats = 0;
              }
	    }
	  } else {
            if(combCallStats != 0) {
	      PRINTCSTT(ios, (*combCallStats), delimString);
              delete combCallStats;
              combCallStats = 0;
            }
	    PRINTCSTT(ios, cs, delimString);
	  }
        }
        ClearBlock(dBlock);
      }
    }

  } else {
    for(int idb(0); idb < dataBlocks.size(); ++idb) {
      DataBlock &dBlock = dataBlocks[idb];
      if(dBlock.proc == whichProc) {
        ReadBlock(dBlock, false, true);  // read only the traces
        Vector<BLProfiler::CallStats> &vCallTrace = dBlock.vCallStats;

        for(int i(0); i < vCallTrace.size(); ++i) {
          BLProfiler::CallStats &cs = vCallTrace[i];
	  if(cs.callStackDepth < 0 || cs.csFNameNumber < 0) {
	    continue;
	  }
	  PRINTCSTT(ios, cs, delimString);
        }
        ClearBlock(dBlock);
      }
    }
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::InitDataFileNames(const Vector<std::string> &hfn) {
  for(int i(0); i < hfn.size(); ++i) {
    std::string dFileName(hfn[i]);
    dFileName.replace(dFileName.find("_H_"), 3, "_D_");
    regDataFileNames.insert(std::pair<std::string, int>(dFileName, i));
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::InitCStatsDataBlock(int proc, long nrss, long ntracestats,
                                           const std::string &filename, long seekpos)
{
  int streamindex;
  std::map<std::string, int>::iterator it =  regDataFileNames.find(filename);
  if(it == regDataFileNames.end()) {
    streamindex = regDataFileNames.size();
    regDataFileNames.insert(std::pair<std::string, int>(filename, streamindex));
  } else {
    streamindex = it->second;
  }
  if(bInitDataBlocks) {
    currentProc = proc;
    dataBlocks.push_back(DataBlock(proc, nrss, ntracestats, filename, seekpos, streamindex));
    currentDataBlock = dataBlocks.size() - 1;
  }
}


// ----------------------------------------------------------------------
void RegionsProfStats::OpenAllStreams(const std::string &dirname) {
  BL_PROFILE_VAR("RegionsProfStats::OpenAllStreams", regsopenallstreams);

  regDataStreams.resize(regDataFileNames.size());
  int dsIndex(0);
#define BL_CYCLEOPENS
#ifdef BL_CYCLEOPENS
  int  myProc(ParallelDescriptor::MyProc());
  int nNames(regDataFileNames.size());
  Vector<std::string> aFullFileNames(nNames);
  for(std::map<std::string, int>::iterator it = regDataFileNames.begin();
      it != regDataFileNames.end(); ++it)
  {
    aFullFileNames[dsIndex] = (dirname + '/' + it->first);
    ++dsIndex;
  }
  int index;
  for(int s(0); s < nNames; ++s) {
    index = (s + myProc) % nNames;
    regDataStreams[index] = new std::ifstream(aFullFileNames[index].c_str());

    if (regDataStreams[index]->fail())
    {
      cout << "****regDataStreams failed. Continuing without persistent streams." << std::endl;
      persistentStreams = false;
      CloseAllStreams();
      break;
    }
 }
#else
  for(std::map<std::string, int>::iterator it = regDataFileNames.begin();
      it != regDataFileNames.end(); ++it)
  {
    std::string fullFileName(dirname + '/' + it->first);
    regDataStreams[dsIndex] = new std::ifstream(fullFileName.c_str());

    if (regDataStreams[dsIndex]->fail())
    {
      cout << "****regDataStreams failed. Continuing without persistent streams." << std::endl;
      persistentStreams = false;
      CloseAllStreams();
      break;
    }
    ++dsIndex;
 }
#endif
  BL_PROFILE_VAR_STOP(regsopenallstreams);
}


// ----------------------------------------------------------------------
void RegionsProfStats::CloseAllStreams() {
  BL_PROFILE_VAR("RegionProfStats::CloseAllStreams", regsclosellstreams);
  for(int i(0); i < regDataStreams.size(); ++i) {
    if (regDataStreams[i] != nullptr)
    {
      if (regDataStreams[i]->is_open())
      {
        regDataStreams[i]->close();
      }
      delete regDataStreams[i];
      regDataStreams[i] = nullptr;
    }
  }
  BL_PROFILE_VAR_STOP(regsclosellstreams);
}


// ----------------------------------------------------------------------
void RegionsProfStats::ReadBlock(DataBlock &dBlock, bool readRSS,
                                 bool readTraces)
{
  BL_PROFILE("RegionsProfStats::ReadBlock()");
  if(dBlock.rStartStop.size() != dBlock.nRSS) {
    dBlock.rStartStop.resize(dBlock.nRSS);
  }
  if(dBlock.vCallStats.size() != dBlock.nTraceStats) {
    dBlock.vCallStats.resize(dBlock.nTraceStats);
  }

  std::string fullFileName(dirName + '/' + dBlock.fileName);
  BL_PROFILE_VAR("OpenStream", openstream);
  std::ifstream instr(fullFileName.c_str());
  BL_PROFILE_VAR_STOP(openstream);

  if(readRSS) {  // ---- read the rstarts and rstops
    long dataSize(dBlock.rStartStop.size() * sizeof(BLProfiler::RStartStop));
    instr.seekg(dBlock.seekpos);
    instr.read(reinterpret_cast<char *>(dBlock.rStartStop.dataPtr()), dataSize);
    for(int i(0); i < dBlock.rStartStop.size(); ++i) {
      minRegionTime = std::min(minRegionTime, dBlock.rStartStop[i].rssTime);
      maxRegionTime = std::max(maxRegionTime, dBlock.rStartStop[i].rssTime);
    }
  } else {
    long dataSize(dBlock.rStartStop.size() * sizeof(BLProfiler::RStartStop));
    instr.seekg(dBlock.seekpos + dataSize);
  }

  if(readTraces) {  // ---- read the trace data
    long dataSize(dBlock.vCallStats.size() * sizeof(BLProfiler::CallStats));
    if(dataSize > 0) {
      instr.read(reinterpret_cast<char *>(dBlock.vCallStats.dataPtr()), dataSize);
    }
  }

  instr.close();
}


// ----------------------------------------------------------------------
void RegionsProfStats::ReadBlockNoOpen(DataBlock &dBlock, bool readRSS,
                                       bool readTraces)
{
  BL_PROFILE("RegionsProfStats::ReadBlockNoOpen()");
  if(dBlock.rStartStop.size() != dBlock.nRSS) {
    dBlock.rStartStop.resize(dBlock.nRSS);
  }
  if(dBlock.vCallStats.size() != dBlock.nTraceStats) {
    dBlock.vCallStats.resize(dBlock.nTraceStats);
  }

  std::ifstream *instr = regDataStreams[dBlock.streamIndex];

  if(readRSS) {  // ---- read the rstarts and rstops
    long dataSize(dBlock.rStartStop.size() * sizeof(BLProfiler::RStartStop));
    instr->seekg(dBlock.seekpos);
    instr->read(reinterpret_cast<char *>(dBlock.rStartStop.dataPtr()), dataSize);
    for(int i(0); i < dBlock.rStartStop.size(); ++i) {
      minRegionTime = std::min(minRegionTime, dBlock.rStartStop[i].rssTime);
      maxRegionTime = std::max(maxRegionTime, dBlock.rStartStop[i].rssTime);
    }
  } else {
    long dataSize(dBlock.rStartStop.size() * sizeof(BLProfiler::RStartStop));
    instr->seekg(dBlock.seekpos + dataSize);
  }

  if(readTraces) {  // ---- read the trace data
    long dataSize(dBlock.vCallStats.size() * sizeof(BLProfiler::CallStats));
    instr->read(reinterpret_cast<char *>(dBlock.vCallStats.dataPtr()), dataSize);
  }
}


// ----------------------------------------------------------------------
bool RegionsProfStats::ReadBlock(DataBlock &/*dBlock*/, const int /*nmessages*/) {
amrex::Abort("not implemented yet.");
return false;
/*
  int leftToRead(dBlock.size - dBlock.readoffset);
  int readSize(std::min(leftToRead, nmessages));
  int readPos(dBlock.seekpos + dBlock.readoffset * csSize);
  //cout << "************" << endl;
  if(dBlock.vCommStats.size() != readSize) {
    dBlock.vCommStats.resize(readSize);
  }
  std::string fullFileName(dirName + '/' + dBlock.fileName);

  std::ifstream instr(fullFileName.c_str());
  int dataSize(readSize * csSize);
  instr.seekg(readPos);
  instr->read(reinterpret_cast<char *>(dBlock.vCommStats.dataPtr()), dataSize);
  instr.close();

  dBlock.readoffset += readSize;
  //SHOWVAL(dBlock.readoffset);
  //cout << "************" << endl;
  return(dBlock.readoffset < dBlock.size);
*/
}


// ----------------------------------------------------------------------
void RegionsProfStats::ClearBlock(DataBlock &dBlock) {
  dBlock.rStartStop.clear();
  Vector<BLProfiler::RStartStop>().swap(dBlock.rStartStop);
  dBlock.vCallStats.clear();
  Vector<BLProfiler::CallStats>().swap(dBlock.vCallStats);
}


// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
