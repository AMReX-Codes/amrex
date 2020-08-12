// ----------------------------------------------------------------------
//  BLProfStats.cpp
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

using std::cout;
using std::endl;
using std::flush;
using std::string;
using std::ifstream;
using std::map;
using std::unordered_map;
using std::unordered_multimap;
using std::vector;
using std::pair;

#include <AMReX_BLProfStats.H>
#include <AMReX_RegionsProfStats.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfUtilities.H>

using namespace amrex;

#ifdef _OPENMP
#include <omp.h>
#endif

#define SHOWVAL(val) { cout << #val << " = " << val << endl; }

int BLProfStats::verbose(-1);
int BLProfStats::blProfVersion(-1);
int BLProfStats::dataNProcs(-1);
int BLProfStats::nOutFiles(-1);
std::string BLProfStats::dirName;
bool BLProfStats::bInitDataBlocks(true);
std::map<std::string, int> BLProfStats::blpDataFileNames;
Real BLProfStats::gPercent(0.10);
Vector<std::ifstream *> BLProfStats::blpDataStreams;
bool BLProfStats::bTimeRangeInitialized(false);

extern std::string SanitizeName(const std::string &s);
extern void PrintTimeRangeList(const std::list<RegionsProfStats::TimeRange> &trList);
extern long amrex::FileSize(const std::string &filename);
extern void amrex::MakeFuncPctTimesMF(const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                               const Vector<std::string> &blpFNames,
			       const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
			       Real runTime, int dataNProcs);
extern void amrex::CollectMProfStats(std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                              const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                              const Vector<std::string> &fNames,
                              Real runTime, int whichProc);
void amrex::GraphTopPct(const std::map<std::string, BLProfiler::ProfStats> &mProfStats,
                 const Vector<Vector<BLProfStats::FuncStat> > &funcStats,
                 const Vector<std::string> &fNames,
                 Real runTime, int dataNProcs, Real gPercent);



// ----------------------------------------------------------------------
bool BLProfStats::AddPiece(std::list<TimeRange> &addToHere,
                           const TimeRange &pieceToAdd)
{
  // ---- add the new piece, sort, then remove overlap
  addToHere.push_back(pieceToAdd);
  addToHere.sort(BLProfStats::TimeRangeCompare());

  std::list<TimeRange>::iterator it, itnext;
  std::list<std::list<TimeRange>::iterator> overlaps;
  for(it = addToHere.begin(); it != addToHere.end(); ++it) {
    itnext = it;
    ++itnext;
    if(itnext != addToHere.end()) {
      TimeRange &tRange = *it;
      TimeRange &tRangeNext = *itnext;
      if(tRange.stopTime >= tRangeNext.startTime) {
        tRangeNext.startTime = tRange.startTime;
        tRangeNext.stopTime  = std::max(tRange.stopTime, tRangeNext.stopTime);
	overlaps.push_back(it);
      }
    }
  }
  std::list<std::list<TimeRange>::iterator >::iterator oit;
  for(oit = overlaps.begin(); oit != overlaps.end(); ++oit) {
    addToHere.erase(*oit);
  }
  return true;
}


// ----------------------------------------------------------------------
std::list<BLProfStats::TimeRange> BLProfStats::RangeIntersection(
                  std::list<TimeRange> &rangeList,
                  const TimeRange &pieceToIntersect)
{
  std::list<TimeRange> intersectList, intersectListComp;
  intersectList.push_back(pieceToIntersect);
  intersectListComp.push_back(pieceToIntersect);

  std::list<TimeRange>::iterator it;
  for(it = rangeList.begin(); it != rangeList.end(); ++it) {
    RemovePiece(intersectListComp, *it);
  }
  for(it = intersectListComp.begin(); it != intersectListComp.end(); ++it) {
    RemovePiece(intersectList, *it);
  }
  return intersectList;
}


// ----------------------------------------------------------------------
bool BLProfStats::RemovePiece(std::list<TimeRange> &removeFromHere,
                              const TimeRange &pieceToRemove)
{
  

  bool piecesRemoved(false);
  std::list<TimeRange>::iterator it;
  std::list<std::list<TimeRange>::iterator> eraseThese;

  for(it = removeFromHere.begin(); it != removeFromHere.end(); ++it) {
    TimeRange &tRangeFrom = *it;
    bool bothLow(pieceToRemove.stopTime  < tRangeFrom.startTime);
    bool loInRange((pieceToRemove.startTime > tRangeFrom.startTime) &&
                   (pieceToRemove.startTime < tRangeFrom.stopTime));
    bool hiInRange((pieceToRemove.stopTime  > tRangeFrom.startTime)  &&
                   (pieceToRemove.stopTime  < tRangeFrom.stopTime));
    bool covered((pieceToRemove.startTime <= tRangeFrom.startTime) &&
                 (pieceToRemove.stopTime  >= tRangeFrom.stopTime));
    bool bothHigh(pieceToRemove.startTime > tRangeFrom.stopTime);
    // Warning! covered test may require an epsilon in future
    //    implementations. Works fine now given fixed nature of
    //    presently used data.

    // ---- there are six cases
    if(bothLow) {                            // ---- do nothing
    } else if( ! loInRange && hiInRange) {   // ---- remove low end piece
      tRangeFrom.startTime = pieceToRemove.stopTime;
      piecesRemoved = true;

    } else if(loInRange && hiInRange) {      // ---- remove middle piece
      removeFromHere.insert(it,
        TimeRange(tRangeFrom.startTime, pieceToRemove.startTime));
      tRangeFrom.startTime = pieceToRemove.stopTime;  // reuse this part
      piecesRemoved = true;

    } else if(loInRange && ! hiInRange) {    // ---- remove high end piece
      tRangeFrom.stopTime = pieceToRemove.startTime;
      piecesRemoved = true;

    } else if(covered) {                     // ---- remove the whole range
      eraseThese.push_back(it);
      piecesRemoved = true;

    } else if(bothHigh) {                    // ---- do nothing
    }
    
  }
  std::list<std::list<TimeRange>::iterator >::iterator eit;
  for(eit = eraseThese.begin(); eit != eraseThese.end(); ++eit) {
    removeFromHere.erase(*eit);
  }

  return piecesRemoved;
}


// ----------------------------------------------------------------------
BLProfStats::BLProfStats() {
  maxRNumber = -1;
  currentProc = -1;
  minRegionTime =  std::numeric_limits<Real>::max();
  maxRegionTime = -std::numeric_limits<Real>::max();
  bDefaultInclude = true;
  currentDataBlock = 0;
}


// ----------------------------------------------------------------------
BLProfStats::~BLProfStats() {
}


// ----------------------------------------------------------------------
void BLProfStats::AddRegionName(const std::string &rname, int rnumber) {
  regionNames.insert(std::pair<std::string, int>(rname, rnumber));
  regionNumbers.insert(std::pair<int, std::string>(rnumber, rname));
  maxRNumber = std::max(maxRNumber, rnumber);
}


// ----------------------------------------------------------------------
std::set<int> BLProfStats::WhichRegions(int proc, Real t) {
  std::set<int> whichRegions;
  if(proc < 0 || proc >= regionTimeRanges.size()) {
    return whichRegions;
  }
  Vector<Vector<TimeRange> > &rtr = regionTimeRanges[proc];
  for(int rnum(0); rnum < rtr.size(); ++rnum) {
    for(int range(0); range < rtr[rnum].size(); ++range) {
      TimeRange &tr = rtr[rnum][range];
      if(tr.Contains(t)) {
        whichRegions.insert(rnum);
      }
    }
  }
  return whichRegions;
}


// ----------------------------------------------------------------------
void BLProfStats::MakeFilterFile(const std::string &ffname) {
  std::ofstream filterFile(ffname.c_str(), std::ios::out | std::ios::trunc);
  if( ! filterFile.good()) {
    cout << "**** Error in BLProfStats::MakeFilterFile:  cannot open file:  "
         << ffname << endl;
  } else {
    filterFile << "__IncludeAll__" << '\n';
    filterFile << "//__IncludeNone__" << '\n';
    std::map<int, std::string>::iterator it;
    for(it = regionNumbers.begin(); it != regionNumbers.end(); ++it) {
      filterFile << ": " << it->second << ' ' << it->first << '\n';
    }
    filterFile.close();
  }
}


// ----------------------------------------------------------------------
void BLProfStats::SetFilter(BLProfStats::FilterStatus fs,
                            const std::string &/*ffname*/, int rnumber)
{
  if(fs == BLProfStats::FilterStatus::ON) {
    includeSet.insert(rnumber);
  } else if(fs == BLProfStats::FilterStatus::OFF) {
    excludeSet.insert(rnumber);
  }
}


// ----------------------------------------------------------------------
void BLProfStats::SetFilter(BLProfStats::FilterStatus fs) {
  if(fs == BLProfStats::FilterStatus::INCLUDEALL) {
    bDefaultInclude = true;
  } else if(fs == BLProfStats::FilterStatus::INCLUDENONE) {
    bDefaultInclude = false;
  } else {
    cout << "**** Error in BLProfStats::SetFilter(fs):  bad fs:   " << fs << endl;
  }
}


// ----------------------------------------------------------------------
void BLProfStats::SetFilterTimeRanges(const Vector<std::list<TimeRange> > &ftr)
{
  filterTimeRanges = ftr;
}


// ----------------------------------------------------------------------
void BLProfStats::InitFilterTimeRanges() {
  BL_PROFILE("BLProfStats::InitFilterTimeRanges()");

  // includeSet  excludeSet
  filterTimeRanges.resize(dataNProcs);
  for(int p(0); p < filterTimeRanges.size(); ++p) {
    filterTimeRanges[p].clear();
  }

  if(bDefaultInclude) {  // ---- include all
    for(int p(0); p < filterTimeRanges.size(); ++p) {
      filterTimeRanges[p].push_back(TimeRange(minRegionTime, maxRegionTime));
    }
  } else {               // ---- include none
  }

  std::set<int>::iterator iit;
  for(iit = includeSet.begin(); iit != includeSet.end(); ++iit) {
    int regNum(*iit);
    for(int p(0); p < regionTimeRanges.size(); ++p) {
      for(int r(0); r < regionTimeRanges[p][regNum].size(); ++r) {
        TimeRange &trange = regionTimeRanges[p][regNum][r];
	//cout << "_here include:  " << trange << endl;
        AddPiece(filterTimeRanges[p], trange);
      }
    }
  }

  for(iit = excludeSet.begin(); iit != excludeSet.end(); ++iit) {
    int regNum(*iit);
    for(int p(0); p < regionTimeRanges.size(); ++p) {
      for(int r(0); r < regionTimeRanges[p][regNum].size(); ++r) {
        TimeRange &trange = regionTimeRanges[p][regNum][r];
	//cout << "_here exclude:  " << trange << endl;
        RemovePiece(filterTimeRanges[p], trange);
      }
    }
  }
  bTimeRangeInitialized = true;
}


// ----------------------------------------------------------------------
bool BLProfStats::InTimeRange(int proc, Real calltime) {
//  if( ! bTimeRangeInitialized) {
//    return true;
//  }
  if(filterTimeRanges.empty()) {
#if 0
  static int count(0);
  if(count++ < 4) {
    cout << "**** BLProfStats::InTimeRange:  init true but range empty." << endl;
  }
#endif
    return true;
  }
  std::list<TimeRange>::iterator iit;
  for(iit = filterTimeRanges[proc].begin(); iit != filterTimeRanges[proc].end(); ++iit) {
    TimeRange &trange = (*iit);
    if(trange.startTime <= calltime && trange.stopTime >= calltime) {
      return true;
    }
  }
  return false;
}


// ----------------------------------------------------------------------
std::ostream &operator<< (std::ostream &os, const BLProfStats::TimeRange &tr)
{
  os << '[' << tr.startTime << ",  " << tr.stopTime << ']';
  return os;
}


// ----------------------------------------------------------------------
void BLProfStats::AddFunctionName(const std::string &fname) {
  std::size_t found;
  std::string fnameNQ(fname.substr(1, fname.length() - 2));  // ---- remove quotes
  while((found = fnameNQ.find(" ")) != std::string::npos) {  // ---- replace spaces
    fnameNQ.replace(found, 1, "_");
  }
  blpFNames.push_back(fnameNQ);
}


// ----------------------------------------------------------------------
void BLProfStats::InitBLProfDataBlock(const int proc, const std::string &filename,
                                      const long seekpos)
{
  int streamindex;
  std::map<std::string, int>::iterator it =  blpDataFileNames.find(filename);
  if(it == blpDataFileNames.end()) {
    streamindex = blpDataFileNames.size();
    blpDataFileNames.insert(std::pair<std::string, int>(filename, streamindex));
  } else {
    streamindex = it->second;
  }
  if(bInitDataBlocks) {
    currentProc = proc;
    blpDataBlocks.push_back(BLPDataBlock(proc, filename, seekpos, streamindex));
    currentDataBlock = blpDataBlocks.size() - 1;
  }
}


// ----------------------------------------------------------------------
void BLProfStats::OpenAllStreams(const std::string &dirname) {
  BL_PROFILE_VAR("BLProfStats::OpenAllStreams", blpsopenallstreams);
  blpDataStreams.resize(blpDataFileNames.size());
  int dsIndex(0);
  for(std::map<std::string, int>::iterator it = blpDataFileNames.begin();
      it != blpDataFileNames.end(); ++it)
  {
    std::string fullFileName(dirname + '/' + it->first);
    blpDataStreams[dsIndex] = new std::ifstream(fullFileName.c_str());
    ++dsIndex;
  }
  BL_PROFILE_VAR_STOP(blpsopenallstreams);
}


// ----------------------------------------------------------------------
void BLProfStats::CloseAllStreams() {
  BL_PROFILE_VAR("BLProfStats::CloseAllStreams", blpsclosellstreams);
  for(int i(0); i < blpDataStreams.size(); ++i) {
    blpDataStreams[i]->close();
    delete blpDataStreams[i];
  }
  BL_PROFILE_VAR_STOP(blpsclosellstreams);
}


// ----------------------------------------------------------------------
void BLProfStats::CheckData() {
  if( ! ParallelDescriptor::IOProcessor()) {
    return;
  }
  int nFuncs(blpFNames.size());
  long dataSizeNC(nFuncs * sizeof(long));
  long dataSizeTT(nFuncs * sizeof(Real));
  long dataSize(dataSizeNC + dataSizeTT);
  std::map<std::string, long> maxSeekPerFile;

  for(int idb(0); idb < blpDataBlocks.size(); ++idb) {
    BLPDataBlock &dBlock = blpDataBlocks[idb];
    long dbSP(dBlock.seekpos);
    std::string fullFileName(dirName + '/' + dBlock.fileName);
    maxSeekPerFile[fullFileName] = std::max(maxSeekPerFile[fullFileName], dbSP);
  }

  bool bad(false);
  long diff(-1);
  std::map<std::string, long>::iterator it;
  for(it = maxSeekPerFile.begin(); it != maxSeekPerFile.end(); ++it) {
    diff = it->second + dataSize - FileSize(it->first); 
    if(diff != 0) {
      bad = true;
      cout << "fName maxSeek fSize diff = " << it->first << "  " << it->second << "  "
           << FileSize(it->first) << "  " << diff << endl;
    }
  }
  if(bad) {
    cerr << "**** Error in BLProfStats::CheckData:  data and file size mismatch." << endl;
  }

}


// ----------------------------------------------------------------------
void BLProfStats::ReadBlock(BLPDataBlock &dBlock) {
  int nFuncs(blpFNames.size());
  if(dBlock.nCalls.size() != nFuncs) {
    dBlock.nCalls.resize(nFuncs);
  }
  if(dBlock.totalTime.size() != nFuncs) {
    dBlock.totalTime.resize(nFuncs);
  }

  std::string fullFileName(dirName + '/' + dBlock.fileName);
  std::ifstream instr(fullFileName.c_str());

  long dataSizeNC(nFuncs * sizeof(long));
  long dataSizeTT(nFuncs * sizeof(Real));

  instr.seekg(dBlock.seekpos);

  instr.read(reinterpret_cast<char *>(dBlock.nCalls.dataPtr()), dataSizeNC);
  instr.read(reinterpret_cast<char *>(dBlock.totalTime.dataPtr()), dataSizeTT);

  if(instr.fail()) {
    cout << "******** Error::ReadBlock:  instr failed: fsize seekpos = "
         << FileSize(fullFileName) << "  " << instr.tellg() << endl;
  }

  instr.close();
}


// ----------------------------------------------------------------------
void BLProfStats::ReadBlockNoOpen(BLPDataBlock &/*dBlock*/)
{
amrex::Abort("not implemented yet.");
}


// ----------------------------------------------------------------------
void BLProfStats::ClearBlock(BLPDataBlock &dBlock) {
  dBlock.nCalls.clear();
  Vector<long>().swap(dBlock.nCalls);
  dBlock.totalTime.clear();
  Vector<Real>().swap(dBlock.totalTime);
}


// ----------------------------------------------------------------------
void BLProfStats::CollectFuncStats(Vector<Vector<FuncStat> > &funcStats)
{
  funcStats.resize(blpFNames.size());  // [fnum][proc]
  for(int n(0); n < funcStats.size(); ++n) {
    funcStats[n].resize(dataNProcs);
  }

  for(int idb(0); idb < blpDataBlocks.size(); ++idb) {
    BLPDataBlock &dBlock = blpDataBlocks[idb];
    ReadBlock(dBlock);
    Vector<long> &nc = dBlock.nCalls;
    Vector<Real> &tt = dBlock.totalTime;

    for(int fnum(0); fnum < blpFNames.size(); ++fnum) {
      FuncStat &fs = funcStats[fnum][dBlock.proc];
      fs.nCalls = nc[fnum];
      fs.totalTime = tt[fnum];
    }
    ClearBlock(dBlock);
  }
}


// ----------------------------------------------------------------------
void BLProfStats::WriteSummary(std::ostream &ios, bool /*bwriteavg*/,
                               int whichProc, bool graphTopPct)
{
  if( ! ParallelDescriptor::IOProcessor()) {
    return;
  }

  Vector<Vector<FuncStat> > funcStats;
  CollectFuncStats(funcStats);

  Real calcRunTime(calcEndTime);
  BLProfiler::SetRunTime(calcRunTime);

  std::map<std::string, BLProfiler::ProfStats> mProfStats;  // [fname, pstats]
  CollectMProfStats(mProfStats, funcStats, blpFNames, calcRunTime, whichProc);

  if(graphTopPct) {
    GraphTopPct(mProfStats, funcStats, blpFNames, calcRunTime, dataNProcs, gPercent);
  }

#ifdef BL_TRACE_PROFILING
  MakeFuncPctTimesMF(funcStats, blpFNames, mProfStats, calcRunTime, dataNProcs);
#endif

  if(ParallelDescriptor::IOProcessor()) {
    std::map<std::string, int> fnameNumbers;
    Vector<BLProfiler::CallStats> vCallStatsAllOneProc;
    bool writeAvg(true), writeInclusive(false);
    BLProfilerUtils::WriteStats(ios, mProfStats, fnameNumbers,
                                vCallStatsAllOneProc, writeAvg, writeInclusive);
  }
}


// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
