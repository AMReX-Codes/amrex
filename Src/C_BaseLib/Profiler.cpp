#ifdef BL_PROFILING

#include <Profiler.H>
#include <REAL.H>
#include <Utility.H>
#include <ParallelDescriptor.H>
#include <Array.H>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

bool Profiler::bWriteAll = true;
bool Profiler::bWriteFabs = true;
bool Profiler::bWriteBLT = false;
int Profiler::currentStep = 0;
Real Profiler::pctTimeLimit = 5.0;
Real Profiler::calcRunTime = 0.0;
std::stack<Real> Profiler::nestedTimeStack;
std::map<int, Real> Profiler::mStepMap;
std::map<std::string, Profiler::ProfStats> Profiler::mProfStats;
std::map<Real, std::string, std::greater<Real> > Profiler::mTimersTotalsSorted;
std::vector<Profiler::CommStats> Profiler::vCommStats;
std::map<std::string, Profiler::CommFuncType> Profiler::CommStats::cftNames;
std::set<Profiler::CommFuncType> Profiler::CommStats::cftExclude;
int Profiler::CommStats::iBarrierNumber = 0;
std::vector<std::string> Profiler::CommStats::vBarrierNames;


Profiler::Profiler(const std::string &funcname)
    : bltstart(0.0), bltelapsed(0.0)
    , fname(funcname)
    , bRunning(false)
{
    start();
}


Profiler::~Profiler() {
  if(bRunning) {
    stop();
  }
}


void Profiler::Initialize() {
    CommStats::cftNames["InvalidCFT"]     = InvalidCFT;
    CommStats::cftNames["AsendTsii"]      = AsendTsii;
    CommStats::cftNames["AsendTsiiM"]     = AsendTsiiM;
    CommStats::cftNames["AsendvTii"]      = AsendvTii;
    CommStats::cftNames["SendTsii"]       = SendTsii;
    CommStats::cftNames["SendvTii"]       = SendvTii;
    CommStats::cftNames["ArecvTsiiM"]     = ArecvTsiiM;
    CommStats::cftNames["ArecvTii"]       = ArecvTii;
    CommStats::cftNames["ArecvvTii"]      = ArecvvTii;
    CommStats::cftNames["RecvTsii"]       = RecvTsii;
    CommStats::cftNames["RecvvTii"]       = RecvvTii;
    CommStats::cftNames["ReduceT"]        = ReduceT; 
    CommStats::cftNames["BCastTsi"]       = BCastTsi; 
    CommStats::cftNames["GatherTsT1Si"]   = GatherTsT1Si; 
    CommStats::cftNames["GatherTi"]       = GatherTi; 
    CommStats::cftNames["ScatterTsT1si"]  = ScatterTsT1si; 
    CommStats::cftNames["Barrier"]        = Barrier;

  // check for exclude file
  std::string exFile("CommFuncExclude.txt");
  std::vector<CommFuncType> vEx;


  Array<char> fileCharPtr;
  bool bExitOnError(false);  // in case the file does not exist
  ParallelDescriptor::ReadAndBcastFile(exFile, fileCharPtr, bExitOnError);

  if(fileCharPtr.size() > 0) {
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream cfex(fileCharPtrString, std::istringstream::in);

    while( ! cfex.eof()) {
        std::string cft;
        cfex >> cft;
        if( ! cfex.eof()) {
	  vEx.push_back(CommStats::StringToCFT(cft));
	}
    }
    for(int i(0); i < vEx.size(); ++i) {
      CommStats::cftExclude.insert(vEx[i]);
    }

  }
}


void Profiler::start() {
  ++mProfStats[fname].nCalls;
  if(ParallelDescriptor::IOProcessor() && bWriteBLT) {
    std::cout << "BLTStart " << fname << std::endl;
  }
  bRunning = true;
  bltstart = ParallelDescriptor::second();
  nestedTimeStack.push(0.0);
}

  
void Profiler::stop() {
  bltelapsed += ParallelDescriptor::second() - bltstart;
  bRunning = false;
  if(ParallelDescriptor::IOProcessor() && bWriteBLT) {
    std::cout << "BLTEnd " << fname << " time = " << bltelapsed << std::endl;
  }
  Real thisFuncTime(bltelapsed);
  if( ! nestedTimeStack.empty()) {
    thisFuncTime -= nestedTimeStack.top();
    nestedTimeStack.pop();
  }
  if( ! nestedTimeStack.empty()) {
    nestedTimeStack.top() += bltelapsed;
  }
  mProfStats[fname].totalTime += thisFuncTime;
}


void Profiler::InitParams(const Real ptl, const bool writeall, const bool writefabs,
                          const bool writeblt)
{
  pctTimeLimit = ptl;
  bWriteAll = writeall;
  bWriteFabs = writefabs;
  bWriteBLT = writeblt;
}


void Profiler::AddStep(const int snum) {
  currentStep = snum;
  mStepMap.insert(std::map<int, Real>::value_type(currentStep,
                                                  ParallelDescriptor::second()));
}


void Profiler::Finalize() {
  // --------------------------------------- gather global stats
  //
  // need to check that the set of profiled functions is the same
  // on all processors
  //
  Real finalizeStart = ParallelDescriptor::second();  // time the timer
  const int nProcs(ParallelDescriptor::NProcs());
  const int myProc(ParallelDescriptor::MyProc());
  const int iopNum(ParallelDescriptor::IOProcessorNumber());

  // filter out profiler communications.
  CommStats::cftExclude.insert(BCastTsi);
  CommStats::cftExclude.insert(GatherTsT1Si);

  int maxlen(0);
  Array<Real> gtimes(1);
  Array<long> ncalls(1);
  if(ParallelDescriptor::IOProcessor()) {
    gtimes.resize(nProcs);
    ncalls.resize(nProcs);
  }
			      
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    std::string profName(it->first);
    int pnLen(profName.size());
    maxlen = std::max(maxlen, pnLen);
    char cpn[profName.size() + 1];
    strcpy(cpn, profName.c_str());
    cpn[profName.size()] = '\0';

    ParallelDescriptor::Bcast(cpn, profName.size() + 1);
    ProfStats &pstats = mProfStats[profName];
    if(nProcs == 1) {
      gtimes[0] = pstats.totalTime;
      ncalls[0] = pstats.nCalls;
    } else {
      ParallelDescriptor::Gather(&pstats.totalTime, 1, gtimes.dataPtr(), 1, iopNum);
      ParallelDescriptor::Gather(&pstats.nCalls, 1, ncalls.dataPtr(), 1, iopNum);
    }
    Real tsum(0.0), tmin(gtimes[0]), tmax(gtimes[0]), tavg(0.0);
    long ncsum(0);
    if(ParallelDescriptor::IOProcessor()) {
      for(int i(0); i < gtimes.size(); ++i) {
        tsum += gtimes[i];
        tmin = std::min(tmin, gtimes[i]);
        tmax = std::max(tmax, gtimes[i]);
      }
      tavg = tsum / static_cast<Real> (gtimes.size());
      pstats.minTime = tmin;
      pstats.maxTime = tmax;
      pstats.avgTime = tavg;

      for(int i(0); i < ncalls.size(); ++i) {
        ncsum += ncalls[i];
      }
      // uncomment for reporting total calls summed over all procs
      //pstats.nCalls = ncsum;
    }
  }

  // --------------------------------------- print global stats to cout
  if(ParallelDescriptor::IOProcessor()) {
    bool bWriteAvg(true);
    if(nProcs == 1) {
      bWriteAvg = false;
    }
    WriteStats(std::cout, bWriteAvg);
  }


  // --------------------------------------- print all procs stats to a file
  if(bWriteAll) {
    std::string outfile("bl_prof.txt");

    // need to do nfiles here
    for(int iproc(0); iproc < nProcs; ++iproc) {  // serialize
      if(myProc == iproc) {
        std::ofstream outfilestr;
        if(iproc == 0) {
          outfilestr.open(outfile.c_str(), std::ios::out|std::ios::trunc);
        } else {
          outfilestr.open(outfile.c_str(), std::ios::out|std::ios::app);
          outfilestr.seekp(0, std::ios::end);
        }
        if( ! outfilestr.good()) {
          BoxLib::FileOpenFailed(outfile);
        }
        WriteStats(outfilestr);

        outfilestr.flush();
        outfilestr.close();

        int iBuff(0), wakeUpPID(myProc + 1), tag(0);
        if(wakeUpPID < nProcs) {
          ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
      }
      if(myProc == (iproc + 1)) {  // next proc waits
        int iBuff(0), waitForPID(myProc - 1), tag(0);
        ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
      }
    }
    ParallelDescriptor::Barrier("Profiler::Finalize");
  }
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Profiler::Finalize():  time:  "   // time the timer
              << ParallelDescriptor::second() - finalizeStart << std::endl;
  }
}


void Profiler::WriteStats(std::ostream &ios, bool bwriteavg) {
  //const int nProcs(ParallelDescriptor::NProcs());
  const int myProc(ParallelDescriptor::MyProc());
  //const int iopNum(ParallelDescriptor::IOProcessorNumber());
  const int colWidth(10);

  mTimersTotalsSorted.clear();

  Real totalTimers(0.0), percent(0.0);
  int maxlen(0);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    std::string profName(it->first);
    int pnLen(profName.size());
    maxlen = std::max(maxlen, pnLen);

    if(bwriteavg) {
      totalTimers += it->second.avgTime;
    } else {
      totalTimers += it->second.totalTime;
    }
  }
  Real pTimeTotal(totalTimers);
  if(calcRunTime > 0.0 && bwriteavg == false) {
    pTimeTotal = calcRunTime;
  }

  ios << std::endl;
  ios << std::endl;
  if( ! bwriteavg) {
    ios << std::setfill('*')
        << std::setw(maxlen + 2 + 3 * (colWidth + 2) - (colWidth+12)) << "";
    ios << std::setfill(' ');
    ios << "  Processor:  " << std::setw(colWidth) << myProc << std::endl;
  }

  // -------- write timers sorted by name
  WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    if(bwriteavg) {
      percent = 100.0 * (it->second.avgTime / pTimeTotal);
    } else {
      percent = 100.0 * (it->second.totalTime / pTimeTotal);
    }
    std::string fname(it->first);
    const ProfStats &pstats = it->second;
    WriteRow(ios, fname, pstats, percent, colWidth, maxlen, bwriteavg);
  }
  ios << std::endl;
  ios << "Total Timers     = " << std::setw(colWidth) << totalTimers
      << " seconds." << std::endl;
  if(calcRunTime > 0.0) {
    percent = 100.0 * totalTimers / calcRunTime;
    ios << "Calc Run Time    = " << std::setw(colWidth) << calcRunTime
        << " seconds." << std::endl;
    ios << "Percent Coverage = " << std::setw(colWidth) << percent << " %" << std::endl;
  }

  // -------- write timers sorted by percent
  ios << std::endl;
  ios << std::endl;
  WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    double dsec;
    if(bwriteavg) {
      dsec = it->second.avgTime;
    } else {
      dsec = it->second.totalTime;
    }
    std::string sfir(it->first);
    mTimersTotalsSorted.insert(std::make_pair(dsec, sfir));
  }

  for(std::map<double, std::string>::const_iterator it = mTimersTotalsSorted.begin();
      it != mTimersTotalsSorted.end(); ++it)
  {
    percent = 100.0 * (it->first / pTimeTotal);
    std::string fname(it->second);
    const ProfStats &pstats = mProfStats[fname];
    WriteRow(ios, fname, pstats, percent, colWidth, maxlen, bwriteavg);
  }
  if(bwriteavg) {
    ios << std::setfill('=') << std::setw(maxlen+4 + 5 * (colWidth+2)) << ""
        << std::endl;
  } else {
    ios << std::setfill('=') << std::setw(maxlen+4 + 3 * (colWidth+2)) << ""
        << std::endl;
  }
  ios << std::setfill(' ');
  ios << std::endl;



  ios << "%%%%%%%%%%%%%%%%%%%%%" << std::endl;
  ios << "vCommStats.size() = " << vCommStats.size() << std::endl;
  ios << "sizeof(vCommStats[0]) = " << sizeof(vCommStats[0]) << std::endl;
  for(int i(0); i < vCommStats.size(); ++i) {
    CommStats &cs = vCommStats[i];
    if(cs.cfType == Barrier) {
      ios << cs.timeStamp << "  " << CommStats::CFTToString(cs.cfType)
          << "  iBarrierNumber = " << cs.size
	  << " vBarrierName = " << CommStats::vBarrierNames[cs.size] << std::endl;
    } else {
      ios << cs.timeStamp << "  " << CommStats::CFTToString(cs.cfType)
          << "  " << cs.dest << "  " << cs.size << std::endl;
    }
  }

}


void Profiler::WriteHeader(std::ostream &ios, const int colWidth,
                           const Real maxlen, const bool bwriteavg)
{
  if(bwriteavg) {
    ios << std::setfill('-') << std::setw(maxlen+4 + 5 * (colWidth+2))
        << std::left << "Total times " << std::endl;
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Min"
        << std::setw(colWidth + 2) << "Avg"
        << std::setw(colWidth + 2) << "Max"
        << std::setw(colWidth + 4) << "Percent %"
        << std::endl;
  } else {
    ios << std::setfill('-') << std::setw(maxlen+4 + 3 * (colWidth+2))
        << std::left << "Total times " << std::endl;
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Time"
        << std::setw(colWidth + 4) << "Percent %"
        << std::endl;
  }
}

void Profiler::WriteRow(std::ostream &ios, const std::string &fname,
                        const ProfStats &pstats, const Real percent,
			const int colWidth, const Real maxlen,
			const bool bwriteavg)
{
    int numPrec(4), pctPrec(2);
    if(bwriteavg) {
      ios << std::right;
      ios << std::setw(maxlen + 2) << fname << "  "
          << std::setw(colWidth) << pstats.nCalls << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.minTime << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.avgTime << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.maxTime << "  "
          << std::setprecision(pctPrec) << std::fixed << std::setw(colWidth)
	  << percent << " %" << std::endl;
    } else {
      ios << std::setw(maxlen + 2) << fname << "  "
          << std::setw(colWidth) << pstats.nCalls << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.totalTime << "  "
          << std::setprecision(pctPrec) << std::fixed << std::setw(colWidth)
	  << percent << " %" << std::endl;
    }
}


void Profiler::AddCommStat(CommFuncType cft, int dest, int size) {
  std::set<CommFuncType>::iterator cfti = CommStats::cftExclude.find(cft);
  if(cfti == CommStats::cftExclude.end()) {
    CommStats cs(cft, ParallelDescriptor::MyProc(), dest, size,
                ParallelDescriptor::second());
    vCommStats.push_back(cs);
  }
}


void Profiler::AddBarrier(CommFuncType cft, std::string message) {
  CommStats cs(cft, ParallelDescriptor::MyProc(), 0, Profiler::CommStats::iBarrierNumber,
              ParallelDescriptor::second());
  vCommStats.push_back(cs);
  CommStats::vBarrierNames.resize(CommStats::iBarrierNumber + 1);
  CommStats::vBarrierNames[CommStats::iBarrierNumber] = message;
  ++CommStats::iBarrierNumber;
}


std::string Profiler::CommStats::CFTToString(CommFuncType cft) {
  switch(cft) {
    case InvalidCFT:     return "InvalidCFT";
    case AsendTsii:      return "AsendTsii";
    case AsendTsiiM:     return "AsendTsiiM";
    case AsendvTii:      return "AsendvTii";
    case SendTsii:       return "SendTsii";
    case SendvTii:       return "SendvTii";
    case ArecvTsiiM:     return "ArecvTsiiM";
    case ArecvTii:       return "ArecvTii";
    case ArecvvTii:      return "ArecvvTii";
    case RecvTsii:       return "RecvTsii";
    case RecvvTii:       return "RecvvTii";
    case ReduceT:        return "ReduceT"; 
    case BCastTsi:       return "BCastTsi"; 
    case GatherTsT1Si:   return "GatherTsT1Si"; 
    case GatherTi:       return "GatherTi"; 
    case ScatterTsT1si:  return "ScatterTsT1si"; 
    case Barrier:        return "Barrier";
  }
  return "*** Error: Bad CommFuncType.";
}


Profiler::CommFuncType Profiler::CommStats::StringToCFT(const std::string &s) {
  return CommStats::cftNames[s];
}



#else

#endif
