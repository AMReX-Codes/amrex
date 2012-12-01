#ifdef BL_PROFILING

#include <Profiler.H>
#include <REAL.H>
#include <Utility.H>
#include <iostream>
#include <iomanip>
#include <fstream>
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
    ParallelDescriptor::Barrier();
  }
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Profiler::Finalize():  time:  "   // time the timer
              << ParallelDescriptor::second() - finalizeStart << std::endl;
  }
}


void Profiler::WriteStats(std::ostream &ios, bool bwriteavg) {
  const int nProcs(ParallelDescriptor::NProcs());
  const int myProc(ParallelDescriptor::MyProc());
  const int iopNum(ParallelDescriptor::IOProcessorNumber());
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
  if(calcRunTime > 0.0) {
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

  WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    percent = 100.0 * (it->second.totalTime / pTimeTotal);
    std::string fname(it->first);
    ProfStats pstats = it->second;
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

  ios << std::endl;
  ios << std::endl;
  WriteHeader(ios, colWidth, maxlen, bwriteavg);
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    double dsec = it->second.totalTime;
    std::string sfir = it->first;
    mTimersTotalsSorted.insert(std::make_pair(dsec, sfir));
  }

  for(std::map<double, std::string>::const_iterator it = mTimersTotalsSorted.begin();
      it != mTimersTotalsSorted.end(); ++it)
  {
    percent = 100.0 * (it->first / pTimeTotal);
    std::string fname(it->second);
    ProfStats pstats(mProfStats[fname]);
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


#else

#endif
