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
bool Profiler::bFirstCommWriteH = true;  // header
bool Profiler::bFirstCommWriteD = true;  // data
bool Profiler::bInitialized = false;

int Profiler::currentStep = 0;
int Profiler::csFlushSize = 2000000;
int Profiler::nProfFiles  = 64;

Real Profiler::pctTimeLimit = 5.0;
Real Profiler::calcRunTime  = 0.0;
Real Profiler::startTime    = 0.0;
Real Profiler::timerTime    = 0.0;

std::stack<Real> Profiler::nestedTimeStack;
std::map<int, Real> Profiler::mStepMap;
std::map<std::string, Profiler::ProfStats> Profiler::mProfStats;
std::map<Real, std::string, std::greater<Real> > Profiler::mTimersTotalsSorted;
std::vector<Profiler::CommStats> Profiler::vCommStats;
std::map<std::string, Profiler::CommFuncType> Profiler::CommStats::cftNames;
std::set<Profiler::CommFuncType> Profiler::CommStats::cftExclude;
int Profiler::CommStats::barrierNumber = 0;
int Profiler::CommStats::reductionNumber = 0;
std::vector<std::pair<std::string,int> > Profiler::CommStats::barrierNames;
std::vector<std::pair<int,int> > Profiler::CommStats::nameTags;
std::vector<std::string> Profiler::CommStats::nameTagNames;
std::vector<int> Profiler::CommStats::reductions;



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
  if(bInitialized) {
    return;
  }

  Real t0, t1;
  int nTimerTimes(1000);
  for(int i(0); i < nTimerTimes; ++i) {  // ---- time the timer
    t0 = ParallelDescriptor::second();
    t1 = ParallelDescriptor::second();
    timerTime += t1 - t0;
  }
  timerTime /= static_cast<Real> (nTimerTimes);

  startTime = ParallelDescriptor::second();

#ifdef BL_COMM_PROFILING
  vCommStats.reserve(csFlushSize);
  // can probably bypass the rest of this function
#endif

  CommStats::cftExclude.insert(AllCFTypes);  // temporarily

  CommStats::cftNames["InvalidCFT"]     = InvalidCFT;
  CommStats::cftNames["AllReduceT"]     = AllReduceT; 
  CommStats::cftNames["AllReduceR"]     = AllReduceR; 
  CommStats::cftNames["AllReduceL"]     = AllReduceL; 
  CommStats::cftNames["AllReduceI"]     = AllReduceI; 
  CommStats::cftNames["AsendTsii"]      = AsendTsii;
  CommStats::cftNames["AsendTsiiM"]     = AsendTsiiM;
  CommStats::cftNames["AsendvTii"]      = AsendvTii;
  CommStats::cftNames["SendTsii"]       = SendTsii;
  CommStats::cftNames["SendvTii"]       = SendvTii;
  CommStats::cftNames["ArecvTsii"]      = ArecvTsii;
  CommStats::cftNames["ArecvTsiiM"]     = ArecvTsiiM;
  CommStats::cftNames["ArecvTii"]       = ArecvTii;
  CommStats::cftNames["ArecvvTii"]      = ArecvvTii;
  CommStats::cftNames["RecvTsii"]       = RecvTsii;
  CommStats::cftNames["RecvvTii"]       = RecvvTii;
  CommStats::cftNames["ReduceT"]        = ReduceT; 
  CommStats::cftNames["ReduceR"]        = ReduceR; 
  CommStats::cftNames["ReduceL"]        = ReduceL; 
  CommStats::cftNames["ReduceI"]        = ReduceI; 
  CommStats::cftNames["BCastTsi"]       = BCastTsi; 
  CommStats::cftNames["GatherTsT1Si"]   = GatherTsT1Si; 
  CommStats::cftNames["GatherTi"]       = GatherTi; 
  CommStats::cftNames["GatherRiRi"]     = GatherRiRi; 
  CommStats::cftNames["ScatterTsT1si"]  = ScatterTsT1si; 
  CommStats::cftNames["Barrier"]        = Barrier;
  CommStats::cftNames["Waitsome"]       = Waitsome;
  CommStats::cftNames["NameTag"]        = NameTag;
  CommStats::cftNames["AllCFTypes"]     = AllCFTypes;
  CommStats::cftNames["NoCFTypes"]      = NoCFTypes;

  // check for exclude file
  std::string exFile("CommFuncExclude.txt");
  std::vector<CommFuncType> vEx;

  Array<char> fileCharPtr;
  bool bExitOnError(false);  // in case the file does not exist
  ParallelDescriptor::ReadAndBcastFile(exFile, fileCharPtr, bExitOnError);

  CommStats::cftExclude.erase(AllCFTypes);

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
  bInitialized = true;
}


void Profiler::start() {
#ifdef _OPENMP
#pragma omp master
#endif
{
  ++mProfStats[fname].nCalls;
  bRunning = true;
  bltstart = ParallelDescriptor::second();
  nestedTimeStack.push(0.0);
}
}

  
void Profiler::stop() {
#ifdef _OPENMP
#pragma omp master
#endif
{
  bltelapsed += ParallelDescriptor::second() - bltstart;
  bRunning = false;
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
}


void Profiler::InitParams(const Real ptl, const bool writeall, const bool writefabs) {
  pctTimeLimit = ptl;
  bWriteAll = writeall;
  bWriteFabs = writefabs;
}


void Profiler::AddStep(const int snum) {
  currentStep = snum;
  mStepMap.insert(std::map<int, Real>::value_type(currentStep,
                                                  ParallelDescriptor::second()));
}


void Profiler::Finalize() {
  if( ! bInitialized) {
    return;
  }
  // --------------------------------------- gather global stats
  Real finalizeStart = ParallelDescriptor::second();  // time the timer
  const int nProcs(ParallelDescriptor::NProcs());
  const int myProc(ParallelDescriptor::MyProc());
  const int iopNum(ParallelDescriptor::IOProcessorNumber());

  // filter out profiler communications.
  CommStats::cftExclude.insert(AllCFTypes);

  int maxlen(0);
  Array<Real> gtimes(1);
  Array<long> ncalls(1);
  if(ParallelDescriptor::IOProcessor()) {
    gtimes.resize(nProcs);
    ncalls.resize(nProcs);
  }


  // -------- make sure the set of profiled functions is the same on all processors
  int pfStringsSize(0);
  std::ostringstream pfStrings;
  if(ParallelDescriptor::IOProcessor()) {
    for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
        it != mProfStats.end(); ++it)
    {
      pfStrings << it->first << '\n';
    }
    pfStrings << std::ends;
    pfStringsSize = pfStrings.str().size();
  }
  ParallelDescriptor::Bcast(&pfStringsSize, 1);

  char *pfChar = new char[pfStringsSize];
  if(ParallelDescriptor::IOProcessor()) {
    std::strcpy(pfChar, pfStrings.str().c_str());
  }
  ParallelDescriptor::Bcast(pfChar, pfStringsSize);

  if( ! ParallelDescriptor::IOProcessor()) {
    std::istringstream pfIn(pfChar);
    std::string pfName;
    while( ! pfIn.eof()) {
      pfIn >> pfName;
      if( ! pfIn.eof()) {
        std::map<std::string, ProfStats>::const_iterator it = mProfStats.find(pfName);
        if(it == mProfStats.end()) {
	  ProfStats ps;
          mProfStats.insert(std::pair<std::string, ProfStats>(pfName, ps));
	  //std::cout << myProc << ":  #### ProfName not found, inserting:  "
	            //<< pfName << std::endl;
        }
      }
    }
  }

  // ------- we really need to send names that are not on the ioproc
  // ------- to the ioproc but for now we will punt
  // ------- make a copy in case there are names not on the ioproc
  std::map<std::string, ProfStats> mProfStatsCopy;
  std::istringstream pfIn(pfChar);
  std::string pfName;
  while( ! pfIn.eof()) {
    pfIn >> pfName;
    if( ! pfIn.eof()) {
      std::map<std::string, ProfStats>::const_iterator it = mProfStats.find(pfName);
      if(it != mProfStats.end()) {
        mProfStatsCopy.insert(std::pair<std::string, ProfStats>(it->first, it->second));
      } else {
	std::cout << myProc << ":  #### Unknown ProfName:  "
	          << it->first << std::endl;
      }
    }
  }
  delete [] pfChar;

  // ------- now check for names not on the ioproc
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStats.begin();
      it != mProfStats.end(); ++it)
  {
    std::map<std::string, ProfStats>::const_iterator itcopy =
                                         mProfStatsCopy.find(it->first);
    //if(itcopy == mProfStatsCopy.end()) {
      //std::cout << myProc << ":  #### ProfName not on ioproc:  "
                //<< it->first << std::endl;
    //}
  }


  // ---------------------------------- now collect global data onto the ioproc
  for(std::map<std::string, ProfStats>::const_iterator it = mProfStatsCopy.begin();
      it != mProfStatsCopy.end(); ++it)
  {
    std::string profName(it->first);
    int pnLen(profName.size());
    maxlen = std::max(maxlen, pnLen);
    char cpn[profName.size() + 1];
    strcpy(cpn, profName.c_str());
    cpn[profName.size()] = '\0';

    ParallelDescriptor::Bcast(cpn, profName.size() + 1);
    ProfStats &pstats = mProfStatsCopy[profName];
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
      ProfStats &pstats = mProfStats[profName];  // not the copy on ioproc
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
    // --------------------- start nfiles block
    std::string cdir("bl_prof");
    if(ParallelDescriptor::IOProcessor()) {
      if( ! BoxLib::UtilCreateDirectory(cdir, 0755)) {
        BoxLib::CreateDirectoryFailed(cdir);
      }
    }
    // Force other processors to wait until directory is built.
    ParallelDescriptor::Barrier("Profiler::Finalize::waitfordir");

    const int   myProc    = ParallelDescriptor::MyProc();
    const int   nProcs    = ParallelDescriptor::NProcs();
    const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
    const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
    const int   mySet     = myProc/nOutFiles;
    std::string cFileName(cdir + '/' + cdir + "_D_");
    std::string FullName  = BoxLib::Concatenate(cFileName, myProc % nOutFiles, 4);

    for(int iSet = 0; iSet < nSets; ++iSet) {
      if(mySet == iSet) {
        {  // scope
          std::ofstream csFile;

          if(iSet == 0) {   // First set.
            csFile.open(FullName.c_str(),
                        std::ios::out|std::ios::trunc|std::ios::binary);
          } else {
            csFile.open(FullName.c_str(),
                        std::ios::out|std::ios::app|std::ios::binary);
            csFile.seekp(0, std::ios::end);   // set to eof
          }
          if( ! csFile.good()) {
            BoxLib::FileOpenFailed(FullName);
          }

          // ----------------------------- write to file here
          WriteStats(csFile);
          // ----------------------------- end write to file here

          csFile.flush();
          csFile.close();
        }  // end scope

        int iBuff     = 0;
        int wakeUpPID = (myProc + nOutFiles);
        int tag       = (myProc % nOutFiles);
        if(wakeUpPID < nProcs) {
          ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
        }
      }
      if(mySet == (iSet + 1)) {   // Next set waits.
        int iBuff;
        int waitForPID = (myProc - nOutFiles);
        int tag        = (myProc % nOutFiles);
        ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
      }
    }
    // --------------------- end nfiles block


    ParallelDescriptor::Barrier("Profiler::Finalize");
  }
  if(ParallelDescriptor::IOProcessor()) {
    std::cout << "Profiler::Finalize():  time:  "   // time the timer
              << ParallelDescriptor::second() - finalizeStart << std::endl;
  }


  // --------------------------------------- print communication stats
#ifdef BL_COMM_PROFILING
  WriteCommStats();
#endif

  bInitialized = false;
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

  ios << '\n' << '\n';
  if( ! bwriteavg) {
    ios << std::setfill('*')
        << std::setw(maxlen + 2 + 3 * (colWidth + 2) - (colWidth+12)) << "";
    ios << std::setfill(' ');
    ios << "  Processor:  " << std::setw(colWidth) << myProc << '\n';
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
  ios << '\n';
  ios << "Total Timers     = " << std::setw(colWidth) << totalTimers
      << " seconds." << '\n';
  if(calcRunTime > 0.0) {
    percent = 100.0 * totalTimers / calcRunTime;
    ios << "Calc Run Time    = " << std::setw(colWidth) << calcRunTime
        << " seconds." << '\n';
    ios << "Percent Coverage = " << std::setw(colWidth) << percent << " %" << '\n';
  }

  // -------- write timers sorted by percent
  ios << '\n' << '\n';
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
        << '\n';
  } else {
    ios << std::setfill('=') << std::setw(maxlen+4 + 3 * (colWidth+2)) << ""
        << '\n';
  }
  ios << std::setfill(' ');
  ios << std::endl;
}


void Profiler::WriteCommStats(const bool bFlushing) {

  bool bAllCFTypesExcluded(OnExcludeList(AllCFTypes));
  if( ! bAllCFTypesExcluded) {
    CommStats::cftExclude.insert(AllCFTypes);  // temporarily
  }

  if(bFlushing) {
    int nCS(vCommStats.size());
    ParallelDescriptor::ReduceIntMax(nCS);
    if(nCS < csFlushSize) {
      if( ! bAllCFTypesExcluded) {
        CommStats::cftExclude.erase(AllCFTypes);
      }
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Bypassing flush, nCS < csFlushSize:  " << nCS
	          << "  " << csFlushSize << std::endl;
      }
      return;
    } else {
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "Flushing commstats:  nCSmax = " << nCS << std::endl;
      }
    }
  }

  std::string cdir("bl_comm_prof");
  if(ParallelDescriptor::IOProcessor()) {
    if( ! BoxLib::UtilCreateDirectory(cdir, 0755)) {
      BoxLib::CreateDirectoryFailed(cdir);
    }
  }
  // Force other processors to wait until directory is built.
  ParallelDescriptor::Barrier("Profiler::WriteCommStats::waitfordir");

  bool bUseRelativeTimeStamp(true);
  if(bUseRelativeTimeStamp) {
    for(int ics(0); ics < Profiler::vCommStats.size(); ++ics) {
      CommStats &cs = Profiler::vCommStats[ics];
      cs.timeStamp -= startTime;
    }
  }

  //std::ostringstream csHeader;

  // --------------------- start nfiles block
  const int   myProc    = ParallelDescriptor::MyProc();
  const int   nProcs    = ParallelDescriptor::NProcs();
  const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
  const int   nSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
  const int   mySet     = myProc/nOutFiles;

  std::string shortDFileName(cdir + "_D_");
  shortDFileName = BoxLib::Concatenate(shortDFileName, myProc % nOutFiles, 4);
  std::string longDFileName  = cdir + '/' + shortDFileName;

  std::string shortHeaderFileName(cdir + "_H_");
  shortHeaderFileName = BoxLib::Concatenate(shortHeaderFileName, myProc % nOutFiles, 4);
  std::string longHeaderFileName  = cdir + '/' + shortHeaderFileName;

  for(int iSet = 0; iSet < nSets; ++iSet) {
    if(mySet == iSet) {
      {  // scope
        std::ofstream csDFile, csHeaderFile;

        if(iSet == 0 && bFirstCommWriteH) {   // First set.
	  bFirstCommWriteH = false;
          csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::trunc);
        } else {
          csHeaderFile.open(longHeaderFileName.c_str(), std::ios::out | std::ios::app);
          csHeaderFile.seekp(0, std::ios::end);   // set to eof
        }
        if( ! csHeaderFile.good()) {
          BoxLib::FileOpenFailed(longHeaderFileName);
        }

        if(iSet == 0 && bFirstCommWriteD) {   // First set.
	  bFirstCommWriteD = false;
          csDFile.open(longDFileName.c_str(),
                      std::ios::out | std::ios::trunc | std::ios::binary);
        } else {
          csDFile.open(longDFileName.c_str(),
                      std::ios::out | std::ios::app | std::ios::binary);
          csDFile.seekp(0, std::ios::end);   // set to eof
        }
        if( ! csDFile.good()) {
          BoxLib::FileOpenFailed(longDFileName);
        }

        // ----------------------------- write to file here
        csHeaderFile << "NProcs  " << nProcs << '\n';
        csHeaderFile << "CommStatsSize  " << sizeof(CommStats) << '\n';
        csHeaderFile << "CommProfProc  " << myProc
                     << "  nCommStats  " << vCommStats.size()
                     << "  datafile  " << shortDFileName
	             << "  seekpos  " << csDFile.tellp() << '\n';
        for(int ib(0); ib < CommStats::barrierNames.size(); ++ib) {
          int seekindex(CommStats::barrierNames[ib].second);
          CommStats &cs = vCommStats[seekindex];
          csHeaderFile << "bNum  " << cs.tag  // tag is used for barrier number
                       << "  " << '"' << CommStats::barrierNames[ib].first << '"'
                       << "  " << seekindex << '\n';
        }
        for(int ib(0); ib < CommStats::nameTags.size(); ++ib) {
          int seekindex(CommStats::nameTags[ib].second);
          csHeaderFile << "nTag  " << CommStats::nameTags[ib].first
                       << seekindex << '\n';
        }
        for(int ib(0); ib < CommStats::reductions.size(); ++ib) {
          int seekindex(CommStats::reductions[ib]);
          CommStats &cs = vCommStats[seekindex];
          csHeaderFile << "red  " << cs.tag  // tag is used for reduction number
	               << "  " << seekindex << '\n';
        }
	if(vCommStats.size() > 0) {
	  csHeaderFile << std::setprecision(16)
	               << "timeMinMax  " << vCommStats[0].timeStamp << "  "
	               << vCommStats[vCommStats.size()-1].timeStamp << '\n';
	  csHeaderFile << std::setprecision(16)
	               << "timerTime  " << timerTime << '\n';
	} else {
	  csHeaderFile << "timeMinMax  0.0  0.0" << '\n';
	}
        for(int i(0); i < CommStats::nameTagNames.size(); ++i) {
          csHeaderFile << "nameTagNames  " << '"' << CommStats::nameTagNames[i]
                       << '"' << '\n';
        }
	csHeaderFile.flush();
        csHeaderFile.close();

	csDFile.write((char *) &vCommStats[0], vCommStats.size() * sizeof(CommStats));

        csDFile.flush();
        csDFile.close();
        // ----------------------------- end write to file here
      }  // end scope

      int iBuff     = 0;
      int wakeUpPID = (myProc + nOutFiles);
      int tag       = (myProc % nOutFiles);
      if(wakeUpPID < nProcs) {
        ParallelDescriptor::Send(&iBuff, 1, wakeUpPID, tag);
      }
    }
    if(mySet == (iSet + 1)) {   // Next set waits.
      int iBuff;
      int waitForPID = (myProc - nOutFiles);
      int tag        = (myProc % nOutFiles);
      ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
    }
  }
  // --------------------- end nfiles block


  // --------------------- flush the data
  vCommStats.clear();
  CommStats::barrierNames.clear();
  CommStats::nameTags.clear();
  CommStats::reductions.clear();
  if( ! bAllCFTypesExcluded) {
    CommStats::cftExclude.erase(AllCFTypes);
  }
}



void Profiler::WriteHeader(std::ostream &ios, const int colWidth,
                           const Real maxlen, const bool bwriteavg)
{
  if(bwriteavg) {
    ios << std::setfill('-') << std::setw(maxlen+4 + 5 * (colWidth+2))
        << std::left << "Total times " << '\n';
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Min"
        << std::setw(colWidth + 2) << "Avg"
        << std::setw(colWidth + 2) << "Max"
        << std::setw(colWidth + 4) << "Percent %"
        << '\n';
  } else {
    ios << std::setfill('-') << std::setw(maxlen+4 + 3 * (colWidth+2))
        << std::left << "Total times " << '\n';
    ios << std::right << std::setfill(' ');
    ios << std::setw(maxlen + 2) << "Function Name"
        << std::setw(colWidth + 2) << "NCalls"
        << std::setw(colWidth + 2) << "Time"
        << std::setw(colWidth + 4) << "Percent %"
        << '\n';
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
	  << percent << " %" << '\n';
    } else {
      ios << std::setw(maxlen + 2) << fname << "  "
          << std::setw(colWidth) << pstats.nCalls << "  "
          << std::setprecision(numPrec) << std::fixed << std::setw(colWidth)
	  << pstats.totalTime << "  "
          << std::setprecision(pctPrec) << std::fixed << std::setw(colWidth)
	  << percent << " %" << '\n';
    }
}


bool Profiler::OnExcludeList(CommFuncType cft) {
  // 
  // the idea for NoCFTypes is to allow local filtering/unfiltering
  // while preserving the users exclude list
  // possibly use a filter stack instead
  // might need caching if performance is a problem
  // what to do if both all and none are on the list?
  // 
  if(CommStats::cftExclude.empty()) {  // nothing on the exclude list
    return false;
  }
  std::set<CommFuncType>::iterator cfti;
  cfti = CommStats::cftExclude.find(NoCFTypes);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude nothing
    return false;
  }
  cfti = CommStats::cftExclude.find(cft);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude this cft
    return true;
  }
  cfti = CommStats::cftExclude.find(AllCFTypes);
  if(cfti != CommStats::cftExclude.end()) {  // found, exclude all types
    return true;
  }
  return false;
}


void Profiler::AddCommStat(const CommFuncType cft, const int size,
                           const int pid, const int tag)
{
  if(OnExcludeList(cft)) {
    return;
  }
  vCommStats.push_back(CommStats(cft, size, pid, tag, ParallelDescriptor::second()));
}


void Profiler::AddBarrier(const std::string &message, const bool beforecall) {
  const CommFuncType cft(Profiler::Barrier);
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    int tag(CommStats::barrierNumber);
    vCommStats.push_back(CommStats(cft, 0, BeforeCall(), tag,
                                   ParallelDescriptor::second()));
    CommStats::barrierNames.push_back(std::make_pair(message, vCommStats.size() - 1));
    ++CommStats::barrierNumber;
  } else {
    int tag(CommStats::barrierNumber - 1);  // it was incremented before the call
    vCommStats.push_back(CommStats(cft, AfterCall(), AfterCall(), tag,
                                   ParallelDescriptor::second()));
  }
}


void Profiler::AddAllReduce(const CommFuncType cft, const int size,
                            const bool beforecall)
{
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    int tag(CommStats::reductionNumber);
    vCommStats.push_back(CommStats(cft, size, BeforeCall(), tag,
                                   ParallelDescriptor::second()));
    CommStats::reductions.push_back(vCommStats.size() - 1);
    ++CommStats::reductionNumber;
  } else {
    int tag(CommStats::reductionNumber - 1);
    vCommStats.push_back(CommStats(cft, size, AfterCall(), tag,
                                   ParallelDescriptor::second()));
  }
}


void Profiler::AddWaitsome(const CommFuncType cft, const Array<MPI_Request> &reqs,
                           const int completed, const Array<int> &indx,
			   const Array<MPI_Status> &status, const bool beforecall)
{
#ifdef BL_USE_MPI
  if(OnExcludeList(cft)) {
    return;
  }
  if(beforecall) {
    vCommStats.push_back(CommStats(cft, BeforeCall(), BeforeCall(), NoTag(),
                         ParallelDescriptor::second()));
  } else {
    for(int i(0); i < completed; ++i) {
      MPI_Status stat(status[i]);
      int c;
      BL_MPI_REQUIRE( MPI_Get_count(&stat, MPI_UNSIGNED_CHAR, &c) );
      vCommStats.push_back(CommStats(cft, c, stat.MPI_SOURCE, stat.MPI_TAG,
                           ParallelDescriptor::second()));
    }
  }
#endif
}


int Profiler::NameTagNameIndex(const std::string &name) {
  for(int i(0); i < CommStats::nameTagNames.size(); ++i) {
    if(CommStats::nameTagNames[i] == name) {
      return i;
    }
  }
  CommStats::nameTagNames.push_back(name);
  return CommStats::nameTagNames.size() - 1;
}


void Profiler::AddNameTag(const std::string &name) {
  const CommFuncType cft(Profiler::NameTag);
  if(OnExcludeList(cft)) {
    return;
  }
  int tag(NameTagNameIndex(name));
  int index(CommStats::nameTags.size());
  vCommStats.push_back(CommStats(cft, index,  vCommStats.size(), tag,
                       ParallelDescriptor::second()));
  CommStats::nameTags.push_back(std::make_pair(tag, vCommStats.size() - 1));
}


Profiler::CommFuncType Profiler::CommStats::StringToCFT(const std::string &s) {
  return CommStats::cftNames[s];
}


void Profiler::CommStats::Filter(CommFuncType cft) {
  if( ! OnExcludeList(cft)) {
    CommStats::cftExclude.insert(cft);
  }
}


void Profiler::CommStats::UnFilter(CommFuncType cft) {
  if(OnExcludeList(cft)) {
    CommStats::cftExclude.erase(cft);
  }
}


#else

#endif
