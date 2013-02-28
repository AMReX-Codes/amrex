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
bool Profiler::bFirstCommWriteH = true;  // header
bool Profiler::bFirstCommWriteD = true;  // data
bool Profiler::bInitialized = false;

int Profiler::currentStep = 0;
int Profiler::csFlushSize = 1000000;
int Profiler::nProfFiles  = 3;

Real Profiler::pctTimeLimit = 5.0;
Real Profiler::calcRunTime  = 0.0;
Real Profiler::startTime    = 0.0;

std::stack<Real> Profiler::nestedTimeStack;
std::map<int, Real> Profiler::mStepMap;
std::map<std::string, Profiler::ProfStats> Profiler::mProfStats;
std::map<Real, std::string, std::greater<Real> > Profiler::mTimersTotalsSorted;
std::vector<Profiler::CommStats> Profiler::vCommStats;
std::map<std::string, Profiler::CommFuncType> Profiler::CommStats::cftNames;
std::set<Profiler::CommFuncType> Profiler::CommStats::cftExclude;
int Profiler::CommStats::barrierNumber = 0;
std::vector<std::pair<std::string,int> > Profiler::CommStats::barrierNames;
std::vector<std::pair<std::string,int> > Profiler::CommStats::nameTags;



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

  startTime = ParallelDescriptor::second();

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
  CommStats::cftNames["ScatterTsT1si"]  = ScatterTsT1si; 
  CommStats::cftNames["Barrier"]        = Barrier;
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
	  std::cout << myProc << ":  #### ProfName not found, inserting:  "
	            << pfName << std::endl;
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
    if(itcopy == mProfStatsCopy.end()) {
      std::cout << myProc << ":  #### ProfName not on ioproc:  "
                << it->first << std::endl;
    }
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
    // Force other processors to wait till directory is built.
    ParallelDescriptor::Barrier();

    const int   myProc    = ParallelDescriptor::MyProc();
    const int   nProcs    = ParallelDescriptor::NProcs();
    const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
    const int   NSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
    const int   MySet     = myProc/nOutFiles;
    std::string cFileName(cdir + '/' + cdir + "_D_");
    std::string FullName  = BoxLib::Concatenate(cFileName, myProc % nOutFiles, 4);

    for(int iSet = 0; iSet < NSets; ++iSet) {
      if(MySet == iSet) {
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
      if(MySet == (iSet + 1)) {   // Next set waits.
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
  WriteCommStats();

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

}


void Profiler::WriteCommStats() {

  bool bBarrierExcluded(OnExcludeList(Barrier));
  if( ! bBarrierExcluded) {
    CommStats::cftExclude.insert(Barrier);  // temporarily
  }

  std::string cdir("bl_comm_prof");
  if(ParallelDescriptor::IOProcessor()) {
    if( ! BoxLib::UtilCreateDirectory(cdir, 0755)) {
      BoxLib::CreateDirectoryFailed(cdir);
    }
  }
  // Force other processors to wait till directory is built.
  ParallelDescriptor::Barrier();

  std::ostringstream csHeader;

  bool bUseRelativeTimeStamp(true);
  if(bUseRelativeTimeStamp) {
    for(int ics(0); ics < Profiler::vCommStats.size(); ++ics) {
      CommStats &cs = Profiler::vCommStats[ics];
      cs.timeStamp -= startTime;
    }
  }

  // --------------------- start nfiles block
  const int   myProc    = ParallelDescriptor::MyProc();
  const int   nProcs    = ParallelDescriptor::NProcs();
  const int   nOutFiles = std::max(1, std::min(nProcs, nProfFiles));
  const int   NSets     = (nProcs + (nOutFiles - 1)) / nOutFiles;
  const int   MySet     = myProc/nOutFiles;
  std::string cShortFileName(cdir + "_D_");
  cShortFileName = BoxLib::Concatenate(cShortFileName, myProc % nOutFiles, 4);
  std::string FullName  = cdir + '/' + cShortFileName;

  for(int iSet = 0; iSet < NSets; ++iSet) {
    if(MySet == iSet) {
      {  // scope
        std::ofstream csFile;

        if(iSet == 0 && bFirstCommWriteD) {   // First set.
	  bFirstCommWriteD = false;
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
        csHeader << "CommProfProc  " << myProc
                 << "  nCommStats  " << vCommStats.size()
                 << "  datafile  " << cShortFileName
	         << "  seekpos  " << csFile.tellp() << std::endl;
        for(int ib(0); ib < CommStats::barrierNames.size(); ++ib) {
          int index(CommStats::barrierNames[ib].second);
          CommStats &cs = vCommStats[index];
          csHeader << "barrierNumber  " << cs.size  // size is used for barrier number
                   << "  name  " << '"' << CommStats::barrierNames[ib].first << '"'
                   << "  index  " << index
                   << std::endl;
        }
        for(int ib(0); ib < CommStats::nameTags.size(); ++ib) {
          int index(CommStats::nameTags[ib].second);
          CommStats &cs = vCommStats[index];
          csHeader << "nameTag  "
                   << "  name  " << '"' << CommStats::nameTags[ib].first << '"'
                   << "  index  " << index
                   << std::endl;
        }

	csFile.write((char *) &vCommStats[0], vCommStats.size() * sizeof(CommStats));
        // ----------------------------- end write to file here
	//==================
	/*
	std::cout << "====================" << std::endl;
	std::cout << myProc << ":  FullName = " << FullName << std::endl;
	int nouts(std::min(16, (int) vCommStats.size()));
	for(int i(0); i < nouts; ++i) {
          CommStats &cs = vCommStats[i];
	  std::cout << Profiler::CommStats::CFTToString(cs.cfType) << "  "
	            << cs.commpid << "  " << cs.size << "  " << cs.tag << "  "
		    << cs.timeStamp << std::endl;
	}
	std::cout << "====================" << std::endl;
	*/
	//==================

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
    if(MySet == (iSet + 1)) {   // Next set waits.
      int iBuff;
      int waitForPID = (myProc - nOutFiles);
      int tag        = (myProc % nOutFiles);
      ParallelDescriptor::Recv(&iBuff, 1, waitForPID, tag);
    }
  }
  // --------------------- end nfiles block

  // --------------------- gather all header data to the ioproc
  int ioProcNum(ParallelDescriptor::IOProcessorNumber());
  int headerSize(csHeader.str().size());

  Array<int> vHeaderSizes(nProcs, 0);
  ParallelDescriptor::Gather(&headerSize, 1, vHeaderSizes.dataPtr(), 1, ioProcNum);

  Array<int> offset(nProcs, 0);
  int totalHeaderSize(0);   // pad(0);

  if(ParallelDescriptor::IOProcessor()) {
    for(int i(1); i < offset.size(); ++i) {
      offset[i] = offset[i-1] + vHeaderSizes[i-1];
    }
    for(int i(0); i < offset.size(); ++i) {
      totalHeaderSize += vHeaderSizes[i];
    }
    //pad = 8 - (totalHeaderSize % 8);
    //pad = 0;
    //totalHeaderSize += pad;
  }

  Array<char> recvdata(totalHeaderSize + 1, '\0');
  if(recvdata.empty()) {
    recvdata.resize(1);
  }

  Array<char> senddata(csHeader.str().size() + 1, '\0');
  if(senddata.empty()) {
    senddata.resize(1);
  }

  BL_ASSERT((senddata.size() - 1) == csHeader.str().size());
  std::strcpy(senddata.dataPtr(), csHeader.str().c_str());

  BL_MPI_REQUIRE( MPI_Gatherv(senddata.dataPtr(),
			      csHeader.str().size(),
                              ParallelDescriptor::Mpi_typemap<char>::type(),
			      recvdata.dataPtr(),
                              vHeaderSizes.dataPtr(),
                              offset.dataPtr(),
                              ParallelDescriptor::Mpi_typemap<char>::type(),
                              ioProcNum,
                              ParallelDescriptor::Communicator()) );


  // --------------------- ioproc writes the header
  if(ParallelDescriptor::IOProcessor()) {
    std::ofstream csHeaderFile;
    std::string cHeaderName(cdir + '/' + cdir + "_H");
    if(bFirstCommWriteH) {
      csHeaderFile.open(cHeaderName.c_str(), std::ios::out | std::ios::trunc);
      bFirstCommWriteH = false;
    } else {
      csHeaderFile.open(cHeaderName.c_str(), std::ios::out | std::ios::app);
    }

    csHeaderFile << "NProcs  " << nProcs << std::endl;
    csHeaderFile << "CommStatsSize  " << sizeof(CommStats) << std::endl;
    csHeaderFile << recvdata.dataPtr();
    csHeaderFile.flush();
    csHeaderFile.close();
  }

  // --------------------- flush the data
  vCommStats.clear();
  CommStats::barrierNames.clear();
  CommStats::nameTags.clear();
  if( ! bBarrierExcluded) {
    CommStats::cftExclude.erase(Barrier);
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


void Profiler::AddCommStat(const CommFuncType cft, const int pid, const int size) {
  if(OnExcludeList(cft)) {
    return;
  }
  CommStats cs(cft, ParallelDescriptor::MyProc(), pid, size,
               ParallelDescriptor::second());
  vCommStats.push_back(cs);
}


void Profiler::AddCommStat(const CommFuncType cft, const int pid,
                           const int size, const int tag)
{
  if(OnExcludeList(cft)) {
    return;
  }
  CommStats cs(cft, pid, size, tag, ParallelDescriptor::second());
  vCommStats.push_back(cs);
}


void Profiler::AddBarrier(const std::string &message, const bool beforecall) {
  const CommFuncType cft(Profiler::Barrier);
  if(OnExcludeList(cft)) {
    return;
  }
  CommStats cs(cft, 0, Profiler::CommStats::barrierNumber, ParallelDescriptor::second());
  vCommStats.push_back(cs);
  CommStats::barrierNames.push_back(std::make_pair(message, vCommStats.size() - 1));
  ++CommStats::barrierNumber;
std::cout << ParallelDescriptor::MyProc() << "::::: vCommStats.size() = " << vCommStats.size() << std::endl;
  if(vCommStats.size() > csFlushSize) {
std::cout << ParallelDescriptor::MyProc() << "*********** writing commstats." << std::endl;
    WriteCommStats();
  }
}


void Profiler::AddNameTag(const std::string &name) {
  const CommFuncType cft(Profiler::NameTag);
  if(OnExcludeList(cft)) {
    return;
  }
  CommStats cs(cft, 0, 0, ParallelDescriptor::second());
  vCommStats.push_back(cs);
  CommStats::nameTags.push_back(std::make_pair(name, vCommStats.size() - 1));
}


void Profiler::AddAllReduce(CommFuncType cft, int size) {
  if(OnExcludeList(cft)) {
    return;
  }
  CommStats cs(cft, 0, 0, size, ParallelDescriptor::second());
  vCommStats.push_back(cs);
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
