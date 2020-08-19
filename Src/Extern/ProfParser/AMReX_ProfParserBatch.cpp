// ----------------------------------------------------------------------
//  AMReX_ProfParserBatch.cpp
// ----------------------------------------------------------------------
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <list>
#include <typeinfo>
#include <limits>
#include <algorithm>
#include <iomanip>

using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::string;

#include <AMReX.H>
# include <AMReX_DataServices.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

extern void PrintTimeRangeList(const std::list<RegionsProfStats::TimeRange> &trList);
extern void amrex::RedistFiles();

namespace {
#define SHOWVAL(val) { cout << #val << " = " << val << endl; }
  const int NTIMESLOTS(25600);
}



// ----------------------------------------------------------------------
void PrintProfParserBatchUsage(std::ostream &os) {
      os << "   [-actpf f] output a plotfile for all call times for func f.\n";
      os << "                f is a quoted string.\n";
      os << "   [-check]   data integrity check.\n";
      os << "   [-dispatch] use the dispatch interface.\n";
      os << "   [-gl]      process only grdlog.\n";
      os << "   [-gpct]    set percent threshold for xgraphs.  range [0, 100]" << '\n';
      os << "   [-html]    write html." << '\n';
      os << "   [-htmlnc]  write html showing ncalls." << '\n';
      os << "   [-msil n]  sets maxSmallImageLength." << '\n';
      os << "   [-mff]     make filter file." << '\n';
      os << "   [-nocomb]  do not combine adjacent call traces." << '\n';
      os << "   [-nts  n]  sets number of time slots (default:  " << NTIMESLOTS << ")." << '\n';
      os << "   [-of  fn]  sets output file name." << '\n';
      os << "   [-proc  n] sets processor number for single processor queries (default:  0).\n";
      os << "   [-prof]    profile the parser." << '\n';
      os << "   [-proxmap] remap ranks to proximity ranks." << '\n';
      os << "   [-pff]     parse filter file." << '\n';
      os << "   [-redist]  redistribute files." << '\n';
      os << "   [-rplt]    make region plot file." << '\n';
      os << "   [-rra  n]  sets refRatioAll." << '\n';
      os << "   [-sendspf] output a sends plotfile." << '\n';
      os << "   [-spd]     process sync point data." << '\n';
      os << "   [-sr]      process sends and receives." << '\n';
      os << "   [-srlist]  list sends and receives." << '\n';
      os << "   [-stats]   print database statistics." << '\n';
      os << "   [-tce]     process only topolcoords for edison." << '\n';
      os << "   [-timelinepf] output a timeline plotfile." << '\n';
      os << "   [-ttrace]  write text call trace." << '\n';
      os << "   [-v n]     verbose:  n can be 0, 1, or 2." << '\n';
      os << "   [-ws]      write summary." << '\n';
      os << "   [-wts]     write trace summary." << '\n';
}



// ----------------------------------------------------------------------
bool ProfParserBatchFunctions(int argc, char *argv[], bool runDefault,
                              bool &bParserProf)
{

  BL_PROFILE_VAR("ProfParserBatchFunctions()", ppbf);

  bool bIOP(ParallelDescriptor::IOProcessor());

  int maxSmallImageLength(800), refRatioAll(4), nTimeSlots(NTIMESLOTS);
  std::map<int, string> mpiFuncNames;

  amrex::BLProfiler::SetBlProfDirName("bl_profprof");

  int verbose(-1), whichProc(0);
  bool runCheck(false), runSendRecv(false), runSendRecvList(false), runSyncPointData(false);
  bool runSendsPF(false), runTimelinePF(false), glOnly(false);
  bool tcEdisonOnly(false), runStats(false), runRedist(false);
  bool statsCollected(false), filenameSet(false), proxMap(false);
  bool bMakeFilterFile(false), bParseFilterFile(true);
  bool bWriteSummary(false), bWriteTraceSummary(false);
  bool bMakeRegionPlt(false), simpleCombine(true);
  bool bWriteHTML(false), bWriteHTMLNC(false), bWriteTextTrace(false);
  bool bRunACTPF(false), bUseDispatch(false);
  string outfileName, delimString("\t");
  Vector<string> actFNames;

  amrex::ignore_unused(bParseFilterFile);

  bParserProf = false;

  if(argc > 2) {  // parse the command line
    int ia(1);
    while(ia < argc-1) {
      if(bIOP) cout << "argv[" << ia << "] = " << argv[ia] << endl;
      if(strcmp(argv[ia], "-v") == 0) {
	if(ia < argc-2) {
	  verbose = atoi(argv[ia+1]);
	}
        if(bIOP) cout << "*** verbose = " << verbose << endl;
	++ia;
      } else if(strcmp(argv[ia], "-ws") == 0) {
        bWriteSummary = true;
      } else if(strcmp(argv[ia], "-check") == 0) {          // ---- commprof options
        if(bIOP) cout << "*** data integrity check." << endl;
        runCheck = true;
      } else if(strcmp(argv[ia], "-stats") == 0) {
        if(bIOP) cout << "*** print database statistics." << endl;
        runStats = true;
      } else if(strcmp(argv[ia], "-timelinepf") == 0) {
        if(bIOP) cout << "*** output a timeline plotfile." << endl;
        runTimelinePF = true;
        runStats = true;
      } else if(strcmp(argv[ia], "-actpf") == 0) {
        bRunACTPF = true;
	if(ia < argc-2) {
          actFNames.push_back(argv[ia+1]);
	}
        if(bIOP) cout << "*** output a plotfile for all call times for function"
	              << actFNames[actFNames.size() - 1] << endl;
	++ia;
      } else if(strcmp(argv[ia], "-sr") == 0) {
        if(bIOP) cout << "*** send receive pairing." << endl;
        runSendRecv = true;
      } else if(strcmp(argv[ia], "-srlist") == 0) {
        if(bIOP) cout << "*** send receive pairing list." << endl;
        runSendRecvList = true;
      } else if(strcmp(argv[ia], "-sendspf") == 0) {
        if(bIOP) cout << "*** sendspf." << endl;
        runSendsPF = true;
      } else if(strcmp(argv[ia], "-gl") == 0) {
        if(bIOP) cout << "*** grdlog." << endl;
        glOnly = true;
      } else if(strcmp(argv[ia], "-tce") == 0) {
        if(bIOP) cout << "*** topolcoords for edison." << endl;
        tcEdisonOnly = true;
      } else if(strcmp(argv[ia], "-spd") == 0) {
        if(bIOP) cout << "*** sync point data." << endl;
        runSyncPointData = true;
      } else if(strcmp(argv[ia], "-redist") == 0) {
        if(bIOP) cout << "*** redist." << endl;
        runRedist = true;
      } else if(strcmp(argv[ia], "-msil") == 0) {
	if(ia < argc-2) {
          maxSmallImageLength = atoi(argv[ia+1]);
	}
        if(bIOP) cout << "*** msil = " << maxSmallImageLength << endl;
	++ia;
      } else if(strcmp(argv[ia], "-rra") == 0) {
	if(ia < argc-2) {
          refRatioAll = atoi(argv[ia+1]);
	}
        if(bIOP) cout << "*** rra = " << refRatioAll << endl;
	++ia;
      } else if(strcmp(argv[ia], "-nts") == 0) {
	if(ia < argc-2) {
          nTimeSlots = atoi(argv[ia+1]);
	}
        if(bIOP) cout << "*** nts = " << nTimeSlots << endl;
	++ia;
      } else if(strcmp(argv[ia], "-proc") == 0) {
	if(ia < argc-2) {
          whichProc = atoi(argv[ia+1]);
	}
        if(bIOP) cout << "*** whichProc = " << whichProc << endl;
	++ia;
      } else if(strcmp(argv[ia], "-of") == 0) {
	if(ia < argc-2) {
          outfileName = argv[ia+1];
	  filenameSet = true;
	}
        if(bIOP) cout << "*** outfileName = " << outfileName << endl;
	++ia;
      } else if(strcmp(argv[ia], "-proxmap") == 0) {
        if(bIOP) cout << "*** proxmap." << endl;
        proxMap = true;

      } else if(strcmp(argv[ia], "-mff") == 0) {     // ---- region and trace options
        bMakeFilterFile = true;
      } else if(strcmp(argv[ia], "-pff") == 0) {
        bParseFilterFile = true;
      } else if(strcmp(argv[ia], "-npff") == 0) {
        bParseFilterFile = false;
      } else if(strcmp(argv[ia], "-rplt") == 0) {
        bMakeRegionPlt = true;
      } else if(strcmp(argv[ia], "-wts") == 0) {
        bWriteTraceSummary = true;
      } else if(strcmp(argv[ia], "-html") == 0) {
        bWriteHTML = true;
      } else if(strcmp(argv[ia], "-htmlnc") == 0) {
        bWriteHTMLNC = true;
      } else if(strcmp(argv[ia], "-ttrace") == 0) {
        bWriteTextTrace = true;
      } else if(strcmp(argv[ia], "-gpct") == 0) {
	if(ia < argc-2) {
          Real gpct(atof(argv[ia+1]));
	  if(gpct >= 0.0 && gpct <= 100.0) {
	    RegionsProfStats::SetGPercent(gpct);
            if(bIOP) cout << "*** gpct = " << gpct << endl;
	  } else {
            if(bIOP) cout << "*** gpct must be in range [0.0, 100.0]" << endl;
	  }
	}
	++ia;
      } else if(strcmp(argv[ia], "-nocomb") == 0) {
        simpleCombine = false;
      } else if(strcmp(argv[ia], "-prof") == 0) {
        bParserProf = true;
      } else if(strcmp(argv[ia], "-dispatch") == 0) {
        if(bIOP) cout << "*** using dispatch interface." << endl;
        bUseDispatch = true;
      } else {
        if(bIOP) cerr << "*** Error:  bad command line arg:  " << argv[ia] << endl;
      }
      ++ia;
    }
  }
  if(argc == 2 && runDefault) {  // ---- set default to write the summary
    bWriteSummary = true;
  }


  // ---------------------------------------
  BLProfStats::SetVerbose(verbose);
  std::string dirName(argv[argc - 1]);

  Amrvis::FileType fileType(Amrvis::PROFDATA);
  DataServices pdServices(dirName, fileType);

  int dataNProcs(BLProfStats::GetNProcs());
  if(whichProc < 0 || whichProc > dataNProcs - 1) {
    if(bIOP) {
      cout << "**** Error:  whichProc out of range:  "
           << whichProc << " [0," << dataNProcs - 1 << "]." << endl;
    }
    whichProc = 0;
  }

  if(bWriteSummary) {
    bool writeAverage(true), useTrace(false), graphTopPct(true);
    if(bIOP) { cout << "Writing summary." << endl; }
    if(bUseDispatch) {
      if(bIOP) {
        DataServices::Dispatch(DataServices::WriteSummaryRequest, &pdServices,
                                   (void *) &(cout),
				   &writeAverage, whichProc, &useTrace,
				   &graphTopPct);
      }
    } else {
      pdServices.WriteSummary(cout, writeAverage, whichProc, useTrace,
                              graphTopPct);
    }
  }


  if(bWriteTraceSummary) {
    bool writeAverage(false), useTrace(true), graphTopPct(true);
    if(bIOP) { cout << "Writing trace summary." << endl; }
    if(bUseDispatch) {
      if(bIOP) {
        DataServices::Dispatch(DataServices::InitTimeRanges, &pdServices);
        DataServices::Dispatch(DataServices::WriteSummaryRequest, &pdServices,
                                   (void *) &(cout),
				   &writeAverage, whichProc, &useTrace,
				   &graphTopPct);
      }
    } else {
      pdServices.InitRegionTimeRanges();
      pdServices.WriteSummary(cout, writeAverage, whichProc, useTrace,
                              graphTopPct);
    }
  }


  if(proxMap) {
    pdServices.InitProxMap();
  }


  if(runCheck) {
    if(bUseDispatch) {
      if(bIOP) {
        DataServices::Dispatch(DataServices::CheckProfDataRequest, &pdServices);
      }
    } else {
      pdServices.CheckProfData();
    }
  }


  if(runStats) {
    if(bIOP) { cout << endl << "---------------- runStats." << endl; }
    if(bUseDispatch) {
      if(bIOP) {
        DataServices::Dispatch(DataServices::RunStatsRequest, &pdServices,
                                   (void *) &(mpiFuncNames),
				   &statsCollected);
      }
    } else {
      pdServices.RunStats(mpiFuncNames, statsCollected);
    }
  }


  if(runSendsPF) {
    std::string plotfileName("pltTSP2P");
    if(filenameSet) {
      plotfileName = outfileName;
    }
    if(bUseDispatch) {
      if(bIOP) {
        DataServices::Dispatch(DataServices::RunSendsPFRequest, &pdServices,
                                   (void *) &(plotfileName),
				   maxSmallImageLength,
				   &proxMap,
				   refRatioAll);
      }
    } else {
      pdServices.RunSendsPF(plotfileName, maxSmallImageLength,
                            proxMap, refRatioAll);
    }
  }


  if(runTimelinePF) {
    std::string plotfileName("pltTimeline");
    if(filenameSet) {
      plotfileName = outfileName;
    }
    if(bUseDispatch) {
      if(bIOP) {
        cout << "Dispatched batch timelinepf currently unavailable." << endl;
/*        DataServices::Dispatch(DataServices::RunTimelinePFRequest, &pdServices,
                                   (void *) &(mpiFuncNames),
                                   (void *) &(plotfileName),
				   maxSmallImageLength,
				   refRatioAll,
				   nTimeSlots,
				   &statsCollected);*/
      }
    } else {
      BLProfStats::TimeRange subTimeRange = pdServices.FindCalcTimeRange();
      pdServices.RunTimelinePF(mpiFuncNames, plotfileName, subTimeRange, 
			       maxSmallImageLength, refRatioAll,
			       nTimeSlots, statsCollected);
    }
  }


  // ---- process the grid log
  if(glOnly) {
    std::string gridlogFileName("grdlog");
    if(filenameSet) {
      // ---- outfileName is used as an infile here
      // ---- dont use this function with other functions using outfileName
      gridlogFileName = outfileName;
    }
    pdServices.ProcessGridLog(gridlogFileName);
  }


  if(tcEdisonOnly) {
    pdServices.TCEdison();
  }


  if(runSyncPointData) {
    pdServices.RunSyncPointData();
  }


  if(runSendRecvList) {
    pdServices.InitRegionTimeRanges();
    pdServices.RunSendRecvList();
  }


  if(runSendRecv) {
    pdServices.RunSendRecv();
  }


  if(runRedist) {
    RedistFiles();
  }


  if(bMakeFilterFile) {
    std::string filterFileName("RegionFilters.txt");
    if(filenameSet) {
      filterFileName = outfileName;
    }
    if(bUseDispatch) {
      if(bIOP) {
        DataServices::Dispatch(DataServices::MakeFilterFileRequest, &pdServices,
                                   (void *) &(filterFileName));
      }
    } else {
      pdServices.MakeFilterFile(filterFileName);
    }
  }


  if(bMakeRegionPlt) {
    std::string plotfileName("pltRegions");
    if(filenameSet) {
      plotfileName = outfileName;
    }
    if(bUseDispatch) {
      if(bIOP) {
        DataServices::Dispatch(DataServices::MakeRegionPltRequest, &pdServices,
                                   (void *) &(plotfileName));
      }
    } else {
      pdServices.MakeRegionPlt(plotfileName);
    }
  }


  if(bWriteHTML) {
    std::string callTraceFileName("CallTrace.html");
    pdServices.WriteHTML(callTraceFileName, simpleCombine, whichProc);
  }

  if(bWriteHTMLNC) {
    std::string callTraceFileName("CallTraceNC.html");
    pdServices.WriteHTMLNC(callTraceFileName, whichProc);
  }

  if(bWriteTextTrace) {
    std::string callTraceFileName("CallTrace.txt");
    pdServices.WriteTextTrace(callTraceFileName, simpleCombine, whichProc, delimString);
  }

  if(bIOP) {
    //PrintTimeRangeList(pdServices.GetRegionsProfStats().GetFilterTimeRanges()[0]);
  }


  // ---- output an all call times plotfile
  if(bRunACTPF) {
    std::string plotfileName;
    if(filenameSet) {
      plotfileName = outfileName;
    }
    pdServices.RunACTPF(plotfileName, maxSmallImageLength, refRatioAll, actFNames);
  }



  bool anyFunctionsRun = runCheck       || runSendRecv     || runSendRecvList || runSyncPointData   ||
                         runSendsPF     || runTimelinePF   || tcEdisonOnly    || runStats           ||
			 runRedist      || bMakeFilterFile || bWriteSummary   || bWriteTraceSummary ||
                         bMakeRegionPlt || bWriteHTML      || bWriteHTMLNC    || bWriteTextTrace    ||
                         glOnly         || bRunACTPF;

  BL_PROFILE_VAR_STOP(ppbf);

  return anyFunctionsRun;
}

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
