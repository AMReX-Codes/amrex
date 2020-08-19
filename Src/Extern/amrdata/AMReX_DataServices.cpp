// ---------------------------------------------------------------
// DataServices.cpp
// ---------------------------------------------------------------
#include <AMReX_AmrvisConstants.H>
#include <AMReX_DataServices.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_XYPlotDataList.H>

#ifdef BL_USE_PROFPARSER
#include <AMReX_BLWritePlotFile.H>
//#include <AMReX_BLProfStats.H>
//#include <AMReX_CommProfStats.H>
//#include <AMReX_RegionsProfStats.H>
//#include <AMReX_BLProfUtilities.cpp>
#include <iomanip>
#include <cstdarg>
#endif

#include <iostream>
#include <fstream>
#include <cstdio>

using std::ios;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;

using namespace amrex;

#ifdef BL_USE_PROFPARSER
extern int yyparse(void *);
extern FILE *yyin;
#endif

namespace amrex {

Vector<DataServices *> DataServices::dsArray;
int DataServices::dsArrayIndexCounter = 0;
int DataServices::dsFabOutSize = 0;
bool DataServices::dsBatchMode = false;
bool DataServices::profiler = false;

#ifdef BL_USE_PROFPARSER
namespace {
  const int XDIR(0);
  const int YDIR(1);
  const int ZDIR(2);
  const int NTIMESLOTS(25600);
}


extern void PrintTimeRangeList(const std::list<RegionsProfStats::TimeRange> &trList);
extern void SimpleRemoveOverlap(BoxArray &ba);
extern void avgDown(MultiFab &S_crse, MultiFab &S_fine, int scomp, int dcomp,
             int ncomp, Vector<int> &ratio);
extern std::string SanitizeName(const std::string &s);
extern void WriteFab(const string &filenameprefix, const int xdim, const int ydim,
                     const double *data);

#define SHOWVAL(val) { cout << #val << " = " << val << endl; }
#endif


// ---------------------------------------------------------------
namespace ParallelDescriptor {
  template <> void Bcast (Box *b, size_t n, int root) {
    const int n3SDim(n * 3 * BL_SPACEDIM);

    Vector<int> tmp(n3SDim);

    int cnt(0);

    for(int j(0); j < n; ++j) {
      for(int i(0); i < BL_SPACEDIM; ++i) {
        tmp[cnt++] = b[j].smallEnd(i);
      }

      for(int i(0); i < BL_SPACEDIM; ++i) {
        tmp[cnt++] = b[j].bigEnd(i);
      }

      IntVect indx = b[j].type();

      for(int i(0); i < BL_SPACEDIM; ++i) {
        tmp[cnt++] = indx[i];
      }
    }

    BL_ASSERT(n3SDim == cnt);

    ParallelDescriptor::Bcast(&tmp[0], n3SDim, root);

    cnt = 0;

    for(int j(0); j < n; ++j) {
      IntVect sm(&tmp[cnt]);
      cnt += BL_SPACEDIM;
      IntVect bg(&tmp[cnt]);
      cnt += BL_SPACEDIM;
      IntVect id(&tmp[cnt]);
      cnt += BL_SPACEDIM;

      b[j] = Box(sm, bg, id);
    }
  }
}  // end namespace


// ---------------------------------------------------------------
DataServices::DataServices(const string &filename, const Amrvis::FileType &filetype)
             : fileName(filename), fileType(filetype), bAmrDataOk(false),
               iWriteToLevel(-1)
{
  numberOfUsers = 0;  // the user must do all incrementing and decrementing
  bAmrDataOk = amrData.ReadData(fileName, fileType);
  profiler = (fileType == Amrvis::PROFDATA);

  if((bAmrDataOk)||(profiler)) {
    dsArrayIndex = DataServices::dsArrayIndexCounter;
    ++DataServices::dsArrayIndexCounter;
    DataServices::dsArray.resize(DataServices::dsArrayIndexCounter);
    DataServices::dsArray[dsArrayIndex] = this;
  }

#ifdef BL_USE_PROFPARSER
  bProfDataAvailable = false;
  bRegionDataAvailable = false;
  bTraceDataAvailable = false;
  bCommDataAvailable = false;

  if(profiler) {
    Init(filename, filetype);
  }
#endif
}


// ---------------------------------------------------------------
DataServices::DataServices() {
  // must call init
  bAmrDataOk = false;
  iWriteToLevel = -1;
}


// ---------------------------------------------------------------
void DataServices::Init(const string &filename, const Amrvis::FileType &filetype) {

  BL_PROFILE("DataServices::Init");

  fileName = filename;
  fileType = filetype;
  bAmrDataOk = false;
  iWriteToLevel = -1;

  numberOfUsers = 0;  // the user must do all incrementing and decrementing
  profiler = (fileType == Amrvis::PROFDATA);

#ifdef BL_USE_PROFPARSER
  #if (BL_SPACEDIM == 2)
  if (profiler)
  {
    bool bIOP(ParallelDescriptor::IOProcessor());
    int  myProc(ParallelDescriptor::MyProc());
    int  nProcs(ParallelDescriptor::NProcs());

    BLProfStats::SetDirName(fileName);

    // -------- parse the main blprof header file.  everyone does this for now
    if(bIOP) { cout << "Parsing main blprof header file." << endl; }
    string blpFileName_H("bl_prof_H");
    string blpFullFileName_H(fileName + '/' + blpFileName_H);
    if( ! (yyin = fopen(blpFullFileName_H.c_str(), "r"))) {
      if(bIOP) {
        cerr << "DataServices::Init:  0:  Cannot open file:  " << blpFullFileName_H << endl;
      }
      bProfDataAvailable = false;
    } else {
      yyparse(&blProfStats_H);
      fclose(yyin);
      bProfDataAvailable = true;
    }

    // -------- parse the main call stats header file.  everyone does this for now
    if(bIOP) { cout << "Parsing main call stats header file." << endl; }
    string regPrefix_H("bl_call_stats_H");
    std::string regFileName_H(fileName + '/' + regPrefix_H);
    if( ! (yyin = fopen(regFileName_H.c_str(), "r"))) {
      if(bIOP) {
        cerr << "DataServices::Init:  1:  Cannot open file:  " << regFileName_H << endl;
      }
      bRegionDataAvailable = false;
      bTraceDataAvailable  = false;
    } else {
      bRegionDataAvailable = true;
      bTraceDataAvailable  = true;
    }
    if(bTraceDataAvailable) {
      yyparse(&regOutputStats_H);
      fclose(yyin);
    } else {
    }

    // ---- make a box for distributing work
    int dataNProcs(BLProfStats::GetNProcs());
    Box procBox(IntVect(0, 0), IntVect(0, dataNProcs - 1));
    IntVect procMaxGrid(1, (dataNProcs / nProcs) + ((dataNProcs % nProcs) > 0 ? 1 : 0));
    BoxArray procBoxArrayTemp(procBox);
    procBoxArrayTemp.maxSize(procMaxGrid);
    // ---- now ensure the boxarray is nprocs long
    Vector<Box> procBoxes;
    int needMoreBoxes(nProcs - procBoxArrayTemp.size());
    for(int ipb(0); ipb < procBoxArrayTemp.size(); ++ipb) {
      Box b(procBoxArrayTemp[ipb]);
      if(needMoreBoxes) {
        Box chopBox(b.chop(YDIR, (b.smallEnd(YDIR) + b.bigEnd(YDIR)) / 2));
        procBoxes.push_back(chopBox);
        --needMoreBoxes;
      }
      procBoxes.push_back(b);
    }
    procBoxArray.resize(procBoxes.size());
    for(int i(0); i < procBoxes.size(); ++i) {
      procBoxArray.set(i, procBoxes[i]);
    }

    if(procBoxArray.size() != nProcs) {
      SHOWVAL(nProcs);
      SHOWVAL(dataNProcs);
      SHOWVAL(procBoxArray.size());
      if(bIOP) cout << "---- procBoxArray = " << procBoxArray << endl;
      amrex::Abort("procBoxArray::Error 0");
    }

    if(bTraceDataAvailable) {
      // -------- parse the data headers.  everyone does this for now
      if(bIOP) { cout << "Parsing data headers." << endl; }
      const Vector<string> &regHeaderFileNames = regOutputStats_H.GetHeaderFileNames();
      for(int i(0); i < regHeaderFileNames.size(); ++i) {
        std::string regFileName_H_nnnn(fileName + '/' + regHeaderFileNames[i]);
        if( ! (yyin = fopen(regFileName_H_nnnn.c_str(), "r"))) {
          if(bIOP) {
            cerr << "DataServices::Init:  2:  Cannot open file:  " << regFileName_H_nnnn
                 << " ... continuing." << endl;
          }
          continue;
        }
        BL_PROFILE_VAR("DataServices::Init(), parsing data headers.", yydheaders);
        yyparse(&regOutputStats_H);
        BL_PROFILE_VAR_STOP(yydheaders);

        fclose(yyin);
      }

      if(regOutputStats_H.TraceDataValid()) {
//        if(bIOP) {
//	  cout << "Calling InitRegionTimeRanges." << endl;
//	}
        RegionsProfStats::OpenAllStreams(fileName);
        Box myBox(procBoxArray[myProc]);
//        bRegionDataAvailable = regOutputStats_H.InitRegionTimeRanges(myBox);
        regOutputStats_H.SyncFNamesAndNumbers();
        RegionsProfStats::CloseAllStreams();
        regOutputStats_H.SetFNames(blProfStats_H.BLPFNames());
//        if(bIOP) {
//	  cout << "Finished InitRegionTimeRanges." << endl;
//	}
      } else {
        bTraceDataAvailable = false;
      }

    }
/*
    bool bParseFilterFile(true);
    if(bParseFilterFile) {
      if(bIOP) cout << "Parsing filter file." << endl;
      ParseFilterFile();
    } 

    if( ! regOutputStats_H.TimeRangeInitialized()) {
      regOutputStats_H.InitFilterTimeRanges();
      if(bIOP) {
        cout << ">>>> timerangelist =" ;
        PrintTimeRangeList(regOutputStats_H.GetFilterTimeRanges()[0]);
      }

    }
*/
    // ----------------------------------------------- comm headers

    // -------- parse the main header file.  everyone does this for now
    if(bIOP) { cout << "Parsing main comm header file." << endl; }
    std::string commPrefix_H("bl_comm_prof_H");
    std::string commFileName_H(fileName + '/' + commPrefix_H);
    if( ! (yyin = fopen(commFileName_H.c_str(), "r"))) {
      if(bIOP) {
        cerr << "DataServices::Init:  3:  Cannot open file:  " << commFileName_H << endl;
      }
      bCommDataAvailable = false;
    } else {
      bCommDataAvailable = true;
    }

    if(bCommDataAvailable) {
      yyparse(&commOutputStats_H);
      fclose(yyin);
    } else {
    }
/*
    if(bRegionDataAvailable) {
      commOutputStats_H.SetRegionTimeRanges(regOutputStats_H.GetRegionTimeRanges());
      commOutputStats_H.SetFilterTimeRanges(regOutputStats_H.GetFilterTimeRanges());
    }
*/
    bAmrDataOk = true;
  } // if (profiler)
  #endif
#endif
}

// ---------------------------------------------------------------
void DataServices::InitRegionTimeRanges() {

  BL_PROFILE("DataServices::InitRegionTimeRanges");

#ifdef BL_USE_PROFPARSER
  #if (BL_SPACEDIM == 2)
  if (profiler)
  {
    bool bIOP(ParallelDescriptor::IOProcessor());
    int  myProc(ParallelDescriptor::MyProc());

    if(bTraceDataAvailable) {
      if(regOutputStats_H.TraceDataValid()) {
        if(bIOP) {
	  cout << "Calling InitRegionTimeRanges." << endl;
	}
        RegionsProfStats::OpenAllStreams(fileName);
        Box myBox(procBoxArray[myProc]);
        bRegionDataAvailable = regOutputStats_H.InitRegionTimeRanges(myBox);
        RegionsProfStats::CloseAllStreams();
        if(bIOP) {
	  cout << "Finished InitRegionTimeRanges." << endl;
	}
      } else {
        bTraceDataAvailable = false;
      }

    }

    if (dsBatchMode) {
      if(bIOP) { cout << "Parsing filter file." << endl; }
      ParseFilterFile();

      if( ! regOutputStats_H.TimeRangeInitialized()) {
        regOutputStats_H.InitFilterTimeRanges();
        if(bIOP) {
          cout << ">>>> timerangelist =" ;
          PrintTimeRangeList(regOutputStats_H.GetFilterTimeRanges()[0]);
        }

      }
    }
    // ----------------------------------------------- comm headers
    if(bRegionDataAvailable) {
      commOutputStats_H.SetRegionTimeRanges(regOutputStats_H.GetRegionTimeRanges());
      commOutputStats_H.SetFilterTimeRanges(regOutputStats_H.GetFilterTimeRanges());
    }

  } // if (profiler)
  #endif
#endif
}


// ---------------------------------------------------------------
DataServices::~DataServices() {
  BL_ASSERT(numberOfUsers == 0);
  if( ! profiler) {
    DataServices::dsArray[dsArrayIndex] = nullptr;
  }
}


// ---------------------------------------------------------------
void DataServices::SetBatchMode() {
    dsBatchMode = true;
}


// ---------------------------------------------------------------
void DataServices::SetFabOutSize(int iSize) {
  if (profiler) // Unused with profiling data
  { 
    if(iSize == 1 || iSize == 8 || iSize == 32) {
      dsFabOutSize = iSize;
    } else {
      cerr << "Warning:  DataServices::SetFabOutSize:  size must be 1, 8 or 32 only."
  	   << "  Defaulting to native." << endl;
      dsFabOutSize = 0;
    }
  }
}


// ---------------------------------------------------------------
void DataServices::Dispatch(DSRequestType requestType, DataServices *ds, ...) {
  bool bContinueLooping(true);
  va_list ap;
  int whichDSIndex;
  int ioProcNumber(ParallelDescriptor::IOProcessorNumber());

 while(bContinueLooping) {
  if(ParallelDescriptor::IOProcessor() || (dsBatchMode && !profiler)) {
    bContinueLooping = false;
  }

  ParallelDescriptor::Barrier();  // nonioprocessors wait here

  {
    int tmp = requestType;
    ParallelDescriptor::Bcast(&tmp, 1, 0);
    requestType = static_cast<DataServices::DSRequestType>(tmp);
  }

  // handle new request
  if(requestType == NewRequest) {

    // ProfDataServices currently skips over this request type
    // Returns to the top of the while loop
    if (profiler) {
      continue;
    }

    // broadcast the fileName and fileType to nonioprocessors
    char *fileNameCharPtr;
    int   fileNameLength(-1), fileNameLengthPadded(-1);
    Amrvis::FileType newFileType(Amrvis::INVALIDTYPE);

    if(ParallelDescriptor::IOProcessor()) {
      fileNameLength = ds->fileName.length();
      newFileType = ds->fileType;
    }

    ParallelDescriptor::Bcast(&fileNameLength, 1, 0);

    {
      int tmp = newFileType;
      ParallelDescriptor::Bcast(&tmp, 1, 0);
      newFileType = static_cast<Amrvis::FileType>(tmp);
    }

    fileNameLengthPadded = fileNameLength + 1;    // for the null
    fileNameLengthPadded += fileNameLengthPadded % 8;  // for alignment on the t3e
    fileNameCharPtr = new char[fileNameLengthPadded];
    if(ParallelDescriptor::IOProcessor()) {
      strcpy(fileNameCharPtr, ds->fileName.c_str());
    }

    ParallelDescriptor::Bcast(fileNameCharPtr, fileNameLengthPadded,0);

    string newFileName(fileNameCharPtr);
    delete [] fileNameCharPtr;

    // make a new DataServices for nonioprocessors
    if( ! ParallelDescriptor::IOProcessor()) {
      ds = new DataServices();
      ds->Init(newFileName, newFileType);

    }

    ds->bAmrDataOk = ds->amrData.ReadData(ds->fileName, ds->fileType);

    if(ds->bAmrDataOk) {
      ds->dsArrayIndex = DataServices::dsArrayIndexCounter;
      ++DataServices::dsArrayIndexCounter;
      DataServices::dsArray.resize(DataServices::dsArrayIndexCounter);
      DataServices::dsArray[ds->dsArrayIndex] = ds;
    } else {
      cerr << "*** Error in DataServices NewRequest:  Bad AmrData." << endl;
      continue;  // go to the top of the while loop
    }

    continue;  // go to the top of the while loop
  }  // end NewRequest


  // handle exit request
  if(requestType == ExitRequest) {                // cleanup memory
    if (!profiler) {
      for(int i(0); i < dsArray.size(); ++i) {
        if(DataServices::dsArray[i] != NULL) {
          BL_ASSERT(DataServices::dsArray[i]->numberOfUsers == 0);
          delete DataServices::dsArray[i];
        }
      }
    }
    ParallelDescriptor::EndParallel();
    exit(0);

  }  // end ExitRequest

  if(ParallelDescriptor::IOProcessor()) {
    va_start(ap, ds);
//    if (!profiler) {
      whichDSIndex = ds->dsArrayIndex;
//    }
  }

//  if (!profiler) {
    ParallelDescriptor::Bcast(&whichDSIndex, 1, 0);
    if( ! ParallelDescriptor::IOProcessor()) {
      ds = DataServices::dsArray[whichDSIndex];
    }
    BL_ASSERT(ds != NULL);
//  }

  ParallelDescriptor::Barrier();

  switch(requestType) {
    case DeleteRequest:
    {
      if (!profiler) {
        bool bDeleteDS(false);
        BL_ASSERT(DataServices::dsArray[whichDSIndex]->numberOfUsers >= 0);
        if(ParallelDescriptor::IOProcessor()) {
	  bDeleteDS = (DataServices::dsArray[whichDSIndex]->numberOfUsers == 0);
        }

        {
          int tmp = bDeleteDS;
          ParallelDescriptor::Bcast(&tmp, 1, 0);
          bDeleteDS = tmp;
        }
        if(bDeleteDS) {
          delete DataServices::dsArray[whichDSIndex];
        }
      }
    }
    break;

    case FillVarOneFab:
    {
      FArrayBox *destFab = NULL;
      Box destBox;
      int fineFillLevel;
      string derivedTemp;
      char *derivedCharPtr;
      int derivedLength, derivedLengthPadded;

      if(ParallelDescriptor::IOProcessor()) {
	destFab = (FArrayBox *) va_arg(ap, void *);
        const Box *boxRef = (const Box *) va_arg(ap, void *);
	destBox = *boxRef;
        fineFillLevel = va_arg(ap, int);
        const string *derivedRef = (const string *) va_arg(ap, void *);
        derivedTemp = *derivedRef;
	derivedLength = derivedTemp.length();
      }

      ParallelDescriptor::Bcast(&destBox, 1, 0);
      ParallelDescriptor::Bcast(&fineFillLevel, 1, 0);
      ParallelDescriptor::Bcast(&derivedLength, 1, 0);

      derivedLengthPadded = derivedLength + 1;
      derivedLengthPadded += derivedLengthPadded % 8;
      derivedCharPtr = new char[derivedLengthPadded];
      if(ParallelDescriptor::IOProcessor()) {
        strcpy(derivedCharPtr, derivedTemp.c_str());
      }

      ParallelDescriptor::Bcast(derivedCharPtr, derivedLengthPadded, 0);

      string derived(derivedCharPtr);
      delete [] derivedCharPtr;

      ds->FillVar(destFab, destBox, fineFillLevel, derived, ioProcNumber);

    }
    break;

    case FillVarArrayOfFabs:
    {
      amrex::Abort("FillVarArrayOfFabs not implemented yet.");
    }
    break;

    case FillVarMultiFab:
    {
      amrex::Abort("FillVarMultiFab not implemented yet.");
    }
    break;

    case WriteFabOneVar:
    {
      // interface: (requestType, dsPtr, fabFileName, box, maxLevel, derivedName)
      Box destBox;
      int fineFillLevel;
      string fabFileName;
      string derivedTemp;
      char *derivedCharPtr;
      int derivedLength, derivedLengthPadded;

      if(ParallelDescriptor::IOProcessor()) {
        const string *fabFileNameRef = (const string *) va_arg(ap, void *);
	fabFileName = *fabFileNameRef;
        const Box *boxRef = (const Box *) va_arg(ap, void *);
	destBox = *boxRef;
        fineFillLevel = va_arg(ap, int);
        const string *derivedRef = (const string *) va_arg(ap, void *);
        derivedTemp = *derivedRef;
	derivedLength = derivedTemp.length();
      }

      ParallelDescriptor::Bcast(&destBox, 1, 0);
      ParallelDescriptor::Bcast(&fineFillLevel, 1, 0);
      ParallelDescriptor::Bcast(&derivedLength, 1, 0);

      derivedLengthPadded = derivedLength + 1;
      derivedLengthPadded += derivedLengthPadded % 8;
      derivedCharPtr = new char[derivedLengthPadded];
      if(ParallelDescriptor::IOProcessor()) {
        strcpy(derivedCharPtr, derivedTemp.c_str());
      }

      ParallelDescriptor::Bcast(derivedCharPtr, derivedLengthPadded, 0);

      string derived(derivedCharPtr);
      delete [] derivedCharPtr;

      ds->WriteFab(fabFileName, destBox, fineFillLevel, derived);

    }
    break;

    case WriteFabAllVars:
    {
      // interface: (requestType, dsPtr, fabFileName, box, maxLevel)
      Box destBox;
      int fineFillLevel;
      string fabFileName;

      if(ParallelDescriptor::IOProcessor()) {
        const string *fabFileNameRef = (const string *) va_arg(ap, void *);
	fabFileName = *fabFileNameRef;
        const Box *boxRef = (const Box *) va_arg(ap, void *);
	destBox = *boxRef;
        fineFillLevel = va_arg(ap, int);
      }

      ParallelDescriptor::Bcast(&destBox, 1, 0);
      ParallelDescriptor::Bcast(&fineFillLevel, 1, 0);

      ds->WriteFab(fabFileName, destBox, fineFillLevel);
    }
    break;

    case DumpSlicePlaneOneVar:
    {
      int slicedir;
      int slicenum;
      string derivedTemp;
      char *derivedCharPtr;
      int derivedLength, derivedLengthPadded;
      if(ParallelDescriptor::IOProcessor()) {
        slicedir = va_arg(ap, int);
        slicenum = va_arg(ap, int);
        const string *derivedRef = (const string *) va_arg(ap, void *);
        derivedTemp = *derivedRef;
	derivedLength = derivedTemp.length();
      }

      ParallelDescriptor::Bcast(&slicedir, 1, 0);
      ParallelDescriptor::Bcast(&slicenum, 1, 0);
      ParallelDescriptor::Bcast(&derivedLength, 1, 0);

      derivedLengthPadded = derivedLength + 1;
      derivedLengthPadded += derivedLengthPadded % 8;
      derivedCharPtr = new char[derivedLengthPadded];
      if(ParallelDescriptor::IOProcessor()) {
        strcpy(derivedCharPtr, derivedTemp.c_str());
      }
      ParallelDescriptor::Bcast(derivedCharPtr, derivedLengthPadded, 0);

      string derived(derivedCharPtr);
      delete [] derivedCharPtr;

      ds->DumpSlice(slicedir, slicenum, derived);

    }
    break;

    case DumpSlicePlaneAllVars:
    {
      int slicedir;
      int slicenum;
      if(ParallelDescriptor::IOProcessor()) {
        slicedir = va_arg(ap, int);
        slicenum = va_arg(ap, int);
      }
      ParallelDescriptor::Bcast(&slicedir, 1, 0);
      ParallelDescriptor::Bcast(&slicenum, 1, 0);

      ds->DumpSlice(slicedir, slicenum);

    }
    break;

    case DumpSliceBoxOneVar:
    {
      Box box;
      string derivedTemp;
      char *derivedCharPtr;
      int derivedLength, derivedLengthPadded;
      if(ParallelDescriptor::IOProcessor()) {
        const Box *boxRef = (const Box *) va_arg(ap, void *);
	box = *boxRef;
        const string *derivedRef = (const string *) va_arg(ap, void *);
	derivedTemp = *derivedRef;
	derivedLength = derivedTemp.length();
      }

      ParallelDescriptor::Bcast(&box, 1, 0);
      ParallelDescriptor::Bcast(&derivedLength, 1, 0);

      derivedLengthPadded = derivedLength + 1;
      derivedLengthPadded += derivedLengthPadded % 8;
      derivedCharPtr = new char[derivedLengthPadded];
      if(ParallelDescriptor::IOProcessor()) {
        strcpy(derivedCharPtr, derivedTemp.c_str());
      }
      ParallelDescriptor::Bcast(derivedCharPtr, derivedLengthPadded, 0);

      string derived(derivedCharPtr);
      delete [] derivedCharPtr;

      ds->DumpSlice(box, derived);
    }
    break;

    case DumpSliceBoxAllVars:
    {
      Box box;
      if(ParallelDescriptor::IOProcessor()) {
        const Box *boxRef = (const Box *) va_arg(ap, void *);
	box = *boxRef;
      }

      ParallelDescriptor::Bcast(&box, 1, 0);

      ds->DumpSlice(box);

    }
    break;

    case MinMaxRequest:
    {
      Box box;
      string derivedTemp;
      char *derivedCharPtr;
      int level;
      int derivedLength, derivedLengthPadded;
      Real dataMin, dataMax;
      bool minMaxValid;
      if(ParallelDescriptor::IOProcessor()) {
        const Box *boxRef = (const Box *) va_arg(ap, void *);
        const string *derivedRef = (const string *) va_arg(ap, void *);
        level = va_arg(ap, int);
	box = *boxRef;
	derivedTemp = *derivedRef;
	derivedLength = derivedTemp.length();
      }

      ParallelDescriptor::Bcast(&box, 1, 0);
      ParallelDescriptor::Bcast(&derivedLength, 1, 0);
      ParallelDescriptor::Bcast(&level, 1, 0);

      derivedLengthPadded = derivedLength + 1;
      derivedLengthPadded += derivedLengthPadded % 8;
      derivedCharPtr = new char[derivedLengthPadded];
      if(ParallelDescriptor::IOProcessor()) {
        strcpy(derivedCharPtr, derivedTemp.c_str());
      }
      ParallelDescriptor::Bcast(derivedCharPtr, derivedLengthPadded, 0);

      string derived(derivedCharPtr);
      delete [] derivedCharPtr;

      ds->MinMax(box, derived, level, dataMin, dataMax, minMaxValid);

      // set the return values
      if(ParallelDescriptor::IOProcessor()) {
        Real *dataMinRef = va_arg(ap, Real *);
        Real *dataMaxRef = va_arg(ap, Real *);
        bool *minMaxValidRef = va_arg(ap, bool *);
	*dataMinRef = dataMin;
	*dataMaxRef = dataMax;
	*minMaxValidRef = minMaxValid;
      }
    }
    break;

    case PointValueRequest:
    {
      // interface: (requestType, dsPtr,
      //             pointBoxArraySize, pointBoxArray *,
      //             derivedName,
      //             coarsestLevelToSearch, finestLevelToSearch,
      //             intersectedLevel,  /* return this value */
      //             intersectedBox,    /* return this value */
      //             dataPointValue,    /* return this value */
      //             bPointIsValid)     /* return this value */

      // need to broadcast pointBoxArraySize, pointBoxArray, derivedName,
      // coarsestLevelToSearch, and finestLevelToSearch

      int pointBoxArraySize;
      Box *pointBoxArrayPtr(NULL), *pointBoxArrayTempPtr(NULL);
      int coarsestLevelToSearch, finestLevelToSearch;

      string derivedTemp;
      char *derivedCharPtr(NULL);
      int derivedLength, derivedLengthPadded;

      if(ParallelDescriptor::IOProcessor()) {
        pointBoxArraySize = va_arg(ap, int);
        pointBoxArrayTempPtr = (Box *) va_arg(ap, void *);
        const string *derivedRef = (const string *) va_arg(ap, void *);
        derivedTemp = *derivedRef;
	derivedLength = derivedTemp.length();
        coarsestLevelToSearch = va_arg(ap, int);
        finestLevelToSearch   = va_arg(ap, int);
      }

      ParallelDescriptor::Bcast(&pointBoxArraySize, 1, 0);
      ParallelDescriptor::Bcast(&derivedLength, 1, 0);
      ParallelDescriptor::Bcast(&coarsestLevelToSearch, 1, 0);
      ParallelDescriptor::Bcast(&finestLevelToSearch, 1, 0);

      pointBoxArrayPtr = new Box[pointBoxArraySize];

      derivedLengthPadded = derivedLength + 1;
      derivedLengthPadded += derivedLengthPadded % 8;
      derivedCharPtr = new char[derivedLengthPadded];
      if(ParallelDescriptor::IOProcessor()) {
        strcpy(derivedCharPtr, derivedTemp.c_str());
	for(int iBox = 0; iBox < pointBoxArraySize; ++iBox) {
	  pointBoxArrayPtr[iBox] = pointBoxArrayTempPtr[iBox];
	}
      }
      ParallelDescriptor::Bcast(derivedCharPtr, derivedLengthPadded, 0);
      ParallelDescriptor::Bcast(pointBoxArrayPtr, pointBoxArraySize, 0);

      string derived(derivedCharPtr);
      delete [] derivedCharPtr;

      // return values
      int intersectedLevel;
      Box intersectedBox;
      Real dataPointValue;
      bool bPointIsValid;

      ds->PointValue(pointBoxArraySize, pointBoxArrayPtr,
		     derived,
		     coarsestLevelToSearch,
		     finestLevelToSearch,
		     intersectedLevel,
		     intersectedBox,
		     dataPointValue,
		     bPointIsValid);

      // set the return values
      if(ParallelDescriptor::IOProcessor()) {
        int *intersectedLevelRef = va_arg(ap, int *);
        Box *intersectedBoxRef   = (Box *) va_arg(ap, void *);
        Real *dataPointValueRef  = va_arg(ap, Real *);
        bool *bPointIsValidRef   = va_arg(ap, bool *);
	*intersectedLevelRef     = intersectedLevel;
	*intersectedBoxRef       = intersectedBox;
	*dataPointValueRef       = dataPointValue;
	*bPointIsValidRef        = bPointIsValid;
      }

      // dont need to broadcast the return values--only the IOProcessor uses them

      delete [] pointBoxArrayPtr;
    }
    break;

    case LineValuesRequest:
    {
      // interface: (requestType, dsPtr,
      //             lineBoxArraySize, lineBoxArray *,
      //             derivedName,
      //             coarsestLevelToSearch, finestLevelToSearch,
      //             dataList,          /* modify this value */
      //             bLineIsValid)      /* return this value */

      // need to broadcast lineBoxArraySize, lineBoxArray, derivedName,
      // coarsestLevelToSearch, and finestLevelToSearch

      int lineBoxArraySize;
      Box *lineBoxArrayPtr(NULL), *lineBoxArrayTempPtr(NULL);
      int coarsestLevelToSearch(-1), finestLevelToSearch(-1), whichDir(-1);
      XYPlotDataList *dataList(NULL);

      string derivedTemp;
      char *derivedCharPtr;
      int derivedLength, derivedLengthPadded;

      if(ParallelDescriptor::IOProcessor()) {
        lineBoxArraySize = va_arg(ap, int);
        lineBoxArrayTempPtr = (Box *) va_arg(ap, void *);
        whichDir = va_arg(ap, int);
        const string *derivedRef = (const string *) va_arg(ap, void *);
        derivedTemp = *derivedRef;
        derivedLength = derivedTemp.length();
        coarsestLevelToSearch = va_arg(ap, int);
        finestLevelToSearch   = va_arg(ap, int);
        dataList = (XYPlotDataList *) va_arg(ap, void *);
      }

      ParallelDescriptor::Bcast(&lineBoxArraySize, 1, 0);
      ParallelDescriptor::Bcast(&derivedLength, 1, 0);
      ParallelDescriptor::Bcast(&coarsestLevelToSearch, 1, 0);
      ParallelDescriptor::Bcast(&finestLevelToSearch, 1, 0);

      lineBoxArrayPtr = new Box[lineBoxArraySize];

      derivedLengthPadded = derivedLength + 1;
      derivedLengthPadded += derivedLengthPadded % 8;
      derivedCharPtr = new char[derivedLengthPadded];
      if(ParallelDescriptor::IOProcessor()) {
        strcpy(derivedCharPtr, derivedTemp.c_str());
        for(int iBox(0); iBox < lineBoxArraySize; ++iBox) {
          lineBoxArrayPtr[iBox] = lineBoxArrayTempPtr[iBox];
        }
      }
      ParallelDescriptor::Bcast(derivedCharPtr, derivedLengthPadded, 0);
      ParallelDescriptor::Bcast(lineBoxArrayPtr, lineBoxArraySize, 0);

      string derived(derivedCharPtr);
      delete [] derivedCharPtr;

      // return values
      bool bLineIsValid;

      ds->LineValues(lineBoxArraySize, lineBoxArrayPtr, whichDir,
                     derived,
                     coarsestLevelToSearch,
                     finestLevelToSearch,
                     dataList,
                     bLineIsValid);

      // set the return values
      if(ParallelDescriptor::IOProcessor()) {
        bool *bLineIsValidRef = va_arg(ap, bool *);
        *bLineIsValidRef      = bLineIsValid;
      }
      // dont need to broadcast the return values--only the IOProcessor uses them

      delete [] lineBoxArrayPtr;

    } 
    break;

    case InvalidRequestType:
    {
      // not used here
    }
    break;

    case ExitRequest:
    {
      // not used here
    }
    break;

    case NewRequest:
    {
      // not used here
    }
    break;

// Profiler Data Requests
#ifdef BL_USE_PROFPARSER
    case InitTimeRanges:
    {
      std::cout << "Dispatch::---- InitTimeRangesDataRequest" << std::endl;
      ds->InitRegionTimeRanges();
    }
    break;

    case CheckProfDataRequest:
    {
      std::cout << "Dispatch::---- CheckProfDataRequest" << std::endl;
      ds->CheckProfData();
    }
    break;

    case WriteSummaryRequest:
    {
      if(ParallelDescriptor::IOProcessor()) { std::cout << "Dispatch::---- WriteSummaryRequest:" << std::endl; }

      std::ofstream *osPtr;
      bool *writeAveragePtr, *useTracePtr, *graphTopPct;
      int whichProc;

      osPtr =  (std::ofstream *) va_arg(ap, void *);
      writeAveragePtr = (bool *) va_arg(ap, bool *);
      whichProc =                va_arg(ap, int);
      useTracePtr =     (bool *) va_arg(ap, bool *);
      graphTopPct =     (bool *) va_arg(ap, bool *);

      ds->WriteSummary(*osPtr, *writeAveragePtr, whichProc, *useTracePtr, *graphTopPct);
    }
    break;

    case RunStatsRequest:
    {
      if(ParallelDescriptor::IOProcessor()) { std::cout << "Dispatch::---- RunStatsRequest:" << std::endl; }

      std::map<int, std::string> *mpiFuncNamesPtr;
      bool *statsCollectedPtr;

      mpiFuncNamesPtr = (std::map<int, std::string> *) va_arg(ap, void *);
      statsCollectedPtr = (bool *) va_arg(ap, void *);

      ds->RunStats(*mpiFuncNamesPtr, *statsCollectedPtr);

      // ---- set return values
      if(ParallelDescriptor::IOProcessor()) {
      }
    }
    break;

    case RunSendsPFRequest:
    {
      if(ParallelDescriptor::IOProcessor()) { std::cout << "Dispatch::---- RunSendsPFRequest:" << std::endl; }

      std::string *plotfileNamePtr;
      int maxSmallImageLength, refRatioAll;
      bool *proxMapPtr;

      if (ParallelDescriptor::IOProcessor()) 
      {
        plotfileNamePtr =  (std::string *) va_arg(ap, void *);
        maxSmallImageLength = va_arg(ap, int);
        proxMapPtr = (bool *) va_arg(ap, bool *);
        refRatioAll = va_arg(ap, int);
      }

       // --- Broadcast Bools
      amrex::BroadcastBool(*proxMapPtr, ParallelDescriptor::MyProc(), ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator());

      // --- Broadcast Ints
      ParallelDescriptor::Bcast(&maxSmallImageLength, 1, ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 
      ParallelDescriptor::Bcast(&refRatioAll, 1, ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 
      
      // --- Broadcast String
      amrex::BroadcastString(*plotfileNamePtr, ParallelDescriptor::MyProc(), ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 

      ds->RunSendsPF(*plotfileNamePtr, maxSmallImageLength, *proxMapPtr, refRatioAll);
    }
    break;

    case RunTimelinePFRequest:
    {
      if(ParallelDescriptor::IOProcessor()) { std::cout << "Dispatch::---- RunTimelinePFRequest:" << std::endl; }

      std::map<int, std::string> mpiFuncNames;
      std::string plotfileName;
      BLProfStats::TimeRange subTimeRange;
      Real start, stop;
      int maxSmallImageLength, refRatioAll, nTimeSlots;
      bool statsCollected;

      // Storage of broken down map.
      Vector<int> mapFirst;
      Vector<std::string> fileNameAndmapSecond;
      Vector<char> serialSecond;

      if (ParallelDescriptor::IOProcessor())
      {
        // Get passed data
        mpiFuncNames = *((std::map<int, std::string> *) va_arg(ap, void *));
        plotfileName = *((std::string *) va_arg(ap, void *));
        subTimeRange = *((BLProfStats::TimeRange *) va_arg(ap, void *)); 
        maxSmallImageLength = va_arg(ap, int);
        refRatioAll = va_arg(ap, int);
        nTimeSlots = va_arg(ap, int);
        statsCollected = *((bool *) va_arg(ap, bool *));

        // Prep data for broadcast
        start = subTimeRange.startTime;
        stop = subTimeRange.stopTime;
        fileNameAndmapSecond.push_back(plotfileName);

        std::map<int, std::string>::iterator it;
        for(it = mpiFuncNames.begin(); it != mpiFuncNames.end(); ++it)
        {
          mapFirst.push_back(it->first);
          fileNameAndmapSecond.push_back(it->second);
        }
        serialSecond = amrex::SerializeStringArray(fileNameAndmapSecond);
      }

      // --- Broadcast Bools
      amrex::BroadcastBool(statsCollected, ParallelDescriptor::MyProc(), ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator());

      // --- Broadcast Ints
      ParallelDescriptor::Bcast(&maxSmallImageLength, 1, ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 
      ParallelDescriptor::Bcast(&refRatioAll, 1, ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 
      ParallelDescriptor::Bcast(&nTimeSlots, 1, ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 

      // --- Broadcast Reals
      ParallelDescriptor::Bcast(&start, 1, ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 
      ParallelDescriptor::Bcast(&stop, 1, ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator()); 

      // --- Broadcast Map as 2 Arrays
      amrex::BroadcastArray(mapFirst, ParallelDescriptor::MyProc(), ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator());
      amrex::BroadcastArray(serialSecond, ParallelDescriptor::MyProc(), ParallelDescriptor::IOProcessorNumber(), ParallelDescriptor::Communicator());

      if(!ParallelDescriptor::IOProcessor())
      {
        subTimeRange = BLProfStats::TimeRange(start, stop);

        fileNameAndmapSecond = amrex::UnSerializeStringArray(serialSecond);
        plotfileName = fileNameAndmapSecond.front();
        for(int i(0); i<mapFirst.size(); ++i)
        {
          mpiFuncNames.insert(std::pair<int, std::string>(mapFirst[i+1], fileNameAndmapSecond[i])); 
        }
      }        

      ds->RunTimelinePF(mpiFuncNames, plotfileName, subTimeRange, maxSmallImageLength,
                         refRatioAll, nTimeSlots, statsCollected);

    }
    break;

    case MakeRegionPltRequest:
    {
      if(ParallelDescriptor::IOProcessor()) { std::cout << "Dispatch::---- MakeRegionPltRequest:" << std::endl; }

      std::string *plotfileNamePtr = (std::string *) va_arg(ap, void *);
      ds->MakeRegionPlt(*plotfileNamePtr);
    }
    break;

    case MakeFilterFileRequest:
    {
      if(ParallelDescriptor::IOProcessor()) { std::cout << "Dispatch::---- MakeFilterFileRequest:" << std::endl; }

      std::string *filterFileNamePtr = (std::string *) va_arg(ap, void *);
      ds->MakeFilterFile(*filterFileNamePtr);
    }
    break;

// May be unnessecary or unwanted, but profiling data has a default.
// Placed here for consistency.

    default:
    break;
#endif

  }  // end switch

  if(ParallelDescriptor::IOProcessor()) {
    va_end(ap);
  }

 }  // end while(bContinueLooping)

  return;
}  // end Dispatch


// ---------------------------------------------------------------
bool DataServices::DumpSlice(int slicedir, int slicenum,
                             const string &varname)
{
  if( ! bAmrDataOk) {
    return false;
  }

  int iWTL(-2);
  if(iWriteToLevel == -1) {
    iWTL = amrData.FinestLevel();
  } else {
    iWTL = iWriteToLevel;
  }

  string sliceFile = fileName;
  sliceFile += ".";
  sliceFile += varname;
  sliceFile += ".";
  if(slicedir == Amrvis::XDIR) {
    sliceFile += "xslice";
  } else if(slicedir == Amrvis::YDIR) {
    sliceFile += "yslice";
  } else if(slicedir == Amrvis::ZDIR) {
    sliceFile += "zslice";
  } else {
    cerr << "bad slicedir = " << slicedir << endl;
    return false;
  }
  sliceFile += ".";
  const int N = 64;
  char slicechar[N];
  if (snprintf(slicechar, N, "%d.Level_%d", slicenum, iWTL) >= N)
    amrex::Abort("DataServices::DumpSlice(1): slicechar buffer too small");
  sliceFile += slicechar;
  sliceFile += ".fab";
  cout << "sliceFile = " << sliceFile << endl;

  Box sliceBox(amrData.ProbDomain()[iWTL]);

  if((BL_SPACEDIM == 2 && slicedir == Amrvis::ZDIR) || (BL_SPACEDIM == 1)) {
    // use probDomain for the sliceBox
  } else {
    // make the box one cell thick in the slice direction
    sliceBox.setSmall(slicedir, slicenum);
    sliceBox.setBig(slicedir, slicenum);
  }

  cout << "sliceBox  = " << sliceBox << endl;
  cout << endl;
  if( ! amrData.ProbDomain()[iWTL].contains(sliceBox)) {
    cerr << "Error:  sliceBox = " << sliceBox << "  slicedir " << slicenum
	 << " on Level " << iWTL
	 << " not in probDomain: " << amrData.ProbDomain()[iWTL]
	 << endl;
    return false;
  }
  bool bWF = WriteFab(sliceFile, sliceBox, iWTL, varname);
  return bWF;
}


// ---------------------------------------------------------------
bool DataServices::DumpSlice(int slicedir, int slicenum) {  // dump all vars
  if( ! bAmrDataOk) {
    return false;
  }

  int iWTL(-2);
  if(iWriteToLevel == -1) {
    iWTL = amrData.FinestLevel();
  } else {
    iWTL = iWriteToLevel;
  }

  string sliceFile = fileName;
  sliceFile += ".";
  if(slicedir == Amrvis::XDIR) {
    sliceFile += "xslice";
  } else if(slicedir == Amrvis::YDIR) {
    sliceFile += "yslice";
  } else if(slicedir == Amrvis::ZDIR) {
    sliceFile += "zslice";
  } else {
    cerr << "bad slicedir = " << slicedir << endl;
    return false;
  }
  sliceFile += ".";
  const int N = 64;
  char slicechar[N];
  if (snprintf(slicechar, N, "%d.Level_%d", slicenum, iWTL) >= N)
    amrex::Abort("DataServices::DumpSlice(2): slicechar buffer too small");
  sliceFile += slicechar;
  sliceFile += ".fab";
  cout << "sliceFile = " << sliceFile << endl;

  Box sliceBox(amrData.ProbDomain()[iWTL]);

  if((BL_SPACEDIM == 2 && slicedir == Amrvis::ZDIR) || (BL_SPACEDIM == 1)) {
    // use probDomain for the sliceBox
  } else {
    // make the box one cell thick in the slice direction
    sliceBox.setSmall(slicedir, slicenum);
    sliceBox.setBig(slicedir, slicenum);
  }

  cout << "sliceBox  = " << sliceBox << endl;
  cout << endl;
  if( ! amrData.ProbDomain()[iWTL].contains(sliceBox)) {
    cerr << "Error:  sliceBox = " << sliceBox << "  slicedir " << slicenum
	 << " on Level " << iWTL
	 << " not in probDomain: " << amrData.ProbDomain()[iWTL]
	 << endl;
    return false;
  }
  bool bWF = WriteFab(sliceFile, sliceBox, iWTL);
  return bWF;
}


// ---------------------------------------------------------------
bool DataServices::DumpSlice(const Box &b, const string &varname) {
  if( ! bAmrDataOk) {
    return false;
  }

  int iWTL(-2);
  if(iWriteToLevel == -1) {
    iWTL = amrData.FinestLevel();
  } else {
    iWTL = iWriteToLevel;
  }

  string sliceFile = fileName;
  sliceFile += ".";
  sliceFile += varname;
  sliceFile += ".";
  const int N = 256;
  char slicechar[N];
# if (BL_SPACEDIM == 1)
  int count = snprintf(slicechar, N, "%d__%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR),
                       b.bigEnd(Amrvis::XDIR), iWTL);
# elif (BL_SPACEDIM == 2)
  int count = snprintf(slicechar, N, "%d_%d__%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR), iWTL);
# else
  int count = snprintf(slicechar, N, "%d_%d_%d__%d_%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR), b.smallEnd(Amrvis::ZDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR),   b.bigEnd(Amrvis::ZDIR), iWTL);
#endif
  if (count >= N) {
    amrex::Abort("DataServices::DumpSlice(3): slicechar buffer too small");      
  }
  sliceFile += slicechar;
  sliceFile += ".fab";
  cout << "sliceFile = " << sliceFile << endl;
  cout << "sliceBox = " << b << endl;
  cout << endl;

  if( ! amrData.ProbDomain()[iWTL].contains(b)) {
    cerr << "Slice box not in probDomain: "
	 << amrData.ProbDomain()[iWTL]
	 << " on Level " << iWTL
	 << endl;
    return false;
  }
  bool bWF = WriteFab(sliceFile, b, iWTL, varname);
  return bWF;
}


// ---------------------------------------------------------------
bool DataServices::DumpSlice(const Box &b) {  // dump all vars
  if( ! bAmrDataOk) {
    return false;
  }

  int iWTL(-2);
  if(iWriteToLevel == -1) {
    iWTL = amrData.FinestLevel();
  } else {
    iWTL = iWriteToLevel;
  }

  string sliceFile = fileName;
  sliceFile += ".";
  const int N = 256;
  char slicechar[N];
# if (BL_SPACEDIM == 1)
  int count = snprintf(slicechar, N, "%d__%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR),
                       b.bigEnd(Amrvis::XDIR), iWTL);
# elif (BL_SPACEDIM == 2)
  int count = snprintf(slicechar, N, "%d_%d__%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR), iWTL);
# else
  int count = snprintf(slicechar, N, "%d_%d_%d__%d_%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR), b.smallEnd(Amrvis::ZDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR),   b.bigEnd(Amrvis::ZDIR), iWTL);
#endif
  if (count >= N)
    amrex::Abort("DataServices::DumpSlice(4): slicechar buffer too small");      
  sliceFile += slicechar;
  sliceFile += ".fab";
  cout << "sliceFile = " << sliceFile << endl;
  cout << "sliceBox = " << b << endl;
  cout << endl;

  if( ! amrData.ProbDomain()[iWTL].contains(b)) {
    cerr << "Slice box not in probDomain: "
	 << amrData.ProbDomain()[iWTL]
	 << " on Level " << iWTL
	 << endl;
    return false;
  }
  bool bWF = WriteFab(sliceFile, b, iWTL);
  return bWF;
}



// ---------------------------------------------------------------
bool DataServices::FillVar(FArrayBox *destFab, const Box &destBox,
			   int finestFillLevel, const string &varname,
			   int procWithFab)
{
  if( ! bAmrDataOk) {
    return false;
  }
  amrData.FillVar(destFab, destBox, finestFillLevel, varname, procWithFab);

  return true;
}  // end FillVar


// ---------------------------------------------------------------
bool DataServices::FillVar(MultiFab &destMultiFab, int finestFillLevel,
			   const string &varname)
{
  if( ! bAmrDataOk) {
    return false;
  }
  amrData.FillVar(destMultiFab, finestFillLevel, varname);

  return true;
}  // end FillVar


// ---------------------------------------------------------------
//
//
// Change this to take an Vector<Box> (or BoxArray?) and create
// a MultiFab and pass the proc number for the multifabs disributionMapping
// to create the FillVared fabs on separate processors
//
//
bool DataServices::WriteFab(const string &fname, const Box &region, int lev,
		            const string &varname)
{
  if( ! bAmrDataOk) {
    return false;
  }
  FArrayBox data;
  if(ParallelDescriptor::IOProcessor()) {
    data.resize(region, 1);
  }

  Vector<FArrayBox *> destFabs(1);
  Vector<Box> destBoxes(1);
  destFabs[0]  = &data;
  destBoxes[0] = region;
  amrData.FillVar(destFabs, destBoxes, lev, varname,
		  ParallelDescriptor::IOProcessorNumber());

  bool bWF(true);
  if(ParallelDescriptor::IOProcessor()) {
    FABio::Format oldFabFormat = FArrayBox::getFormat();
    if(dsFabOutSize == 1) {
      FArrayBox::setFormat(FABio::FAB_ASCII);
    }
    if(dsFabOutSize == 8) {
      FArrayBox::setFormat(FABio::FAB_8BIT);
    }
    if(dsFabOutSize == 32) {
      FArrayBox::setFormat(FABio::FAB_IEEE_32);
    }

    ofstream os;
    os.open(fname.c_str(), ios::out);
    if(os) {
      data.writeOn(os);
      os.close();
    } else {
      cerr << "*** Error:  cannot open file:  " << fname << endl;
      bWF = false;
    }

    FArrayBox::setFormat(oldFabFormat);
  }

  return bWF;
}  // end WriteFab


// ---------------------------------------------------------------
//
//
// Change this to take an Vector<Box> (or BoxArray?) and create
// a MultiFab and pass the proc number for the multifabs disributionMapping
// to create the FillVared fabs on separate processors
//
//
bool DataServices::WriteFab(const string &fname, const Box &region, int lev) {
  if( ! bAmrDataOk) {
    return false;
  }

  // write all fab vars
  FArrayBox tempdata;
  FArrayBox data;
  if(ParallelDescriptor::IOProcessor()) {
    tempdata.resize(region, 1);
    data.resize(region, amrData.NComp());
  }
  for(int ivar = 0; ivar < amrData.NComp(); ++ivar) {
    //amrData.FillVar(tempdata, lev, amrData.PlotVarNames()[ivar]);
    Vector<FArrayBox *> destFabs(1);
    Vector<Box> destBoxes(1);
    destFabs[0]  = &tempdata;
    destBoxes[0] = region;
    amrData.FillVar(destFabs, destBoxes, lev, amrData.PlotVarNames()[ivar],
		    ParallelDescriptor::IOProcessorNumber());
    int srccomp(0);
    int destcomp(ivar);
    int ncomp(1);
    if(ParallelDescriptor::IOProcessor()) {
      data.copy<RunOn::Host>(tempdata, srccomp, destcomp, ncomp);
    }
    amrData.FlushGrids(ivar);
  }

  bool bWF(true);
  if(ParallelDescriptor::IOProcessor()) {
    FABio::Format oldFabFormat = FArrayBox::getFormat();
    if(dsFabOutSize == 1) {
      FArrayBox::setFormat(FABio::FAB_ASCII);
    }
    if(dsFabOutSize == 8) {
      FArrayBox::setFormat(FABio::FAB_8BIT);
    }
    if(dsFabOutSize == 32) {
      FArrayBox::setFormat(FABio::FAB_IEEE_32);
    }

    ofstream os;
    os.open(fname.c_str(), ios::out);
    if(os) {
      data.writeOn(os);
      os.close();
    } else {
      cerr << "*** Error:  cannot open file:  " << fname << endl;
      bWF = false;
    }

    FArrayBox::setFormat(oldFabFormat);
  }

  return bWF;
}  // end WriteFab


// ---------------------------------------------------------------
bool DataServices::CanDerive(const string &name) const {
  if( ! bAmrDataOk) {
    return false;
  }
  return amrData.CanDerive(name);
}


// ---------------------------------------------------------------
bool DataServices::CanDerive(const Vector<string> &names) const {
  if( ! bAmrDataOk) {
    return false;
  }
  return amrData.CanDerive(names);
}


// ---------------------------------------------------------------
// output the list of variables that can be derived
void DataServices::ListDeriveFunc(std::ostream &os) const {
  if( ! bAmrDataOk) {
    return;
  }
  amrData.ListDeriveFunc(os);
}


// ---------------------------------------------------------------
int DataServices::NumDeriveFunc() const {
  return amrData.NumDeriveFunc();
}


// ---------------------------------------------------------------
void DataServices::PointValue(int /*pointBoxArraySize*/, Box *pointBoxArray,
		              const string &currentDerived,
		              int coarsestLevelToSearch,
			      int finestLevelToSearch,
		              int &intersectedLevel,
		              Box &intersectedGrid,
			      Real &dataPointValue,
		              bool &bPointIsValid)
{
  bPointIsValid = false;
  if( ! bAmrDataOk) {
    return;
  }

  intersectedLevel =
	  amrData.FinestContainingLevel(pointBoxArray[finestLevelToSearch],
					finestLevelToSearch);

  if(intersectedLevel < coarsestLevelToSearch) {
    return;
  }

  Box destBox(pointBoxArray[intersectedLevel]);
  if(destBox.numPts() != 1) {
    cout << "Error in DS::PointValue:bad destBox:  " << destBox << endl;
  }
  BL_ASSERT(destBox.numPts() == 1);

  const BoxArray &intersectedBA = amrData.boxArray(intersectedLevel);
  for(int iGrid = 0; iGrid < intersectedBA.size(); ++iGrid) {
    if(destBox.intersects(intersectedBA[iGrid])) {
      intersectedGrid = intersectedBA[iGrid];
      break;
    }
  }

  FArrayBox *destFab = NULL;
  if(ParallelDescriptor::IOProcessor()) {
    destFab = new FArrayBox(destBox, 1);
  }
  amrData.FillVar(destFab, destBox, intersectedLevel, currentDerived,
		  ParallelDescriptor::IOProcessorNumber());

  if(ParallelDescriptor::IOProcessor()) {
    dataPointValue = (destFab->dataPtr())[0];
    delete destFab;
  }
  bPointIsValid = true;

}  // end PointValue


// ---------------------------------------------------------------
void DataServices::LineValues(int /*lineBoxArraySize*/, Box *lineBoxArray, int whichDir,
                              const string &currentDerived,
                              int coarsestLevelToSearch, int finestLevelToSearch,
                              XYPlotDataList *dataList, bool &bLineIsValid) {
  bLineIsValid = false;
  if( ! bAmrDataOk) {
    return;
  }

  for(int lev(coarsestLevelToSearch); lev <= finestLevelToSearch; ++lev) {
    const BoxArray &intersectedBA = amrData.boxArray(lev);
    int numGrids(intersectedBA.size());
    for(int iGrid(0); iGrid != numGrids; ++iGrid) {
      if(lineBoxArray[lev].intersects(intersectedBA[iGrid])) {
        bLineIsValid = true;
        FArrayBox *destFab = NULL;
        Box destFabBox(lineBoxArray[lev] & intersectedBA[iGrid]);
        if(ParallelDescriptor::IOProcessor()) {
          destFab = new FArrayBox(destFabBox, 1);
        }
        amrData.FillVar(destFab, destFabBox,
                        lev, currentDerived,
                        ParallelDescriptor::IOProcessorNumber());
        if(ParallelDescriptor::IOProcessor()) {
          dataList->AddFArrayBox(*destFab, whichDir, lev);
          delete destFab;
        }
      }
    }
  }
}


// ---------------------------------------------------------------
bool DataServices::MinMax(const Box &onBox, const string &derived, int level,
		          Real &dataMin, Real &dataMax, bool &minMaxValid)
{
  minMaxValid =  amrData.MinMax(onBox, derived, level, dataMin, dataMax);
  return minMaxValid;
}


// ---------------------------------------------------------------
#ifdef BL_USE_PROFPARSER
// profiler functions
// ----------------------------------------------------------------------
void DataServices::ParseFilterFile()
{
    bool bIOP(ParallelDescriptor::IOProcessor());
    std::string filterFileName("RegionFilters.txt");

    if( ! (yyin = fopen(filterFileName.c_str(), "r"))) {
      if(bIOP) {
        cerr << "DataServices::ParseFilterFile:  Cannot open file:  " << filterFileName << endl;
      }
    } else {
      yyparse(&regOutputStats_H);
      fclose(yyin);

      regOutputStats_H.InitFilterTimeRanges();
      if(ParallelDescriptor::IOProcessor()) {
        PrintTimeRangeList(regOutputStats_H.GetFilterTimeRanges()[0]);
      }
    }
}


// ----------------------------------------------------------------------
void DataServices::WriteSummary(std::ostream &os, bool bWriteAverage,
                                int whichProc, bool bUseTrace,
				bool graphTopPct)
{
  bool bIOP(ParallelDescriptor::IOProcessor());
  if(bUseTrace) {
    if( ! bTraceDataAvailable) {
      if(bIOP) {
        cout << "ProfDataServices::WriteSummary:  trace data is not available." << std::endl;
      }
      return;
    }
    regOutputStats_H.WriteSummary(os, bWriteAverage, whichProc, graphTopPct);
  } else {
    blProfStats_H.WriteSummary(os, bWriteAverage, whichProc, graphTopPct);
  }
}


// -----------------------------------------------------------------------
void DataServices::CheckProfData()
{
    bool bIOP(ParallelDescriptor::IOProcessor());
    int myProc(ParallelDescriptor::MyProc());
    int nProcs(ParallelDescriptor::NProcs());

    if(bIOP) { cout << endl << "---------------- checking profiling data." << endl; }

    double dstart(ParallelDescriptor::second());

    if(bProfDataAvailable) {
      if(bIOP) { cout << "Checking BLProfStats." << endl; }
      blProfStats_H.CheckData();
    }

    if(bCommDataAvailable) {
      if(bIOP) { cout << "Checking CommProfStats." << endl; }
      Vector<Long> nBMin, nBMax, nRMin, nRMax;
      const Vector<string> &commHeaderFileNames = CommProfStats::GetHeaderFileNames();

      if(myProc < commHeaderFileNames.size()) {
        for(int hfnI(0); hfnI < commHeaderFileNames.size(); ++hfnI) {
          if(myProc == hfnI % nProcs) {
            CommProfStats commOutputStats;
            std::string commFileName_H_nnnn(fileName + '/' + commHeaderFileNames[hfnI]);
            if( ! ( yyin = fopen(commFileName_H_nnnn.c_str(), "r"))) {
              if(bIOP) {
                cerr << "DataServices::CheckProfData:  Cannot open file:  " << commFileName_H_nnnn
                     << "  continuing ...." << endl;
                continue;
              }
            }
            yyparse(&commOutputStats);
            fclose(yyin);
            commOutputStats.CheckCommData(nBMin, nBMax, nRMin, nRMax);
          }
        }
      }
    }

    if(bRegionDataAvailable) {
      if(bIOP) { cout << endl << "---------------- checking regions data." << endl; }

      const Vector<string> &regionsHeaderFileNames = RegionsProfStats::GetHeaderFileNames();
      cout << "# of RegionFiles: " << regionsHeaderFileNames.size() << endl;
      RegionsProfStats regionsOutputStats;

      string regPrefix_H("bl_call_stats_H");
      std::string regFileName_H(fileName + '/' + regPrefix_H);
      if((yyin = fopen(regFileName_H.c_str(), "r"))) {
        yyparse(&regionsOutputStats);
        fclose(yyin);
      } 
      else {
        cerr << "DataServices::CheckProfData: Cannot open file  " << regPrefix_H << endl;
      }

      if(myProc < regionsHeaderFileNames.size()) {
        for(int hfnI(0); hfnI < regionsHeaderFileNames.size(); ++hfnI) {
          if(myProc == hfnI % nProcs) {
            std::string regionsFileName_H_nnnn(fileName + '/' + regionsHeaderFileNames[hfnI]);
            if( ! ( yyin = fopen(regionsFileName_H_nnnn.c_str(), "r"))) {
              if(bIOP) {
                cerr << "DataServices::CheckProfData:  Cannot open file:  " << regionsFileName_H_nnnn
                     << "  continuing ...." << endl;
                continue;
              }
            }
            yyparse(&regionsOutputStats);
            fclose(yyin);
          }
        }
      }
      cout << "Number of regions = " << regionsOutputStats.RegionNumbers().size() << endl;
      regionsOutputStats.CheckRegionsData();
    }
                   
    ParallelDescriptor::Barrier();
    if(bIOP) { cout << "---------------- finished checking profiling data." << endl << endl; }
    if(bIOP) { cout << "Check Time = " << ParallelDescriptor::second() - dstart << " s." << endl; }
}


// ----------------------------------------------------------------------
void DataServices::ProcessGridLog(const std::string &gridlogFileName) {
    if(ParallelDescriptor::IOProcessor()) {
      CommProfStats glOutputStats;
      if( ! ( yyin = fopen(gridlogFileName.c_str(), "r"))) {
        cout << "DataServices::ProcessGridLog:  Cannot open file:  " << gridlogFileName << endl;
      } else {
        cout << "---------------- parsing " << gridlogFileName << endl;
        yyparse(&glOutputStats);
        fclose(yyin);
        cout << endl;

        const std::map<int, Long> &glMap = glOutputStats.GLMap();
        std::map<int, Long>::const_iterator it;
        std::ofstream glout("grdlogRankNPoints.xgr");
        for(it = glMap.begin(); it != glMap.end(); ++it) {
          glout << it->first << ' ' << it->second << '\n';
        }
        glout.close();

        const std::map<int, int> &glSizeMap = glOutputStats.GLSizeMap();
        std::map<int, int>::const_iterator its;
        std::string gridGraphName("grdlogSizeNGrids.xgr");
        std::ofstream glsizeout(gridGraphName);
        cout << "---------------- writing " << gridGraphName << endl;
        for(its = glSizeMap.begin(); its != glSizeMap.end(); ++its) {
          glsizeout << its->first << ' ' << its->second << '\n';
        }
        glsizeout.close();
      }
      cout << "---------------- finished processing " << gridlogFileName << endl;
    }
}


// ----------------------------------------------------------------------
void DataServices::PrintCommStats(std::ostream &os,
                                  bool printHeaderNames)
{
  if( ! bCommDataAvailable) {
    os << "ProfDataServices::PrintCommStats:  comm data is not available." << std::endl;
    return;
  }

  os << "CommProfVersion = " << commOutputStats_H.GetCPVersion() << '\n';
  os << "DataNProcs      = " << commOutputStats_H.GetNProcs() << '\n';
  os << "CSSize          = " << commOutputStats_H.GetCSSize() << '\n';
  if(commOutputStats_H.GetFinestLevel() >= 0) {
    os << "FinestLevel     = " << commOutputStats_H.GetFinestLevel() << '\n';
  }
  if(commOutputStats_H.GetMaxLevel() >= 0) {
    os << "MaxLevel        = " << commOutputStats_H.GetMaxLevel() << '\n';
  }

  if(commOutputStats_H.GetRefRatio().size() > 0) {
    os << "RefRatios:" << '\n';
    for(int i(0); i < commOutputStats_H.GetRefRatio().size(); ++i) {
      os << "  Level[" << i << "] = " << commOutputStats_H.GetRefRatio()[i] << '\n';
    }
  }
  if(commOutputStats_H.GetProbDomain().size() > 0) {
    os << "ProbDomain:" << '\n';
    for(int i(0); i < commOutputStats_H.GetProbDomain().size(); ++i) {
      os << "  Level[" << i << "] = " << commOutputStats_H.GetProbDomain()[i] << '\n';
    }
  }
  if(printHeaderNames) {
    const Vector<string> &headerFileNames = commOutputStats_H.GetHeaderFileNames();
    if(headerFileNames.size() > 0) {
      os << "headerFileNames:" << '\n';
      for(int i(0); i < headerFileNames.size(); ++i) {
        os << "  " << headerFileNames[i] << '\n';
      }
    }
  }
  os << std::flush;
}


// ----------------------------------------------------------------------
void DataServices::RunStats(std::map<int, string> &mpiFuncNames,
                                bool &statsCollected)
{
    bool bIOP(ParallelDescriptor::IOProcessor());
    if( ! bCommDataAvailable) {
      if(bIOP) { cout << "DataServices::RunStats:  comm data is not available." << std::endl; }
      return;
    }

    double dstart(ParallelDescriptor::second());
    int  myProc(ParallelDescriptor::MyProc());
    int  nProcs(ParallelDescriptor::NProcs());
    int nMsgSizes(10000), bytesPerSlot(100);
    int minMsgSize(std::numeric_limits<int>::max());
    int maxMsgSize(std::numeric_limits<int>::min());
    Vector<Long> msgSizes(nMsgSizes, 0);
    Vector<Long> totalFunctionCalls(BLProfiler::NUMBER_OF_CFTS, 0);
    Real timeMin(std::numeric_limits<Real>::max());
    Real timeMax(-std::numeric_limits<Real>::max());
    Real timerTime(0.0);
    Long totalNCommStats(0), totalSentData(0);
    int dataNProcs(BLProfStats::GetNProcs());
    Vector<int> rankNodeNumbers(dataNProcs, 0);
    const Vector<string> &commHeaderFileNames = CommProfStats::GetHeaderFileNames();

    CommProfStats::SetInitDataBlocks(true);
    CommProfStats::InitDataFileNames(commHeaderFileNames);
    CommProfStats::OpenAllStreams(fileName);

    if(myProc < commHeaderFileNames.size()) {
      for(int hfnI(0); hfnI < commHeaderFileNames.size(); ++hfnI) {
        if(myProc == hfnI % nProcs) {
          CommProfStats commOutputStats;
          if(bRegionDataAvailable) {
            commOutputStats.SetRegionTimeRanges(commOutputStats_H.GetRegionTimeRanges());
            commOutputStats.SetFilterTimeRanges(commOutputStats_H.GetFilterTimeRanges());
          }

          cout << endl << myProc << ":  commHeaderFileNames[" << hfnI << "] = "
               << commHeaderFileNames[hfnI] << endl;
          std::string commDataHeaderFileName(fileName + '/' + commHeaderFileNames[hfnI]);

          if( ! ( yyin = fopen(commDataHeaderFileName.c_str(), "r"))) {
            if(bIOP) {
              cerr << "DataServices::RunStats:  Cannot open file:  " << commDataHeaderFileName
                   << "  continuing ...." << endl;
            }
            continue;
          }

          yyparse(&commOutputStats);
          fclose(yyin);
          commOutputStats.ReportStats(totalSentData, totalNCommStats,
                                      totalFunctionCalls, bytesPerSlot,
                                      msgSizes, minMsgSize, maxMsgSize,
                                      timeMin, timeMax, timerTime, rankNodeNumbers);
        }
      }
    }

    CommProfStats::CloseAllStreams();

    ParallelDescriptor::ReduceLongSum(totalSentData);
    ParallelDescriptor::ReduceLongSum(totalNCommStats);
    ParallelDescriptor::ReduceLongSum(&totalFunctionCalls[0], totalFunctionCalls.size());
    ParallelDescriptor::ReduceLongSum(&msgSizes[0], msgSizes.size());
    ParallelDescriptor::ReduceIntMin(minMsgSize);
    ParallelDescriptor::ReduceIntMax(maxMsgSize);
    ParallelDescriptor::ReduceRealMin(timeMin);
    ParallelDescriptor::ReduceRealMax(timeMax);
    ParallelDescriptor::ReduceRealSum(timerTime);
    ParallelDescriptor::ReduceIntSum(&rankNodeNumbers[0], rankNodeNumbers.size());
    timerTime /= nProcs;

    if(bIOP) {
      const Vector<std::list<BLProfStats::TimeRange> > &ftRange =
                                    commOutputStats_H.GetFilterTimeRanges();
      if( ! ftRange.empty()) {
        if( ! ftRange[0].empty()) {
          timeMin = std::max(timeMin, ftRange[0].front().startTime);
          timeMax = std::min(timeMax, ftRange[0].back().stopTime);
        }
      }
      cout << endl;
      cout << "---------------- Communication data." << endl;

      PrintCommStats(cout, false);
      const std::locale oldLoc(cout.std::ios_base::getloc());
      cout.imbue(std::locale(""));
      cout << endl;
      cout << "Total Sent Data for all procs = "
           << totalSentData / 1.0e+09 << " GB" << endl;
      cout << "MsgSizeMin MsgSizeMax   = " << minMsgSize << "  "
           << maxMsgSize << endl;
      cout << "Total Comm Stats = " << totalNCommStats << endl;
      cout << "Total Function Calls = " << endl;
      for(int i(0); i < totalFunctionCalls.size(); ++i) {
        if(totalFunctionCalls[i] > 0) {
          const BLProfiler::CommFuncType cft = static_cast<BLProfiler::CommFuncType> (i);
          string fName(BLProfiler::CommStats::CFTToString(cft));
          cout << "  " << fName << "  " << totalFunctionCalls[i] << endl;
          mpiFuncNames.insert(std::make_pair(i, fName));
        }
      }
      cout << "TimeMin TimeMax   = " << timeMin << "  " << timeMax << endl;
      cout << "Average TimerTime = " << timerTime << endl;

      cout << "Stats Time = " << ParallelDescriptor::second() - dstart << " s." << endl;
      cout << "---------------- End Communication data." << endl;
      cout << endl;
      cout.imbue(oldLoc);

      std::string msgSFileName("msgSizes.xgr");
      std::ofstream msgsizesmout(msgSFileName);
      cout << "---------------- Writing:  " << msgSFileName << endl;
      for(int i(0); i < msgSizes.size(); ++i) {
        msgsizesmout << i << " " << msgSizes[i] << endl;
      }
      msgsizesmout.close();

      std::string pnnFName("rankNodeNumbers.xgr");
      std::ofstream pnnout(pnnFName);
      cout << "---------------- Writing:  " << pnnFName << endl;
      for(int ip(0); ip < dataNProcs; ++ip) {
        pnnout << ip << " " << rankNodeNumbers[ip] << '\n';
      }
      pnnout.close();

    }
    statsCollected = true;
}


// ----------------------------------------------------------------------
void DataServices::RunSendsPF(std::string &plotfileName,
                                  int maxSmallImageLength,
                                  bool proxMap, int refRatioAll)
{
#if (BL_SPACEDIM != 2)
  cout << "**** Error:  DataServices::RunSendsPF is only supported for 2D" << endl;
#else
    bool bIOP(ParallelDescriptor::IOProcessor());
    //int  myProc(ParallelDescriptor::MyProc());
    int  nProcs(ParallelDescriptor::NProcs());
    if( ! bCommDataAvailable) {
      if(bIOP) {
        cout << "DataServices::RunSendsPF:  comm data is not available." << endl;
      }
      return;
    }

    BL_PROFILE_VAR("runSendsPF_ALL", runsendspfall)
    BL_PROFILE_VAR("runSendsPF_TOP", runsendspftop)
    double dstart(ParallelDescriptor::second());
    if(bIOP) { cout << endl << "---------------- Process Sends PF." << endl; }

    int dataNProcs(BLProfStats::GetNProcs());
    Long totalSends(0), totalSentData(0);
    Vector<Long> totalSendsPerProc(dataNProcs, 0);
    Vector<Long> totalSentDataPerProc(dataNProcs, 0);

    IntVect maxGrid((dataNProcs / nProcs) + refRatioAll, dataNProcs);
    int nLevels(1), finestLevel(0), numState(2), nGrow(0);
    Box dnpBox(IntVect(0, 0), IntVect(dataNProcs - 1, dataNProcs - 1));

    Box dnpBoxBlocked(dnpBox);
    if(bIOP) { cout << ")))) dnpBoxBlocked = " << dnpBoxBlocked << endl; }
    while(dnpBoxBlocked.length(XDIR) > maxSmallImageLength) {
      dnpBoxBlocked.coarsen(refRatioAll);
      if(bIOP) { cout << ")))) coarsened dnpBoxBlocked = " << dnpBoxBlocked << endl; }
      ++nLevels;
    }
    finestLevel = nLevels - 1;
    Vector<Box> probDomain(nLevels);
    for(int i(0); i < nLevels; ++i) {
      probDomain[i] = dnpBoxBlocked;
      if(bIOP) { cout << ")))) probDomain[" << i << "] =  " << probDomain[i] << endl; }
      dnpBoxBlocked.refine(refRatioAll);
    }

    BoxArray dnpBoxArray(probDomain[finestLevel]);
    dnpBoxArray.maxSize(maxGrid);
    if(dnpBoxArray.size() != nProcs) {
      if(bIOP) {
        cout << "---->>>> nGrids ! = nProcs:  " << nProcs << "  "
             << dnpBoxArray.size() << endl;
      }
    }

    if(bIOP) {
      cout << "Filling finest level sendPF:  nGrids = "
           << dnpBoxArray.size() << endl;
    }

    Vector<MultiFab> state(nLevels);
    const DistributionMapping stateDM(dnpBoxArray);
    state[finestLevel].define(dnpBoxArray, stateDM, numState, nGrow);
    MultiFab &sendMF = state[finestLevel];
    sendMF.setVal(0.0);
    if(bIOP) cout << "DistMap:  " << sendMF.DistributionMap() << endl;
    BL_PROFILE_VAR_STOP(runsendspftop)

    BL_PROFILE_VAR("runSendsPF_PP", runsendspfpp)
    const Vector<string> &commHeaderFileNames = CommProfStats::GetHeaderFileNames();
    CommProfStats::SetInitDataBlocks(true);
    CommProfStats::InitDataFileNames(commHeaderFileNames);
    CommProfStats::OpenAllStreams(fileName);
    BL_PROFILE_VAR_STOP(runsendspfpp)

    for(int hfnI_I(0); hfnI_I < commHeaderFileNames.size(); ++hfnI_I) {
          // everyone reads all headers
          //int hfnI((hfnI_I + myProc) % commHeaderFileNames.size());  // cycle reads
          int hfnI(hfnI_I);
          //cout << myProc << ":  parsing file # " << hfnI << endl;

          CommProfStats commOutputStats;

          std::string commDataHeaderFileName(fileName + '/' + commHeaderFileNames[hfnI]);

          if( ! ( yyin = fopen(commDataHeaderFileName.c_str(), "r"))) {
            if(bIOP) {
              cerr << "DataServices::RunSendsPF:  Cannot open file:  " << commDataHeaderFileName
                   << " ... continuing." << endl;
            }
            continue;
          }

          yyparse(&commOutputStats);
          fclose(yyin);

          if(bRegionDataAvailable) {
            commOutputStats.SetRegionTimeRanges(commOutputStats_H.GetRegionTimeRanges());
            commOutputStats.SetFilterTimeRanges(commOutputStats_H.GetFilterTimeRanges());
          }

          for(MFIter mfi(sendMF); mfi.isValid(); ++mfi) {
            FArrayBox &sendFAB  = sendMF[mfi.index()];

            commOutputStats.FillSendFAB(totalSends, totalSentData,
                                    totalSendsPerProc, totalSentDataPerProc,
                                    sendFAB, proxMap);
          }
    }

    CommProfStats::CloseAllStreams();

    ParallelDescriptor::Barrier();
    ParallelDescriptor::ReduceLongSum(totalSends);
    ParallelDescriptor::ReduceLongSum(totalSentData);
    ParallelDescriptor::ReduceLongSum(&totalSendsPerProc[0], totalSendsPerProc.size());
    ParallelDescriptor::ReduceLongSum(&totalSentDataPerProc[0], totalSentDataPerProc.size());

    if(bIOP) cout << "Finished filling finest level sendMF." << endl;

    // ---- now the real data is in the multifab, avg down to a reasonable size
    Vector<Vector<int> > adRefRatio(state.size() - 1);
    for(int i(0); i < finestLevel; ++i) {
      adRefRatio[i].resize(BL_SPACEDIM, refRatioAll);
    }

    for(int cLev(finestLevel - 1); cLev >= 0; --cLev) {
      if(bIOP) cout << "Averaging down level " << cLev << endl;
      BoxArray ba(BoxArray(state[cLev + 1].boxArray()).coarsen(adRefRatio[cLev][0]));
      // ---- call uniqify, otherwise ba is just a reference to the
      // ---- original boxarray with a coarsening factor
      ba.uniqify();
      if( ! ba.isDisjoint()) {
        if(bIOP) cout << "BA:  Coarsened BoxArray not disjoint:  " << ba << endl;
        SimpleRemoveOverlap(ba);
        if(bIOP) cout << "BA:  Coarsened BoxArray after removeOverlap:  " << ba << endl;
        if( ! ba.isDisjoint()) {
          if(bIOP) cout << "BA:  Coarsened BoxArray still not disjoint:  " << ba << endl;
        }
      }
      if(bIOP) cout << "Coarsened BoxArray size = " << ba.size() << endl;
      state[cLev].define(ba, stateDM, numState, nGrow);
      state[cLev].setVal(0.0);

      avgDown(state[cLev], state[cLev+1], 0, 0, numState, adRefRatio[cLev]);
    }
    ParallelDescriptor::Barrier();

    bool bCreateSwapPairs(false);
    if(bCreateSwapPairs) {   // ---- create a list of swap pairs

    Real mfMin(sendMF.min(0));
    Real mfMax(sendMF.max(0));
    Real pctHighVals(0.10);
    Real highVal(mfMax - (pctHighVals * (mfMax - mfMin)));
    if(bIOP) cout << "mfMin mfMax highVal = " << mfMin << "  " << mfMax << "  " << highVal << endl;
    if(bIOP) cout << "dnpBox = " << dnpBox << endl;
    Vector<int> dmapArray(1);
    dmapArray[0] = ParallelDescriptor::IOProcessorNumber();
    DistributionMapping dmapOneProc(dmapArray);
    Vector<int> swapPairs;
    std::set<int> availableIndicies;
    for(int iY(0); iY < dnpBox.length(YDIR); ++iY) {
      availableIndicies.insert(iY);
    }
    for(int iX(dnpBox.smallEnd(XDIR)); iX <= dnpBox.bigEnd(XDIR); ++iX) {
      std::multimap<Real, int> nCallMap;
      Box b(dnpBox);
      b.setSmall(XDIR, iX);
      b.setBig(XDIR, iX);
      BoxArray bA(b);
      MultiFab mfLine(bA, dmapOneProc, 1, 0);
      mfLine.copy(sendMF);
      Real value(0.0);
      if(bIOP) {
        for(int iLine(0); iLine < b.length(YDIR); ++iLine) {
          value = mfLine[0].dataPtr(0)[iLine];
          if(iLine != iX) {
            nCallMap.insert(std::pair<Real, int>(value, iLine));
          }
          if(value >= highVal) {
            if(bIOP) cout << "]]]] mfLine[" << iLine << "] = " << value << endl;
          }
        }
        cout << "----------------" << endl;

        std::multimap<Real,int>::iterator ubMIter, mIter;
        for(ubMIter = nCallMap.upper_bound(highVal); ubMIter != nCallMap.end(); ++ubMIter) {
          int highSwap(ubMIter->second);
          std::set<int>::iterator aIIH = availableIndicies.find(highSwap);
          if(aIIH != availableIndicies.end()) {  // ---- find a proc to swap with
            for(mIter = nCallMap.begin(); mIter != nCallMap.end(); ++mIter) {
              int lowSwap(mIter->second);
              std::set<int>::iterator aIIL = availableIndicies.find(lowSwap);
              if(aIIL != availableIndicies.end() && mIter->first < highVal) {
                swapPairs.push_back(highSwap);
                cout << "pushing:  " << highSwap << endl;
                swapPairs.push_back(lowSwap);
                cout << "pushing:  " << lowSwap << endl;
                availableIndicies.erase(highSwap);
                availableIndicies.erase(lowSwap);
                break;
              }
            }
          }
        }
      }
    }
    if(bIOP) {
      cout << "++++ swapPairs.size() = " << swapPairs.size() << endl;
      for(int isp(0); isp < swapPairs.size(); isp += 2) {
        cout << "---- " << swapPairs[isp] << " <<-->> " << swapPairs[isp+1] << endl;
      }
    }
    if(bIOP) {
      std::ofstream spo("SwapPairs.txt");
      spo << swapPairs.size() << '\n';
      for(int isp(0); isp < swapPairs.size(); isp += 2) {
        spo << swapPairs[isp] << ' ' << swapPairs[isp+1] << '\n';
      }
      spo.close();
    }

    }

    // ---- copy data into a more ddio favorable configuration
    int sqMG(128);
    Vector<MultiFab> sqState(nLevels);
    for(int i(0); i < state.size(); ++i) {
      BoxArray sqBA(state[i].boxArray().minimalBox());
      sqBA.maxSize(sqMG);
      const DistributionMapping sqDM(sqBA);
      sqState[i].define(sqBA, sqDM, numState, nGrow);
      sqState[i].setVal(0.0);
      sqState[i].copy(state[i]);
    }


    // ---- write the data as a plot file
    BL_PROFILE_VAR("WritePlotfile", writeplotfile);
    if(bIOP) { cout << "Writing plotfile:  " << plotfileName << endl; }

    std::string plotFileVersion("NavierStokes-V1.1");
    Real time(0.0);
    Vector<Real> probLo(BL_SPACEDIM, 0.0);
    Vector<Real> probHi(BL_SPACEDIM, 1.0);
    Vector<int>  refRatioPerLevel(sqState.size() - 1);
    for(int i(0); i < refRatioPerLevel.size(); ++i) {
      refRatioPerLevel[i] = refRatioAll;
    }
    Vector<Vector<Real> > dxLevel(sqState.size());
    for(int i(0); i < sqState.size(); ++i) {
      dxLevel[i].resize(BL_SPACEDIM);
      dxLevel[i][0] = 1.0;
      dxLevel[i][1] = 1.0;
    }
    int coordSys(0);
    Vector<std::string> inVarNames(numState);
    inVarNames[0] = "totalSendsP2P";
    inVarNames[1] = "totalSentDataP2P";
    bool verb(false);
    FABio::Format oldFabFormat(FArrayBox::getFormat());
    FArrayBox::setFormat(FABio::FAB_NATIVE_32);

    WritePlotfile(plotFileVersion, sqState, time,
                  probLo, probHi, refRatioPerLevel,
                  probDomain, dxLevel, coordSys,
                  plotfileName, inVarNames, verb);

    FArrayBox::setFormat(oldFabFormat);
    BL_PROFILE_VAR_STOP(writeplotfile);

    if(bIOP) {
      cout << "*******************************************************" << endl;
      PrintCommStats(cout, false);
      const std::locale oldLoc(cout.std::ios_base::getloc());
      cout.imbue(std::locale(""));
      cout << endl;
      cout << "Total Sends for all procs = " << totalSends << endl;
      cout << "Total Sent Data for all procs = "
           << totalSentData / 1.0e+09 << " GB" << endl;
      cout.imbue(oldLoc);

      std::ofstream tsppout("totalSendsPerProc.xgr");
      for(int i(0); i < totalSendsPerProc.size(); ++i) {
        tsppout << i << " " << totalSendsPerProc[i] << '\n';
      }
      tsppout.close();

      std::ofstream tsdppout("totalSentDataPerProc.xgr");
      for(int i(0); i < totalSentDataPerProc.size(); ++i) {
        tsdppout << i << " " << totalSentDataPerProc[i] << '\n';
      }
      tsdppout.close();

      cout << "Process Sends Time = " << ParallelDescriptor::second() - dstart << " s." << endl;
      cout << "---------------- End Process Sends MF." << endl << endl;
    }
    BL_PROFILE_VAR_STOP(runsendspfall)
#endif
}

// ----------------------------------------------------------------------
BLProfStats::TimeRange DataServices::FindCalcTimeRange()
{
    BL_PROFILE("DataServices::FindCalcTimeRange()");

    bool bIOP(ParallelDescriptor::IOProcessor());
    int  myProc(ParallelDescriptor::MyProc());
    int  nProcs(ParallelDescriptor::NProcs());
 
    const Vector<string> &commHeaderFileNames = CommProfStats::GetHeaderFileNames();
    BLProfStats::TimeRange calcTimeRange(std::numeric_limits<Real>::max(), -std::numeric_limits<Real>::max());

    // find the calc's min and max times.  the user could set these, too.
    if(myProc < commHeaderFileNames.size()) {
      for(int hfnI(0); hfnI < commHeaderFileNames.size(); ++hfnI) {
        if(myProc == hfnI % nProcs) {
          CommProfStats commOutputStats;
          std::string commDataHeaderFileName(fileName + '/' + commHeaderFileNames[hfnI]);

          if( ! ( yyin = fopen(commDataHeaderFileName.c_str(), "r"))) {
            if(bIOP) {
              cerr << "DataServices::RunTimelinePF:  1:  Cannot open file:  " << commDataHeaderFileName
                   << " ... continuing." << endl;
            }
            continue;
          }
          yyparse(&commOutputStats);
          fclose(yyin);
          commOutputStats.FindTimeRange(calcTimeRange);
        }
      }
    }
    ParallelDescriptor::ReduceRealMin(calcTimeRange.startTime);
    ParallelDescriptor::ReduceRealMax(calcTimeRange.stopTime);
    if(bIOP) {
      cout << "++++ calcTimeRange = " << calcTimeRange << endl;
    }

    return calcTimeRange;
}

// ----------------------------------------------------------------------
void DataServices::RunTimelinePF(std::map<int, string> &mpiFuncNames,
                                     std::string &plotfileName,
                                     BLProfStats::TimeRange &subTimeRange,
                                     int maxSmallImageLength,
                                     int refRatioAll, int nTimeSlots,
	                             bool &statsCollected)
{
    BL_PROFILE("DataServices::RunTimelinePF()");

#if (BL_SPACEDIM != 2)
  cout << "**** Error:  DataServices::RunTimelinePF is only supported for 2D" << endl;
#else
    bool bIOP(ParallelDescriptor::IOProcessor());
    int  nProcs(ParallelDescriptor::NProcs());
    int dataNProcs(BLProfStats::GetNProcs());

    if(BL_SPACEDIM != 2) {
      if(bIOP) { cout << "DataServices::RunTimelinePF only supported for 2D." << endl; }
      return;
    }

    if( ! bCommDataAvailable) {
      if(bIOP) {
        cout << "DataServices::RunTimelinePF:  comm data is not available." << std::endl;
      }
      return;
    }

    if(bIOP) { cout << endl << "---------------- Timeline." << endl; }

    const Vector<string> &commHeaderFileNames = CommProfStats::GetHeaderFileNames();
    if( ! statsCollected) {
      RunStats(mpiFuncNames, statsCollected);
    }

    double dstart(ParallelDescriptor::second());

    int rankMin(0), rankMax(dataNProcs - 1), rankStride(1);

    IntVect maxGrid(nTimeSlots, (dataNProcs / nProcs) + 1);
    int nLevels(1), finestLevel(0), numState(2), nGrow(0);
    Box dnpBox(IntVect(0, 0), IntVect(nTimeSlots - 1, dataNProcs - 1));
    Box dnpBoxBlocked(dnpBox);
    if(bIOP) cout << ")))) dnpBoxBlocked = " << dnpBoxBlocked << endl;

    while(dnpBoxBlocked.length(XDIR) > maxSmallImageLength) {
      dnpBoxBlocked.coarsen(refRatioAll);
      if(bIOP) cout << ")))) coarsened dnpBoxBlocked = " << dnpBoxBlocked << endl;
      ++nLevels;
    }
    finestLevel = nLevels - 1;
    Vector<Box> probDomain(nLevels);
    for(int i(0); i < nLevels; ++i) {
      probDomain[i] = dnpBoxBlocked;
      if(bIOP) cout << ")))) probDomain[" << i << "] =  " << probDomain[i] << endl;
      dnpBoxBlocked.refine(refRatioAll);
    }

    CommProfStats::SetInitDataBlocks(true);
    CommProfStats::InitDataFileNames(commHeaderFileNames);
    CommProfStats::OpenAllStreams(fileName);

    Vector<MultiFab> state(nLevels);
    IntVect cRR(1, 1);
    Vector<string> nameTagNames, barrierNames;
    Real ntnMultiplier(0.0), bnMultiplier(0.0);
    Vector<Real> ntnNumbers, bnNumbers;

    for(int iLevel(finestLevel); iLevel >= 0; --iLevel) {
      BoxArray dnpBoxArray(probDomain[iLevel]);
      IntVect levelMaxGrid(maxGrid);
      levelMaxGrid[YDIR] /= cRR[YDIR];
      levelMaxGrid[YDIR] += 1;
      //dnpBoxArray.maxSize(maxGrid / cRR);
      dnpBoxArray.maxSize(levelMaxGrid);
      if(bIOP) cout << "---->>>> iLevel nGrids = " << iLevel << "  "
                    << dnpBoxArray.size() << endl;
      if(dnpBoxArray.size() != nProcs) {
        if(bIOP) cout << "---->>>> nGrids ! = nProcs:  " << nProcs << "  "
                      << dnpBoxArray.size() << endl;
      }

      if(bIOP) cout << "************ dnpBoxArray = " << dnpBoxArray << endl;
      const DistributionMapping dnpDM(dnpBoxArray);
      state[iLevel].define(dnpBoxArray, dnpDM, numState, nGrow);

      MultiFab &timelineMF = state[iLevel];
      timelineMF.setVal(-1.0, 0, 1);       // Timeline initialization
      timelineMF.setVal(0.0, 1, 1);        // MPI count initialization

      for(int hfnI_I(0); hfnI_I < commHeaderFileNames.size(); ++hfnI_I) {
        // everyone reads all headers
        //int hfnI((hfnI_I + myProc) % commHeaderFileNames.size());  // cycle reads
        int hfnI(hfnI_I);

        CommProfStats commOutputStats;
        if(bRegionDataAvailable) {
//          commOutputStats.SetRegionTimeRanges(commOutputStats_H.GetRegionTimeRanges());
          commOutputStats.SetFilterTimeRanges(commOutputStats_H.GetFilterTimeRanges());
        }

        std::string commDataHeaderFileName(fileName + '/' + commHeaderFileNames[hfnI]);

        if( ! ( yyin = fopen(commDataHeaderFileName.c_str(), "r"))) {
          if(bIOP) {
            cerr << "DataServices::RunTimelinePF:  2:  Cannot open file:  " << commDataHeaderFileName
                 << " ... continuing." << endl;
          }
          continue;
        }

        yyparse(&commOutputStats);
        fclose(yyin);
 
        if(hfnI_I == 0) {  // this assumes all headers have the same nametag and barrier names
          // ---- this part encodes the name tag name into the NameTag cfType value
          nameTagNames = commOutputStats.NameTagNames();
          int ntnSize(nameTagNames.size());
          ntnNumbers.resize(ntnSize, 0.0);
          if(ntnSize > 0) {
            ntnMultiplier = pow(10, static_cast<int>( 1 + log10(ntnSize)));
          }
          for(int i(0); i < ntnSize; ++i) {
            if(ntnMultiplier > 0.0) {
              ntnNumbers[i] = i / ntnMultiplier;
            }
          }

          // ---- this part encodes the barrier name into the Barrier cfType value
          barrierNames = commOutputStats.BarrierNames();
          int bnSize(barrierNames.size());
          bnNumbers.resize(bnSize, 0.0);
          if(bnSize > 0) {
            bnMultiplier = pow(10, static_cast<int>( 1 + log10(bnSize)));
          }
          for(int i(0); i < bnSize; ++i) {
            if(bnMultiplier > 0.0) {
              bnNumbers[i] = i / bnMultiplier;
            }
          }
        }

        for(MFIter mfi(timelineMF); mfi.isValid(); ++mfi) {
          FArrayBox &timelineFAB = timelineMF[mfi.index()];
          commOutputStats.TimelineFAB(timelineFAB, probDomain[iLevel],
                                  subTimeRange, rankMin, rankMax,
                                  rankStride * cRR[YDIR],
                                  ntnMultiplier, ntnNumbers, bnMultiplier, bnNumbers);
        }
      }
      cRR *= refRatioAll;
      if(bIOP) { cout << "Finished filling level " << iLevel << " timelineMF." << endl; }
    }
    CommProfStats::CloseAllStreams();
    ParallelDescriptor::Barrier();

    // ---- copy data into a more ddio favorable configuration
    int sqMG(128);
    Vector<MultiFab> sqState(nLevels);
    for(int i(0); i < state.size(); ++i) {
      BoxArray sqBA(state[i].boxArray().minimalBox());
      sqBA.maxSize(sqMG);
      const DistributionMapping sqDM(sqBA);
      sqState[i].define(sqBA, sqDM, numState, nGrow);
      sqState[i].setVal(0.0);
      sqState[i].copy(state[i]);
    }


    // ---- write the data as a plot file
    BL_PROFILE_VAR("WritePlotfile", writeplotfile);
    if(bIOP) { cout << "Writing plotfile:  " << plotfileName << endl; }

    std::string plotFileVersion("CommProfTimeline-V1.0");
    Real time(subTimeRange.stopTime);
    Vector<Real> probLo(BL_SPACEDIM, 0.0);
    Vector<Real> probHi(BL_SPACEDIM, 1.0);
    Vector<int>  refRatioPerLevel(sqState.size() - 1);
    for(int i(0); i < refRatioPerLevel.size(); ++i) {
      refRatioPerLevel[i] = refRatioAll;
    }
    Vector<Vector<Real> > dxLevel(sqState.size());
    for(int i(0); i < sqState.size(); ++i) {
      dxLevel[i].resize(BL_SPACEDIM);
      dxLevel[i][0] = 1.0;
      dxLevel[i][1] = 1.0;
    }
    int coordSys(0);
    Vector<std::string> inVarNames(numState);
    inVarNames[0] = "timeline";
    inVarNames[1] = "mpiCount";
    bool verb(false);
    FABio::Format oldFabFormat(FArrayBox::getFormat());
    FArrayBox::setFormat(FABio::FAB_NATIVE_32);

    WritePlotfile(plotFileVersion, sqState, time,
                  probLo, probHi, refRatioPerLevel,
                  probDomain, dxLevel, coordSys,
                  plotfileName, inVarNames, verb);

    FArrayBox::setFormat(oldFabFormat);

    if(bIOP) {
      string fnoutFileName(plotfileName + "/MPIFuncNames.txt");
      std::ofstream fnout(fnoutFileName.c_str());
      for(std::map<int, std::string>::const_iterator it = mpiFuncNames.begin();
          it != mpiFuncNames.end(); ++it)
      {
        fnout << it->first << ' ' << it->second << '\n';
      }
      fnout.close();
    }
    if(bIOP) {
      string fnoutFileName(plotfileName + "/NameTagNames.txt");
      std::ofstream fnout(fnoutFileName.c_str());
      int ntnSize(nameTagNames.size());
      fnout << ntnSize << ' ' << ntnMultiplier << '\n';
      for(int i(0); i < ntnSize; ++i) {
        //fnout << ntnNumbers[i] << ' ' << nameTagNames[i] << '\n';
        fnout << nameTagNames[i] << '\n';
      }
      fnout.close();
    }
    if(bIOP) {
      string fnoutFileName(plotfileName + "/BarrierNames.txt");
      std::ofstream fnout(fnoutFileName.c_str());
      int bnSize(barrierNames.size());
      fnout << bnSize << ' ' << bnMultiplier << '\n';
      for(int i(0); i < bnSize; ++i) {
        //fnout << bnNumbers[i] << ' ' << barrierNames[i] << '\n';
        fnout << barrierNames[i] << '\n';
      }
      fnout.close();
    }
    if(bIOP) {
      string fnoutFileName(plotfileName + "/CallTrace.txt");
      WriteTextTrace(fnoutFileName);
    }
    if(bIOP) {
      string fnoutFileName(plotfileName + "/TimeRange.txt");
      std::ofstream fnout(fnoutFileName.c_str());
      fnout << std::setprecision(16)
            << subTimeRange.startTime << ' ' << subTimeRange.stopTime << '\n';
      fnout.close();
    }

    BL_PROFILE_VAR_STOP(writeplotfile);

    if(bIOP) {
      cout << "Timeline Time = " << ParallelDescriptor::second() - dstart << " s." << endl;
      cout << "------------------------------------ End timeline." << endl;
      cout << endl;
    }
#endif
}


// ----------------------------------------------------------------------
void DataServices::MakeFilterFile(std::string &fFileName)
{
  bool bIOP(ParallelDescriptor::IOProcessor());
  if(bIOP) {
    blProfStats_H.MakeFilterFile(fFileName);
  }
}


// ----------------------------------------------------------------------
void DataServices::MakeRegionPlt(std::string &plotfileName)
{
    bool bIOP(ParallelDescriptor::IOProcessor());
    //int  myProc(ParallelDescriptor::MyProc());
    //int  nProcs(ParallelDescriptor::NProcs());
    //int dataNProcs(BLProfStats::GetNProcs());

    if( ! bRegionDataAvailable) {
      if(bIOP) {
        cout << "DataServices::MakeRegionPlt:  region data is not available." << std::endl;
      }
      return;
    }

    FArrayBox rFab;
    int noRegionNumber(0);
    std::string rname("\"__NoRegion__\"");
    const std::map<std::string, int> &regionNamesH = regOutputStats_H.RegionNames();
    std::map<std::string, int>::const_iterator itn = regionNamesH.find(rname);
    if(itn == regionNamesH.end()) {
      cout << "rname not found:  " << rname << endl;
    } else {
      cout << "found rname:  " << rname << "  rnum = " << itn->second << endl;
      noRegionNumber = itn->second;
    }

    Vector<Vector<Box>> regionBoxArray;
    regOutputStats_H.MakeRegionPlt(rFab, noRegionNumber, 320, 12, regionBoxArray);

    // ---- write the data as a plot file
    int finestLevel(0), nLevels(1), numState(1), nGrow(0);
    Vector<Box> probDomain(nLevels);
    BoxArray dnpBoxArray(rFab.box());
    probDomain[0] = rFab.box();
    Vector<MultiFab> state(nLevels);
    const DistributionMapping dnpDM(dnpBoxArray);
    state[finestLevel].define(dnpBoxArray, dnpDM, numState, nGrow);
    MultiFab &regionsMF = state[finestLevel];
    if(bIOP) {
      regionsMF[0].copy<RunOn::Host>(rFab);
    }
    if(bIOP) { cout << "Writing plotfile:  " << plotfileName << endl; }

    std::string plotFileVersion("RegionsProf-V1.0");
    Real time(0.0);
    Vector<Real> probLo(BL_SPACEDIM, 0.0);
    Vector<Real> probHi(BL_SPACEDIM, 1.0);
    Vector<int>  refRatioPerLevel(state.size() - 1);
    int refRatioAll(2);
    for(int i(0); i < refRatioPerLevel.size(); ++i) {
      refRatioPerLevel[i] = refRatioAll;
    }
    Vector<Vector<Real> > dxLevel(state.size());
    for(int i(0); i < state.size(); ++i) {
      dxLevel[i].resize(BL_SPACEDIM);
      dxLevel[i][0] = 1.0;
      dxLevel[i][1] = 1.0;
    }
    int coordSys(0);
    Vector<std::string> inVarNames(numState);
    inVarNames[0] = "regions";
    bool verb(false);
    FABio::Format oldFabFormat(FArrayBox::getFormat());
    FArrayBox::setFormat(FABio::FAB_NATIVE_32);

    WritePlotfile(plotFileVersion, state, time,
                  probLo, probHi, refRatioPerLevel,
                  probDomain, dxLevel, coordSys,
                  plotfileName, inVarNames, verb);

    FArrayBox::setFormat(oldFabFormat);

    std::string rFileName(plotfileName + "/RegionNames.txt");
    std::ofstream rnames(rFileName.c_str());
    std::map<int, std::string>::const_iterator it;
    const std::map<int, std::string> &regionNumbersH = regOutputStats_H.RegionNumbers();
    for(it = regionNumbersH.begin(); it != regionNumbersH.end(); ++it) {
      rnames << it->second << ' ' << it->first << '\n';
    }
    rnames.close();
}


// ----------------------------------------------------------------------
void DataServices::RunACTPF(std::string &plotfileName,
                                int maxSmallImageLength, int refRatioAll,
	                        const Vector<string> &actFNames)
{
#if (BL_SPACEDIM != 2)
  cout << "**** Error:  DataServices::RunACTPF is only supported for 2D" << endl;
#else
    bool bIOP(ParallelDescriptor::IOProcessor());
    int  myProc(ParallelDescriptor::MyProc());
    int  nProcs(ParallelDescriptor::NProcs());
    int dataNProcs(BLProfStats::GetNProcs());

    if(BL_SPACEDIM != 2) {
      if(bIOP) {
        cout << "DataServices::RunACTPF:  bRunACTPF only supported for 2D."
	     << endl; 
      }
      return;
    }

    if( ! bTraceDataAvailable) {
      if(bIOP) {
        cout  << "DataServices::RunACTPF:  trace data is not available." << endl;
      }
      return;
    }

    RegionsProfStats::OpenAllStreams(fileName);
    bool bTestFab(false);
    if(bTestFab) {
      FArrayBox actFab;
      regOutputStats_H.AllCallTimesFAB(actFab, actFNames[0]);
    }
    BL_PROFILE_VAR("RunACTPF", RunACTPF)
    if(bIOP) cout << endl << "------------------------------------ Process ACT PF." << endl;
/*
    Box procBox(IntVect(0, 0), IntVect(0, dataNProcs - 1));
    IntVect procMaxGrid(1, (dataNProcs / nProcs) + ((dataNProcs % nProcs) > 0 ? 1 : 0));
    BoxArray procBoxArrayTemp(procBox);
    procBoxArrayTemp.maxSize(procMaxGrid);
    // ---- now ensure the boxarray is nprocs long
    Vector<Box> procBoxes;
    int needMoreBoxes(nProcs - procBoxArrayTemp.size());
    for(int ipb(0); ipb < procBoxArrayTemp.size(); ++ipb) {
      Box b(procBoxArrayTemp[ipb]);
      if(needMoreBoxes) {
        Box chopBox(b.chop(YDIR, (b.smallEnd(YDIR) + b.bigEnd(YDIR)) / 2));
        procBoxes.push_back(chopBox);
        --needMoreBoxes;
      }
      procBoxes.push_back(b);
    }
    BoxArray procBoxArray(procBoxes.dataPtr(), procBoxes.size());
    if(procBoxArray.size() != nProcs) {
      SHOWVAL(nProcs);
      SHOWVAL(dataNProcs);
      SHOWVAL(procBoxArray.size());
      if(bIOP) cout << "---- procBoxArray = " << procBoxArray << endl;
      amrex::Abort("bRunACTPF::Error 0");
    }
*/
    //const BoxArray &procBoxArray = pdServices.ProcBoxArray();
    Box myBox(procBoxArray[myProc]);

    if(bIOP) cout << "---- procBoxArray = " << procBoxArray << endl;

    Vector<Vector<Real> > whichFuncAllTimes(dataNProcs);  // [proc][functime]

    RegionsProfStats::SetInitDataBlocks(true);
    std::string whichFuncName(actFNames[0]);

    int whichFuncNameInt(-1);
    for(int i(0); i < regOutputStats_H.NumbersToFName().size(); ++i) {
      if(regOutputStats_H.NumbersToFName()[i] == whichFuncName) {
        whichFuncNameInt = i;
      }
    }
    if(bIOP) {
      cout << "**** whichFuncName whichFuncNameInt = " << whichFuncName
           << "  " <<  whichFuncNameInt << endl;
    }


    RegionsProfStats::OpenAllStreams(fileName);
    regOutputStats_H.FillAllCallTimes(whichFuncAllTimes, whichFuncName,
                                      whichFuncNameInt, myBox);

    RegionsProfStats::CloseAllStreams();

    int smallY(myBox.smallEnd(YDIR)), bigY(myBox.bigEnd(YDIR));
    int whichFuncNCalls(whichFuncAllTimes[smallY].size());
    bool bSameNCalls(true);
    for(int p(smallY); p <= bigY; ++p) {
      if(whichFuncAllTimes[p].size() != whichFuncNCalls) {
        cout << "==== bSameNCalls = false" << endl;
        bSameNCalls = false;
      }
    }
    int ncMin(whichFuncNCalls), ncMax(whichFuncNCalls);
    ParallelDescriptor::ReduceIntMin(ncMin);
    ParallelDescriptor::ReduceIntMax(ncMax);
    if( ! bSameNCalls || ncMin != ncMax) {
      if(bIOP) {
        cout << "**** bSameNCalls == false for:  " << whichFuncName
	     << " ::: unsupported." << endl;
        SHOWVAL(ncMin);
        SHOWVAL(ncMax);
        SHOWVAL(whichFuncNCalls);
      }
    } else {

    // ---- fill the local fab with function call times
    Box myWFNBox(myBox);
    myWFNBox.setBig(XDIR, whichFuncNCalls - 1);
    FArrayBox myFab(myWFNBox, 1);
    myFab.setVal(0.0);
    Real *dptr = myFab.dataPtr(0);
    int nX(myWFNBox.length(XDIR)), nY(myWFNBox.length(YDIR));
    if(bIOP) { SHOWVAL(nX); }
    if(bIOP) { SHOWVAL(nY); }
    for(int p(0); p < nY; ++p) {
      for(int cnum(0); cnum < nX; ++cnum) {
        int index((p * nX) + cnum);
        dptr[index] = whichFuncAllTimes[p + smallY][cnum];
      }
    }

    BoxArray wfnBoxArray(procBoxArray);
    for(int pba(0); pba < wfnBoxArray.size(); ++pba) {
      Box b(wfnBoxArray[pba]);
      b.setBig(XDIR, whichFuncNCalls - 1);
      wfnBoxArray.set(pba, b);
    }
    if(bIOP) { SHOWVAL(wfnBoxArray); }

    Vector<int> myMap(nProcs);
    for(int i(0); i < myMap.size() - 1; ++i) {
      myMap[i] = i;
    }

    DistributionMapping myDMap(myMap);
    MultiFab mfWFN(wfnBoxArray, myDMap, 1, 0);
    for(MFIter mfi(mfWFN); mfi.isValid(); ++mfi) {
      mfWFN[mfi.index()].copy<RunOn::Host>(myFab);
    }


    int nLevels(1), finestLevel(0), numState(1), nGrow(0);

    Box dnpBox(IntVect(0, 0), IntVect(whichFuncNCalls - 1, dataNProcs - 1));
    Box dnpBoxBlocked(dnpBox);
    if(bIOP) cout << ")))) dnpBoxBlocked = " << dnpBoxBlocked << endl;
    while(dnpBoxBlocked.length(YDIR) > maxSmallImageLength) {
      dnpBoxBlocked.coarsen(refRatioAll);
      if(bIOP) cout << ")))) coarsened dnpBoxBlocked = " << dnpBoxBlocked << endl;
      ++nLevels;
    }
    finestLevel = nLevels - 1;
    Vector<Box> probDomain(nLevels);
    for(int i(0); i < nLevels; ++i) {
      probDomain[i] = dnpBoxBlocked;
      if(bIOP) { cout << ")))) probDomain[" << i << "] =  " << probDomain[i] << endl; }
      dnpBoxBlocked.refine(refRatioAll);
    }

    BoxArray dnpBoxArray(probDomain[finestLevel]);
    IntVect maxGrid((dataNProcs / nProcs) + refRatioAll, dataNProcs);
    dnpBoxArray.maxSize(maxGrid);

    Vector<MultiFab> state(nLevels);
    const DistributionMapping dnpDM(dnpBoxArray);
    state[finestLevel].define(dnpBoxArray, dnpDM, numState, nGrow);
    MultiFab &fLMF = state[finestLevel];
    fLMF.setVal(0.0);
    fLMF.copy(mfWFN);

    // ---- make an xgraph of coefficient of variation for each call
    Vector<Real> coeffVar(whichFuncNCalls, 0.0);
    if(bIOP) cout << ")))) fLMF.boxArray = " << fLMF.boxArray() << endl;
    for(MFIter mfi(fLMF); mfi.isValid(); ++mfi) {
      const FArrayBox &fab = fLMF[mfi];
      const Box &b = fab.box();
      int smallX(b.smallEnd(XDIR)), bigX(b.bigEnd(XDIR));
      for(int c(smallX); c <= bigX; ++c) {
        Box cBox(b);
        cBox.setSmall(XDIR, c);
        cBox.setBig(XDIR, c);
        FArrayBox cFab(cBox, 1);
        int np(cFab.box().numPts());
        cFab.copy<RunOn::Host>(fab);
        Real cAvg(cFab.sum(0) / static_cast<Real>(np));
        Real variance(0.0);
        Real *ptr = cFab.dataPtr(0);
        for(int ip(0); ip < np; ++ip) {
          Real r(ptr[ip]);
          variance += (r - cAvg) * (r - cAvg);
        }
        variance /= static_cast<Real>(np);  // ---- np - 1 for sample
        if(variance > 0.0 && cAvg > 0.0) {
          coeffVar[c] = 100.0 * (std::sqrt(variance) / cAvg);  // ---- percent
        }
      }
    }
    ParallelDescriptor::ReduceRealSum(coeffVar.dataPtr(), coeffVar.size(),
                                      ParallelDescriptor::IOProcessorNumber());
    if(bIOP) {
      std::string cvfileNameUS("CV_");
      cvfileNameUS += whichFuncName;
      cvfileNameUS += ".xgr";
      std::string cvfileName(SanitizeName(cvfileNameUS));
      std::ofstream cvarout(cvfileName.c_str());
      for(int i(0); i < coeffVar.size(); ++i) {
        cvarout << i << " " << coeffVar[i] << endl;
      }
      cvarout.close();
    }


    // ---- now the real data is in the multifab, avg down to a reasonable size
    Vector<Vector<int> > adRefRatio(state.size() - 1);
    for(int i(0); i < finestLevel; ++i) {
      adRefRatio[i].resize(BL_SPACEDIM, refRatioAll);
    }

    for(int cLev(finestLevel - 1); cLev >= 0; --cLev) {
      if(bIOP) { cout << "Averaging down level " << cLev << endl; }
      BoxArray ba(BoxArray(state[cLev + 1].boxArray()).coarsen(adRefRatio[cLev][0]));
      // ---- call uniqify, otherwise ba is just a reference to the
      // ---- original boxarray with a coarsening factor
      ba.uniqify();
      if( ! ba.isDisjoint()) {
        if(bIOP) { cout << "BA:  Coarsened BoxArray not disjoint:  " << ba << endl; }
        SimpleRemoveOverlap(ba);
        if(bIOP) { cout << "BA:  Coarsened BoxArray after removeOverlap:  " << ba << endl; }
        if( ! ba.isDisjoint()) {
          if(bIOP) { cout << "BA:  Coarsened BoxArray still not disjoint:  " << ba << endl; }
        }
      }
      if(bIOP) { cout << "Coarsened BoxArray size = " << ba.size() << endl; }
      const DistributionMapping baDM(ba);
      state[cLev].define(ba, baDM, numState, nGrow);
      state[cLev].setVal(0.0);

      avgDown(state[cLev], state[cLev+1], 0, 0, numState, adRefRatio[cLev]);
    }
    ParallelDescriptor::Barrier();

    // ---- copy data into a more ddio favorable configuration
    int sqMG(128);
    Vector<MultiFab> sqState(nLevels);
    for(int i(0); i < state.size(); ++i) {
      BoxArray sqBA(state[i].boxArray().minimalBox());
      const DistributionMapping sqDM(sqBA);
      sqBA.maxSize(sqMG);
      sqState[i].define(sqBA, sqDM, numState, nGrow);
      sqState[i].setVal(0.0);
      sqState[i].copy(state[i]);
    }


    // ---- write the data as a plot file
    std::string plotfileNameOut;
    if(plotfileName.empty()) {
      std::string plotfileNameUS("plt");
      plotfileNameUS += whichFuncName;
      plotfileNameOut = SanitizeName(plotfileNameUS);
    } else {
      plotfileNameOut = plotfileName;
    }
    if(bIOP) { cout << "Writing plotfile:  " << plotfileNameOut << endl; }

    std::string plotFileVersion("NavierStokes-V1.1");
    Real time(0.0);
    Vector<Real> probLo(BL_SPACEDIM, 0.0);
    Vector<Real> probHi(BL_SPACEDIM, 1.0);
    Vector<int>  refRatioPerLevel(sqState.size() - 1);
    for(int i(0); i < refRatioPerLevel.size(); ++i) {
      refRatioPerLevel[i] = refRatioAll;
    }
    Vector<Vector<Real> > dxLevel(sqState.size());
    for(int i(0); i < sqState.size(); ++i) {
      dxLevel[i].resize(BL_SPACEDIM);
      dxLevel[i][0] = 1.0;
      dxLevel[i][1] = 1.0;
    }
    int coordSys(0);
    Vector<std::string> inVarNames(numState);
    inVarNames[0] = whichFuncName;
    bool verb(false);
    FABio::Format oldFabFormat(FArrayBox::getFormat());
    FArrayBox::setFormat(FABio::FAB_NATIVE_32);

    WritePlotfile(plotFileVersion, sqState, time,
                  probLo, probHi, refRatioPerLevel,
                  probDomain, dxLevel, coordSys,
                  plotfileNameOut, inVarNames, verb);

    FArrayBox::setFormat(oldFabFormat);

    }

    BL_PROFILE_VAR_STOP(RunACTPF)
#endif
}


// ----------------------------------------------------------------------
void DataServices::RunSyncPointData()
{
   bool bIOP(ParallelDescriptor::IOProcessor());
    if( ! bCommDataAvailable) {
      if(bIOP) {
        cout << "DataServices::RunSyncPointData:  comm data is not available." << std::endl;
      }
      return;
    }

    int  myProc(ParallelDescriptor::MyProc());
    int  nProcs(ParallelDescriptor::NProcs());
    int dataNProcs(BLProfStats::GetNProcs());

    bool bDoReductions(true);
    Vector<Vector<Real> > barrierExitTimes(dataNProcs);  // [proc, bnum]
    Vector<Vector<Real> > barrierWaitTimes(dataNProcs);
    Vector<Vector<Real> > barrierSkewTimes(dataNProcs);
    Vector<Vector<Real> > reductionWaitTimes(dataNProcs);
    Long nBMax(0), nRMax(0);

    const Vector<string> &commHeaderFileNames = CommProfStats::GetHeaderFileNames();
    CommProfStats::SetInitDataBlocks(true);
    CommProfStats::InitDataFileNames(commHeaderFileNames);
    CommProfStats::OpenAllStreams(fileName);

    Vector<CommProfStats> commOutputStats(commHeaderFileNames.size());
    for(int hfnI(0); hfnI < commHeaderFileNames.size(); ++hfnI) {
      if(bRegionDataAvailable) {
        commOutputStats[hfnI].SetRegionTimeRanges(commOutputStats_H.GetRegionTimeRanges());
        commOutputStats[hfnI].SetFilterTimeRanges(commOutputStats_H.GetFilterTimeRanges());
      }

      std::string commDataHeaderFileName(fileName + '/' + commHeaderFileNames[hfnI]);
      if( ! ( yyin = fopen(commDataHeaderFileName.c_str(), "r"))) {
        if(bIOP) {
          cerr << "DataServices::RunSyncPointData:  Cannot open file:  " << commDataHeaderFileName
               << " ... continuing." << endl;
        }
        continue;
      }

      yyparse(&commOutputStats[hfnI]);
      fclose(yyin);
      commOutputStats[hfnI].ReportSyncPointDataSetup(nBMax, nRMax);
    }

    ParallelDescriptor::ReduceLongMax(nBMax);
    ParallelDescriptor::ReduceLongMax(nRMax);

    if(nRMax > 2048) {
      bDoReductions = false;
    }

    for(int i(0); i < dataNProcs; ++i) {
      barrierExitTimes[i].resize(nBMax + 1, 0.0);
      barrierWaitTimes[i].resize(nBMax + 1, 0.0);
      //barrierSkewTimes[i].resize(nBMax + 1, 0.0);
      if(bDoReductions) {
        reductionWaitTimes[i].resize(nRMax + 1, 0.0);
      } else {
        reductionWaitTimes[i].resize(1, 0.0);
      }
    }

    for(int hfnI(0); hfnI < commHeaderFileNames.size(); ++hfnI) {
      if(myProc == hfnI % nProcs) {
        commOutputStats[hfnI].ReportSyncPointData(barrierExitTimes,
                                                  barrierWaitTimes,
                                                  reductionWaitTimes,
                                                  bDoReductions);
      }
    }

    CommProfStats::CloseAllStreams();

    Long nBarriers(nBMax + 1);
    Long nReductions(nRMax + 1);
    Vector<Real> bExitAll(dataNProcs * nBarriers);
    Vector<Real> bWaitAll(dataNProcs * nBarriers);
    Vector<Real> rwAll;
    if(bDoReductions) {
      rwAll.resize(dataNProcs * nReductions);
    }

    // pack
    int bCount(0), rCount(0);
    for(int p(0); p < dataNProcs; ++p) {
      for(int b(0); b < nBarriers; ++b) {
        bExitAll[bCount] = barrierExitTimes[p][b];
        bWaitAll[bCount] = barrierWaitTimes[p][b];
        ++bCount;
      }
      if(bDoReductions) {
        for(int r(0); r < nReductions; ++r) {
          rwAll[rCount] = reductionWaitTimes[p][r];
          ++rCount;
        }
      }
    }

    ParallelDescriptor::ReduceRealSum(bExitAll.dataPtr(), bExitAll.size());
    ParallelDescriptor::ReduceRealSum(bWaitAll.dataPtr(), bWaitAll.size());
    if(bDoReductions) {
      ParallelDescriptor::ReduceRealSum(rwAll.dataPtr(), rwAll.size());
    }

    // unpack
    bCount = 0;
    rCount = 0;
    for(int p(0); p < dataNProcs; ++p) {
      for(int b(0); b < nBarriers; ++b) {
        barrierExitTimes[p][b] = bExitAll[bCount];
        barrierWaitTimes[p][b] = bWaitAll[bCount];
        ++bCount;
      }
      if(bDoReductions) {
        for(int r(0); r < nReductions; ++r) {
          reductionWaitTimes[p][r] = rwAll[rCount];
          ++rCount;
        }
      }
    }

    // ------------------------------------------------ print barrier wait times
    if(bIOP) {
      ::WriteFab("bwaits", nBarriers, dataNProcs, &bWaitAll[0]);
      ::WriteFab("bexits", nBarriers, dataNProcs, &bExitAll[0]);
      if(bDoReductions) {
       ::WriteFab("rwaits", nReductions, dataNProcs, &rwAll[0]);
      }
    }

}


// ----------------------------------------------------------------------
void DataServices::RunSendRecv()
{
    bool bIOP(ParallelDescriptor::IOProcessor());
    if( ! bCommDataAvailable) {
      if(bIOP) {
        cout << "DataServices::RunSendRecv:  comm data is not available." << std::endl;
      }
      return;
    }

    /*
    if(myProc < commHeaderFileNames.size()) {
      for(int hfnI(0); hfnI < commHeaderFileNames.size(); ++hfnI) {
        if(myProc == hfnI % nProcs) {
          CommProfStats commOutputStats;
          std::string commFileName_H_nnnn(dirName + '/' + commHeaderFileNames[hfnI]);
          if( ! ( yyin = fopen(commFileName_H_nnnn.c_str(), "r"))) {
            if(bIOP) {
              cerr << "DataServices::RunSendRecv:  Cannot open file:  " << commFileName_H_nnnn
                   << " ... continuing." << endl;
            }
            continue;
          }
          if(bRegionDataAvailable) {
            commOutputStats.SetRegionTimeRanges(commOutputStats_H.GetRegionTimeRanges());
            commOutputStats.SetFilterTimeRanges(commOutputStats_H.GetFilterTimeRanges());
          }
          yyparse(&commOutputStats);
          fclose(yyin);
          commOutputStats.SendRecvData("SendRecvData", 0.0, -1.0);
        }
      }
    }
    */
}


// ----------------------------------------------------------------------
void DataServices::RunSendRecvList()
{
    // ---- this is not parallelized yet

    bool bIOP(ParallelDescriptor::IOProcessor());
    if( ! bCommDataAvailable) {
      if(bIOP) {
        cout << "DataServices::RunSendRecvList:  comm data is not available." << std::endl;
      }
      return;
    }

//    int  myProc(ParallelDescriptor::MyProc());
//    int  nProcs(ParallelDescriptor::NProcs());

    if(ParallelDescriptor::IOProcessor()) {
    const Vector<string> &commHeaderFileNames = CommProfStats::GetHeaderFileNames();
    std::multimap<Real, CommProfStats::SendRecvPairUnpaired> srMMap;  // [call time, sr]
//    if(myProc < commHeaderFileNames.size()) {
      for(int hfnI(0); hfnI < commHeaderFileNames.size(); ++hfnI) {
//        if(myProc == hfnI % nProcs) {
          CommProfStats commOutputStats;
          std::string commFileName_H_nnnn(fileName + '/' + commHeaderFileNames[hfnI]);
          if( ! ( yyin = fopen(commFileName_H_nnnn.c_str(), "r"))) {
            if(bIOP) {
              cerr << "DataServices::RunSendRecvList:  Cannot open file:  " << commFileName_H_nnnn
                   << " ... continuing." << endl;
            }
            continue;
          }
          if(bRegionDataAvailable) {
            commOutputStats.SetRegionTimeRanges(commOutputStats_H.GetRegionTimeRanges());
            commOutputStats.SetFilterTimeRanges(commOutputStats_H.GetFilterTimeRanges());
          }
          yyparse(&commOutputStats);
          fclose(yyin);
          commOutputStats.SendRecvList(srMMap);
//        }
      }
//    }

      std::ofstream srlOut("SendRecvList.txt");
      srlOut << std::left << std::setw(16) << "#### Time" << '\t' << "Type    " << '\t'
             << "From" << '\t' << "To" << '\t' << "Size" << '\t'
             << "Tag" << '\t' << "Regions"
             << std::endl;
      srlOut << std::setprecision(16) << std::fixed;
      const std::map<int, std::string> &regionNumbersH = regOutputStats_H.RegionNumbers();
      std::map<int, std::string>::const_iterator itn;
      std::list<CommProfStats::SendRecvPairUnpaired> unpairedSR;

      std::multimap<Real, CommProfStats::SendRecvPairUnpaired>::iterator it;
      blProfStats_H.SetRegionTimeRanges(regOutputStats_H.GetRegionTimeRanges());
      blProfStats_H.SetFilterTimeRanges(regOutputStats_H.GetFilterTimeRanges());
      for(it = srMMap.begin(); it != srMMap.end(); ++it) {
        CommProfStats::SendRecvPairUnpaired &srp = it->second;
        if(srp.fromProc < 0) {
          unpairedSR.push_back(srp);
        }
        int timedProc;
        if ( commOutputStats_H.IsSend(srp.unmatchedCFType) )
        { timedProc = srp.fromProc; }
        else
        { timedProc = srp.toProc; }
        std::set<int> whichRegions(blProfStats_H.WhichRegions(timedProc, it->first));
        if(srp.fromProc >= 0) {
          srlOut << it->first << '\t'
                 << BLProfiler::CommStats::CFTToString(srp.unmatchedCFType) << '\t'
                 << srp.fromProc << '\t' << srp.toProc << '\t' << srp.dataSize << '\t'
                 << srp.tag << '\t';
                 for(std::set<int>::iterator itwr = whichRegions.begin();
                     itwr != whichRegions.end(); ++itwr)
                 {
                   itn = regionNumbersH.find(*itwr);
                   srlOut << itn->second << '\t';
                 }
                 srlOut << '\n';
        }
      }
      srlOut.close();
      if( ! unpairedSR.empty()) {
        std::ofstream upSROut("UnpairedSendRecv.txt");
        upSROut << std::setprecision(8) << std::fixed;
        upSROut << std::left << std::setw(10) << "Type    " << '\t'
                << "From" << '\t' << "To" << '\t' << "Size" << '\t'
                << "Tag" << '\t' << "sendTime" << '\t' << "recvTime" << '\t'
                << "totalTime" << '\n';
        for(auto uit = unpairedSR.begin(); uit != unpairedSR.end(); ++uit) {
          CommProfStats::SendRecvPairUnpaired &srp = *uit;
          upSROut << BLProfiler::CommStats::CFTToString(srp.unmatchedCFType) << '\t'
                  << srp.fromProc << '\t' << srp.toProc << '\t' << srp.dataSize << '\t'
                  << srp.tag << '\t'
                  << srp.sendTime << '\t' << srp.recvTime << '\t' << srp.totalTime << '\n';
        }
        upSROut.close();
      }
    }
//    ParallelDescriptor::Barrier();
}


// ----------------------------------------------------------------------
void DataServices::InitProxMap()
{
  CommProfStats::InitProxMap();
}


// ----------------------------------------------------------------------
void DataServices::TCEdison()
{
  bool bIOP(ParallelDescriptor::IOProcessor());
  if(bIOP) {
    commOutputStats_H.InitEdisonTopoMF();
    std::string topoFileName("edisontopo.out");
    if( ! ( yyin = fopen(topoFileName.c_str(), "r"))) {
      cout << "DataServices::TCEdison:  Cannot open file:  " << topoFileName << endl;
    } else {
      yyparse(&commOutputStats_H);
      fclose(yyin);
      commOutputStats_H.WriteEdisonTopoMF();
    }
  }
}


// ----------------------------------------------------------------------
void DataServices::WriteHTML(const std::string &hFileName,
                                 bool simpleCombine, int whichProc)
{
  bool bIOP(ParallelDescriptor::IOProcessor());
  if( ! bTraceDataAvailable) {
    if(bIOP) {
      cout << "DataServices::WriteHTML:  trace data is not available." << std::endl;
    }
    return;
  }

  if(bIOP) {
    RegionsProfStats::OpenAllStreams(hFileName);

    std::ofstream outStream(hFileName.c_str(), std::ios::out | std::ios::trunc);
    if( ! outStream.good()) {
      cerr << "**** Error in DataServices::WriteHTML:  could not open "
           << hFileName << endl;
    } else {
      regOutputStats_H.WriteHTML(outStream, simpleCombine, whichProc);
      outStream.close();
    }

    RegionsProfStats::CloseAllStreams();
  }
}


// ----------------------------------------------------------------------
void DataServices::WriteHTMLNC(const std::string &hFileName, int whichProc)
{
  bool bIOP(ParallelDescriptor::IOProcessor());
  if( ! bTraceDataAvailable) {
    if(bIOP) {
      cout << "DataServices::WriteHTMLNC:  trace data is not available." << std::endl;
    }
    return;
  }

  if(bIOP) {
    RegionsProfStats::OpenAllStreams(hFileName);

    std::ofstream outStream(hFileName.c_str(), std::ios::out | std::ios::trunc);
    if( ! outStream.good()) {
      cerr << "**** Error in DataServices::WriteHTML:  could not open "
           << hFileName << endl;
    } else {
      regOutputStats_H.WriteHTMLNC(outStream, whichProc);
      outStream.close();
    }

    RegionsProfStats::CloseAllStreams();
  }
}


// ----------------------------------------------------------------------
void DataServices::WriteTextTrace(const std::string &tFileName, bool simpleCombine,
                                  int whichProc, std::string delimString)
{
  bool bIOP(ParallelDescriptor::IOProcessor());
  if( ! bTraceDataAvailable) {
    if(bIOP) {
      cout << "DataServices::WriteTextTrace:  trace data is not available."
           << std::endl;
    }
    return;
  }

  if(bIOP) {
    RegionsProfStats::OpenAllStreams(tFileName);

    std::ofstream outStream(tFileName.c_str(), std::ios::out | std::ios::trunc);
    if( ! outStream.good()) {
      cerr << "**** Error in DataServices::WriteHTML:  could not open "
           << tFileName << endl;
    } else {
      regOutputStats_H.WriteTextTrace(outStream, simpleCombine, whichProc, delimString);
      outStream.close();
    }

    RegionsProfStats::CloseAllStreams();
  }
}
#endif

// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
}  // namespace amrex
