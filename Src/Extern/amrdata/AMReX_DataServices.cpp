
// ---------------------------------------------------------------
// DataServices.cpp
// ---------------------------------------------------------------
#include <winstd.H>

#include <AmrvisConstants.H>
#include <DataServices.H>
#include <ParallelDescriptor.H>

#ifndef BL_NOLINEVALUES
# include <XYPlotDataList.H>
#endif

#include <iostream>
#include <fstream>
#include <cstdio>
using std::ios;
using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;

Array<DataServices *> DataServices::dsArray;
int DataServices::dsArrayIndexCounter = 0;
int DataServices::dsFabOutSize = 0;
bool DataServices::dsBatchMode = false;


// ---------------------------------------------------------------
namespace ParallelDescriptor {
  template <> void Bcast (Box *b, size_t n, int root) {
    const int n3SDim(n * 3 * BL_SPACEDIM);

    Array<int> tmp(n3SDim);

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

  if(bAmrDataOk) {
    dsArrayIndex = DataServices::dsArrayIndexCounter;
    ++DataServices::dsArrayIndexCounter;
    DataServices::dsArray.resize(DataServices::dsArrayIndexCounter);
    DataServices::dsArray[dsArrayIndex] = this;
  }
}


// ---------------------------------------------------------------
DataServices::DataServices() {
  // must call init
  bAmrDataOk = false;
  iWriteToLevel = -1;
}


// ---------------------------------------------------------------
void DataServices::Init(const string &filename, const Amrvis::FileType &filetype) {
  fileName = filename;
  fileType = filetype;
  bAmrDataOk = false;
  iWriteToLevel = -1;

  numberOfUsers = 0;  // the user must do all incrementing and decrementing
}


// ---------------------------------------------------------------
DataServices::~DataServices() {
  BL_ASSERT(numberOfUsers == 0);
  DataServices::dsArray[dsArrayIndex] = NULL;
}


// ---------------------------------------------------------------
void DataServices::SetBatchMode() {
  dsBatchMode = true;
}


// ---------------------------------------------------------------
void DataServices::SetFabOutSize(int iSize) {
  if(iSize == 1 || iSize == 8 || iSize == 32) {
    dsFabOutSize = iSize;
  } else {
    cerr << "Warning:  DataServices::SetFabOutSize:  size must be 1, 8 or 32 only."
	 << "  Defaulting to native." << endl;
    dsFabOutSize = 0;
  }
}


// ---------------------------------------------------------------
void DataServices::Dispatch(DSRequestType requestType, DataServices *ds, ...) {
  bool bContinueLooping(true);
  va_list ap;
  int whichDSIndex;
  int ioProcNumber(ParallelDescriptor::IOProcessorNumber());

 while(bContinueLooping) {
  if(ParallelDescriptor::IOProcessor() || dsBatchMode) {
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
    for(int i(0); i < dsArray.size(); ++i) {
      if(DataServices::dsArray[i] != NULL) {
	BL_ASSERT(DataServices::dsArray[i]->numberOfUsers == 0);
	delete DataServices::dsArray[i];
      }
    }
    ParallelDescriptor::EndParallel();
    exit(0);
  }  // end ExitRequest

  if(ParallelDescriptor::IOProcessor()) {
    va_start(ap, ds);
    whichDSIndex = ds->dsArrayIndex;
  }

  ParallelDescriptor::Bcast(&whichDSIndex, 1, 0);
  if( ! ParallelDescriptor::IOProcessor()) {
    ds = DataServices::dsArray[whichDSIndex];
  }
  BL_ASSERT(ds != NULL);

  switch(requestType) {
    case DeleteRequest:
    {
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
      BoxLib::Abort("FillVarArrayOfFabs not implemented yet.");
    }
    break;

    case FillVarMultiFab:
    {
      BoxLib::Abort("FillVarMultiFab not implemented yet.");
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

#ifndef BL_NOLINEVALUES
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
#endif

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
    BoxLib::Abort("DataServices::DumpSlice(1): slicechar buffer too small");
  sliceFile += slicechar;
  sliceFile += ".fab";
  cout << "sliceFile = " << sliceFile << endl;

  Box sliceBox(amrData.ProbDomain()[iWTL]);

  if(BL_SPACEDIM == 2 && slicedir == Amrvis::ZDIR) {
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
    BoxLib::Abort("DataServices::DumpSlice(2): slicechar buffer too small");
  sliceFile += slicechar;
  sliceFile += ".fab";
  cout << "sliceFile = " << sliceFile << endl;

  Box sliceBox(amrData.ProbDomain()[iWTL]);

  if(BL_SPACEDIM == 2 && slicedir == Amrvis::ZDIR) {
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
# if (BL_SPACEDIM == 2)
  int count = snprintf(slicechar, N, "%d_%d__%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR), iWTL);
# else
  int count = snprintf(slicechar, N, "%d_%d_%d__%d_%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR), b.smallEnd(Amrvis::ZDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR),   b.bigEnd(Amrvis::ZDIR), iWTL);
#endif
  if (count >= N)
    BoxLib::Abort("DataServices::DumpSlice(3): slicechar buffer too small");      
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
# if (BL_SPACEDIM == 2)
  int count = snprintf(slicechar, N, "%d_%d__%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR), iWTL);
# else
  int count = snprintf(slicechar, N, "%d_%d_%d__%d_%d_%d.Level_%d",
                       b.smallEnd(Amrvis::XDIR), b.smallEnd(Amrvis::YDIR), b.smallEnd(Amrvis::ZDIR),
                       b.bigEnd(Amrvis::XDIR),   b.bigEnd(Amrvis::YDIR),   b.bigEnd(Amrvis::ZDIR), iWTL);
#endif
  if (count >= N)
    BoxLib::Abort("DataServices::DumpSlice(4): slicechar buffer too small");      
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
// Change this to take an Array<Box> (or BoxArray?) and create
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

  Array<FArrayBox *> destFabs(1);
  Array<Box> destBoxes(1);
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
// Change this to take an Array<Box> (or BoxArray?) and create
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
    Array<FArrayBox *> destFabs(1);
    Array<Box> destBoxes(1);
    destFabs[0]  = &tempdata;
    destBoxes[0] = region;
    amrData.FillVar(destFabs, destBoxes, lev, amrData.PlotVarNames()[ivar],
		    ParallelDescriptor::IOProcessorNumber());
    int srccomp(0);
    int destcomp(ivar);
    int ncomp(1);
    if(ParallelDescriptor::IOProcessor()) {
      data.copy(tempdata, srccomp, destcomp, ncomp);
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
bool DataServices::CanDerive(const Array<string> &names) const {
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
void DataServices::PointValue(int pointBoxArraySize, Box *pointBoxArray,
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


#ifndef BL_NOLINEVALUES
// ---------------------------------------------------------------
void DataServices::LineValues(int lineBoxArraySize, Box *lineBoxArray, int whichDir,
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
#endif


// ---------------------------------------------------------------
bool DataServices::MinMax(const Box &onBox, const string &derived, int level,
		          Real &dataMin, Real &dataMax, bool &minMaxValid)
{
  minMaxValid =  amrData.MinMax(onBox, derived, level, dataMin, dataMax);
  return minMaxValid;
}  // end MinMax

// ---------------------------------------------------------------
// ---------------------------------------------------------------
