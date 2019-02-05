// ---------------------------------------------------------------
// DataServices.H
// ---------------------------------------------------------------
#ifndef _DATASERVICES_H_
#define _DATASERVICES_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_AmrData.H>

#include <cstdarg>
#include <string>

#ifdef BL_USE_PROFPARSER
#include <map>
#include <AMReX_BLProfUtilities.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_BLProfStats.H>
#include <BLProfParser.tab.H>
#include <AMReX_CommProfStats.H>
#include <AMReX_RegionsProfStats.H>
#endif


using std::string;

namespace amrex {

class XYPlotDataList;

class DataServices {

  public:

    enum DSRequestType {
      InvalidRequestType = 0,
      ExitRequest,
      NewRequest,
      DeleteRequest,
      FillVarOneFab,
      FillVarArrayOfFabs,
      FillVarMultiFab,
      WriteFabOneVar,
      WriteFabAllVars,
      DumpSlicePlaneOneVar,
      DumpSlicePlaneAllVars,
      DumpSliceBoxOneVar,
      DumpSliceBoxAllVars,
      MinMaxRequest,
      PointValueRequest,
      LineValuesRequest,
#ifdef BL_USE_PROFPARSER
      InitTimeRanges,
      WriteSummaryRequest,
      CheckProfDataRequest,
      RunStatsRequest,
      RunSendsPFRequest,
      RunTimelinePFRequest,
      MakeRegionPltRequest,
      MakeFilterFileRequest,
#endif
    };

    DataServices(const string &filename, const Amrvis::FileType &filetype);
    DataServices();
    ~DataServices();

    void Init(const string &filename, const Amrvis::FileType &filetype);
    void InitRegionTimeRanges();

    static void SetBatchMode();
    static void Dispatch(DSRequestType requestType, DataServices *ds, ...);

    static void SetFabOutSize(int iSize);

  private:

    // Dispatchable functions

      //! FillVarOneFab
      bool FillVar(FArrayBox *destFab, const Box &destBox,
                   int finestFillLevel, const string &varname,
                   int procWithFab);

      // FillVarArrayOfFabs
      //bool FillVar(Vector<FArrayBox *> &destFabs, const Vector<Box> &destBoxes,
                   //int finestFillLevel, const string &varname,
                   //int procWithFab);

      //! FillVarMultiFab
      bool FillVar(MultiFab &destMultiFab, int finestFillLevel,
                   const string &varname);

      //! WriteFabOneVar:  write a single variable to one fab
      bool WriteFab(const string &fabfilename, const Box &region, int lev,
                    const string &varname);

      //! WriteFabAllVars:  write all vars in one fab
      bool WriteFab(const string &fabfilename, const Box &region, int lev);

      //! DumpSlicePlaneOneVar
      bool DumpSlice(int slicedir, int slicenum, const string &varname);

      //! DumpSlicePlaneAllVars
      bool DumpSlice(int slicedir, int slicenum);

      //! DumpSliceBoxOneVar
      bool DumpSlice(const Box &slicebox, const string &varname);

      //! DumpSliceBoxAllVars
      bool DumpSlice(const Box &slicebox);

      //! MinMaxRequest
      bool MinMax(const Box &onBox, const string &derived, int level,
                  Real &dataMin, Real &dataMax, bool &minMaxValid);

      //! PointValueRequest
      void PointValue(int pointBoxArraySize, Box *pointBoxArray,
		      const string &currentDerived,
		      int coarsestLevelToSearch, int finestLevelToSearch,
		      int &intersectedLevel, Box &intersectedGrid,
		      Real &dataPointValue, bool &bPointIsValid);

      //! LineValuesRequest
      void LineValues(int lineBoxArraySize, Box *lineBoxArray, int whichDir,
                      const string &currentDerived,
                      int coarsestLevelToSearch, int finestLevelToSearch,
                      XYPlotDataList *dataList, bool &lineIsValid);

  public:

    // information functions
    const string &GetFileName() const { return fileName; }
    Amrvis::FileType GetFileType()      const { return fileType; }
    bool CanDerive(const string &name) const;
    bool CanDerive(const Vector<string> &names) const;
    int  NumDeriveFunc() const;
    void ListDeriveFunc(std::ostream &os) const;
    const Vector<string> &PlotVarNames() const { return amrData.PlotVarNames(); }
    AmrData &AmrDataRef() { return amrData; }
    bool AmrDataOk() { return bAmrDataOk; }
    void SetWriteToLevel(const int iwtlev)  { iWriteToLevel = iwtlev; }

    // other functions
    int GetNumberOfUsers() const  { return numberOfUsers; }
    void IncrementNumberOfUsers() { ++numberOfUsers; }
    void DecrementNumberOfUsers() { --numberOfUsers; }

    // profiling functions
#ifdef BL_USE_PROFPARSER
    bool ProfDataAvailable()   const { return bProfDataAvailable;   }
    bool RegionDataAvailable() const { return bRegionDataAvailable; }
    bool TraceDataAvailable()  const { return bTraceDataAvailable;  }
    bool CommDataAvailable()   const { return bCommDataAvailable;   }

    BLProfStats &GetBLProfStats() { return blProfStats_H; }
    RegionsProfStats &GetRegionsProfStats() { return regOutputStats_H; }
    CommProfStats &GetCommOutputStats() { return commOutputStats_H; }
    const amrex::BoxArray &ProcBoxArray() const { return procBoxArray; }

    void CheckProfData();

    void RunStats(std::map<int, string> &mpiFuncNames,
                  bool &statsCollected);

    void RunSendsPF(std::string &plotfileName,
                    int  maxSmallImageLength,
                    bool proxMap,
                    int refRatioAll);

    void RunTimelinePF(std::map<int, string> &mpiFuncNames,
                       std::string &plotfileName,
                       BLProfStats::TimeRange& subTimeRange,
                       int maxSmallImageLength,
                       int refRatioAll,
                       int nTimeSlots,
                       bool &statsCollected);

    void RunACTPF(std::string &plotfileName,
                  int maxSmallImageLength,
                  int refRatioAll,
                  const amrex::Vector<string> &actFNames);

    void RunSyncPointData();
    void RunSendRecvList();
    void RunSendRecv();

    BLProfStats::TimeRange FindCalcTimeRange();
    void MakeRegionPlt(std::string &plotfileName);
    void MakeFilterFile(std::string &fFileName);

    void ProcessGridLog(const std::string &gridlogFileName);
    void ParseFilterFile();

    void PrintCommStats(std::ostream &os, bool printHeaderNames);
    void WriteSummary(std::ostream &os, bool bWriteAverage, int whichProc,
                      bool bUseTrace, bool graphTopPct);
    void InitProxMap();
    void TCEdison();
    void WriteHTML(const std::string &fFileName, bool simpleCombine, int whichProc);
    void WriteHTMLNC(const std::string &fFileName, int whichProc);
    void WriteTextTrace(const std::string &fFileName, bool simpleCombine = true,
                        int whichProc = 0, std::string delimString = "\t");
#endif

  private:
    string fileName;
    Amrvis::FileType fileType;
    AmrData amrData;
    bool bAmrDataOk;
    int dsArrayIndex;
    int numberOfUsers;
    int iWriteToLevel;

    static Vector<DataServices *> dsArray;
    static int dsArrayIndexCounter;
    static int dsFabOutSize;
    static bool dsBatchMode;
    static bool profiler;
#ifdef BL_USE_PROFPARSER
    bool bProfDataAvailable;
    bool bRegionDataAvailable;
    bool bTraceDataAvailable;
    bool bCommDataAvailable;

    BLProfStats blProfStats_H;
    RegionsProfStats regOutputStats_H;
    CommProfStats commOutputStats_H;

    amrex::BoxArray procBoxArray;
#endif
};

}

#endif
