
// ---------------------------------------------------------------
// DataServices.H
// ---------------------------------------------------------------
#ifndef _DATASERVICES_H_
#define _DATASERVICES_H_

#include <Array.H>
#include <AmrData.H>

#include <cstdarg>
#include <string>
using std::string;

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
#ifndef BL_NOLINEVALUES
      LineValuesRequest
#endif
    };

    DataServices(const string &filename, const Amrvis::FileType &filetype);
    DataServices();
    void Init(const string &filename, const Amrvis::FileType &filetype);
    ~DataServices();

    static void SetBatchMode();
    static void Dispatch(DSRequestType requestType, DataServices *ds, ...);

    static void SetFabOutSize(int iSize);

  private:

    // Dispatchable functions

      // FillVarOneFab
      bool FillVar(FArrayBox *destFab, const Box &destBox,
                   int finestFillLevel, const string &varname,
                   int procWithFab);

      // FillVarArrayOfFabs
      //bool FillVar(Array<FArrayBox *> &destFabs, const Array<Box> &destBoxes,
                   //int finestFillLevel, const string &varname,
                   //int procWithFab);

      // FillVarMultiFab
      bool FillVar(MultiFab &destMultiFab, int finestFillLevel,
                   const string &varname);

      // WriteFabOneVar:  write a single variable to one fab
      bool WriteFab(const string &fabfilename, const Box &region, int lev,
                    const string &varname);

      // WriteFabAllVars:  write all vars in one fab
      bool WriteFab(const string &fabfilename, const Box &region, int lev);

      // DumpSlicePlaneOneVar
      bool DumpSlice(int slicedir, int slicenum, const string &varname);

      // DumpSlicePlaneAllVars
      bool DumpSlice(int slicedir, int slicenum);

      // DumpSliceBoxOneVar
      bool DumpSlice(const Box &slicebox, const string &varname);

      // DumpSliceBoxAllVars
      bool DumpSlice(const Box &slicebox);

      // MinMaxRequest
      bool MinMax(const Box &onBox, const string &derived, int level,
                  Real &dataMin, Real &dataMax, bool &minMaxValid);

      // PointValueRequest
      void PointValue(int pointBoxArraySize, Box *pointBoxArray,
		      const string &currentDerived,
		      int coarsestLevelToSearch, int finestLevelToSearch,
		      int &intersectedLevel, Box &intersectedGrid,
		      Real &dataPointValue, bool &bPointIsValid);

#ifndef BL_NOLINEVALUES
      // LineValuesRequest
      void LineValues(int lineBoxArraySize, Box *lineBoxArray, int whichDir,
                      const string &currentDerived,
                      int coarsestLevelToSearch, int finestLevelToSearch,
                      XYPlotDataList *dataList, bool &lineIsValid);
#endif

  public:

    // information functions
    const string &GetFileName() const { return fileName; }
    Amrvis::FileType GetFileType()      const { return fileType; }
    bool CanDerive(const string &name) const;
    bool CanDerive(const Array<string> &names) const;
    int  NumDeriveFunc() const;
    void ListDeriveFunc(std::ostream &os) const;
    const Array<string> &PlotVarNames() const { return amrData.PlotVarNames(); }
    AmrData &AmrDataRef() { return amrData; }
    bool AmrDataOk() { return bAmrDataOk; }  
    void SetWriteToLevel(const int iwtlev)  { iWriteToLevel = iwtlev; }

    // other functions
    int GetNumberOfUsers() const  { return numberOfUsers; }
    void IncrementNumberOfUsers() { ++numberOfUsers; }
    void DecrementNumberOfUsers() { --numberOfUsers; }


  private:
    string fileName;
    Amrvis::FileType fileType;
    AmrData amrData;
    bool bAmrDataOk;
    int dsArrayIndex;
    int numberOfUsers;
    int iWriteToLevel;

    static Array<DataServices *> dsArray;
    static int dsArrayIndexCounter;
    static int dsFabOutSize;
    static bool dsBatchMode;
};

#endif
