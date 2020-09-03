// ----------------------------------------------------------------------
//  RegionsProfStats.H
// ----------------------------------------------------------------------
#ifndef BL_REGIONSPROFSTATS_H
#define BL_REGIONSPROFSTATS_H

#include <AMReX_BLProfStats.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <stack>
#include <sstream>
#include <stdlib.h>
using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::multimap;

using amrex::Real;


// ----------------------------------------------------------------------
class RegionsProfStats : public BLProfStats {
  public:

    struct DataBlock {
      DataBlock()
        :  proc(0), nRSS(0), nTraceStats(0), fileName(""), seekpos(0), readoffset(0),
           streamIndex(-1)
        { }
      DataBlock(int p, long nrss, long nts, const std::string &fn, long sp, int si)
        :  proc(p), nRSS(nrss), nTraceStats(nts), fileName(fn), seekpos(sp),
           readoffset(0), streamIndex(si)
        { }

      int proc;
      long nRSS, nTraceStats;
      std::string fileName;
      long seekpos;
      long readoffset;
      int streamIndex;

      double timeMin, timeMax;
      amrex::Vector<amrex::BLProfiler::RStartStop> rStartStop;
      amrex::Vector<amrex::BLProfiler::CallStats> vCallStats;
    };


  public:
    RegionsProfStats();
    ~RegionsProfStats();

    virtual void AddCStatsHeaderFileName(const std::string &hfn);

    virtual void SetCSVersion(int csv)  { csVersion = csv; }
    int  GetCSVersion() const  { return csVersion; }

    void SetMaxFNumber(int n)  { maxFNumber = n; }
    int  GetMaxFNumber() const  { return maxFNumber; }

    void AddFunctionName(const std::string &/*fname*/)
        { amrex::Abort("Wrong AddFunctionName for RegionsProfStats. Use (string, int)."); }
    virtual void AddFunctionName(const std::string &fname, int fnumber);
    virtual void AddTimeMinMax(double tmin, double tmax);

    TimeRange MakeRegionPlt(amrex::FArrayBox &rFab, int noregionnumber,
                            int width, int height,
			    amrex::Vector<amrex::Vector<amrex::Box>> &regionBoxes);
    virtual void CollectFuncStats(amrex::Vector<amrex::Vector<FuncStat> > &funcStats);
    virtual void WriteSummary(std::ostream &os, bool bwriteavg = false, int whichProc = 0,
                              bool graphTopPct = true);
    bool AllCallTimesFAB(amrex::FArrayBox &actFab, const std::string &whichFuncName);
    void FillAllCallTimes(amrex::Vector<amrex::Vector<Real> > &allCallTimes,
                          const std::string whichFuncName,
                          int whichFuncNumber, const amrex::Box &procBox);
    void CheckRegionsData();

    virtual void InitCStatsDataBlock(int proc, long nrss, long ntracestats,
                                     const std::string &filename, long seekpos);

    static void SetInitDataBlocks(bool b) { bInitDataBlocks = b;    }
    static bool InitDBlocks()             { return bInitDataBlocks; }
    static void InitDataFileNames(const amrex::Vector<std::string> &hfn);
    bool TraceDataValid() const { return mFNameNumbersPerProc.size() > 0; }
		     
    void SyncFNamesAndNumbers();
    bool InitRegionTimeRanges(const amrex::Box &procBox);

    void WriteHTML(std::ostream &os, bool simpleCombine = true, int whichProc = 0);
    void PrintCallTreeNode(CallTreeNode &callTree, amrex::Vector<std::string> &fNumberNames);
    void CreateVCallStats(CallTreeNode &callTree,
                          amrex::Vector<amrex::BLProfiler::CallStats> &vCallStatsNC);
    void WriteHTMLNC(std::ostream &os, int whichProc = 0);
    void WriteTextTrace(std::ostream &os, bool simpleCombine = true, int whichProc = 0,
                        std::string delimString = "\t");

    static const amrex::Vector<std::string> &GetHeaderFileNames() { return regHeaderFileNames; }
    const amrex::Vector<std::string> &NumbersToFName() const { return numbersToFName; }

    static void OpenAllStreams(const std::string &dirname);
    static void CloseAllStreams();

    void FillRegionTimeRanges(amrex::Vector<amrex::Vector<TimeRange>> &rtr,  // ---- [region][range]
                              int whichProc);

  private:
    int currentDataBlock;
    int csVersion;
    int maxFNumber;

    amrex::Vector<std::map<std::string, int> > mFNameNumbersPerProc;  // [proc][fname, fnumber]
    amrex::Vector<amrex::Vector<int> > fnameRemap;  // [dataNProc][remappedfnumber]
    std::map<std::string, int> mFNameNumbersLocal;  // [fname, fnumber]
    amrex::Vector<std::string> numbersToFName;  // [fname]
    static amrex::Vector<std::string> regHeaderFileNames;
    static amrex::Vector<std::ifstream *> regDataStreams;
    static std::map<std::string, int> regDataFileNames;  // [filename, stream index]

    static std::map<Real, std::string, std::greater<Real> > mTimersTotalsSorted;

    bool Include(const FuncStat &fs);

    amrex::Vector<DataBlock> dataBlocks;
    static bool bInitDataBlocks, persistentStreams;

    void ReadBlock(DataBlock &dBlock, bool readRSS = true,
                   bool readTraces = true);  // reads whole block
    void ReadBlockNoOpen(DataBlock &dBlock, bool readRSS = true,
                   bool readTraces = true);  // files must be open already
    bool ReadBlock(DataBlock &dBlock, const int nmessages);  // reads nmessages
    void ClearBlock(DataBlock &dBlock);

    friend int yyparse(void *);

};

#endif
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
