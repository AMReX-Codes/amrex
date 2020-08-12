// ----------------------------------------------------------------------
//  BLProfStats.H
// ----------------------------------------------------------------------
#ifndef BL_BLPROFSTATS_H
#define BL_BLPROFSTATS_H

#include <AMReX_BLProfiler.H>
#include <AMReX_SPACE.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <map>
#include <sstream>

namespace {
  static std::stringstream hs;
}

using amrex::Real;


// ----------------------------------------------------------------------
class BLProfStats {
  public:

    enum FilterStatus { UNDEFINED, ON, OFF, INCLUDEALL, INCLUDENONE };

    struct BLPDataBlock {
      BLPDataBlock()
        :  proc(0), fileName(""), seekpos(0), readoffset(0),
           streamIndex(-1)
        { }
      BLPDataBlock(int p, const std::string &fn, long sp, int si)
        :  proc(p), fileName(fn), seekpos(sp),
           readoffset(0), streamIndex(si)
        { }

      int proc;
      std::string fileName;
      long seekpos;
      long readoffset;
      int streamIndex;

      amrex::Vector<long> nCalls;
      amrex::Vector<Real> totalTime;
    };

    struct FuncStat {
      FuncStat()
        : nCalls(0), totalTime(0.0)
      { }

      FuncStat(long ncalls, Real tt)
        : nCalls(ncalls), totalTime(tt)
      { }

      long nCalls;
      Real totalTime;
    };

    struct CallTreeNode {
      CallTreeNode()
        : fnameNumber(-1), stackDepth(-1), nCalls(0), totalTime(0.0), stackTime(0.0)
      { }

      int fnameNumber, stackDepth;
      long nCalls;
      Real totalTime, stackTime;
      std::map<int, CallTreeNode> calledFunctions;  // [fnn, ctn]
    };

    struct TimeRange {
      TimeRange() : startTime(0.0), stopTime(0.0)  { }
      TimeRange(Real startT, Real stopT) : startTime(startT), stopTime(stopT)  { }

      Real startTime, stopTime;
      bool Contains(Real t) const {
        return(t >= startTime && t <= stopTime);
      }
      bool operator==(const TimeRange &rhs) const {
        return((startTime == rhs.startTime) && (stopTime == rhs.stopTime));
      }
      bool operator!=(const TimeRange &rhs) const {
        return((startTime != rhs.startTime) || (stopTime != rhs.stopTime));
      }
      TimeRange operator=(const TimeRange &rhs) {
        if (rhs != *this)
        {
          startTime = rhs.startTime;
          stopTime = rhs.stopTime;
        }
        return *this;
      }
      std::ostream& operator<<(std::ostream& os) {
        os << startTime << " - " << stopTime;
        if (os.fail()) { amrex::Error("output for TimeRange failed"); }
        return os;
      }
    };

    struct TimeRangeCompare {
      bool operator()(const TimeRange &lhs, const TimeRange &rhs) const {
        return lhs.startTime < rhs.startTime;
      }
    };

    BLProfStats();
    ~BLProfStats();
    static int  Verbose() { return verbose; }
    static void SetVerbose(int vlevel = 0) { verbose = vlevel; }
    static void SetDirName(const std::string &dirname)  { dirName = dirname; }
    static void SetNOutFiles(const int nof)  { nOutFiles = nof; }
    static int  GetNOutFiles() { return nOutFiles; }

    void AddRegionName(const std::string &rname, int rnumber);
    std::set<int> WhichRegions(int proc, Real t);

    void SetProc(int p)  { currentProc = p; }
    int  GetProc() const  { return currentProc; }

    static void SetNProcs(int p)  { dataNProcs = p; }
    static int  GetNProcs() { return dataNProcs; }

    void SetMaxRNumber(int n)  { maxRNumber = n; }
    int  GetMaxRNumber() const  { return maxRNumber; }

    void InitFilterTimeRanges();
    void SetFilterTimeRanges(const amrex::Vector<std::list<TimeRange> > &ftr);
    const amrex::Vector<std::list<TimeRange> > &GetFilterTimeRanges() const
                                                { return filterTimeRanges; }

    void SetRegionTimeRanges(const amrex::Vector<amrex::Vector<amrex::Vector<TimeRange> > > &rtr)
                                             { regionTimeRanges = rtr; }
    void SetRegionTimeRanges(const TimeRange &tr, const int n, const int r, const int t)
                                             { regionTimeRanges[n][r][t] = tr; }
    const amrex::Vector<amrex::Vector<amrex::Vector<TimeRange> > > &GetRegionTimeRanges()
                                             { return regionTimeRanges; }

    const std::map<std::string, int> &RegionNames() { return regionNames; }
    const std::map<int, std::string> &RegionNumbers() { return regionNumbers; }

    //const std::map<std::string, int> &BarrierNamesToNumbers() { return barrierNamesToNumbers; }
    const std::map<int, std::string> &BarrierNumbersToNames() { return barrierNumbersToNames; }

    void MakeFilterFile(const std::string &ffname);
    void SetFilter(FilterStatus fs, const std::string &rname, int rnumber);
    void SetFilter(FilterStatus fs);

    static bool AddPiece(std::list<TimeRange> &addToHere,
                         const TimeRange &pieceToAdd);
    static std::list<TimeRange> RangeIntersection(std::list<TimeRange> &rangeList,
                                                  const TimeRange &pieceToIntersect);
    static bool RemovePiece(std::list<TimeRange> &removeFromHere,
                            const TimeRange &pieceToRemove);

    virtual void SetBLPVersion(const int blpv) { blProfVersion = blpv; }
    virtual void AddFunctionName(const std::string &fname);
    virtual void InitBLProfDataBlock(const int proc, const std::string &filename,
                                     const long seekpos);
    virtual void CollectFuncStats(amrex::Vector<amrex::Vector<FuncStat> > &funcStats);
    virtual void WriteSummary(std::ostream &ios, bool bwriteavg = false, int whichProc = 0,
                              bool graphTopPct = true);
    virtual void AddCalcEndTime(double cet) { calcEndTime = cet; }

    virtual void SetCPVersion(const int /*cpv*/) { }
    virtual void SetCSSize(const int /*css*/) { }
    virtual void InitCommDataBlock(const int /*proc*/, const long /*ncommstats*/,
                               const std::string &/*filename*/, const long /*seekpos*/,
			       const std::string &/*nodename*/ = "",
			       const int /*nodenumber*/ = -1) { }
    virtual void AddBarrier(long /*bnum*/, const std::string &/*bname*/, long /*index*/) { }
    virtual void AddReduction(const long /*rnum*/, const long /*index*/) { }
    virtual void AddTimeMinMax(const double /*tmin*/, const double /*tmax*/) { }
    virtual void AddTimerTime(const double /*tt*/) { }
    virtual void AddNameTag(const long /*ntnindex*/, const long /*seekindex*/) { }
    virtual void AddNameTagName(const std::string &/*name*/) { }
    virtual void AddTagRange(const long /*tmin*/, const long /*tmax*/) { }
    virtual void AddGridLevel(const int /*level*/, const int /*ngrids*/) { }
    virtual void AddGrid3D(int /*level*/, int /*xlo*/, int /*ylo*/, int /*zlo*/,
                           int /*xhi*/,   int /*yhi*/, int /*zhi*/,
                           int /*xc*/,    int /*yc*/,  int /*zc*/,
                           int /*xn*/,    int /*yn*/,  int /*zn*/, int /*proc*/) { }
    virtual void AddFinestLevel(const int /*lev*/) { }
    virtual void AddMaxLevel(const int /*lev*/) { }
    virtual void AddRefRatio(const int /*lev*/, const amrex::IntVect &/*rr*/) { }
    virtual void AddTopoCoord(const int /*nid*/, const int /*node*/,
                      const int /*tx*/, const int /*ty*/, const int /*tz*/,
                      const bool /*servicenode*/ = false) { }

    virtual void AddCommHeaderFileName(const std::string &/*hfn*/) { }
    virtual void SetCSVersion(int /*csv*/)  { }
    virtual void AddCStatsHeaderFileName(const std::string &/*hfn*/) { }
    virtual void AddFunctionName(const std::string &/*fname*/, int /*fnumber*/) { }
    virtual void InitCStatsDataBlock(int /*proc*/, long /*nrss*/, long /*ntracestats*/,
                                     const std::string &/*filename*/, long /*seekpos*/) { }
    static void SetInitDataBlocks(bool b) { bInitDataBlocks = b;    }
    static bool InitDBlocks()             { return bInitDataBlocks; }

    static void SetGPercent(Real p)  { gPercent = p/100.0; }

    virtual void AddProbDomain(const int /*lev*/, const amrex::Box &/*pd*/) { }
    virtual TimeRange MakeRegionPlt(amrex::FArrayBox &/*rFab*/, int /*noregionnumber*/,
                               int /*width*/, int /*height*/,
			       amrex::Vector<amrex::Vector<amrex::Box>> &/*regionBoxes*/)
			       { return TimeRange(); }

    static void OpenAllStreams(const std::string &dirname);
    static void CloseAllStreams();

    virtual void AddEdisonPID(int /*X*/,   int /*Y*/, int /*Z*/,
                              int /*col*/, int /*row*/, int /*cage*/, int /*slot*/,
                              int /*cpu*/, int /*pid*/) { }

    static bool TimeRangeInitialized() { return bTimeRangeInitialized; }

    void CheckData();

    const amrex::Vector<std::string> &BLPFNames() { return blpFNames; }
    void SetFNames(const amrex::Vector<std::string> &fNames) { blpFNames = fNames; }

  protected:
    static int verbose;
    static int blProfVersion, dataNProcs, nOutFiles;
    int maxRNumber, currentProc, currentDataBlock;
    Real calcEndTime;
    Real minRegionTime, maxRegionTime;
    static std::string dirName;
    static bool bInitDataBlocks;

    static bool bTimeRangeInitialized;
    bool InTimeRange(int proc, Real calltime);

    std::map<std::string, int> regionNames;    // [rname, rnumber]
    std::map<int, std::string> regionNumbers;  // [rnumber, rname]

    std::map<int, std::string> barrierNumbersToNames;  // [bnumber, bname]

    std::set<int> includeSet, excludeSet;
    bool bDefaultInclude;

    amrex::Vector<amrex::Vector<amrex::Vector<TimeRange> > > regionTimeRanges;  // [proc][rnum][range]
    amrex::Vector<std::list<TimeRange> > filterTimeRanges;  // [proc]

    static std::map<std::string, int> blpDataFileNames;  // [filename, stream index]

    amrex::Vector<std::string> blpFNames;
    amrex::Vector<BLPDataBlock> blpDataBlocks;

    static Real gPercent;

    static amrex::Vector<std::ifstream *> blpDataStreams;

  private:

    void ReadBlock(BLPDataBlock &dBlock);  // reads whole block
    void ReadBlockNoOpen(BLPDataBlock &dBlock);  // files must be open already
    void ClearBlock(BLPDataBlock &dBlock);

    friend int yyparse(void *);
};

std::ostream &operator<< (std::ostream &os, const BLProfStats::TimeRange &tr);

#endif
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
