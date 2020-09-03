// ----------------------------------------------------------------------
//  CommProfStats.H
// ----------------------------------------------------------------------
#ifndef BL_COMMPROFSTATS_H
#define BL_COMMPROFSTATS_H

#include <AMReX_BLProfStats.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_SPACE.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
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
class CommProfStats : public BLProfStats {
  public:
    struct BarrierEntry {
      BarrierEntry()
        :  number(-1), name("__uninitialized_Barrier_Entry_name__"),
	   seekIndex(-1), bNameNumber(-1)  { }
      BarrierEntry(long n, const std::string &nm, long i, int bnn)
        :  number(n), name(nm), seekIndex(i), bNameNumber(bnn)  { }

      BarrierEntry &operator=(const BarrierEntry &be);

      long number;
      std::string name;
      long seekIndex;
      int bNameNumber;
    };

    struct ReductionEntry {
      ReductionEntry(const long n, const long i)
        :  number(n), seekIndex(i)  { }

      long number, seekIndex;
    };

    struct NameTagEntry {
      NameTagEntry(const long ntindex, const long seekindex)
        :  ntIndex(ntindex), seekIndex(seekindex)  { }

      long ntIndex, seekIndex;
    };

    struct DataBlock {
      DataBlock()
        :  proc(0), size(0), fileName(""), seekpos(0), readoffset(0), nodeName(""),
	   nodeNumber(-1), nTimerTimes(0), timerTime(0.0), streamIndex(-1)
	{ }
      DataBlock(const int p, const long s, const std::string &fn, const long sp,
                const std::string &pn, const int nn, const int si)
        :  proc(p), size(s), fileName(fn), seekpos(sp), readoffset(0), nodeName(pn),
	   nodeNumber(nn), nTimerTimes(0), timerTime(0.0), streamIndex(si)
	{ }

      int proc;
      long size;
      std::string fileName;
      long seekpos;
      long readoffset;
      std::string nodeName;
      int nodeNumber;
      int nTimerTimes;
      double timerTime;
      int streamIndex;

      double timeMin, timeMax;
      amrex::Vector<BarrierEntry> barriers;
      amrex::Vector<ReductionEntry> reductions;
      amrex::Vector<NameTagEntry> nameTags;
      amrex::Vector<amrex::BLProfiler::CommStats> vCommStats;
    };

    struct SendRecvPair {
      SendRecvPair(const int sp,
                   const int rp, const int ds, const int tg,
                   const Real st, const Real rt, const Real tt)
        : fromProc(sp), toProc(rp), dataSize(ds),
	  tag(tg), sendTime(st), recvTime(rt), totalTime(tt)  { }

      int fromProc, toProc, dataSize, tag;
      Real sendTime, recvTime, totalTime;
      //BLProfiler::CommFuncType unmatchedCFType;

      bool Matches(int fromproc, int toproc, int size, int mtag) {
        return(fromproc == fromProc && toproc == toProc &&
	       size == dataSize && mtag == tag);
      }
      bool Matches(int fromproc, int toproc, int size) {
        return(fromproc == fromProc && toproc == toProc &&
	       size == dataSize);
      }
      long HashLong() {
        std::stringstream hss;
	hss << fromProc << toProc << dataSize << tag;
	return atol(hss.str().c_str());
      }
      inline static long HashLong(int fp, int tp, int t) {
	return(fp + tp + 64*t);
      }
      inline static long HashLong(int fp, int tp, int /*ds*/, int t) {
	//if(fp < 0 || tp < 0 || ds < 0 || t < 0) {
	  //cout << ":::: HashLong:  negative values:  fp tp ds t = "
	       //<< fp << "  " << tp << "  " << ds << "  " << t << endl;
	//}
        //std::stringstream hss;
	//// hss.seekp(0);
	//// hss << fp << tp << ds << t;
	//// return atol(hss.str().c_str());
	//return(fp + 100*tp + 1000*ds + 10000*t);
	//return(fp + tp + 10*ds + 100*t);
	//return(fp + tp + ds + 10*t);
	//return(fp + tp + 32*t);
	return(fp + tp + 64*t);
	//return(fp + tp + ds + t);
      }
    };

    struct SendRecvPairUnpaired : public SendRecvPair {
      SendRecvPairUnpaired(const amrex::BLProfiler::CommFuncType cft, const int sp,
                           const int rp, const int ds, const int tg,
                           const Real st, const Real rt, const Real tt)
        : SendRecvPair(sp, rp, ds, tg, st, rt, tt), unmatchedCFType(cft)
	{ }

      amrex::BLProfiler::CommFuncType unmatchedCFType;
    };

  public:
    CommProfStats();
    ~CommProfStats();
    void CheckCommData(amrex::Vector<long> &nBMin, amrex::Vector<long> &nBMax,
                       amrex::Vector<long> &nRMin, amrex::Vector<long> &nRMax);
    void FillSendFAB(long &totalSends, long &totalSentData,
		     amrex::Vector<long> &totalSendsPerProc,
		     amrex::Vector<long> &totalSentDataPerProc,
		     amrex::FArrayBox &sendFAB, bool proxmap);
    void ReportSyncPointDataSetup(long &nBMax, long &nRMax);
    void ReportSyncPointData(amrex::Vector<amrex::Vector<Real> > &barrierExitTimes,
                             amrex::Vector<amrex::Vector<Real> > &barrierWaitTimes,
                             amrex::Vector<amrex::Vector<Real> > &reductionWaitTimes,
			     bool bDoReductions = false);
    void ReportStats(long &totalSentData, long &totalNCommStats,
                     amrex::Vector<long> &totalFuncCalls,
		     int bytesPerSlot, amrex::Vector<long> &msgSizes,
		     int &minMsgSize, int &maxMsgSize,
                     Real &timeMin, Real &timeMax, Real &timerTime,
		     amrex::Vector<int> &rankNodeNumbers);
    void WriteTopoFab();
    void FindTimeRange(BLProfStats::TimeRange& tr);
    void TimelineFAB(amrex::FArrayBox &timelineFAB, const amrex::Box &probDomain,
                     const BLProfStats::TimeRange tr,
		     const int rankMin, const int rankMax,
		     const int rankStride,
		     const Real ntnMultiplier, const amrex::Vector<Real> &ntnNumbers,
		     const Real bnMultiplier, const amrex::Vector<Real> &bnNumbers);
    void SendRecvData(const std::string &filename,
                       const double tlo, const double thi);
    void SendRecvList(std::multimap<Real, SendRecvPairUnpaired> &srMMap);
    virtual void InitCommDataBlock(const int proc, const long ncommstats,
                       const std::string &filename, const long seekpos,
		       const std::string &nodename = "", const int nodenumber = -1);

    int  AfterBarrier(const int proc, const double t);

    virtual void AddBarrier(long bnum, const std::string &bname, long index);
    virtual void AddReduction(const long rnum, const long index);
    virtual void AddTimeMinMax(const double tmin, const double tmax);
    virtual void AddTimerTime(const double tt);
    virtual void AddNameTag(const long ntnindex, const long seekindex);
    virtual void AddNameTagName(const std::string &name);
    const amrex::Vector<std::string> &NameTagNames() { return nameTagNames; }
    const amrex::Vector<std::string> &BarrierNames() { return barrierNames; }
    virtual void AddTagRange(const long tmin, const long tmax);
    virtual void AddGridLevel(const int level, const int ngrids);
    virtual void AddGrid3D(int level, int xlo, int ylo, int zlo,
                           int xhi, int yhi, int zhi,
                           int xc,  int yc,  int zc,
                           int xn,  int yn,  int zn, int proc);

    virtual void SetCPVersion(const int cpv)  { cpVersion = cpv; }
    static int  GetCPVersion() { return cpVersion; }
    virtual void SetCSSize(const int css)  { csSize = css; }
    static int  GetCSSize() { return csSize; }
    virtual void AddFinestLevel(const int lev) { finestLevel = lev; }
    static int  GetFinestLevel() { return finestLevel; }
    virtual void AddMaxLevel(const int lev) { maxLevel = lev; }
    static int  GetMaxLevel() { return maxLevel; }
    virtual void AddRefRatio(const int lev, const amrex::IntVect &rr);
    void SetRefRatio(const amrex::Vector<amrex::IntVect> &rr) { calcRefRatios = rr; }
    const amrex::Vector<amrex::IntVect> &GetRefRatio() const  { return calcRefRatios; }
    virtual void AddTopoCoord(const int nid, const int node,
                      const int tx, const int ty, const int tz,
                      const bool servicenode = false);
    virtual void AddProbDomain(const int lev, const amrex::Box &pd);
    void SetProbDomain(const amrex::Vector<amrex::Box> &pd) { probDomain = pd; }
    const amrex::Vector<amrex::Box> &GetProbDomain() const  { return probDomain; }
    virtual void AddCommHeaderFileName(const std::string &hfn);
    const std::map<int, long> &GLMap() const  { return glMap; }
    const std::map<int, int> &GLSizeMap() const  { return glSizeMap; }

    static void SetInitDataBlocks(bool b) { bInitDataBlocks = b;    }
    static bool InitDBlocks()             { return bInitDataBlocks; }
    static void InitDataFileNames(const amrex::Vector<std::string> &hfn);
    static void InitProxMap();
    static const amrex::Vector<std::string> &GetHeaderFileNames() { return commHeaderFileNames; }

    static void OpenAllStreams(const std::string &dirname);
    static void CloseAllStreams();

    void InitEdisonTopoMF();
    void WriteEdisonTopoMF();
    virtual void AddEdisonPID(int X, int Y, int Z,
                              int cab, int row, int cage, int slot,
                              int cpu, int pid);

    bool IsSend(const amrex::BLProfiler::CommFuncType &cft);
    bool IsRecv(const amrex::BLProfiler::CommFuncType &cft);
    bool IsBlockingRecv(const amrex::BLProfiler::CommFuncType &cft);

  private:
    int currentDataBlock;
    static int cpVersion, csSize, finestLevel, maxLevel;
    static bool persistentStreams;
    int tagMin, tagMax;
    std::stack<double> nestedTimeStack;
    amrex::Vector<DataBlock> dataBlocks;
    amrex::Vector<std::string> nameTagNames;
    amrex::Vector<std::string> barrierNames;
    static amrex::Vector<amrex::IntVect> calcRefRatios;
    amrex::Vector<std::map<int, amrex::IntVect> > TopoMap;
    static amrex::Vector<amrex::Box> probDomain;
    std::map<int, long> glMap;     // [proc, npoints]
    std::map<int, int> glSizeMap;  // [grid size, ngrids of that size]
    static bool bInitDataBlocks;
    static amrex::Vector<int> rankFromProx;  // [prox]
    static amrex::Vector<int> proxFromRank;  // [rank]
    static bool bProxMapOK;
    static amrex::Vector<std::string> commHeaderFileNames;
    static std::map<std::string, int> commDataFileNames;  // [filename, stream index]
    static amrex::Vector<std::ifstream *> commDataStreams;

    void ReadCommStats(DataBlock &dBlock);  // reads whole block
    void ReadCommStatsNoOpen(DataBlock &dBlock);  // files must be open already
    bool ReadCommStats(DataBlock &dBlock, const int nmessages);  // reads nmessages
    void ClearCommStats(DataBlock &dBlock);

    // ---- Edison support
    amrex::Vector<Real> edisonNodeFab;  // ---- not a fab yet
    amrex::Vector<Real> edisonCPUFab;   // ---- not a fab yet
    amrex::Vector<amrex::Vector<amrex::Vector<int> > > edisonNodeXYZ;   // ---- [X][Y][Z]

    friend int yyparse(void *);
};


std::ostream &operator<<(std::ostream &os, const CommProfStats::BarrierEntry &be);

#endif
// ----------------------------------------------------------------------
// ----------------------------------------------------------------------
