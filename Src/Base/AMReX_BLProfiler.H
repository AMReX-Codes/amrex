#ifndef BL_PROFILER_H_
#define BL_PROFILER_H_

#ifdef BL_PROFILING

#include <AMReX_REAL.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#ifndef BL_AMRPROF
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#endif

#include <iostream>
#include <ostream>
#include <string>
#include <stack>
#include <set>
#include <map>
#include <list>
#include <utility>
#include <typeinfo>

namespace amrex {

class BLProfiler 
{
  public:
    struct ProfStats {
      ProfStats() : nCalls(0), totalTime(0.0), minTime(0.0),
                    maxTime(0.0), avgTime(0.0), variance(0.0) { }
      long nCalls;
      Real totalTime, minTime, maxTime, avgTime, variance;
    };

    struct CallStats {
      CallStats() : callStackDepth(-2), csFNameNumber(-2),
                    nCSCalls(0), totalTime(0.0),
		    stackTime(0.0), callTime(0.0)  { }
      CallStats(int depth, int fnamenumber) : callStackDepth(depth),
                                              csFNameNumber(fnamenumber),
                                              nCSCalls(0), totalTime(0.0),
					      stackTime(0.0), callTime(0.0)  { }
      CallStats(int depth, int fnamenumber, int ncalls, Real totaltime,
                Real stacktime, Real calltime) : callStackDepth(depth),
                                                 csFNameNumber(fnamenumber),
                                                 nCSCalls(ncalls), totalTime(totaltime),
					         stackTime(stacktime), callTime(calltime)  { }
      int callStackDepth, csFNameNumber;
      long nCSCalls;
      Real totalTime, stackTime, callTime;

      static int cstatsVersion;
      static Real minCallTime, maxCallTime;
    };

    struct CallStatsStack {
      CallStatsStack()        : bFlushed(false), index(-1)   { }
      CallStatsStack(int idx) : bFlushed(false), index(idx)  { }
      bool bFlushed;
      int index;
    };

    struct CallStatsPatch {
      CallStatsPatch() : seekPos(-1), callStats(), fileName("")        { }
      CallStatsPatch(long spos, const CallStats &cs, const std::string &fname)
                       : seekPos(spos), callStats(cs), fileName(fname) { }
      long seekPos;
      CallStats callStats;
      std::string fileName;
    };

    struct RStartStop {
      RStartStop()
        : rssStart(false), rssRNumber(-2), rssTime(0.0)  { }
      RStartStop(bool s, int r, Real t)
        : rssStart(s), rssRNumber(r), rssTime(t)  { }
      bool rssStart;
      int  rssRNumber;
      Real rssTime;
    };

    enum CommFuncType {
      InvalidCFT = 0,  //  0
      AllReduceT,      //  1
      AllReduceR,      //  2
      AllReduceL,      //  3
      AllReduceI,      //  4
      AsendTsii,       //  5
      AsendTsiiM,      //  6
      AsendvTii,       //  7
      SendTsii,        //  8
      SendvTii,        //  9
      ArecvTsii,       // 10
      ArecvTsiiM,      // 11
      ArecvTii,        // 12
      ArecvvTii,       // 13
      RecvTsii,        // 14
      RecvvTii,        // 15
      ReduceT,         // 16
      ReduceR,         // 17
      ReduceL,         // 18
      ReduceI,         // 19
      BCastTsi,        // 20
      GatherTsT1Si,    // 21
      GatherTi,        // 22
      GatherRiRi,      // 23
      ScatterTsT1si,   // 24
      Barrier,         // 25
      Waitsome,        // 26
      NameTag,         // 27
      AllCFTypes,      // 28
      NoCFTypes,       // 29
      IOStart,         // 30
      IOEnd,           // 31
      TagWrap,         // 32
      Allgather,       // 33
      Alltoall,        // 34
      Alltoallv,       // 35
      Gatherv,         // 36
      Get_count,       // 37
      Iprobe,          // 38
      Test,            // 39
      Wait,            // 40
      Waitall,         // 41
      NUMBER_OF_CFTS   // 42
    };

    struct CommStats {
      CommStats() : cfType(InvalidCFT), size(-2), commpid(-2), tag(-2),
                    timeStamp(-2.0) { }
      CommStats(CommFuncType cft, int sz, int cpid, int tg, Real ts)
                  : cfType(cft), size(sz), commpid(cpid), tag(tg),
                    timeStamp(ts) { }

      static std::string CFTToString(CommFuncType cft);
      static CommFuncType StringToCFT(const std::string &s);
      static void Filter(CommFuncType cft);
      static void UnFilter(CommFuncType cft);

      static std::map<std::string, CommFuncType> cftNames;
      static std::set<CommFuncType> cftExclude;
      static int barrierNumber;
      static int reductionNumber;
      static int tagWrapNumber;
      static int tagMin;
      static int tagMax;
      static Array<std::pair<std::string,int> > barrierNames;  // [name, seek]
      static Array<std::pair<int,int> > nameTags;              // [nameindex, seek]
      static Array<std::string> nameTagNames;                  // [name]
      static Array<int> tagWraps;                              // [index]
      static int csVersion;

      CommFuncType cfType;
      int size, commpid, tag;
      Real timeStamp;
    };

    static std::map<std::string, BLProfiler *> mFortProfs;  // [fname, fortfunc]
    static Array<std::string> mFortProfsErrors;     // [error string]
    static Array<BLProfiler *> mFortProfsInt;     // [fortfuncindex]
    static Array<std::string> mFortProfsIntNames;     // [fortfuncindex name]

    explicit BLProfiler(const std::string &funcname);
    BLProfiler(const std::string &funcname, bool bstart);

    ~BLProfiler();

    static void Initialize();
    static void Finalize();

    static void WriteCallTrace(bool bFlushing = false);
    static void WriteCommStats(bool bFlushing = false);
    static void WriteFortProfErrors();

    static void AddCommStat(const CommFuncType cft, const int size,
                            const int pid, const int tag);
    static void AddWaitsome(const CommFuncType cft, const Array<MPI_Request> &reqs,
                            const int completed, const Array<int> &indx,
			    const Array<MPI_Status> &status, const bool bc);
    static void AddBarrier(const std::string &message, const bool bc);
    static void AddNameTag(const std::string &name);
    static void AddAllReduce(const CommFuncType cft, const int size, const bool bc);
    static void TagRange(const int tagmin, const int tagmax);
    static void AddTagWrap();
#ifndef BL_AMRPROF
    static void InitAMR(const int flev, const int mlev, const Array<IntVect> &rr,
                        const Array<Box> pd);
#endif

    void start();
    void stop();
    void PStart();
    void PStop();
    static void ChangeFortIntName(const std::string &fname, int intname);

    static void InitParams(const Real ptl, const bool writeall,
                           const bool writefabs);
    static void AddStep(const int snum);
    static Real GetRunTime() { return calcRunTime; }
    static void SetRunTime(Real rtime) { calcRunTime = rtime; }

    static void RegionStart(const std::string &rname);
    static void RegionStop(const std::string &rname);

    static inline int NoTag()      { return -3; }
    static inline int BeforeCall() { return -5; }
    static inline int AfterCall()  { return -7; }

    static inline int ProcNumber() { return procNumber; }
    static void SetBlProfDirName(const std::string &name) { blProfDirName = name; }
    static void SetNoOutput() { bNoOutput = true; }

    static int GetCSFlushSize() { return csFlushSize; }
    static void SetCSFlushSize(int fsize) { csFlushSize = fsize; }
    static int GetTraceFlushSize() { return traceFlushSize; }
    static void SetTraceFlushSize(int fsize) { traceFlushSize = fsize; }

    static void SetNFiles(int nfiles) { nProfFiles = nfiles; }

  private:
    Real bltstart, bltelapsed;
    std::string fname;
    bool bRunning;

    static bool bWriteAll, bWriteFabs, groupSets;
    static bool bFirstCommWrite;
    static bool bInitialized, bNoOutput;
    static int  currentStep, csFlushSize, traceFlushSize, nProfFiles;
    static int  finestLevel, maxLevel;
    static Real pctTimeLimit;
    static Real calcRunTime;
    static Real startTime;
    static Real timerTime;
#ifndef BL_AMRPROF
    static Array<IntVect> refRatio;
    static Array<Box> probDomain;
#endif
    static std::stack<Real> nestedTimeStack;
    static std::map<int, Real> mStepMap;  // [step, time]
    static std::map<std::string, ProfStats> mProfStats;  // [fname, pstats]
    static Array<CommStats> vCommStats;
    static std::string procName;
    static int procNumber;
    static bool blProfDirCreated;
    static std::string blProfDirName;
    static int BLProfVersion;

    static bool OnExcludeList(CommFuncType cft);
    static int  NameTagNameIndex(const std::string &name);

    static std::map<std::string, int> mFNameNumbers;  // [fname, fnamenumber]
    static Array<CallStats> vCallTrace;

    // region support
    static std::map<std::string, int> mRegionNameNumbers;  // [rname, rnamenumber]
    static int inNRegions;
    static Array<RStartStop> rStartStop;
    static const std::string noRegionName;

    static bool bFirstTraceWrite;

    static Array<CallStatsStack> callIndexStack;  // need Array for iterator
    static Array<CallStatsPatch> callIndexPatch;

#ifdef BL_TRACE_PROFILING
    static int callStackDepth;
    static int prevCallStackDepth;

  public:
    typedef std::pair<Real, int> RIpair;
    struct fTTComp {
      bool operator()(const RIpair &lhs, const RIpair &rhs) const {
        return lhs.first > rhs.first;
      }
    };
#endif

};


namespace BLProfilerUtils {
  void WriteHeader(std::ostream &os, const int colWidth,
                   const Real maxlen, const bool bwriteavg);

  void WriteRow(std::ostream &os, const std::string &fname,
                const BLProfiler::ProfStats &pstats, const Real percent,
	        const int colWidth, const Real maxlen,
	        const bool bwriteavg);

  void WriteStats(std::ostream &os,
                  const std::map<std::string, BLProfiler::ProfStats> &mpStats,
		  const std::map<std::string, int> &fnameNumbers,
		  const Array<BLProfiler::CallStats> &callTraces,
		  bool bwriteavg = false, bool bwriteinclusivetimes = false);
}


inline std::ostream &operator<< (std::ostream &os, const BLProfiler::CommStats &cs) {
  os << BLProfiler::CommStats::CFTToString(cs.cfType) << "   " << cs.size
     << "  " << cs.commpid << "  " << cs.tag << "  " << cs.timeStamp;
  return os;
}



inline std::string BLProfiler::CommStats::CFTToString(CommFuncType cft) {
  switch(cft) {
    case InvalidCFT:     return "InvalidCFT";
    case AllReduceT:     return "AllReduceT";
    case AllReduceR:     return "AllReduceR";
    case AllReduceL:     return "AllReduceL";
    case AllReduceI:     return "AllReduceI";
    case AsendTsii:      return "AsendTsii";
    case AsendTsiiM:     return "AsendTsiiM";
    case AsendvTii:      return "AsendvTii";
    case SendTsii:       return "SendTsii";
    case SendvTii:       return "SendvTii";
    case ArecvTsii:      return "ArecvTsii";
    case ArecvTsiiM:     return "ArecvTsiiM";
    case ArecvTii:       return "ArecvTii";
    case ArecvvTii:      return "ArecvvTii";
    case RecvTsii:       return "RecvTsii";
    case RecvvTii:       return "RecvvTii";
    case ReduceT:        return "ReduceT";
    case ReduceR:        return "ReduceR";
    case ReduceL:        return "ReduceL";
    case ReduceI:        return "ReduceI";
    case BCastTsi:       return "BCastTsi";
    case GatherTsT1Si:   return "GatherTsT1Si";
    case GatherTi:       return "GatherTi";
    case GatherRiRi:     return "GatherRiRi";
    case ScatterTsT1si:  return "ScatterTsT1si";
    case Barrier:        return "Barrier";
    case Waitsome:       return "Waitsome";
    case NameTag:        return "NameTag";
    case AllCFTypes:     return "AllCFTypes";
    case NoCFTypes:      return "NoCFTypes";
    case IOStart:        return "IOStart";
    case IOEnd:          return "IOEnd";
    case TagWrap:        return "TagWrap";
    case Allgather:      return "Allgather";
    case Alltoall:       return "Alltoall";
    case Alltoallv:      return "Alltoallv";
    case Gatherv:        return "Gatherv";
    case Get_count:      return "Get_count";
    case Iprobe:         return "Iprobe";
    case Test:           return "Test";
    case Wait:           return "Wait";
    case Waitall:        return "Waitall";
    case NUMBER_OF_CFTS: return "NUMBER_OF_CFTS";
  }
  return "*** Error: Bad CommFuncType.";
}

}

#define BL_PROFILE_INITIALIZE()  amrex::BLProfiler::Initialize();
#define BL_PROFILE_FINALIZE()    amrex::BLProfiler::Finalize();

#define BL_PROFILE(fname) amrex::BLProfiler bl_profiler__((fname));
#define BL_PROFILE_T(fname, T) amrex::BLProfiler bl_profiler__((std::string(fname) + typeid(T).name()));                
#ifdef BL_PROFILING_SPECIAL
#define BL_PROFILE_S(fname) amrex::BLProfiler bl_profiler__((fname));
#define BL_PROFILE_T_S(fname, T) amrex::BLProfiler bl_profiler__((std::string(fname) + typeid(T).name()));                
#else
#define BL_PROFILE_S(fname)
#define BL_PROFILE_T_S(fname, T)
#endif
 
#define BL_PROFILE_VAR(fname, vname) amrex::BLProfiler bl_profiler__##vname((fname));
#define BL_PROFILE_VAR_NS(fname, vname) amrex::BLProfiler bl_profiler__##vname(fname, false);
#define BL_PROFILE_VAR_START(vname) bl_profiler__##vname.start();
#define BL_PROFILE_VAR_STOP(vname) bl_profiler__##vname.stop();

#define BL_PROFILE_INIT_PARAMS(ptl,wall,wfabs)  \
    amrex::BLProfiler::InitParams(ptl,wall, wfabs);
#define BL_PROFILE_ADD_STEP(snum)  amrex::BLProfiler::AddStep(snum);
#define BL_PROFILE_SET_RUN_TIME(rtime)  amrex::BLProfiler::SetRunTime(rtime);

#define BL_PROFILE_REGION_START(rname) amrex::BLProfiler::RegionStart(rname);
#define BL_PROFILE_REGION_STOP(rname)  amrex::BLProfiler::RegionStop(rname);

// these combine regions with profile variables
#define BL_PROFILE_REGION_VAR(fname, rvname) amrex::BLProfiler::RegionStart(fname); \
    amrex::BLProfiler bl_profiler__##rvname((fname));
#define BL_PROFILE_REGION_VAR_START(fname, rvname) amrex::BLProfiler::RegionStart(fname); \
                                            bl_profiler__##rvname.start();
#define BL_PROFILE_REGION_VAR_STOP(fname, rvname)  bl_profiler__##rvname.stop();   \
    amrex::BLProfiler::RegionStop(fname);

#define BL_TRACE_PROFILE_FLUSH() { amrex::BLProfiler::WriteCallTrace(true); }
#define BL_TRACE_PROFILE_SETFLUSHSIZE(fsize) { amrex::BLProfiler::SetTraceFlushSize(fsize); }

#define BL_PROFILE_CHANGE_FORT_INT_NAME(fname, intname) { amrex::BLProfiler::ChangeFortIntName(fname, intname); }

#ifdef BL_COMM_PROFILING

#define BL_COMM_PROFILE(cft, size, pid, tag) {                       \
        amrex::BLProfiler::AddCommStat(cft, size, pid, tag);         \
}
#define BL_COMM_PROFILE_BARRIER(message, bc) { amrex::BLProfiler::AddBarrier(message, bc); }
#define BL_COMM_PROFILE_ALLREDUCE(cft, size, bc) { amrex::BLProfiler::AddAllReduce(cft, size, bc); }
#define BL_COMM_PROFILE_REDUCE(cft, size, pid) {   \
        amrex::BLProfiler::AddCommStat(cft, size, pid, amrex::BLProfiler::NoTag()); }
#define BL_COMM_PROFILE_WAITSOME(cft, reqs, completed, indx, status, bc) {                \
        amrex::BLProfiler::AddWaitsome(cft, reqs, completed, indx, status, bc); \
}
#define BL_COMM_PROFILE_NAMETAG(message) { amrex::BLProfiler::AddNameTag(message); }
#define BL_COMM_PROFILE_FILTER(cft) { amrex::BLProfiler::CommStats::Filter(cft); }
#define BL_COMM_PROFILE_UNFILTER(cft) { amrex::BLProfiler::CommStats::UnFilter(cft); }
#define BL_COMM_PROFILE_FLUSH() { amrex::BLProfiler::WriteCommStats(true); }
#define BL_COMM_PROFILE_SETFLUSHSIZE(fsize) { amrex::BLProfiler::SetCSFlushSize(fsize); }
#define BL_COMM_PROFILE_TAGRANGE(tagmin, tagmax) {  \
        amrex::BLProfiler::TagRange(tagmin, tagmax); }
#define BL_COMM_PROFILE_TAGWRAP() { amrex::BLProfiler::AddTagWrap(); }
#define BL_COMM_PROFILE_INITAMR(flev, mlev, rr, pd) {                \
        amrex::BLProfiler::InitAMR(flev, mlev, rr, pd);              \
}

#endif


// --------------------------------------------
#elif defined(BL_TINY_PROFILING)

#include <AMReX_TinyProfiler.H>

#define BL_PROFILE_INITIALIZE()   amrex::TinyProfiler::Initialize();
#define BL_PROFILE_FINALIZE()     amrex::TinyProfiler::Finalize();
#define BL_PROFILE(fname)         amrex::TinyProfiler tiny_profiler__((fname));
#define BL_PROFILE_T(a, T)
#define BL_PROFILE_S(fname)
#define BL_PROFILE_T_S(fname, T)

#define BL_PROFILE_VAR(fname, vname)      amrex::TinyProfiler tiny_profiler__##vname((fname));
#define BL_PROFILE_VAR_NS(fname, vname)   amrex::TinyProfiler tiny_profiler__##vname(fname, false);
#define BL_PROFILE_VAR_START(vname)       tiny_profiler__##vname.start();
#define BL_PROFILE_VAR_STOP(vname)        tiny_profiler__##vname.stop();
#define BL_PROFILE_INIT_PARAMS(ptl,wall,wfabs)
#define BL_PROFILE_ADD_STEP(snum)
#define BL_PROFILE_SET_RUN_TIME(rtime)
#define BL_PROFILE_REGION_START(rname)
#define BL_PROFILE_REGION_STOP(rname)
#define BL_PROFILE_REGION_VAR(fname, rvname)
#define BL_PROFILE_REGION_VAR_START(fname, rvname)
#define BL_PROFILE_REGION_VAR_STOP(fname, rvname)
#define BL_TRACE_PROFILE_FLUSH()
#define BL_TRACE_PROFILE_SETFLUSHSIZE(fsize)
#define BL_PROFILE_CHANGE_FORT_INT_NAME(fname, intname)

#else

#include <string>
#include <AMReX_REAL.H>

namespace amrex {

class BLProfiler
{
  public:
    BLProfiler(const std::string &funcname) { }
    static void Initialize() { }
    static void Finalize() { }
    static void WriteStats(std::ostream &os) { }
    static void WriteCommStats() { }
    void start() { }
    void stop() { }
    static void InitParams(const Real ptl, const bool writeall,
                           const bool writefabs) { }
    static void AddStep(const int snum) { }
};

}

#define BL_PROFILE_INITIALIZE()
#define BL_PROFILE_FINALIZE()
#define BL_PROFILE(a)
#define BL_PROFILE_T(a, T)
#define BL_PROFILE_S(fname)
#define BL_PROFILE_T_S(fname, T)
#define BL_PROFILE_TIMER(var, a)
#define BL_PROFILE_START(var)
#define BL_PROFILE_STOP(var)
#define BL_PROFILE_THIS_NAME()

#define BL_PROFILE_VAR(fname, vname)
#define BL_PROFILE_VAR_NS(fname, vname)
#define BL_PROFILE_VAR_START(vname)
#define BL_PROFILE_VAR_STOP(vname)
#define BL_PROFILE_INIT_PARAMS(ptl,wall,wfabs)
#define BL_PROFILE_ADD_STEP(snum)
#define BL_PROFILE_SET_RUN_TIME(rtime)
#define BL_PROFILE_REGION_START(rname)
#define BL_PROFILE_REGION_STOP(rname)
#define BL_PROFILE_REGION_VAR(fname, rvname)
#define BL_PROFILE_REGION_VAR_START(fname, rvname)
#define BL_PROFILE_REGION_VAR_STOP(fname, rvname)
#define BL_TRACE_PROFILE_FLUSH()
#define BL_TRACE_PROFILE_SETFLUSHSIZE(fsize)
#define BL_PROFILE_CHANGE_FORT_INT_NAME(fname, intname)


#endif

#ifndef BL_COMM_PROFILING

#define BL_COMM_PROFILE(cft, size, pid, tag)
#define BL_COMM_PROFILE_BARRIER(message, bc)
#define BL_COMM_PROFILE_ALLREDUCE(cft, size, bc)
#define BL_COMM_PROFILE_REDUCE(cft, size, pid)
#define BL_COMM_PROFILE_WAITSOME(cft, reqs, completed, indx, status, bc)
#define BL_COMM_PROFILE_NAMETAG(message)
#define BL_COMM_PROFILE_FILTER(cft)
#define BL_COMM_PROFILE_UNFILTER(cft)
#define BL_COMM_PROFILE_FLUSH()
#define BL_COMM_PROFILE_SETFLUSHSIZE(fsize)
#define BL_COMM_PROFILE_TAGRANGE(tagmin, tagmax)
#define BL_COMM_PROFILE_TAGWRAP()
#define BL_COMM_PROFILE_INITAMR(flev, mlev, rr, pd)

#endif

#endif
