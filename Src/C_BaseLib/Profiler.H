#ifndef _BL_PROFILER_H_
#define _BL_PROFILER_H_


#ifdef BL_PROFILING

#include <REAL.H>

#include <iostream>
#include <string>
#include <stack>
#include <set>
#include <map>
#include <vector>
#include <utility>
#include <typeinfo>


class Profiler 
{
  public:
    struct ProfStats {
      ProfStats() : nCalls(0), totalTime(0.0), minTime(0.0),
                    maxTime(0.0), avgTime(0.0) { }
      long nCalls;
      Real totalTime, minTime, maxTime, avgTime;
    };

    enum CommFuncType {
      InvalidCFT = 0,
      AllReduceT,
      AllReduceR,
      AllReduceL,
      AllReduceI,
      AsendTsii,
      AsendTsiiM,
      AsendvTii,
      SendTsii,
      SendvTii,
      ArecvTsii,
      ArecvTsiiM,
      ArecvTii,
      ArecvvTii,
      RecvTsii,
      RecvvTii,
      ReduceT,
      ReduceR,
      ReduceL,
      ReduceI,
      BCastTsi,
      GatherTsT1Si,
      GatherTi,
      ScatterTsT1si,
      Barrier,
      NameTag,
      AllCFTypes,
      NoCFTypes
    };
                 
    struct CommStats {
      CommStats() : cfType(InvalidCFT), commpid(-1), size(-1), tag(-1),
                    timeStamp(-1.0) { }
      CommStats(CommFuncType cft, int commpid, int size, Real ts)
                  : cfType(cft), commpid(commpid), size(size), tag(-1),
                    timeStamp(ts) { }
      CommStats(CommFuncType cft, int commpid, int size, int tag, Real ts)
                  : cfType(cft), commpid(commpid), size(size), tag(tag),
                    timeStamp(ts) { }

      static std::string CFTToString(CommFuncType cft);
      static CommFuncType StringToCFT(const std::string &s);
      static void Filter(CommFuncType cft);
      static void UnFilter(CommFuncType cft);

      static std::map<std::string, CommFuncType> cftNames;
      static std::set<CommFuncType> cftExclude;
      static int barrierNumber;
      static int reductionNumber;
      static std::vector<std::pair<std::string,int> > barrierNames;  // [name,index]
      static std::vector<std::pair<std::string,int> > nameTags;      // [name,index]
      static std::vector<int> reductions;      // [index]

      CommFuncType cfType;
      int commpid, size, tag;
      Real timeStamp;
    };


    explicit Profiler(const std::string &funcname);

    ~Profiler();

    static void Initialize();
    static void Finalize();
    static void WriteStats(std::ostream &os, bool bwriteavg = false);
    static void WriteCommStats();

    static void AddCommStat(const CommFuncType cft, const int pid, const int size);
    static void AddCommStat(const CommFuncType cft, const int pid, const int size, const int tag);
    static void AddBarrier(const std::string &message, const bool bc);
    static void AddNameTag(const std::string &name);
    static void AddAllReduce(const CommFuncType cft, const int size, const bool bc);

    void start();
    void stop();

    static void InitParams(const Real ptl, const bool writeall,
                           const bool writefabs, const bool writeblt);
    static void AddStep(const int snum);
    static void SetRunTime(const Real rtime) { calcRunTime = rtime; }

  private:
    Real bltstart, bltelapsed;
    std::string fname;
    bool bRunning;

    static bool bWriteAll, bWriteFabs, bWriteBLT;
    static bool bFirstCommWriteH, bFirstCommWriteD;
    static bool bInitialized;
    static int currentStep, csFlushSize, nProfFiles;
    static Real pctTimeLimit;
    static Real calcRunTime;
    static Real startTime;
    static std::stack<Real> nestedTimeStack;
    static std::map<int, Real> mStepMap;  // [step, time]
    static std::map<std::string, ProfStats> mProfStats;  // [fname, pstats]
    static std::map<Real, std::string, std::greater<Real> > mTimersTotalsSorted;
    static std::vector<CommStats> vCommStats;

    static void WriteHeader(std::ostream &os, const int colWidth,
                            const Real maxlen, const bool bwriteavg);
    static void WriteRow(std::ostream &os, const std::string &fname,
                         const ProfStats &pstats, const Real percent,
			 const int colWidth, const Real maxlen,
			 const bool bwriteavg);

    static bool OnExcludeList(CommFuncType cft);
};



inline std::string Profiler::CommStats::CFTToString(CommFuncType cft) {
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
    case ScatterTsT1si:  return "ScatterTsT1si";
    case Barrier:        return "Barrier";
    case NameTag:        return "NameTag";
    case AllCFTypes:     return "AllCFTypes";
    case NoCFTypes:      return "NoCFTypes";
  }
  return "*** Error: Bad CommFuncType.";
}



#define BL_PROFILE_INITIALIZE()  Profiler::Initialize();
#define BL_PROFILE(fname) Profiler bl_profiler__((fname));
#define BL_PROFILE_T(fname, T) {                                    \
		     std::string pname(fname);                      \
		     pname += typeid(T).name();                     \
                     Profiler bl_profiler__((pname));               \
}
#define BL_PROFILE_VAR(fname, vname) Profiler bl_profiler__##vname((fname));
#define BL_PROFILE_VAR_START(vname) bl_profiler__##vname.start();
#define BL_PROFILE_VAR_STOP(vname) bl_profiler__##vname.stop();
#define BL_PROFILE_INIT_PARAMS(ptl,wall,wfabs,wblt)  \
                                  Profiler::InitParams(ptl,wall, wfabs,wblt);
#define BL_PROFILE_ADD_STEP(snum)  Profiler::AddStep(snum);
#define BL_PROFILE_SET_RUN_TIME(rtime)  Profiler::SetRunTime(rtime);
#define BL_PROFILE_FINALIZE()  Profiler::Finalize();

#define BL_COMM_PROFILE_T(cft, T, n, pid) {                                \
                          Profiler::AddCommStat(cft, pid, n * sizeof(T));  \
}
#define BL_COMM_PROFILE_TTAG(cft, T, n, pid, tag) {                             \
                          Profiler::AddCommStat(cft, pid, n * sizeof(T), tag);  \
}
#define BL_COMM_PROFILE_BARRIER(message, bc) { Profiler::AddBarrier(message, bc); }
#define BL_COMM_PROFILE_ALLREDUCE(cft, size, bc) { Profiler::AddAllReduce(cft, size, bc); }
#define BL_COMM_PROFILE_REDUCE(cft, size, pid) {   \
                          Profiler::AddCommStat(cft, pid, size); }
#define BL_COMM_PROFILE_NAMETAG(message) { Profiler::AddNameTag(message); }
#define BL_COMM_PROFILE_FILTER(cft) { Profiler::CommStats::Filter(cft); }
#define BL_COMM_PROFILE_UNFILTER(cft) { Profiler::CommStats::UnFilter(cft); }

#else

#include <string>
#include <REAL.H>

class Profiler
{
  public:
    Profiler(const std::string &funcname) { }
    static void Initialize() { }
    static void Finalize() { }
    static void WriteStats(std::ostream &os) { }
    static void WriteCommStats() { }
    void start() { }
    void stop() { }
    static void InitParams(const Real ptl, const bool writeall,
                           const bool writefabs, const bool writeblt) { }
    static void AddStep(const int snum) { }
};

#define BL_PROFILE_INITIALIZE()
#define BL_PROFILE(a)
#define BL_PROFILE_T(a, T)
#define BL_PROFILE_TIMER(var, a)
#define BL_PROFILE_START(var)
#define BL_PROFILE_STOP(var)
#define BL_PROFILE_THIS_NAME()

#define BL_PROFILE_VAR(fname, vname)
#define BL_PROFILE_VAR_START(vname)
#define BL_PROFILE_VAR_STOP(vname)
#define BL_PROFILE_INIT_PARAMS(ptl,wall,wfabs,wblt)
#define BL_PROFILE_ADD_STEP(snum)
#define BL_PROFILE_SET_RUN_TIME(rtime)
#define BL_PROFILE_FINALIZE()

#define BL_COMM_PROFILE_T(cft, T, n, pid)
#define BL_COMM_PROFILE_TTAG(cft, T, n, pid, tag)
#define BL_COMM_PROFILE_BARRIER(message, bc)
#define BL_COMM_PROFILE_ALLREDUCE(cft, size, bc)
#define BL_COMM_PROFILE_REDUCE(cft, size, pid)
#define BL_COMM_PROFILE_NAMETAG(message)
#define BL_COMM_PROFILE_FILTER(cft)
#define BL_COMM_PROFILE_UNFILTER(cft)

#endif

#endif
