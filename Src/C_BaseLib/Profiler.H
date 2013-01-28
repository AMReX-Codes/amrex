#ifndef _BL_PROFILER_H_
#define _BL_PROFILER_H_


#ifdef BL_PROFILING

#include <REAL.H>

#include <iostream>
#include <string>
#include <stack>
#include <map>
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

    explicit Profiler(const std::string &funcname);

    ~Profiler();

    static void Initialize();
    static void Finalize();
    static void WriteStats(std::ostream &os, bool bwriteavg = false);

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
    static int currentStep;
    static Real pctTimeLimit;
    static Real calcRunTime;
    static std::stack<Real> nestedTimeStack;
    static std::map<int, Real> mStepMap;  // [step, time]
    static std::map<std::string, ProfStats> mProfStats;  // [fname, pstats]
    static std::map<Real, std::string, std::greater<Real> > mTimersTotalsSorted;

    static void WriteHeader(std::ostream &os, const int colWidth,
                            const Real maxlen, const bool bwriteavg);
    static void WriteRow(std::ostream &os, const std::string &fname,
                         const ProfStats &pstats, const Real percent,
			 const int colWidth, const Real maxlen,
			 const bool bwriteavg);
};



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
    void start() { }
    void stop() { }
    static void InitParams(const Real ptl, const bool writeall,
                           const bool writefabs, const bool writeblt) { }
    static void AddStep(const int snum) { }
};

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


#endif

#endif
