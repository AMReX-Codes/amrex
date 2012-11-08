#ifndef _BL_PROFILER_H_
#define _BL_PROFILER_H_

#include <iostream>
#include <string>

#ifdef BL_PROFILING

#include <REAL.H>
#include <ParallelDescriptor.H>

class Profiler 
{
public:
  class Tag
  {
  public:
    explicit Tag(const char* str) {}
  };
  Profiler(const Tag &funcname, bool iverbose)
    : bltstart(0.0), bltelapsed(0.0)
    , verbose(iverbose)
    , bRunning(false)
  {
    if(ParallelDescriptor::IOProcessor()) {
      std::cerr << "Profiler(Tag,...) interface not supported." << std::endl;
    }
  }

  Profiler(const std::string &funcname, int iverbose = 1)
    : bltstart(0.0), bltelapsed(0.0)
    , fname(funcname)
    , verbose(iverbose)
    , bRunning(false)
  {
    start();
  }

  ~Profiler() {
    if(bRunning) {
      stop();
    }
  }

  void Initialize() {
  }

  void start() {
    if(ParallelDescriptor::IOProcessor()) {
      if(verbose > 0) std::cout << "BLTStart " << fname << std::endl;
    }
    bRunning = true;
    bltstart = ParallelDescriptor::second();
  }

  void stop() {
    bltelapsed += ParallelDescriptor::second() - bltstart;
    bRunning = false;
    if(verbose > 0) {
      //ParallelDescriptor::ReduceRealMax(bltelapsed,
                                        //ParallelDescriptor::IOProcessorNumber());
      if(ParallelDescriptor::IOProcessor()) {
        std::cout << "BLTEnd " << fname << " time = " << bltelapsed << std::endl;
      }
    }
  }

  void Finalize() {
  }

  private:
    Real bltstart, bltelapsed;
    std::string fname;
    int verbose;
    bool bRunning;
};

#define BL_PROFILE(a) Profiler bl_profilerP((a));
#define BL_PROFILE_TIMER(var, a)
#define BL_PROFILE_START(var)
#define BL_PROFILE_STOP(var)
#define BL_PROFILE_THIS_NAME()

#else

class Profiler
{
public:
  class Tag
  {
  public:
    explicit Tag(const char* str) {}
  };
  Profiler(const Tag&, bool) {}
  Profiler(const std::string &funcname, int iverbose = 1) {}
  static void Initialize() {}
  static void Finalize() {}
  void start() {}
  void stop() {}
};

#define BL_PROFILE(a)
#define BL_PROFILE_TIMER(var, a)
#define BL_PROFILE_START(var)
#define BL_PROFILE_STOP(var)
#define BL_PROFILE_THIS_NAME()

#endif

#endif
