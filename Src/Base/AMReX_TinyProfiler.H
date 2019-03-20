#ifndef AMREX_TINY_PROFILER_H_
#define AMREX_TINY_PROFILER_H_

#include <string>
#include <stack>
#include <map>
#include <vector>
#include <utility>
#include <limits>

#include <AMReX_REAL.H>

#ifdef AMREX_USE_CUDA
#include "nvToolsExt.h"
#endif

namespace amrex {

//! A simple profiler that returns basic performance information (e.g. min, max, and average running time)
class TinyProfiler
{
public:
    TinyProfiler (std::string funcname) noexcept;
    TinyProfiler (std::string funcname, bool start_) noexcept;
    TinyProfiler (const char* funcname) noexcept;
    TinyProfiler (const char* funcname, bool start_) noexcept;
    ~TinyProfiler ();

    void start () noexcept;
    void stop () noexcept;

    static void Initialize () noexcept;
    static void Finalize (bool bFlushing = false) noexcept;

    static void StartRegion (std::string regname) noexcept;
    static void StopRegion (const std::string& regname) noexcept;

private:
    //! stats on a single process
    struct Stats
    {
	Stats () noexcept : depth(0), n(0L), dtin(0.0), dtex(0.0) { }
	int  depth; //!< recursive depth
	long n;     //!< number of calls
	double dtin;  //!< inclusive dt
	double dtex;  //!< exclusive dt
    };

    //! stats across processes
    struct ProcStats
    {
	ProcStats () : nmin(std::numeric_limits<long>::max()),
		       navg(0L), nmax(0L),
		       dtinmin(std::numeric_limits<double>::max()),
		       dtinavg(0.0), dtinmax(0.0),
		       dtexmin(std::numeric_limits<double>::max()),
		       dtexavg(0.0), dtexmax(0.0)  {}
	long nmin, navg, nmax;
	double dtinmin, dtinavg, dtinmax;
	double dtexmin, dtexavg, dtexmax;
	std::string fname;
	static bool compex (const ProcStats& lhs, const ProcStats& rhs) {
	    return lhs.dtexmax > rhs.dtexmax;
	}
	static bool compin (const ProcStats& lhs, const ProcStats& rhs) {
	    return lhs.dtinmax > rhs.dtinmax;
	}
    };

    std::string fname;
    int global_depth;
    std::vector<Stats*> stats;

    static std::vector<std::string> regionstack;
    static std::stack<std::pair<double,double> > ttstack;
    static std::map<std::string,std::map<std::string, Stats> > statsmap;
    static double t_init;

#ifdef AMREX_USE_CUDA
    nvtxRangeId_t nvtx_id;
#endif

    static void PrintStats(std::map<std::string,Stats>& regstats, double dt_max);
};

class TinyProfileRegion
{
public:
    TinyProfileRegion (std::string a_regname) noexcept;
    TinyProfileRegion (const char* a_regname) noexcept;
    ~TinyProfileRegion ();
private:
    std::string regname;
    TinyProfiler tprof;
};

}

#endif
