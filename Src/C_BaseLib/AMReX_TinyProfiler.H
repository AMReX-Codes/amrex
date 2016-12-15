#ifndef _TINY_PROFILER_H_
#define _TINY_PROFILER_H_

#include <string>
#include <stack>
#include <map>
#include <utility>
#include <limits>

#include <REAL.H>

class TinyProfiler
{
public:
    TinyProfiler (const std::string &funcname);
    TinyProfiler (const std::string &funcname, bool start_);
    ~TinyProfiler ();

    void start ();
    void stop ();

    static void Initialize ();
    static void Finalize ();

private:
    struct Stats   // stats on a single process
    {
	Stats () : depth(0), n(0L), dtin(0.0), dtex(0.0) { } 
	int  depth; // recursive depth
	long n;     // number of calls
	Real dtin;  // inclusive dt
	Real dtex;  // exclusive dt
    };

    struct ProcStats // stats across processes
    {
	ProcStats () : nmin(std::numeric_limits<long>::max()), 
		       navg(0L), nmax(0L), 
		       dtinmin(std::numeric_limits<Real>::max()), 
		       dtinavg(0.0), dtinmax(0.0),
		       dtexmin(std::numeric_limits<Real>::max()), 
		       dtexavg(0.0), dtexmax(0.0)  {}
	long nmin, navg, nmax;
	Real dtinmin, dtinavg, dtinmax;
	Real dtexmin, dtexavg, dtexmax;
	std::string fname;
	static bool compex (const ProcStats& lhs, const ProcStats& rhs) {
	    return lhs.dtexmax > rhs.dtexmax;
	}
	static bool compin (const ProcStats& lhs, const ProcStats& rhs) {
	    return lhs.dtinmax > rhs.dtinmax;
	}
    };

    std::string fname;
    bool running;
    int global_depth;

    static std::stack<std::pair<Real,Real> > ttstack;
    static std::map<std::string, Stats> statsmap;
    static Real t_init;
};

#endif
