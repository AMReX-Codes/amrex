#ifndef AMREX_TINY_PROFILER_H_
#define AMREX_TINY_PROFILER_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_REAL.H>

#ifdef AMREX_USE_CUDA
#include <nvToolsExt.h>
#endif

#if defined(AMREX_USE_HIP) && defined(AMREX_USE_ROCTX)
#include <roctracer/roctx.h>
#endif

#include <deque>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace amrex {

class CArena;

struct MemStat
{
    enum memory_type {
        device = 0,
        managed = 1,
        pinned = 2,
        cpu = 3,
    };

    static constexpr std::array<const char *, 4> memory_name = {
        "Device", "Managed", "Pinned", "Cpu"
    };

    Long nalloc = 0;
    Long nfree = 0;
    Long currentmem = 0;
    Long maxmem = 0;
};

//! A simple profiler that returns basic performance information (e.g. min, max, and average running time)
class TinyProfiler
{
public:
    explicit TinyProfiler (std::string funcname) noexcept;
    TinyProfiler (std::string funcname, bool start_, bool useCUPTI=false) noexcept;
    explicit TinyProfiler (const char* funcname) noexcept;
    TinyProfiler (const char* funcname, bool start_, bool useCUPTI=false) noexcept;
    ~TinyProfiler ();

    void start () noexcept;
    void stop () noexcept;
#ifdef AMREX_USE_CUPTI
    void stop (unsigned boxUintID) noexcept;
#endif // AMREX_USE_CUPTI


    void memory_start () const noexcept;
    void memory_stop () const noexcept;

    static MemStat* memory_alloc (std::size_t nbytes, int memtype) noexcept;
    static void memory_free (std::size_t nbytes, MemStat* stat) noexcept;

    static void Initialize () noexcept;
    static void Finalize (bool bFlushing = false) noexcept;

    static void StartRegion (std::string regname) noexcept;
    static void StopRegion (const std::string& regname) noexcept;

    static void PrintCallStack (std::ostream& os);

private:
    struct Stats
    {
        Stats () noexcept : depth(0), n(0L), dtin(0.0), dtex(0.0),
                            usesCUPTI(false), nk(0) { }
        int  depth;     //!< recursive depth
        Long n;         //!< number of calls
        double dtin;    //!< inclusive dt
        double dtex;    //!< exclusive dt
        bool usesCUPTI; //!< uses CUPTI
        Long nk;        //!< number of kernel calls
    };

    //! stats across processes
    struct ProcStats
    {
        ProcStats () : nmin(std::numeric_limits<Long>::max()),
                       navg(0L), nmax(0L),
                       dtinmin(std::numeric_limits<double>::max()),
                       dtinavg(0.0), dtinmax(0.0),
                       dtexmin(std::numeric_limits<double>::max()),
                       dtexavg(0.0), dtexmax(0.0),
                       usesCUPTI(false) {}
        Long nmin, navg, nmax;
        double dtinmin, dtinavg, dtinmax;
        double dtexmin, dtexavg, dtexmax;
        bool usesCUPTI;
        std::string fname;
        static bool compex (const ProcStats& lhs, const ProcStats& rhs) {
            return lhs.dtexmax > rhs.dtexmax;
        }
        static bool compin (const ProcStats& lhs, const ProcStats& rhs) {
            return lhs.dtinmax > rhs.dtinmax;
        }
    };

    struct MemProcStats
    {
        Long nalloc = 0;
        Long nfree = 0;
        Long maxmem_min = std::numeric_limits<Long>::max();
        Long maxmem_avg = 0;
        Long maxmem_max = 0;
        std::string fname;
        static bool compmem (const MemProcStats& lhs, const MemProcStats& rhs) {
            return lhs.maxmem_max > rhs.maxmem_max;
        }
    };

    std::string fname;
    bool uCUPTI;
    bool in_parallel_region;
    int global_depth;
    std::vector<Stats*> stats;


    static std::deque<const TinyProfiler*> mem_stack;

#ifdef AMREX_USE_OMP
    struct aligned_deque {
        // align thread-private data to cache lines to prevent false sharing
        alignas(64) std::deque<const TinyProfiler*> deque;
    };

    static std::vector<aligned_deque> mem_stack_thread_private;
#endif
    static std::array<std::map<std::string, MemStat>, 4> mem_statsmap;

    static std::vector<std::string> regionstack;
    static std::deque<std::tuple<double,double,std::string*> > ttstack;
    static std::map<std::string,std::map<std::string, Stats> > statsmap;
    static double t_init;
    static int device_synchronize_around_region;
    static int n_print_tabs;
    static int verbose;
    static void PrintStats (std::map<std::string,Stats>& regstats, double dt_max);
    static void PrintMemStats (int mem_type);
};

class TinyProfileRegion
{
public:
    explicit TinyProfileRegion (std::string a_regname) noexcept;
    explicit TinyProfileRegion (const char* a_regname) noexcept;
    ~TinyProfileRegion ();
private:
    std::string regname;
    TinyProfiler tprof;
};

}
#endif
