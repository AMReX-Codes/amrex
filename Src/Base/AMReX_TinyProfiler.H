#ifndef AMREX_TINY_PROFILER_H_
#define AMREX_TINY_PROFILER_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_TinyProfiler.H
 * \brief Lightweight RAII profiler that records inclusive and exclusive timing.
 */

#include <AMReX_INT.H>
#include <AMReX_REAL.H>

#include <array>
#include <deque>
#include <iosfwd>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace amrex {

/**
 * \struct amrex::MemStat
 * \ingroup amrex_utilities
 * \brief Per-region memory accounting counters updated by TinyProfiler.
 */
struct MemStat
{
    Long nalloc = 0;        //!< number of allocations
    Long nfree = 0;         //!< number of frees
    Long currentmem = 0;    //!< amount of currently used memory in bytes
    double avgmem = 0.;     //!< memory used (bytes) times time in use (seconds)
    Long maxmem = 0;        //!< running maximum of currentmem
};

/**
 * \class amrex::TinyProfiler
 * \ingroup amrex_utilities
 * \brief Lightweight RAII profiler that records per-function timing statistics.
 *
 * Construct a TinyProfiler at the entry of a region to be timed;
 * it automatically starts the timer and stops it when destroyed.
 * Timing is accumulated across calls and reported at finalization.
 */
class TinyProfiler
{
public:
    //! Construct and immediately start timing the region named \p funcname.
    explicit TinyProfiler (std::string funcname) noexcept;
    //! Construct; start timing if \p start_ is \c true.
    TinyProfiler (std::string funcname, bool start_) noexcept;
    //! Construct and immediately start timing the region named \p funcname.
    explicit TinyProfiler (const char* funcname) noexcept;
    //! Construct; start timing if \p start_ is \c true.
    TinyProfiler (const char* funcname, bool start_) noexcept;
    ~TinyProfiler ();

    TinyProfiler (TinyProfiler const&) = delete;
    TinyProfiler (TinyProfiler &&) = delete;
    TinyProfiler& operator= (TinyProfiler const&) = delete;
    TinyProfiler& operator= (TinyProfiler &&) = delete;

    //! Start (or resume) the timer for this region.
    void start () noexcept;
    //! Stop (or pause) the timer for this region.
    void stop () noexcept;

    //! Record a memory-tracking start event for this region.
    void memory_start () const noexcept;
    //! Record a memory-tracking stop event for this region.
    void memory_stop () const noexcept;

    //! Record an allocation of \p nbytes to the named memory-stat entry; returns the stat pointer.
    static MemStat* memory_alloc (std::size_t nbytes,
                                  std::map<std::string, MemStat>& memstats) noexcept;
    //! Record a deallocation of \p nbytes from the stat entry \p stat.
    static void memory_free (std::size_t nbytes, MemStat* stat) noexcept;

    //! Initialize the TinyProfiler subsystem.
    static void Initialize ();
    //! Finalize the TinyProfiler subsystem and print a timing report.
    static void Finalize (bool bFlushing = false) noexcept;

    //! Initialize memory-tracking subsystem.
    static void MemoryInitialize ();
    //! Finalize memory-tracking subsystem and print a memory report.
    static void MemoryFinalize (bool bFlushing = false) noexcept;

    //! Register an arena's memory-stat map under \p memory_name; returns \c true on success.
    static bool RegisterArena (const std::string& memory_name,
                               std::map<std::string, MemStat>& memstats) noexcept;

    //! Deregister the arena identified by \p memstats.
    static void DeregisterArena (std::map<std::string, MemStat>& memstats) noexcept;

    //! Open a named timing region \p regname.
    static void StartRegion (std::string regname) noexcept;
    //! Close the named timing region \p regname.
    static void StopRegion (const std::string& regname) noexcept;

    //! Print the current call stack to \p os.
    static void PrintCallStack (std::ostream& os);
    //! Print memory usage statistics to \p os.
    static void PrintMemoryUsage (std::ostream* os, bool only_local) noexcept;

private:
    struct Stats
    {
        Stats () noexcept  = default;
        int  depth{0};     //!< recursive depth
        Long n{0L};         //!< number of calls
        double dtin{0.0};    //!< inclusive dt
        double dtex{0.0};    //!< exclusive dt
    };

    //! stats across processes
    struct ProcStats
    {

        Long nmin{std::numeric_limits<Long>::max()};
        Long navg{0L}, nmax{0L};
        double dtinmin{std::numeric_limits<double>::max()};
        double dtinavg{0.0}, dtinmax{0.0};
        double dtexmin{std::numeric_limits<double>::max()};
        double dtexavg{0.0}, dtexmax{0.0};
        bool do_print{true};
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
        Long curmem_max = 0;
        Long avgmem_min = std::numeric_limits<Long>::max();
        Long avgmem_avg = 0;
        Long avgmem_max = 0;
        Long maxmem_min = std::numeric_limits<Long>::max();
        Long maxmem_avg = 0;
        Long maxmem_max = 0;
        std::string fname;
        static bool compmem (const MemProcStats& lhs, const MemProcStats& rhs) {
            return lhs.maxmem_max > rhs.maxmem_max;
        }
    };

    std::string fname;
    bool in_parallel_region = false;
    int global_depth = -1;
    std::vector<Stats*> stats;

    static std::deque<const TinyProfiler*> mem_stack;

#ifdef AMREX_USE_OMP
    struct aligned_deque {
        // align thread-private data to cache lines to prevent false sharing
        alignas(64) std::deque<const TinyProfiler*> deque;
    };

    static std::vector<aligned_deque> mem_stack_thread_private;
#endif
    static std::vector<std::map<std::string, MemStat>*> all_memstats;
    static std::vector<std::string> all_memnames;

    static std::vector<std::string> regionstack;
    static std::deque<std::tuple<double,double,std::string*> > ttstack;
    static std::map<std::string,std::map<std::string, Stats> > statsmap;
    static double t_init;
    static double t_memory_init;
    static bool device_synchronize_around_region;
    static int n_print_tabs;
    static int verbose;
    static double print_threshold;
    static bool enabled;
    static bool memprof_enabled;
    static std::string output_file;

    static std::string const& get_output_file ();
    static void PrintStats (std::map<std::string,Stats>& regstats, double dt_max,
                            std::ostream* os);
    static void PrintMemStats (std::map<std::string, MemStat>& memstats,
                               std::string const& memname, double dt_max,
                               double t_final, std::ostream* os, bool only_local);
};

class TinyProfileRegion
{
public:
    explicit TinyProfileRegion (std::string a_regname) noexcept;
    explicit TinyProfileRegion (const char* a_regname) noexcept;
    TinyProfileRegion (TinyProfileRegion const&) = delete;
    TinyProfileRegion (TinyProfileRegion &&) = delete;
    TinyProfileRegion& operator= (TinyProfileRegion const&) = delete;
    TinyProfileRegion& operator= (TinyProfileRegion &&) = delete;
    ~TinyProfileRegion ();
private:
    std::string regname;
    TinyProfiler tprof;
};

}
#endif
