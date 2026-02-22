
#ifndef BL_MEMPROFILER_H
#define BL_MEMPROFILER_H
#include <AMReX_Config.H>

#include <AMReX_INT.H>

#include <stack>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <iosfwd>
#include <memory>

/**
 * \file AMReX_MemProfiler.H
 * \brief Lightweight memory usage profiler and reporting hooks.
 */

namespace amrex {

/**
 * \class amrex::MemProfiler
 * \brief Singleton that aggregates memory usage callbacks and prints reports.
 */
class MemProfiler
{
public:

    struct MemInfo {
        Long current_bytes;
        Long hwm_bytes;
    };

    struct NBuildsInfo {
        int  current_builds;
        int  hwm_builds;
    };

    //! Register a callback returning current bytes / HWM stats.
    static void add (const std::string& name, std::function<MemInfo()>&& f);
    //! Register a callback returning current/high-water build counts.
    static void add (const std::string& name, std::function<NBuildsInfo()>&& f);

    //! Emit a memory report to stdout (and optionally to a log file via prefix).
    static void report (const std::string& prefix = std::string());

    //! Destroy the global profiler instance and release callbacks.
    static void Finalize ();

    MemProfiler (const MemProfiler&) = delete;
    MemProfiler& operator= (const MemProfiler&) = delete;

    MemProfiler () {}
    ~MemProfiler () {}

private:

    void report_ (const std::string& prefix, const std::string& memory_log_name) const;

    /// \cond DOXYGEN_IGNORE
    struct Bytes {
        Long mn;
        Long mx;
    };
    /// \endcond
    friend std::ostream& operator<< (std::ostream& os,
                                     const MemProfiler::Bytes& bytes);

    /// \cond DOXYGEN_IGNORE
    struct Builds {
        int mn;
        int mx;
    };
    /// \endcond
    friend std::ostream& operator<< (std::ostream& os,
                                     const MemProfiler::Builds& builds);

    static MemProfiler& getInstance ();

    static std::unique_ptr<MemProfiler> the_instance;

    std::vector<std::string>               the_names;
    std::vector<std::function<MemInfo()> > the_funcs;

    std::vector<std::string>                   the_names_builds;
    std::vector<std::function<NBuildsInfo()> > the_funcs_builds;
};

}

#endif
