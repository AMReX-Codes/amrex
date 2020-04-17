
#ifndef BL_MEMPROFILER_H
#define BL_MEMPROFILER_H

#include <AMReX_INT.H>
#include <stack>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <memory>

namespace amrex {

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

    static void add (const std::string& name, std::function<MemInfo()>&& f);
    static void add (const std::string& name, std::function<NBuildsInfo()>&& f);

    static void report (const std::string& prefix = std::string());

    static void Finalize ();

    MemProfiler (const MemProfiler&) = delete;
    MemProfiler& operator= (const MemProfiler&) = delete;

    MemProfiler () {}
    ~MemProfiler () {}

private:

    void report_ (const std::string& prefix, const std::string& memory_log_name) const;

    struct Bytes {
        Long mn;
        Long mx;
    };
    friend std::ostream& operator<< (std::ostream& os, 
                                     const MemProfiler::Bytes& bytes);

    struct Builds {
        int mn;
        int mx;
    };
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
