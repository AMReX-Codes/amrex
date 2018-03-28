
#ifndef BL_MEMPROFILER_H
#define BL_MEMPROFILER_H

#include <stack>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace amrex {

class MemProfiler
{
public:

    struct MemInfo {
	long current_bytes;
	long hwm_bytes;
    };

    struct NBuildsInfo {
	int  current_builds;
	int  hwm_builds;
    };

    static void add (const std::string& name, std::function<MemInfo()>&& f);
    static void add (const std::string& name, std::function<NBuildsInfo()>&& f);

    static void report (const std::string& prefix = std::string());

    static void Finalize ();

private:

    MemProfiler (const MemProfiler&) = delete;
    MemProfiler& operator= (const MemProfiler&) = delete;

    MemProfiler () {}
    ~MemProfiler () {}

    void report_ (const std::string& prefix, const std::string& memory_log_name) const;

    struct Bytes {
	long mn;
	long mx;
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

    std::vector<std::string>               the_names;
    std::vector<std::function<MemInfo()> > the_funcs;

    std::vector<std::string>                   the_names_builds;
    std::vector<std::function<NBuildsInfo()> > the_funcs_builds;
};

}

#endif
