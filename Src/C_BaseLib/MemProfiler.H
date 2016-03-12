
#ifndef BL_MEMPROFILER_H
#define BL_MEMPROFILER_H

#include <stack>
#include <functional>
#include <string>
#include <vector>
#include <iostream>

class MemProfiler
{
public:

    struct MemInfo {
	long current_bytes;
	long hwm_bytes;
    };

    static void add (const std::string& name, std::function<MemInfo()>&& f);

    static void report (const std::string& prefix = std::string());

private:

    MemProfiler (const MemProfiler&) = delete;
    MemProfiler& operator= (const MemProfiler&) = delete;

    MemProfiler () {}
    ~MemProfiler () {}

    void report_ (const std::string& prefix) const;

    static MemProfiler& getInstance ();

    static std::string convertUnit (long& mn, long& mx);

    std::vector<std::string>               the_names;
    std::vector<std::function<MemInfo()> > the_funcs;
};

#endif
