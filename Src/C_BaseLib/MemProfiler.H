
#ifndef BL_MEMPROFILER_H
#define BL_MEMPROFILER_H

#include <stack>
#include <functional>
#include <string>
#include <vector>
#include <map>
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

    // For gcc, there are extra 32 bytes for each map node.
    // For others, this number may be different.
    static const long gcc_map_node_extra_bytes = 32L;

private:

    MemProfiler (const MemProfiler&) = delete;
    MemProfiler& operator= (const MemProfiler&) = delete;

    MemProfiler () {}
    ~MemProfiler () {}

    void report_ (const std::string& prefix) const;

    struct Bytes {
	long mn;
	long mx;
    };
    friend std::ostream& operator<< (std::ostream& os, 
				     const MemProfiler::Bytes& bytes);

    static MemProfiler& getInstance ();

    std::vector<std::string>               the_names;
    std::vector<std::function<MemInfo()> > the_funcs;
};

#endif
