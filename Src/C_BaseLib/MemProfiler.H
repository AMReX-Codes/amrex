
#ifndef BL_MEMPROFILER_H
#define BL_MEMPROFILER_H

#include <stack>
#include <functional>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include <Box.H>

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

namespace BoxLib
{
    // To avoid erroenous template deduction, we do not
    // define template<typename T> long bytesOf (const T& a);
    extern long bytesOf (int);
    extern long bytesOf (const IntVect&);
    extern long bytesOf (const Box&);
    
    template<typename T>
    long
    bytesOf (const std::vector<T>& v)
    {
	long r = sizeof(v) + (v.capacity()-v.size())*sizeof(T);
	for (auto&& x: v)
	    r += bytesOf(x);
	return r;
    }
    
    template<typename Key, typename T, class Compare, class Alloc>
    long
    bytesOf (const std::map<Key,T,Compare,Alloc>& m)
    {
	long r = sizeof(m);
	for (auto&& x: m)
	    // For glibc, there are extra 32 bytes for each mamp node.
	    // For others, this number may be different.
	    r += 32L + bytesOf(x.first) + bytesOf(x.second);
	return r;
    }
}

#endif
