
#ifndef BL_MEMPROFILER_H
#define BL_MEMPROFILER_H

#include <stack>
#include <functional>
#include <string>
#include <vector>
#include <iostream>

#include <Array.H>

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

template<typename T>
long
bytesOf (const T& a)
{
    return sizeof(T);
}

template<typename T>
long
bytesOf (const std::vector<T>& v)
{
    long r = sizeof(v) + (v.capacity()-v.size())*sizeof(T);
    for (auto&& x: v)
	r += bytesOf(x);
    return r;
}

template<typename T>
long
bytesOf (const Array<T>& v)
{
    return bytesOf(static_cast<std::vector<T> >(v));
}

#endif
