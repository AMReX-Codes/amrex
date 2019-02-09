#ifndef AMREX_BL_BACKTRACE_H_
#define AMREX_BL_BACKTRACE_H_

#include <stack>
#include <string>
#include <utility>
#include <sstream>

#include <cstdio>
#include <cstdlib>

#if defined(__linux__) && !defined(__NEC__)
#include <execinfo.h>
#endif

#include <csignal>
#include <cfenv>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef BL_BACKTRACING
#define BL_PASTE2(x, y) x##y
#define BL_PASTE(x, y) BL_PASTE2(x, y)
#define BL_BACKTRACE_PUSH( S )  amrex::BLBTer BL_PASTE( bl_bter, __COUNTER__ )( S, __FILE__, __LINE__ )
#define BL_BACKTRACE_POP()   amrex::BLBackTrace::bt_stack.pop()
#else
#define BL_BACKTRACE_PUSH(S) ((void)0)
#define BL_BACKTRACE_POP(S) ((void)0)
#endif

namespace amrex {

struct BLBackTrace
{
    static void handler(int s);
    //! Abort backtrace. Used automatically.
    static void print_backtrace_info (FILE* f);
    //! Non-abort backtrace. Prints to specified file and continues.
    static void print_backtrace_info (const std::string& filename);

#ifdef BL_BACKTRACING
    static std::stack<std::pair<std::string, std::string> > bt_stack;
// threadprivate here doesn't work with Cray and Intel
#if defined(_OPENMP) && !defined(_CRAYC) && !defined(__INTEL_COMPILER)
#pragma omp threadprivate(bt_stack)
#endif
#endif // BL_BACKTRACING
};

#ifdef BL_BACKTRACING
class BLBTer
{
public:
    BLBTer(const std::string& s, const char* file, int line);
    ~BLBTer();
private:
    std::string line_file;
    void pop_bt_stack();
};
#endif // BL_BACKTRACING

}

#endif
