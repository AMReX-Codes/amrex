#ifndef AMREX_BL_BACKTRACE_H_
#define AMREX_BL_BACKTRACE_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <stack>
#include <string>
#include <utility>
#include <cstdlib>

#define BL_PASTE2(x, y) x##y
#define BL_PASTE(x, y) BL_PASTE2(x, y)
#define BL_BACKTRACE_PUSH( S )  amrex::BLBTer BL_PASTE( bl_bter, __COUNTER__ )( S, __FILE__, __LINE__ )
#define BL_BACKTRACE_POP()   amrex::BLBackTrace::bt_stack.pop()

namespace amrex {

struct BLBackTrace
{
    static void handler(int s);
    //! Abort backtrace. Used automatically.
    static void print_backtrace_info (FILE* f);
    //! Non-abort backtrace. Prints to specified file and continues.
    static void print_backtrace_info (const std::string& filename);

    static std::stack<std::pair<std::string, std::string> > bt_stack;
// threadprivate here doesn't work with Cray, Intel, and Fujitsu
#if defined(AMREX_USE_OMP) && !defined(_CRAYC) && !defined(__INTEL_COMPILER) && !defined(__PGI) && !defined(__NVCOMPILER) && !defined(__FUJITSU)
#pragma omp threadprivate(bt_stack)
#endif
};

class BLBTer
{
public:
    BLBTer(const std::string& s, const char* file, int line);
    BLBTer (BLBTer const&) = delete;
    BLBTer (BLBTer &&) = delete;
    BLBTer& operator= (BLBTer const&) = delete;
    BLBTer& operator= (BLBTer &&) = delete;
    ~BLBTer();
private:
    std::string line_file;
    void pop_bt_stack();
};

}

#endif
