#ifndef AMREX_BL_BACKTRACE_H_
#define AMREX_BL_BACKTRACE_H_
#include <AMReX_Config.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

/**
 * \file AMReX_BLBackTrace.H
 * \brief Utilities for capturing and printing lightweight backtraces.
 */

#include <stack>
#include <string>
#include <utility>
#include <cstdlib>

/// \cond DOXYGEN_IGNORE
#define BL_PASTE2(x, y) x##y
#define BL_PASTE(x, y) BL_PASTE2(x, y)
/// \endcond
#define BL_BACKTRACE_PUSH( S )  amrex::BLBTer BL_PASTE( bl_bter, __COUNTER__ )( S, __FILE__, __LINE__ )
#define BL_BACKTRACE_POP()   amrex::BLBackTrace::bt_stack.pop()

namespace amrex {

struct BLBackTrace
{
    /** \brief Signal handler that prints and aborts on fatal signals. */
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
    /**
     * \brief Push an entry onto the backtrace stack.
     *
     * \param s Message describing the scope being entered.
     * \param file Source file name.
     * \param line Source line number.
     */
    BLBTer(const std::string& s, const char* file, int line);
    BLBTer (BLBTer const&) = delete;
    BLBTer (BLBTer &&) = delete;
    BLBTer& operator= (BLBTer const&) = delete;
    BLBTer& operator= (BLBTer &&) = delete;
    ~BLBTer();
private:
    std::string line_file;
    /**
     * \brief Pop the most recent backtrace entry.
     *
     * Invoked automatically when the guard object falls out of scope.
     */
    void pop_bt_stack();
};

}

#endif
