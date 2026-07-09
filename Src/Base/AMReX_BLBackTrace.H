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
/**
 * \def BL_BACKTRACE_PUSH(S)
 * \brief Push a named scope onto the backtrace stack via an RAII guard.
 *
 * Place at the top of a function to record its entry in backtrace output.
 * The guard is automatically popped when it falls out of scope.
 *
 * \param S String literal describing the scope (e.g., function name).
 */
#define BL_BACKTRACE_PUSH( S )  amrex::BLBTer BL_PASTE( bl_bter, __COUNTER__ )( S, __FILE__, __LINE__ )
/**
 * \def BL_BACKTRACE_POP()
 * \brief Manually pop the top entry from the backtrace stack.
 *
 * Normally not needed; the RAII guard from BL_BACKTRACE_PUSH handles this.
 * Use only when manual stack control is required.
 */
#define BL_BACKTRACE_POP()   amrex::BLBackTrace::bt_stack.pop()

namespace amrex {

struct BLBackTrace
{
    /**
     * \brief Signal handler that prints and aborts on fatal signals.
     * \param s Signal number provided by the runtime.
     */
    static void handler(int s);
    /**
     * \brief Print backtrace details to \p f for abort handling.
     * \param f Destination file handle (for example, stderr).
     */
    static void print_backtrace_info (FILE* f);
    /**
     * \brief Non-aborting backtrace that prints to a file and continues.
     * \param filename Output file that receives the textual stack trace.
     */
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
     * \brief Remove this guard's entry from the backtrace stack when it is on top.
     *
     * Pops only when the stack top's scope identifier matches this guard's
     * \c line_file. Invoked from the destructor when the guard falls out of scope.
     */
    void pop_bt_stack();
};

}

#endif
