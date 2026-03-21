#ifndef AMREX_BLFORT_H_
#define AMREX_BLFORT_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

/**
 * \file AMReX_BLFort.H
 * \brief Macros for interoperating between AMReX C++ and legacy Fortran kernels.
 */

#ifndef BL_LANG_FORT

/**
 * \def AMREX_FORT_NAME(A,B)
 * \brief Helper that expands to the correct Fortran symbol decoration.
 *
 * When Fortran linking is enabled (i.e., `BL_NO_FORT` is not defined), exactly
 * one of `BL_FORT_USE_UPPERCASE`, `BL_FORT_USE_LOWERCASE`, or
 * `BL_FORT_USE_UNDERSCORE` must be defined at build time (see `AMReX_Utility.H`).
 * The macro expands to `A` for uppercase linking, `B` for lowercase linking, and
 * `B_` when a trailing underscore is required.  When `BL_NO_FORT` is defined,
 * none of the `BL_FORT_USE_*` macros are needed.
 *
 * \param A Uppercase variant of the symbol name.
 * \param B Lowercase variant of the symbol name.
 */
#if   defined(BL_FORT_USE_UPPERCASE)
#define AMREX_FORT_NAME(A,B) A
#elif defined(BL_FORT_USE_LOWERCASE)
#define AMREX_FORT_NAME(A,B) B
#elif defined(BL_FORT_USE_UNDERSCORE)
#define AMREX_FORT_NAME(A,B) B##_
#elif !defined(BL_NO_FORT)
#error "One of BL_FORT_USE_{UPPERCASE,LOWERCASE,UNDERSCORE} must be defined"
#endif

/**
 * \def BL_FORT_PROC_DECL(A,B)
 * \brief Declare a Fortran-callable routine using the configured symbol naming.
 * \param A Uppercase variant of the symbol name.
 * \param B Lowercase variant of the symbol name.
 */
#if defined(__cplusplus) && !defined(BL_NO_FORT)
#define BL_FORT_PROC_DECL(A,B) extern "C" void AMREX_FORT_NAME(A,B)
#else
#define BL_FORT_PROC_DECL(A,B) void AMREX_FORT_NAME(A,B)
#endif
/**
 * \def BL_FORT_PROC_CALL(A,B)
 * \brief Invoke a Fortran routine using the configured symbol decoration.
 */
#define BL_FORT_PROC_CALL(A,B) AMREX_FORT_NAME(A,B)
/**
 * \def BL_FORT_PROC_NAME(A,B)
 * \brief Alias that yields the decorated Fortran symbol for use in arrays or tables.
 */
#define BL_FORT_PROC_NAME(A,B) AMREX_FORT_NAME(A,B)

/**
 * \def BL_FORT_FAB_ARG(A)
 * \brief Expand to the canonical `(data, lo, hi)` argument list for an FArrayBox.
 * \param A Base name used to form the argument identifiers.
 */
#define BL_FORT_FAB_ARG(A) amrex_real* A##_fab, AMREX_ARLIM_P(A##_lo), AMREX_ARLIM_P(A##_hi)
/**
 * \def BL_FORT_IFAB_ARG(A)
 * \brief Expand to the canonical `(data, lo, hi)` argument list for an IArrayBox.
 */
#define BL_FORT_IFAB_ARG(A) int* A##_fab, AMREX_ARLIM_P(A##_lo), AMREX_ARLIM_P(A##_hi)
/**
 * \def BL_FORT_FAB_VAL_3D(A)
 * \brief Expand to `(data, lo, hi)` using raw pointers for explicit 3D kernels.
 */
#define BL_FORT_FAB_VAL_3D(A) A##_fab, A##_lo, A##_hi
/// \brief Expand to `(data, lo, hi)` with explicit pointers (3D).
#define BL_FORT_FAB_ARG_3D(A) amrex_real* A##_fab, const int* A##_lo, const int* A##_hi
/// \brief Expand to `(data, lo, hi)` arguments using ANYD bounds pointers.
#define BL_FORT_FAB_ARG_ANYD(A) amrex_real* A##_fab, const int* A##_lo, const int* A##_hi
/// \brief Expand to `(data, lo, hi)` arguments for integer fabs in 3D.
#define BL_FORT_IFAB_ARG_3D(A) int* A##_fab, const int* A##_lo, const int* A##_hi
/// \brief Expand to `(data, lo, hi)` arguments for integer fabs using ANYD bounds pointers.
#define BL_FORT_IFAB_ARG_ANYD(A) int* A##_fab, const int* A##_lo, const int* A##_hi

#endif

#endif
