#ifndef AMREX_BLFORT_H_
#define AMREX_BLFORT_H_
#include <AMReX_Config.H>

#include <AMReX_INT.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

#ifndef BL_LANG_FORT

#if   defined(BL_FORT_USE_UPPERCASE)
#define AMREX_FORT_NAME(A,B) A
#elif defined(BL_FORT_USE_LOWERCASE)
#define AMREX_FORT_NAME(A,B) B
#elif defined(BL_FORT_USE_UNDERSCORE)
#define AMREX_FORT_NAME(A,B) B##_
#elif !defined(BL_NO_FORT)
#error "One of BL_FORT_USE_{UPPERCASE,LOWERCASE,UNDERSCORE} must be defined"
#endif

#if defined(__cplusplus) && !defined(BL_NO_FORT)
#define BL_FORT_PROC_DECL(A,B) extern "C" void AMREX_FORT_NAME(A,B)
#else
#define BL_FORT_PROC_DECL(A,B) void AMREX_FORT_NAME(A,B)
#endif
#define BL_FORT_PROC_CALL(A,B) AMREX_FORT_NAME(A,B)
#define BL_FORT_PROC_NAME(A,B) AMREX_FORT_NAME(A,B)

#define BL_FORT_FAB_ARG(A) amrex_real* A##_fab, AMREX_ARLIM_P(A##_lo), AMREX_ARLIM_P(A##_hi)
#define BL_FORT_IFAB_ARG(A) int* A##_fab, AMREX_ARLIM_P(A##_lo), AMREX_ARLIM_P(A##_hi)
#define BL_FORT_FAB_VAL_3D(A) A##_fab, A##_lo, A##_hi
#define BL_FORT_FAB_ARG_3D(A) amrex_real* A##_fab, const int* A##_lo, const int* A##_hi
#define BL_FORT_FAB_ARG_ANYD(A) amrex_real* A##_fab, const int* A##_lo, const int* A##_hi
#define BL_FORT_IFAB_ARG_3D(A) int* A##_fab, const int* A##_lo, const int* A##_hi
#define BL_FORT_IFAB_ARG_ANYD(A) int* A##_fab, const int* A##_lo, const int* A##_hi

#endif

#endif
