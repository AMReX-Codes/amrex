#ifndef AMREX_MULITIFAB_UTIL_C_H_
#define AMREX_MULITIFAB_UTIL_C_H_

#include <AMReX_MultiFabUtil_nd_C.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MultiFabUtil_1D_C.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MultiFabUtil_2D_C.H>
#else
#include <AMReX_MultiFabUtil_3D_C.H>
#endif

#endif
