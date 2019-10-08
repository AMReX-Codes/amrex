#ifndef AMREX_EB_MULTIFAB_UTIL_C_H_
#define AMREX_EB_MULTIFAB_UTIL_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>

#if (AMREX_SPACEDIM == 2)
#include <AMReX_EBMultiFabUtil_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EBMultiFabUtil_3D_C.H>
#endif

#endif
