#ifndef AMREX_EB2_C_H_
#define AMREX_EB2_C_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EB2_Graph.H>

#if (AMREX_SPACEDIM == 2)
#include <AMReX_EB2_2D_C.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EB2_3D_C.H>
#endif

#endif
