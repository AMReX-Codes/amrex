#ifndef AMREX_MLTENSOR_K_H_
#define AMREX_MLTENSOR_K_H_

#include <AMReX_FArrayBox.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLTensor_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLTensor_2D_K.H>
#else
#include <AMReX_MLTensor_3D_K.H>
#endif

#endif
