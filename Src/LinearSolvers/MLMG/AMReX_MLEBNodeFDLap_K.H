#ifndef AMREX_MLEBNODEFDLAP_K_H_
#define AMREX_MLEBNODEFDLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

#if (AMREX_SPACEDIM == 2)
#include <AMReX_MLEBNodeFDLap_2D_K.H>
#else
#include <AMReX_MLEBNodeFDLap_3D_K.H>
#endif

#endif
