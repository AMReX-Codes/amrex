#ifndef AMREX_EB_SLOPES_K_H_
#define AMREX_EB_SLOPES_K_H_

#if (AMREX_SPACEDIM == 2)
#include <AMReX_EB_Slopes_2D_K.H>
#elif (AMREX_SPACEDIM == 3)
#include <AMReX_EB_Slopes_3D_K.H>
#endif

#endif
