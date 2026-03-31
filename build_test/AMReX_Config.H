#ifndef AMREX_CONFIG_H_
#define AMREX_CONFIG_H_

#ifndef AMREX_SPACEDIM
#  define AMREX_SPACEDIM 3
#  define BL_SPACEDIM 3
#else
#  define BL_SPACEDIM AMREX_SPACEDIM
#endif

#if AMREX_SPACEDIM==1
#  include "AMReX_Config_1D.H"
#elif AMREX_SPACEDIM==2
#  include "AMReX_Config_2D.H"
#elif AMREX_SPACEDIM==3
#  include "AMReX_Config_3D.H"
#else
#  error "AMREX_SPACEDIM not 1, 2 or 3"
#endif

#endif /* AMREX_CONFIG_H_ */
