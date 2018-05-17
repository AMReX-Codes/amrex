#ifndef AMREX_COORDSYS_F_H_
#define AMREX_COORDSYS_F_H_

#ifdef BL_LANG_FORT
#    define FORT_SETVOL      setvol
#    define FORT_SETVOLPT      setvolpt
#    define FORT_SETAREA     setarea
#    define FORT_SETDLOGA    setdloga
#else
#  ifdef BL_FORT_USE_UPPERCASE
#    define FORT_SETVOL      SETVOL
#    define FORT_SETVOLPT      SETVOLPT
#    define FORT_SETAREA     SETAREA
#    define FORT_SETDLOGA    SETDLOGA
#  else
#  ifdef BL_FORT_USE_LOWERCASE
#    define FORT_SETVOL      setvol
#    define FORT_SETVOLPT      setvolpt
#    define FORT_SETAREA     setarea
#    define FORT_SETDLOGA    setdloga
#  else
#    define FORT_SETVOL      setvol_
#    define FORT_SETVOLPT      setvolpt_
#    define FORT_SETAREA     setarea_
#    define FORT_SETDLOGA    setdloga_
#  endif
#  endif

#include <AMReX_ArrayLim.H>
#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif
     void FORT_SETVOL (ARLIM_P(rlo), ARLIM_P(rhi),
		       amrex_real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                       const amrex_real* offset,
                       const amrex_real* dx, const int* coord);

     void FORT_SETVOLPT (amrex_real* vol, ARLIM_P(vollo), ARLIM_P(volhi),
			 amrex_real* rolo, ARLIM_P(rolo), ARLIM_P(rohi),
			 amrex_real* rilo, ARLIM_P(rilo), ARLIM_P(rihi),
			 const amrex_real* dx, const int* coord);

     void FORT_SETDLOGA (amrex_real* vol, ARLIM_P(vlo), ARLIM_P(vhi),
                         const amrex_real* offset,
                         const amrex_real* dx, const int* dir,
                         const int* coord);

     void FORT_SETAREA (ARLIM_P(rlo), ARLIM_P(rhi),
			amrex_real* area,  ARLIM_P(alo), ARLIM_P(ahi),
                        const amrex_real* offset,
                        const amrex_real* dx, const int* dir,
                        const int* coord);
#ifdef __cplusplus
}
#endif

#endif

#endif /*_COORDSYS_F_H_*/
