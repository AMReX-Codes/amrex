#ifndef _MAKESLICE_F_H_
#define _MAKESLICE_F_H_

#ifdef BL_LANG_FORT

#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice
#    define FORT_MAKEMASKSLICE  makemaskslice
#endif
#else

#  ifdef BL_FORT_USE_UPPERCASE
#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      MAKESLICE
#    define FORT_MAKEMASKSLICE  MAKEMASKSLICE
#endif

#  else
#  ifdef BL_FORT_USE_LOWERCASE
#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice
#    define FORT_MAKEMASKSLICE  makemaskslice
#endif
#  else
#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice_
#    define FORT_MAKEMASKSLICE  makemaskslice_
#endif
#  endif
#  endif

#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C"
{
#endif
#if (BL_SPACEDIM==3)
    void FORT_MAKESLICE(const amrex::Real* data1,ARLIM_P(data1l), ARLIM_P(data1h),
                        amrex::Real* data2,ARLIM_P(data2l), ARLIM_P(data2h),
                        amrex::Real* data_min, amrex::Real* data_max,amrex::Real* slice_val,
                        int* slice,const amrex::Real* dx);

    void FORT_MAKEMASKSLICE(const int* data1,ARLIM_P(data1l), ARLIM_P(data1h),
                            int* data2,ARLIM_P(data2l), ARLIM_P(data2h),
                            amrex::Real* data_min, amrex::Real* data_max,amrex::Real* slice_val,
                            int* slice,const amrex::Real* dx);
#endif
#ifdef __cplusplus
};
#endif
#endif

#endif /*_MAKESLICE_F_H_*/
