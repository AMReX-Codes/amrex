#ifndef _MAKESLICE_F_H_
#define _MAKESLICE_F_H_

/*
** $Id: MAKESLICE_F.H,v 1.6 1999-04-05 17:39:21 car Exp $
*/

#ifdef BL_LANG_FORT

#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice
#    define FORT_MAKEMASKSLICE  makemaskslice
#endif
#endif

#ifdef BL_LANG_CC

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

extern "C"
{
#if (BL_SPACEDIM==3)
    void FORT_MAKESLICE(const Real* data1,ARLIM_P(data1l), ARLIM_P(data1h),
                        Real* data2,ARLIM_P(data2l), ARLIM_P(data2h),
                        Real* data_min, Real* data_max,Real* slice_val,
                        int* slice,const Real* dx);

    void FORT_MAKEMASKSLICE(const int* data1,ARLIM_P(data1l), ARLIM_P(data1h),
                            int* data2,ARLIM_P(data2l), ARLIM_P(data2h),
                            Real* data_min, Real* data_max,Real* slice_val,
                            int* slice,const Real* dx);
#endif
};
#endif

#endif /*_MAKESLICE_F_H_*/
