#ifndef _MAKESLICE_F_H_
#define _MAKESLICE_F_H_

/*
** $Id: MAKESLICE_F.H,v 1.4 1998-10-07 21:18:20 vince Exp $
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
#  ifdef BL_AIX
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
    void FORT_MAKESLICE(const Real* data1,const int* dlo1,const int* dhi1,
                        Real* data2,const int* dlo2,const int* dhi2,
                        Real* data_min, Real* data_max,Real* slice_val,
                        int* slice,const Real* dx);

    void FORT_MAKEMASKSLICE(const int* data1,const int* dlo1,const int* dhi1,
                            int* data2,const int* dlo2,const int* dhi2,
                            Real* data_min, Real* data_max,Real* slice_val,
                            int* slice,const Real* dx);
#endif
};
#endif

#endif /*_MAKESLICE_F_H_*/
