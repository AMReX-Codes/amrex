#ifndef _MAKESLICE_F_H_
#define _MAKESLICE_F_H_  "%W% %G%"

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
#if (BL_SPACEDIM==3)
#    define FORT_MAKESLICE      makeslice_
#    define FORT_MAKEMASKSLICE  makemaskslice_
#endif
#  endif

  extern "C" {
#if (BL_SPACEDIM==3)
    void FORT_MAKESLICE(const REAL* data1,const int* dlo1,const int* dhi1,
     REAL* data2,const int* dlo2,const int* dhi2,REAL* data_min,
     REAL* data_max,REAL* slice_val,int* slice,const REAL* dx);

    void FORT_MAKEMASKSLICE(const int* data1,const int* dlo1,const int* dhi1,
     int* data2,const int* dlo2,const int* dhi2,REAL* data_min,
     REAL* data_max,REAL* slice_val,int* slice,const REAL* dx);
#endif
  };
#endif

#endif



