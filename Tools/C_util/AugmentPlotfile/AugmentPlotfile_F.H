
#ifndef _AUGPLT_F_H_
#define _AUGPLT_F_H_ 

//#ifdef BL_LANG_FORT
//#    define FORT_VORT_3D    vort_3d
//#    define FORT_DIVU_3D    divu_3d
//#    define FORT_COPY_3D    copy_3d
//#else
//#  if  defined(BL_FORT_USE_UPPERCASE)
//#    define FORT_VORT_3D    VORT_3D
//#    define FORT_DIVU_3D    DIVU_3D
//#    define FORT_COPY_3D    COPY_3D
//#  elif defined(BL_FORT_USE_LOWERCASE)
//#    define FORT_VORT_3D    vort_3d
//#    define FORT_DIVU_3D    divu_3d
//#    define FORT_COPY_3D    copy_3d
//#  elif defined(BL_FORT_USE_UNDERSCORE)
//#    define FORT_VORT_3D    vort_3d_
//#    define FORT_DIVU_3D    divu_3d_
//#    define FORT_COPY_3D    copy_3d_
//#  endif

#define FORT_VORT_3D     vort_3d
#define FORT_DIVU_3D     divu_3d
#define FORT_COPY_3D     copy_3d


#include <AMReX_BLFort.H>

using namespace amrex;

#ifdef __cplusplus
extern "C" {
#endif
    void FORT_VORT_3D(
            const int* lo , const int* hi,
            amrex_real* d, const int* nComps,
            const int* dlo, const int* dhi,
            const int* vel, const int* vort,
            amrex_real* delta
            );

    void FORT_DIVU_3D(
            const int* lo , const int* hi,
            amrex_real* d, const int* nComps,
            const int* dlo, const int* dhi,
            const int* vel, const int* divu,
            amrex_real* delta
            );
            
    void FORT_COPY_3D(
            const int* srclo , const int* srchi,
            const amrex_real* srcd, const int* srcNComps,
            const int* srcdlo, const int* srcdhi,
            const int* dstlo , const int* dsthi,
                  amrex_real* dstd, const int* dstNComps,
            const int* dstdlo, const int* dstdhi,
            const int* srccomp, const int* dstcomp
            );
#ifdef __cplusplus
}
#endif

#endif /*_AUGPLT_F_H_*/
