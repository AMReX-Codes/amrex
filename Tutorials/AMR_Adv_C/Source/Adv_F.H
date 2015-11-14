#ifndef _Adv_F_H_
#define _Adv_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(GET_TAGGING_PARAMS,get_tagging_params)(const int* name, const int* namlen);

BL_FORT_PROC_DECL(NULLFILL,nullfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

BL_FORT_PROC_DECL(ESTDT,estdt)
    (const int* lo, const int* hi, 
     const BL_FORT_FAB_ARG_3D(state),
     const Real* dx, Real* dt);

BL_FORT_PROC_DECL(STATE_ERROR, state_error)
    (int* tag, const int* tag_lo, const int* tag_hi,
     BL_FORT_FAB_ARG_3D(state),
     const int* tagval, const int* clearval,
     const int* lo, const int* hi,
     const Real* dx, const Real* problo,
     const Real* time, const int* level);

BL_FORT_PROC_DECL(INITDATA,initdata)
    (const int& level, const Real& time, 
     const int* lo, const int* hi,
     BL_FORT_FAB_ARG_3D(state),
     const Real* dx, const Real* xlo, const Real* xhi);

#endif
