#ifndef _Adv_F_H_
#define _Adv_F_H_
#include <BLFort.H>

BL_FORT_PROC_DECL(GET_TAGGING_PARAMS,get_tagging_params)(const int* name, const int* namlen);

BL_FORT_PROC_DECL(NULLFILL,nullfill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);

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
     const Real* dx, const Real* problo);

BL_FORT_PROC_DECL(GET_FACE_VELOCITY,get_face_velocity)
    (const int& level, const Real& time, 
     D_DECL(BL_FORT_FAB_ARG(xvel),
            BL_FORT_FAB_ARG(yvel),
            BL_FORT_FAB_ARG(zvel)),
     const Real* dx, const Real* problo);

BL_FORT_PROC_DECL(ADVECT,advect)
    (const Real& time, const int* lo, const int*hi,
     const BL_FORT_FAB_ARG_3D(statein),
     BL_FORT_FAB_ARG_3D(stateout),
     D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
            const BL_FORT_FAB_ARG_3D(yvel),
            const BL_FORT_FAB_ARG_3D(zvel)),
     D_DECL(BL_FORT_FAB_ARG_3D(fx),
            BL_FORT_FAB_ARG_3D(fy),
            BL_FORT_FAB_ARG_3D(fz)),
     const Real* dx, const Real& dt);

#endif
