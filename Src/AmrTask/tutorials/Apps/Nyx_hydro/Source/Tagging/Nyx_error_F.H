#ifndef _Nyx_error_F_H_
#define _Nyx_error_F_H_
#include <AMReX_BLFort.H>

BL_FORT_PROC_DECL(TAG_LAPLAC_ERROR, tag_laplac_error)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(var),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_DENERROR, tag_denerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_VELERROR, tag_velerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(vel),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_TEMPERROR, tag_temperror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(temp),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_PRESSERROR, tag_presserror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_PART_CNT_ERR, tag_part_cnt_err)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(ls),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_MFRACERROR, tag_mfracerror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(press),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_ENTERROR, tag_enterror)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(den),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_STATE_ERROR, tag_state_error)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const amrex::Real* time, const int* level);

BL_FORT_PROC_DECL(TAG_OVERDENSITY, tag_overdensity)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, const int* clearval,
     BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     const int* ncomp,
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const int* level,
     const amrex::Real* avg);

BL_FORT_PROC_DECL(TAG_REGION, tag_region)
    (int* tag, ARLIM_P(tag_lo), ARLIM_P(tag_hi),
     const int* tagval, 
     const int lo[], const int hi[],
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], 
     const amrex::Real problo[], const int* level);

BL_FORT_PROC_DECL(SET_TRIGGER, set_trigger)
    (BL_FORT_FAB_ARG(state),
     const int    lo[], const int    hi[],
     const int domlo[], const int domhi[],
     const amrex::Real dx[], const amrex::Real xlo[], const amrex::Real problo[],
     const int* level, const amrex::Real* avg_den, const int* trigger);

#endif
