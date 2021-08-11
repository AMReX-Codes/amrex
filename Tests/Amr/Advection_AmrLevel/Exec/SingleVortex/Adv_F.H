#ifndef _Adv_F_H_
#define _Adv_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

  void initdata(const int* level, const amrex_real* time,
                const int* lo, const int* hi,
                BL_FORT_FAB_ARG_3D(state),
                const amrex_real* dx, const amrex_real* problo);

#ifdef __cplusplus
}
#endif

#endif
