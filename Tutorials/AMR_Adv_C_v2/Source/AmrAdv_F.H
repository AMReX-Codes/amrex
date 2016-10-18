#ifndef _Adv_F_H_
#define _Adv_F_H_
#include <BLFort.H>

extern "C" 
{
    void initdata(const int& level, const Real& time, 
		  const int* lo, const int* hi,
		  BL_FORT_FAB_ARG_3D(state),
		  const Real* dx, const Real* problo);

    void get_face_velocity(const int& level, const Real& time, 
			   D_DECL(BL_FORT_FAB_ARG(xvel),
				  BL_FORT_FAB_ARG(yvel),
				  BL_FORT_FAB_ARG(zvel)),
			   const Real* dx, const Real* problo);

    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
		     const BL_FORT_FAB_ARG_3D(state),
		     const int* tagval, const int* clearval,
		     const int* lo, const int* hi,
		     const Real* dx, const Real* problo,
		     const Real* time, const Real* phierr);

    void advect(const Real& time, const int* lo, const int*hi,
		const BL_FORT_FAB_ARG_3D(statein),
		BL_FORT_FAB_ARG_3D(stateout),
		D_DECL(const BL_FORT_FAB_ARG_3D(xvel),
		       const BL_FORT_FAB_ARG_3D(yvel),
		       const BL_FORT_FAB_ARG_3D(zvel)),
		D_DECL(BL_FORT_FAB_ARG_3D(fx),
		       BL_FORT_FAB_ARG_3D(fy),
		       BL_FORT_FAB_ARG_3D(fz)),
		const Real* dx, const Real& dt);
}

#endif
