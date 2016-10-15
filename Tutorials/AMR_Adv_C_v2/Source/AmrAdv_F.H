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
    
}

#endif
