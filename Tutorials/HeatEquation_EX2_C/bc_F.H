#ifndef _bc_F_H_ 
#define _bc_F_H_ 
#include <BLFort.H>

extern "C"
{
    void phifill (BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const Real* dx, const Real* glo, 
		  const Real* time, const int* bc);
}
#endif
