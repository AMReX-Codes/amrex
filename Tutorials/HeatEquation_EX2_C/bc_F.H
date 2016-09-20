#include <BLFort.H>

BL_FORT_PROC_DECL(PHIFILL,phifill)
    (BL_FORT_FAB_ARG_3D(state),
     const int* dlo, const int* dhi,
     const Real* dx, const Real* glo, 
     const Real* time, const int* bc);
