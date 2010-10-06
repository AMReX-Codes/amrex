#ifndef _Particles_F_H_
#define _Particles_F_H_
#include <BLFort.H>

#ifdef USE_PARTICLES
BL_FORT_PROC_DECL(CA_ESTPARTDT,ca_estpartdt)
    (const BL_FORT_FAB_ARG(state),
     const int lo[], const int hi[],
     Real* dt);
#endif
#endif
