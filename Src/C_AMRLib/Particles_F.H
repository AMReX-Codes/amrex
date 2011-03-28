#ifndef _Particles_F_H_

#define _Particles_F_H_

#include <BLFort.H>

#ifdef USE_PARTICLES

BL_FORT_PROC_DECL(PART_SUMMASSDOWN,part_summassdown)
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

BL_FORT_PROC_DECL(PART_SUMMASSUP,part_summassup)
    (BL_FORT_FAB_ARG(crse_fab), const int& nc,
     const BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

#endif

#endif /*_Particles_F_H_*/
