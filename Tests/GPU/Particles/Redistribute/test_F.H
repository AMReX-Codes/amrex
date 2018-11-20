#ifndef _TEST_F_H_
#define _TEST_F_H_

#include <AMReX_BLFort.H>

#ifdef IBM
#define NOUNDERSCORE
#endif

#ifdef NOUNDERSCORE
#define FTOC(x) x
#else
#define FTOC(x) x##_
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void FTOC(move_particles)
        (const int np,
         void* particles,
         amrex::Real* vxp,
         amrex::Real* vyp,
         amrex::Real* vzp,
         const amrex::Real* dx);

#ifdef __cplusplus
};
#endif

#endif /*_TEST_F_H_*/
