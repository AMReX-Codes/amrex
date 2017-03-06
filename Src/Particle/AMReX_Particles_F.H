#ifndef AMREX_Particles_F_H_
#define AMREX_Particles_F_H_

#include <AMReX_BLFort.H>

#ifdef USE_PARTICLES

extern "C"
{

    void amrex_particle_set_position (amrex::Real*, int ns, int np,
                                      D_DECL(const amrex::Real* x,
                                             const amrex::Real* y,
                                             const amrex::Real* z));

    void amrex_particle_get_position (const amrex::Real*, int ns, int np,
                                      D_DECL(amrex::Real* x,
                                             amrex::Real* y,
                                             amrex::Real* z));

    void amrex_deposit_cic(const amrex::Real*, int ns, int np,
                           amrex::Real* rho, const int* lo, const int* hi,
                           const amrex::Real* plo, const amrex::Real* dx);


BL_FORT_PROC_DECL(PART_SUMMASSDOWN,part_summassdown)
    (BL_FORT_FAB_ARG(crse_fab),
     const BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

BL_FORT_PROC_DECL(PART_SUMMASSUP,part_summassup)
    (const BL_FORT_FAB_ARG(crse_fab),
     BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

BL_FORT_PROC_DECL(PART_SUMDENSUP,part_sumdensup)
    (const BL_FORT_FAB_ARG(crse_fab),
     BL_FORT_FAB_ARG(fine_fab),
     const int* ovlo, const int* ovhi, const int* rat);

#if (BL_SPACEDIM==2)
    BL_FORT_PROC_DECL(BICUBICINTERP,bicubicinterp)
    (const amrex::Real* x, const amrex::Real* y, const amrex::Real* fab,
     const int* lo, const int* hi,
     const amrex::Real* plo, const amrex::Real* phi, const amrex::Real* dx, amrex::Real* z);
#endif

};

#endif

#endif /*_Particles_F_H_*/
