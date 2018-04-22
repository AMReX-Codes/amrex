#ifndef AMREX_Particles_F_H_
#define AMREX_Particles_F_H_

#include <AMReX_BLFort.H>

#ifdef BL_SINGLE_PRECISION_PARTICLES
typedef float amrex_particle_real;
#else
typedef double amrex_particle_real;
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    void amrex_particle_set_position (amrex_particle_real*, int ns, int np,
                                      AMREX_D_DECL(const amrex_real* x,
                                             const amrex_real* y,
                                             const amrex_real* z));

    void amrex_particle_get_position (const amrex_particle_real*, int ns, int np,
                                      AMREX_D_DECL(amrex_real* x,
                                             amrex_real* y,
                                             amrex_real* z));

    void amrex_build_nl(void*, int ns, int np,
                        int* nbor_list, int* size_nl, const int* lo, const int* hi,
                        const amrex_real* plo, const amrex_real* dx);

    void amrex_deposit_cic(const amrex_particle_real*, int ns, int np, int nc,
                           amrex_real* rho, const int* lo, const int* hi,
                           const amrex_real* plo, const amrex_real* dx);

    void amrex_deposit_particle_dx_cic(const amrex_particle_real*, int ns, int np, int nc,
                                       amrex_real* rho, const int* lo, const int* hi,
                                       const amrex_real* plo, const amrex_real* dx,
                                       const amrex_real* particle_dx);

    void amrex_interpolate_cic(const amrex_particle_real*, int ns, int np,
                               const amrex_real* acc, const int* lo, const int* hi, int ncomp,
                               const amrex_real* plo, const amrex_real* dx);

    void amrex_atomic_accumulate_fab(const amrex_real*, const int*, const int*,
                                     amrex_real*, const int*, const int*, int);

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

#if (AMREX_SPACEDIM==2)
    BL_FORT_PROC_DECL(BICUBICINTERP,bicubicinterp)
    (const amrex_real* x, const amrex_real* y, const amrex_real* fab,
     const int* lo, const int* hi,
     const amrex_real* plo, const amrex_real* phi, const amrex_real* dx, amrex_real* z);
#endif

#ifdef __cplusplus
}
#endif

#endif /*_Particles_F_H_*/
