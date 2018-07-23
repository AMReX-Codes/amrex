#ifndef _DEPOSIT_F_H_
#define _DEPOSIT_F_H_

#include <AMReX_BLFort.H>
#include <AMReX_SPACE.H>

#ifdef USE_PARTICLES

#ifdef __cplusplus
extern "C"
{
#endif
    
    void push_particles(const amrex_real*, int ns, int np);
    
    void deposit_cic(const amrex_real*, int ns, int np,
                     amrex_real* rho, const int* lo, const int* hi,
                     const amrex_real* plo, const amrex_real* dx);
    
    void interpolate_cic(const amrex_real*, int ns, int np,
                         amrex_real* acc, const int* lo, const int* hi,
                         const amrex_real* plo, const amrex_real* dx);
    
    void cuda_push_particles(const amrex_real*, int ns, int np);
    
    void cuda_deposit_cic(const amrex_real*, int ns, int np, int gid,
                          amrex_real* rho, const int* lo, const int* hi,
                          const amrex_real* plo, const amrex_real* dx);
    
    void cuda_interpolate_cic(const amrex_real*, int ns, int np, int gid,
                              amrex_real* acc, const int* lo, const int* hi,
                              const amrex_real* plo, const amrex_real* dx);

#ifdef __cplusplus
};
#endif

#endif

#endif /*_DEPOSIT_F_H_*/
