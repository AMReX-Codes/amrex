#ifndef AMREX_COMPUTE_FORCE_F_H_
#define AMREX_COMPUTE_FORCE_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif

    void sum_fine_to_crse_nodal(const int* lo, const int* hi, const int* lrat,
                                      amrex::Real* crse, const int* clo, const int* chi,
                                const amrex::Real* fine, const int* flo, const int* fhi);

    void zero_out_bndry(const int* lo, const int* hi,
                        amrex::Real* input_data, amrex::Real* bndry_data, 
                        const int* mask);

    void build_mask(const int* lo, const int* hi,
                    const int* tmp_mask, int* mask, const int* ncells);

    void compute_E_nodal (const int* lo, const int* hi,
                          const amrex::Real* phi,
                          amrex::Real* Ex, amrex::Real* Ey, 
#if BL_SPACEDIM == 3                          
                          amrex::Real* Ez, 
#endif
                          const amrex::Real* dx);
    
    void deposit_cic(const amrex::Real* particles, int ns, int np,
                     const amrex::Real* weights, 
                     const amrex::Real* charge,
                     amrex::Real* rho, const int* lo, const int* hi,
                     const amrex::Real* plo, const amrex::Real* dx,
                     const int* ng);
    
    void interpolate_cic_two_levels(const amrex::Real* particles, int ns, int np,
                                    amrex::Real* Ex_p, amrex::Real* Ey_p, 
#if BL_SPACEDIM == 3                                                              
                                    amrex::Real* Ez_p,
#endif
                                    const amrex::Real* Ex, const amrex::Real* Ey, 
#if BL_SPACEDIM == 3                          
                                    const amrex::Real* Ez,
#endif
                                    const int* lo, const int* hi, const amrex::Real* dx,
                                    const amrex::Real* cEx, const amrex::Real* cEy, 
#if BL_SPACEDIM == 3                          
                                    const amrex::Real* cEz,
#endif
                                    const int* mask,
                                    const int* clo, const int* chi, const amrex::Real* cdx,
                                    const amrex::Real* plo, const int* ng, const int* lev);

    void interpolate_cic(const amrex::Real* particles, int ns, int np,
                         amrex::Real* Ex_p, amrex::Real* Ey_p, 
#if BL_SPACEDIM == 3                          
                                    amrex::Real* Ez_p,
#endif
                         const amrex::Real* Ex, const amrex::Real* Ey, 
#if BL_SPACEDIM == 3                          
                         const amrex::Real* Ez,
#endif
                         const int* lo, const int* hi,
                         const amrex::Real* plo, const amrex::Real* dx,
                         const int* ng);

    void push_leapfrog(amrex::Real* particles, int ns, int np, 
                       amrex::Real* vx_p, amrex::Real* vy_p, 
#if BL_SPACEDIM == 3                          
                         amrex::Real* vz_p,
#endif
                       const amrex::Real* Ex_p, const amrex::Real* Ey_p, 
#if BL_SPACEDIM == 3                          
                         const amrex::Real* Ez_p,
#endif
                       const amrex::Real* charge, const amrex::Real* mass, const amrex::Real* dt,
                       const amrex::Real* prob_lo, const amrex::Real* prob_hi);

    void push_leapfrog_positions(amrex::Real* particles, int ns, int np,
                                 amrex::Real* vx_p, amrex::Real* vy_p, 
#if BL_SPACEDIM == 3                          
                                 amrex::Real* vz_p,
#endif
                                 const amrex::Real* dt, const amrex::Real* prob_lo,
                                 const amrex::Real* prob_hi);

#ifdef __cplusplus
}
#endif

#endif /* _COMPUTE_FORCE_F_H_ */
