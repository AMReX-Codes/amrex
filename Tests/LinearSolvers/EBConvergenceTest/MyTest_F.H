#ifndef MyTest_F_H_
#define MyTest_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif 

#if AMREX_SPACEDIM == 2
    void build_rhs_2d(const int* lo, const int* hi, 
                   amrex_real* rhs, const int* rlo, const int* rhi, 
                   const amrex_real* a, const int* alo, const int* ahi, 
                   const amrex_real* bx, const int* bxlo, const int* bxhi, 
                   const amrex_real* by, const int* bylo, const int* byhi, 
                   const amrex_real* dx, const amrex_real* prob_lo, 
                   const amrex_real* prob_hi);

    void build_a_2d(const int* lo, const int* hi,
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi,
                    amrex_real* a, const int* alo, const int* ahi, 
                    const amrex_real* dx);

    void build_b_2d(const int* lo, const int* hi,
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi,
                    amrex_real* bx, const int* bxlo, const int* bxhi, 
                    amrex_real* by, const int* bylo, const int* byhi, 
                    const amrex_real* dx);

#elif AMREX_SPACEDIM ==3
    void build_rhs_3d(const int* lo, const int* hi, 
                   amrex_real* rhs, const int* rlo, const int* rhi, 
                   const amrex_real* a, const int* alo, const int* ahi, 
                   const amrex_real* bx, const int* bxlo, const int* bxhi, 
                   const amrex_real* by, const int* bylo, const int* byhi,
                   const amrex_real* bz, const int* bzlo, const int* bzhi, 
                   const amrex_real* dx, const amrex_real* prob_lo, 
                   const amrex_real* prob_hi);

    void build_a_3d(const int* lo, const int* hi,
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi,
                    amrex_real* a, const int* alo, const int* ahi, 
                    const amrex_real* dx);

    void build_b_3d(const int* lo, const int* hi,
                    const amrex_real* prob_lo, 
                    const amrex_real* prob_hi,
                    amrex_real* bx, const int* bxlo, const int* bxhi, 
                    amrex_real* by, const int* bylo, const int* byhi, 
                    amrex_real* bz, const int* bzlo, const int* bzhi, 
                    const amrex_real* dx);

#endif
    void apply_bc(const int* lo, const int* hi, 
                  amrex_real* phi, const int* phlo, const int* phhi, 
                  const amrex_real* dx, const amrex_real* prob_lo, 
                  const amrex_real* prob_hi); 
#ifdef __cplusplus
}
#endif
#endif 
