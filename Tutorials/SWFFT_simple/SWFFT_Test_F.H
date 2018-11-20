#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void fort_init_phi_spatial (const int* lo, const int* hi, 
                        amrex_real* phi_spatial, const int* rlo, const int* rhi,
                        const int* dom_lo, const int* dom_hi, 
                        const amrex_real* problo, const amrex_real* probhi,
                        const amrex_real* dx, int* prob_type);

    void fort_comp_asol (const int* lo, const int* hi, 
                         const amrex_real* phi_dft, const int* slo, const int* shi,
                         const int* dom_lo, const int* dom_hi, 
                         const amrex_real* problo, const amrex_real* probhi,
                         const amrex_real* dx);

#ifdef __cplusplus
}
#endif
