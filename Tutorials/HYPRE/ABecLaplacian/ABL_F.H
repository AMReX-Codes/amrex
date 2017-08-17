#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C" {
#endif

    void fort_init_rhs (const int* lo, const int* hi, amrex_real* rhs, const int* rlo, const int* rhi,
                        const amrex_real* dx, const amrex_real* a, const amrex_real* b,
                        const amrex_real* sigma, const amrex_real* w, const int* ibnd);

    void fort_init_cc_coef (const int* lo, const int* hi, amrex_real* beta, const int* blo, const int*bhi,
                            const amrex_real* dx, const amrex_real* sigma, const amrex_real* w);

    void fort_comp_asol (const amrex_real* soln, const int* slo, const int* shi,
                         const amrex_real* dx, const int* ibnd);

#ifdef __cplusplus
}
#endif
