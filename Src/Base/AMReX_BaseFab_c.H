#include <AMReX_REAL.H>

#ifdef __cplusplus
extern "C" {
#endif

    amrex_real amrex_c_fab_dot (const int* lo, const int* hi,
                                const amrex_real* x, const int* xlo, const int* xhi,
                                const amrex_real* y, const int* ylo, const int* yhi, const int* yblo,
                                const int* ncomp);


#ifdef __cplusplus
}
#endif
