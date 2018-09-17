#include <AMReX_BaseFab_c.H>

#define DEFINE_STRIDES(q, qlo, qhi) \
    const int jstride_##q = qhi[0]-qlo[0]+1; \
    const int kstride_##q = jstride_##q*(qhi[1]-qlo[1]+1); \
    const int nstride_##q = kstride_##q*(qhi[2]-qlo[2]+1);

#define INDEX(q,i,j,k) \
    i-q##lo[0] + (j-q##lo[1])*jstride_##q + (k-q##lo[2])*kstride_##q
#define INDEX_N(q,i,j,k,n) \
    INDEX(q,i,j,k) + n*nstride_##q

#define GET_VALUE(q,i,j,k) q[INDEX(q,i,j,k)]
#define GET_VALUE_N(q,i,j,k,n) q[INDEX_N(q,i,j,k,n)]

#define GET_PTR(q,i,j,k) q+INDEX(q,i,j,k)
#define GET_PTR_N(q,i,j,k,n) q+INDEX_N(q,i,j,k,n)

amrex_real
amrex_c_fab_dot (const int* lo, const int* hi,
                 const amrex_real* restrict x, const int* xlo, const int* xhi,
                 const amrex_real* restrict y, const int* ylo, const int* yhi, const int* yblo,
                 const int* ncomp)
{
    amrex_real r = 0.0;

    DEFINE_STRIDES(x, xlo, xhi);
    DEFINE_STRIDES(y, ylo, yhi);

    for (int n = 0; n < *ncomp; ++n) {
        for (int k = lo[2]; k <= hi[2]; ++k) {
            for (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    r = r + GET_VALUE_N(x,i,j,k,n) * GET_VALUE_N(y,i+yblo[0],j+yblo[1],k+yblo[2],n);
                }
            }
        }
    }

    return r;
}

