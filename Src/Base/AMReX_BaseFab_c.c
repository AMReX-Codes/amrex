#include <AMReX_BaseFab_c.H>
#include <math.h>

#define DEFINE_STRIDES(q, qlo, qhi) \
    const int jstride_##q = qhi[0]-qlo[0]+1; \
    const int kstride_##q = jstride_##q*(qhi[1]-qlo[1]+1); \
    const int nstride_##q = kstride_##q*(qhi[2]-qlo[2]+1); \
    const int xlbound_##q = qlo[0]; \
    const int ylbound_##q = qlo[1]; \
    const int zlbound_##q = qlo[2];

#define INDEX(q,i,j,k) \
    i-xlbound_##q + (j-ylbound_##q)*jstride_##q + (k-zlbound_##q)*kstride_##q
#define INDEX_N(q,i,j,k,n) \
    INDEX(q,i,j,k) + n*nstride_##q

#define GET_VALUE(q,i,j,k) q[INDEX(q,i,j,k)]
#define GET_VALUE_N(q,i,j,k,n) q[INDEX_N(q,i,j,k,n)]

#define GET_PTR(q,i,j,k) q+INDEX(q,i,j,k)
#define GET_PTR_N(q,i,j,k,n) q+INDEX_N(q,i,j,k,n)

/*****************************************************************************/

void
amrex_c_fab_copy (const int* lo, const int* hi,
                  amrex_real* restrict dst, const int* dlo, const int* dhi,
                  const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                  const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) = GET_VALUE_N(src,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }
}

long
amrex_c_fab_copytomem (const int* lo, const int* hi,
                       amrex_real* restrict dst,
                       const amrex_real* restrict src, const int* slo, const int* shi,
                       const int* ncomp)
{
    DEFINE_STRIDES(src,slo,shi);
    long xlength = hi[0]-lo[0]+1;
    long offset = -lo[0];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    dst[offset+i] = GET_VALUE_N(src,i,j,k,n);
                }
                offset += xlength;
            }
        }
    }
    return offset + lo[0];
}

long
amrex_c_fab_copyfrommem (const int* lo, const int* hi,
                         amrex_real* restrict dst, const int* dlo, const int* dhi, const int* ncomp,
                         const amrex_real* restrict src)
{
    DEFINE_STRIDES(dst,dlo,dhi);
    long xlength = hi[0]-lo[0]+1;
    long offset = -lo[0];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) = src[offset+i];
                }
                offset += xlength;
            }
        }
    }
    return offset + lo[0];
}


void
amrex_c_fab_setval (const int* lo, const int* hi, 
                    amrex_real* restrict dst, const int* dlo, const int* dhi, const int* ncomp,
                    const amrex_real* restrict val)
{
    DEFINE_STRIDES(dst,dlo,dhi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) = *val;
                }
            }
        }
    }
}
    
amrex_real
amrex_c_fab_norminfmask (const int* lo, const int* hi,
                         const int* restrict msk, const int* mlo, const int* mhi,
                         const amrex_real* restrict src, const int* slo, const int* shi, const int* ncomp)
{
    amrex_real nrm = 0.0;
    DEFINE_STRIDES(msk,mlo,mhi);
    DEFINE_STRIDES(src,slo,shi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    if (1 == GET_VALUE(msk,i,j,k)) {
                        nrm = fmax(nrm,GET_VALUE_N(src,i,j,k,n));
                    }
                }
            }
        }
    }
    return nrm;
}

amrex_real
amrex_c_fab_norm (const int* lo, const int* hi,
                  const amrex_real* restrict src, const int* slo, const int* shi, const int* ncomp,
                  const int* p)
{
    amrex_real nrm = 0.0;
    DEFINE_STRIDES(src,slo,shi);
    if (0 == *p) {
        for (int n = 0; n < *ncomp; ++n) {
            for         (int k = lo[2]; k <= hi[2]; ++k) {
                for     (int j = lo[1]; j <= hi[1]; ++j) {
                    for (int i = lo[0]; i <= hi[0]; ++i) {
                        nrm = fmax(nrm,fabs(GET_VALUE_N(src,i,j,k,n)));
                    }
                }
            }
        }
    } else if (1 == *p) {
        for (int n = 0; n < *ncomp; ++n) {
            for         (int k = lo[2]; k <= hi[2]; ++k) {
                for     (int j = lo[1]; j <= hi[1]; ++j) {
                    for (int i = lo[0]; i <= hi[0]; ++i) {
                        nrm += fabs(GET_VALUE_N(src,i,j,k,n));
                    }
                }
            }
        }
    }
    return nrm;    
}

amrex_real
amrex_c_fab_sum (const int* lo, const int* hi,
                 const amrex_real* restrict src, const int* slo, const int* shi, const int* ncomp)
{
    amrex_real sm = 0.0;
    DEFINE_STRIDES(src,slo,shi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    sm += GET_VALUE_N(src,i,j,k,n);
                }
            }
        }
    }
    return sm;    
}
    
void
amrex_c_fab_plus (const int* lo, const int* hi,
                  amrex_real* restrict dst, const int* dlo, const int* dhi,
                  const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                  const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) += GET_VALUE_N(src,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }    
}

void
amrex_c_fab_minus (const int* lo, const int* hi,
                   amrex_real* restrict dst, const int* dlo, const int* dhi,
                   const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                   const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) -= GET_VALUE_N(src,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }    
}
    
void
amrex_c_fab_mult (const int* lo, const int* hi,
                  amrex_real* restrict dst, const int* dlo, const int* dhi,
                  const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                  const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) *= GET_VALUE_N(src,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }    
}

void
amrex_c_fab_divide (const int* lo, const int* hi,
                    amrex_real* restrict dst, const int* dlo, const int* dhi,
                    const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                    const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) /= GET_VALUE_N(src,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }    
}

void
amrex_c_fab_protdivide (const int* lo, const int* hi,
                        amrex_real* restrict dst, const int* dlo, const int* dhi,
                        const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                        const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    if (0. != GET_VALUE_N(src,i+ioff,j+joff,k+koff,n)) {
                        GET_VALUE_N(dst,i,j,k,n) /= GET_VALUE_N(src,i+ioff,j+joff,k+koff,n);
                    }
                }
            }
        }
    }    
}
    
void
amrex_c_fab_invert (const int* lo, const int* hi,
                    amrex_real* restrict dst, const int* dlo, const int* dhi, const int* ncomp,
                    const amrex_real* restrict a)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) = (*a) / GET_VALUE_N(dst,i,j,k,n);
                }
            }
        }
    }    
}

void
amrex_c_fab_saxpy (const int* lo, const int* hi,
                   amrex_real* restrict dst, const int* dlo, const int* dhi,
                   const amrex_real* restrict a,
                   const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                   const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) += (*a) * GET_VALUE_N(src,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }
}

    
void amrex_c_fab_xpay (const int* lo, const int* hi,
                       amrex_real* restrict dst, const int* dlo, const int* dhi,
                       const amrex_real* restrict a,
                       const amrex_real* restrict src, const int* slo, const int* shi, const int* sblo,
                       const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(src, slo, shi);
    const int ioff = sblo[0]-lo[0];
    const int joff = sblo[1]-lo[1];
    const int koff = sblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) = GET_VALUE_N(src,i+ioff,j+joff,k+koff,n) + (*a) * GET_VALUE_N(dst,i,j,k,n);
                }
            }
        }
    }
}

void
amrex_c_fab_lincomb (const int* lo, const int* hi,
                     amrex_real* restrict dst, const int* dlo, const int* dhi,
                     const amrex_real* restrict a, const amrex_real* restrict x, const int* xlo, const int* xhi, const int* xblo,
                     const amrex_real* restrict b, const amrex_real* restrict y, const int* ylo, const int* yhi, const int* yblo,
                     const int* ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(x, xlo, xhi);
    DEFINE_STRIDES(y, ylo, yhi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
//                    GET_VALUE_N(dst,i,j,k,n) =
//                        (*a) * GET_VALUE_N(x,i+
                }
            }
        }
    }
}

    
    void amrex_c_fab_addproduct (const int* lo, const int* hi,
                                 amrex_real* restrict dst, const int* dlo, const int* dhi,
                                 const amrex_real* restrict src1, const int* s1lo, const int* s1hi,
                                 const amrex_real* restrict src2, const int* s2lo, const int* s2hi,
                                 const int* ncomp);

amrex_real
amrex_c_fab_dot (const int* lo, const int* hi,
                 const amrex_real* restrict x, const int* xlo, const int* xhi,
                 const amrex_real* restrict y, const int* ylo, const int* yhi, const int* yblo,
                 const int* ncomp)
{
    amrex_real r = 0.0;
    DEFINE_STRIDES(x, xlo, xhi);
    DEFINE_STRIDES(y, ylo, yhi);
    const int ioff = yblo[0]-lo[0];
    const int joff = yblo[1]-lo[1];
    const int koff = yblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    r = r + GET_VALUE_N(x,i,j,k,n) * GET_VALUE_N(y,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }
    return r;
}
    
    amrex_real amrex_c_fab_dot_mask (const int* lo, const int* hi,
                                     const amrex_real* restrict x, const int* xlo, const int* xhi,
                                     const amrex_real* restrict y, const int* ylo, const int* yhi, const int* yblo,
                                     const int* m, const int* mlo, const int* mhi,
                                     const int* ncomp);

    void amrex_c_ifab_copy (const int* lo, const int* hi,
                            int* dst, const int* dlo, const int* dhi,
                            const int* src, const int* slo, const int* shi, const int* sblo,
                            const int* ncomp);

    long amrex_c_ifab_copytomem (const int* lo, const int* hi,
                                 int* dst,
                                 const int* src, const int* slo, const int* shi,
                                 const int* ncomp);

    long amrex_c_ifab_copyfrommem (const int* lo, const int* hi,
                                   const int* dst, const int* dlo, const int* dhi, const int* ncomp,
                                   const int* src);
    
    void amrex_c_ifab_setval (const int* lo, const int* hi, 
                              const int* dst, const int* dlo, const int* dhi, const int* ncomp,
                              const int* val);

    void amrex_c_ifab_plus (const int* lo, const int* hi,
                            int* dst, const int* dlo, const int* dhi,
                            const int* src, const int* slo, const int* shi, const int* sblo,
                            const int* ncomp);    
    
    void amrex_c_ifab_minus (const int* lo, const int* hi,
                             int* dst, const int* dlo, const int* dhi,
                             const int* src, const int* slo, const int* shi, const int* sblo,
                             const int* ncomp);    
    
    void amrex_c_fab_setval_ifnot (const int* lo, const int* hi,
                                   amrex_real* restrict dst, const int* dlo, const int* dhi, const int* ncomp,
                                   const int* msk, const int* mlo, const int* mhi,
                                   const int val);



