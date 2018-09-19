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
amrex_c_fab_copy (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                  amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                  const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                  const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_copytomem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi, amrex_real* AMREX_RESTRICT dst,
                       const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                       const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_copyfrommem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                         amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                         const int* AMREX_RESTRICT ncomp, const amrex_real* AMREX_RESTRICT src)
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
amrex_c_fab_setval (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi, 
                    amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                    const int* AMREX_RESTRICT ncomp, const amrex_real* AMREX_RESTRICT val)
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

void
amrex_c_fab_setval_ifnot (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                          amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                          const int* AMREX_RESTRICT ncomp,
                          const int* AMREX_RESTRICT msk, const int* AMREX_RESTRICT mlo, const int* AMREX_RESTRICT mhi,
                          const amrex_real* AMREX_RESTRICT val)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(msk, mlo, mhi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    if (0 == GET_VALUE(msk,i,j,k)) {
                        GET_VALUE_N(dst,i,j,k,n) = *val;
                    }
                }
            }
        }
    }
}

amrex_real
amrex_c_fab_norminfmask (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                         const int* AMREX_RESTRICT msk, const int* AMREX_RESTRICT mlo, const int* AMREX_RESTRICT mhi,
                         const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo,
                         const int* AMREX_RESTRICT shi, const int* AMREX_RESTRICT ncomp)
{
    amrex_real nrm = 0.0;
    DEFINE_STRIDES(msk,mlo,mhi);
    DEFINE_STRIDES(src,slo,shi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    if (1 == GET_VALUE(msk,i,j,k)) {
                        nrm = fmax(nrm,fabs(GET_VALUE_N(src,i,j,k,n)));
                    }
                }
            }
        }
    }
    return nrm;
}

amrex_real
amrex_c_fab_norm (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                  const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                  const int* AMREX_RESTRICT ncomp, const int* AMREX_RESTRICT p)
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
amrex_c_fab_sum (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                 const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                 const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_plus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                  amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                  const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                  const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_minus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                   amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                   const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                   const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_mult (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                  amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                  const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                  const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_divide (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                    amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                    const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                    const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_protdivide (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                        amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                        const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                        const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_fab_invert (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                    amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                    const int* AMREX_RESTRICT ncomp, const amrex_real* AMREX_RESTRICT a)
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
amrex_c_fab_saxpy (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                   amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                   const amrex_real* AMREX_RESTRICT a,
                   const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                   const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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

    
void amrex_c_fab_xpay (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                       amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                       const amrex_real* AMREX_RESTRICT a,
                       const amrex_real* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                       const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
                    GET_VALUE_N(dst,i,j,k,n) = GET_VALUE_N(src,i+ioff,j+joff,k+koff,n)
                        + (*a) * GET_VALUE_N(dst,i,j,k,n);
                }
            }
        }
    }
}

void
amrex_c_fab_lincomb (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                     amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                     const amrex_real* AMREX_RESTRICT a, const amrex_real* AMREX_RESTRICT x,
                     const int* AMREX_RESTRICT xlo, const int* AMREX_RESTRICT xhi, const int* AMREX_RESTRICT xblo,
                     const amrex_real* AMREX_RESTRICT b, const amrex_real* AMREX_RESTRICT y,
                     const int* AMREX_RESTRICT ylo, const int* AMREX_RESTRICT yhi, const int* AMREX_RESTRICT yblo,
                     const int* AMREX_RESTRICT ncomp)
{
    DEFINE_STRIDES(dst, dlo, dhi);
    DEFINE_STRIDES(x, xlo, xhi);
    DEFINE_STRIDES(y, ylo, yhi);
    const int ixoff = xblo[0]-lo[0];
    const int jxoff = xblo[1]-lo[1];
    const int kxoff = xblo[2]-lo[2];
    const int iyoff = yblo[0]-lo[0];
    const int jyoff = yblo[1]-lo[1];
    const int kyoff = yblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) =
                        (*a) * GET_VALUE_N(x,i+ixoff,j+jxoff,k+kxoff,n) +
                        (*b) * GET_VALUE_N(y,i+iyoff,j+jyoff,k+kyoff,n);
                }
            }
        }
    }
}

    
void
amrex_c_fab_addproduct (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                        amrex_real* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                        const amrex_real* AMREX_RESTRICT src1, const int* AMREX_RESTRICT s1lo, const int* AMREX_RESTRICT s1hi,
                        const amrex_real* AMREX_RESTRICT src2, const int* AMREX_RESTRICT s2lo, const int* AMREX_RESTRICT s2hi,
                        const int* AMREX_RESTRICT ncomp)
{
    DEFINE_STRIDES(dst,dlo,dhi);
    DEFINE_STRIDES(src1,s1lo,s1hi);
    DEFINE_STRIDES(src2,s2lo,s2hi);
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    GET_VALUE_N(dst,i,j,k,n) += GET_VALUE_N(src1,i,j,k,n) * GET_VALUE_N(src2,i,j,k,n);
                }
            }
        }
    }
}

amrex_real
amrex_c_fab_dot (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                 const amrex_real* AMREX_RESTRICT x, const int* AMREX_RESTRICT xlo, const int* AMREX_RESTRICT xhi,
                 const amrex_real* AMREX_RESTRICT y, const int* AMREX_RESTRICT ylo, const int* AMREX_RESTRICT yhi,
                 const int* AMREX_RESTRICT yblo, const int* AMREX_RESTRICT ncomp)
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
                    r += GET_VALUE_N(x,i,j,k,n) * GET_VALUE_N(y,i+ioff,j+joff,k+koff,n);
                }
            }
        }
    }
    return r;
}
    
amrex_real
amrex_c_fab_dot_mask (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                      const amrex_real* AMREX_RESTRICT x, const int* AMREX_RESTRICT xlo, const int* AMREX_RESTRICT xhi,
                      const amrex_real* AMREX_RESTRICT y, const int* AMREX_RESTRICT ylo, const int* AMREX_RESTRICT yhi,
                      const int* AMREX_RESTRICT yblo,
                      const int* AMREX_RESTRICT m, const int* AMREX_RESTRICT mlo, const int* AMREX_RESTRICT mhi,
                      const int* AMREX_RESTRICT ncomp)
{
    amrex_real r = 0.0;
    DEFINE_STRIDES(x, xlo, xhi);
    DEFINE_STRIDES(y, ylo, yhi);
    DEFINE_STRIDES(m, mlo, mhi);
    const int ioff = yblo[0]-lo[0];
    const int joff = yblo[1]-lo[1];
    const int koff = yblo[2]-lo[2];
    for (int n = 0; n < *ncomp; ++n) {
        for         (int k = lo[2]; k <= hi[2]; ++k) {
            for     (int j = lo[1]; j <= hi[1]; ++j) {
                for (int i = lo[0]; i <= hi[0]; ++i) {
                    r += GET_VALUE_N(x,i,j,k,n) * GET_VALUE_N(y,i+ioff,j+joff,k+koff,n)
                        * GET_VALUE(m,i,j,k);
                }
            }
        }
    }
    return r;
}

void
amrex_c_ifab_copy (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                   int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                   const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                   const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_ifab_copytomem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi, int* AMREX_RESTRICT dst,
                        const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                        const int* AMREX_RESTRICT ncomp)
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
amrex_c_ifab_copyfrommem (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                          int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                          const int* AMREX_RESTRICT ncomp, const int* AMREX_RESTRICT src)
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
amrex_c_ifab_setval (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi, 
                     int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                     const int* AMREX_RESTRICT ncomp, const int* AMREX_RESTRICT val)
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

void
amrex_c_ifab_plus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                   int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                   const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                   const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
amrex_c_ifab_minus (const int* AMREX_RESTRICT lo, const int* AMREX_RESTRICT hi,
                    int* AMREX_RESTRICT dst, const int* AMREX_RESTRICT dlo, const int* AMREX_RESTRICT dhi,
                    const int* AMREX_RESTRICT src, const int* AMREX_RESTRICT slo, const int* AMREX_RESTRICT shi,
                    const int* AMREX_RESTRICT sblo, const int* AMREX_RESTRICT ncomp)
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
    
