#ifndef BASEFAB_F_H_
#define BASEFAB_F_H_

#include <REAL.H>
#include <ArrayLim.H>

extern "C"
{
    void fort_fab_copy (const int* lo, const int* hi,
			Real* dst, const int* dlo, const int* dhi,
			const Real* src, const int* slo, const int* shi, const int* sblo,
			const int* ncomp);

    void fort_fab_copytomem (const int* lo, const int* hi,
			     Real* dst,
			     const Real* src, const int* slo, const int* shi,
			     const int* ncomp);

    void fort_fab_copyfrommem (const int* lo, const int* hi,
			       const Real* dst, const int* dlo, const int* dhi, const int* ncomp,
			       const Real* src);

    void fort_fab_setval (const int* lo, const int* hi, 
			  const Real* dst, const int* dlo, const int* dhi, const int* ncomp,
			  const Real* val);

    Real fort_fab_norm (const int* lo, const int* hi,
			const Real* src, const int* slo, const int* shi, const int* ncomp,
			const int* p);

    Real fort_fab_sum (const int* lo, const int* hi,
		       const Real* src, const int* slo, const int* shi, const int* ncomp);

    void fort_fab_plus (const int* lo, const int* hi,
			Real* dst, const int* dlo, const int* dhi,
			const Real* src, const int* slo, const int* shi, const int* sblo,
			const int* ncomp);    

    void fort_fab_minus (const int* lo, const int* hi,
			 Real* dst, const int* dlo, const int* dhi,
			 const Real* src, const int* slo, const int* shi, const int* sblo,
			 const int* ncomp);    

    void fort_fab_mult (const int* lo, const int* hi,
			Real* dst, const int* dlo, const int* dhi,
			const Real* src, const int* slo, const int* shi, const int* sblo,
			const int* ncomp);    

    void fort_fab_divide (const int* lo, const int* hi,
			  Real* dst, const int* dlo, const int* dhi,
			  const Real* src, const int* slo, const int* shi, const int* sblo,
			  const int* ncomp);    

    void fort_fab_protdivide (const int* lo, const int* hi,
			      Real* dst, const int* dlo, const int* dhi,
			      const Real* src, const int* slo, const int* shi, const int* sblo,
			      const int* ncomp);        

    void fort_fab_invert (const int* lo, const int* hi,
			  Real* dst, const int* dlo, const int* dhi, const int* ncomp,
			  const Real* a);

    void fort_fab_saxpy (const int* lo, const int* hi,
			 Real* dst, const int* dlo, const int* dhi,
			 const Real* a,
			 const Real* src, const int* slo, const int* shi, const int* sblo,
			 const int* ncomp);

    void fort_fab_xpay (const int* lo, const int* hi,
			Real* dst, const int* dlo, const int* dhi,
			const Real* a,
			const Real* src, const int* slo, const int* shi, const int* sblo,
			const int* ncomp);

    void fort_fab_lincomb (const int* lo, const int* hi,
			   Real* dst, const int* dlo, const int* dhi,
			   const Real* a, const Real* x, const int* xlo, const int* xhi, const int* xblo,
			   const Real* b, const Real* y, const int* ylo, const int* yhi, const int* yblo,
			   const int* ncomp);

    void fort_fab_addproduct (const int* lo, const int* hi,
			      Real* dst, const int* dlo, const int* dhi,
			      const Real* src1, const int* s1lo, const int* s1hi,
			      const Real* src2, const int* s2lo, const int* s2hi,
			      const int* ncomp);

    Real fort_fab_dot (const int* lo, const int* hi,
		       const Real* x, const int* xlo, const int* xhi,
		       const Real* y, const int* ylo, const int* yhi, const int* yblo,
		       const int* ncomp);
}

#endif
