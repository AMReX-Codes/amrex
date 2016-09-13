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

    void fort_fab_saxpy (const int* lo, const int* hi,
			 Real* dst, const int* dlo, const int* dhi,
			 const Real* a,
			 const Real* src, const int* slo, const int* shi, const int* sblo,
			 const int* ncomp);
}

#endif
