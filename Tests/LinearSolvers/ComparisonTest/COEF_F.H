#ifndef _COEFF_H_
#define _COEFF_H_

#include <REAL.H>

#if        defined(BL_LANG_FORT)
#define FORT_SET_COEF        set_coef
#define FORT_COEF_TO_EDGES   coef_to_edges
#else

#if    defined(BL_FORT_USE_UPPERCASE)
#define FORT_SET_COEF        SET_COEF
#define FORT_COEF_TO_EDGES   COEF_TO_EDGES
#elif  defined(BL_FORT_USE_LOWERCASE)
#define FORT_SET_COEF        set_coef
#define FORT_COEF_TO_EDGES   coef_to_edges
#elif  defined(BL_FORT_USE_UNDERSCORE)
#define FORT_SET_COEF        set_coef_
#define FORT_COEF_TO_EDGES   coef_to_edges_
#endif

#include <ArrayLim.H>

extern "C"
{
  void FORT_SET_COEF(Real* exac, Real* alph, Real* beta, Real* rhs,
		     const int* lo, const int* hi, const Real* prob_lo, const Real* prob_hi, 
		     const Real* dx, const Real& a, const Real& b, 
		     const Real& sigma, const Real& w, int& ibnd);

  void FORT_COEF_TO_EDGES (int * dir, 
			   Real* edge, ARLIM_P(edge_lo), ARLIM_P(edge_hi),
			   const Real* coef, ARLIM_P(coef_lo), ARLIM_P(coef_hi),
			   const int* lo, const int* hi);
}
#endif

#endif /*_COEF_F_H_*/
