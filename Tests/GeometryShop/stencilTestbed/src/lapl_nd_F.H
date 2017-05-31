#ifndef _lapl_nd_F_H_
#define _lapl_nd_F_H_

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

#ifdef __cplusplus
extern "C"
{
#endif
  void fort_lapl_simple(amrex::Real*       lph, const int* lolph, const int* hilph,
                        const amrex::Real* phi, const int* lophi, const int* hiphi,
                        const int* reglo, const int* reghi, const amrex::Real* dx); 
#ifdef __cplusplus
}
#endif

#endif
