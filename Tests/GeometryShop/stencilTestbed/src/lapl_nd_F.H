#ifndef _lapl_nd_F_H_
#define _lapl_nd_F_H_

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>

namespace amrex
{
extern "C"
{
  void fort_lapl_simple(Real*       lph, const int* lolph, const int* hilph,
                        const Real* phi, const int* lophi, const int* hiphi,
                        const int* reglo, const int* reghi, const Real* dx); 
}
}
#endif
