#ifndef _EBEllipticFriggingFortran_F_H___
#define _EBEllipticFriggingFortran_F_H___

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>


#ifdef __cplusplus
extern "C"
{
#endif
  void ebefnd_decrinvrelcoefebco(Real*       relco, const int* relcolo, const int* relcohi, const int* incorelco,
                                 const Real* bcoef, const int* bcoeflo, const int* bcoefhi, const int* incobcoef,
                                 const int*  gridlo, const int* gridhi, 
                                 const Real*  beta, const Real* dx, const int* idir);


  void ebefnd_invertlambdaebco   (Real*       relco, const int* relcolo, const int* relcohi, const int* incorelco,
                                  const int*  gridlo, const int* gridhi, const Real& safety);

    
#ifdef __cplusplus
}
#endif

#endif
