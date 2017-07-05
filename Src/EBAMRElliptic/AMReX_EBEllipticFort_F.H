#ifndef _EBEllipticFriggingFortran_F_H___
#define _EBEllipticFriggingFortran_F_H___

#include <AMReX_REAL.H>
#include <AMReX_ArrayLim.H>


#ifdef __cplusplus
extern "C"
{
#endif
  void ebefnd_decrinvrelcoefebco(amrex::Real*       relco, const int* relcolo, const int* relcohi, const int* incorelco,
                                 const amrex::Real* bcoef, const int* bcoeflo, const int* bcoefhi, const int* incobcoef,
                                 const int*  gridlo, const int* gridhi, 
                                 const amrex::Real*  beta, const amrex::Real* dx, const int* idir);

  void ebefnd_invertlambdaebco(amrex::Real*       relco, const int* relcolo, const int* relcohi, const int* incorelco,
                               const int*  gridlo, const int* gridhi, const amrex::Real& safety);

  void ebefnd_applyop_ebcond_nobcs(amrex::Real*       lph , const int* lphlo , const int* lphhi , const int* incolph,
                                   const amrex::Real* phi , const int* philo , const int* phihi , const int* incophi,
                                   const amrex::Real* aco , const int* acolo , const int* acohi , const int* incoaco,
                                   const amrex::Real* bco0, const int* bco0lo, const int* bco0hi, const int* incobco0,
                                   const amrex::Real* bco1, const int* bco1lo, const int* bco1hi, const int* incobco1,
                                   const amrex::Real* bco2, const int* bco2lo, const int* bco2hi, const int* incobco2,
                                   const int*  gridlo, const int* gridhi,
                                   const amrex::Real* dx,  const amrex::Real* alpha, const amrex::Real*  beta);


  void ebefnd_gscolor_ebcond(amrex::Real*       phi , const int* philo , const int* phihi , const int* incophi,
                             const amrex::Real* rhs , const int* rhslo , const int* rhshi , const int* incorhs,
                             const amrex::Real* lam , const int* lamlo , const int* lamhi , const int* incolam,
                             const amrex::Real* aco , const int* acolo , const int* acohi , const int* incoaco,
                             const amrex::Real* bco0, const int* bco0lo, const int* bco0hi, const int* incobco0,
                             const amrex::Real* bco1, const int* bco1lo, const int* bco1hi, const int* incobco1,
                             const amrex::Real* bco2, const int* bco2lo, const int* bco2hi, const int* incobco2,
                             const int*  gridlo, const int* gridhi,
                             const amrex::Real* dx,  const amrex::Real* alpha, const amrex::Real*  beta);
    

  void ebefnd_getflux_ebco(amrex::Real*       flux, const int* fluxlo, const int* fluxhi, const int* incoflux,
                           const amrex::Real* bcoef, const int* bcoeflo, const int* bcoefhi, const int* incobcoef,
                           const int*  gridlo, const int* gridhi, 
                           const amrex::Real*  beta, const amrex::Real* dx, const int* idir);
#ifdef __cplusplus
}
#endif

#endif
