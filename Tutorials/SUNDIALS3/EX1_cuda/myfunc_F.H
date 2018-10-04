#ifndef MYFUNC_F_H_
#define MYFUNC_F_H_

#ifdef __cplusplus
extern "C"
{
#endif


  void FSetIC(amrex_real* phi,
                     const int* lo,
                     const int* hi,
	             double* dptr);

  void FSetIC_mfab(amrex_real* phi,
                     const int* lo,
		   const int* hi);

  void FSetSol(amrex_real* mf,
                     const int* lo,
                     const int* hi,
	             double* dptr);

  void RhsFn(realtype t, N_Vector u, N_Vector udot, void *user_data);

#ifdef __cplusplus
}
#endif

#endif
