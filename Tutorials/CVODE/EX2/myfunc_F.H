#ifndef MYFUNC_F_H__
#define MUFUNC_F_H__

#ifdef __cplusplus
extern "C"
{
#endif


  void integrate_ode_no_jac(amrex_real* phi,
                            const int* lo,
                            const int* hi,
                            const int* Ncomp,
                            const int* cvode_meth,
                            const int* cvode_itmeth);

  void integrate_ode_with_jac(amrex_real* phi,
                              const int* lo,
                              const int* hi,
                              const int* Ncomp,
                              const int* cvode_meth,
                              const int* cvode_itmeth);

#ifdef __cplusplus
}
#endif

#endif
