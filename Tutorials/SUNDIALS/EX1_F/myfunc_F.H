#ifndef MYFUNC_F_H_
#define MYFUNC_F_H_

#ifdef __cplusplus
extern "C"
{
#endif


  void integrate_ode(amrex_real* phi,
                     const int* lo,
                     const int* hi,
                     const int* cvode_meth,
                     const int* cvode_itmeth);

#ifdef __cplusplus
}
#endif

#endif
