#ifndef TEST_REACT_H
#define TEST_REACT_H

#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */

#ifdef CUDA_CVODE

#include <nvector/nvector_cuda.h>         /* access to serial N_Vector                */

#ifdef CVODE_CUSOLVER

#include <cusolver/cvode_cusolver_spqr.h>    /* access to cuSolver interface             */

#else

#include <sunlinsol/sunlinsol_spgmr.h>    /* access to SPGMR SUNLinearSolver          */
#include <cvode/cvode_spils.h>            /* access to CVSpils interface              */

#endif

#else

#include <nvector/nvector_serial.h>         /* access to serial N_Vector                */
#include <sunmatrix/sunmatrix_dense.h>        /* access to dense SUNMatrix                */
#include <sunlinsol/sunlinsol_dense.h>        /* access to dense SUNLinearSolver          */
#include <cvode/cvode_direct.h>           /* access to CVDls interface                */

#endif

#include <sundials/sundials_types.h>         /* definition of realtype                   */
#include <sundials/sundials_math.h>          /* contains the macros ABS, SUNSQR, and EXP */

typedef struct CVodeUserData {
  int num_cells;
  int num_eqs_per_cell;
}* UserData;

void main_main();

void do_react(const int* lo, const int* hi,
              amrex::Real* state, const int* s_lo, const int* s_hi,
              const int ncomp, const amrex::Real dt);

static int fun_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#ifdef CUDA_CVODE

__global__ static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
				      void *user_data);

#ifdef CVODE_CUSOLVER

int fun_csr_jac(realtype t, N_Vector y, N_Vector fy,
		CV_cuSolver_csr_sys csr_sys, void* user_data);

__global__ static void fun_csr_jac_kernel(realtype t, realtype* y, realtype* fy,
					  realtype* csr_jac, void* user_data,
					  const int size, const int nnz, const int nbatched);

#else

static int fun_jac_times_vec(N_Vector v, N_Vector Jv, realtype t,
			     N_Vector y, N_Vector fy,
			     void *user_data, N_Vector tmp);

__global__ static void fun_jtv_kernel(realtype* v, realtype* Jv, realtype t,
				      realtype* u, realtype* fu,
				      void* user_data, realtype* tmp);

#endif

static void set_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

static void get_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

#else

static int fun_jac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                   void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

#endif

#endif
