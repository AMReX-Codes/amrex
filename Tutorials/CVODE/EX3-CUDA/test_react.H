#ifndef TEST_REACT_H
#define TEST_REACT_H

#include "cvode/cvode.h"                  /* prototypes for CVODE fcts., consts.      */
#include "nvector/nvector_cuda.h"         /* access to serial N_Vector                */
#include "sunmatrix/sunmatrix_dense.h"        /* access to dense SUNMatrix                */
#include "sunlinsol/sunlinsol_dense.h"        /* access to dense SUNLinearSolver          */
#include "cvode/cvode_direct.h"           /* access to CVDls interface                */
#include "sundials/sundials_types.h"         /* definition of realtype                   */
#include "sundials/sundials_math.h"          /* contains the macros ABS, SUNSQR, and EXP */

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

static int fun_jac_times_vec(N_Vector v, N_Vector Jv, realtype t,
			     N_Vector y, N_Vector fy,
			     void *user_data, N_Vector tmp);

__global__ static void fun_jtv_kernel(realtype* v, realtype* Jv, realtype t,
				      realtype* u, realtype* fu,
				      void* user_data, realtype* tmp);

static void set_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

static void get_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

#else

static int fun_jac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                   void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

#endif

#endif
