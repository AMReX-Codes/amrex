#include <AMReX_MultiFab.H>
#include "test_react.H"
#include "test_react_F.H"
#include <iostream>

using namespace amrex;

void do_react(const int* lo, const int* hi,
	      amrex::Real* state, const int* s_lo, const int* s_hi,
	      const int ncomp, const amrex::Real dt)
{
  const int size_x = hi[0]-lo[0]+1;
  const int size_y = hi[1]-lo[1]+1;
  const int size_z = hi[2]-lo[2]+1;
  const int size_state = size_x * size_y * size_z;
  const int neqs = 3;
  const int size_flat = neqs * size_state;

  UserData user_data;
  cudaMallocManaged(&user_data, sizeof(struct CVodeUserData));
  user_data->num_cells = size_state;
  user_data->num_eqs_per_cell = neqs;

  realtype reltol=1.0e-4, time=0.0e0, tout;

  realtype abstol_values[size_flat];
  realtype state_y[size_flat];  

  N_Vector y = NULL, yout=NULL;
  N_Vector abstol = NULL;
  SUNLinearSolver Linsol = NULL;
  void* cvode_mem = NULL;
  int flag;

  // Create NVectors
  y = N_VNew_Cuda(size_flat);
  yout = N_VNew_Cuda(size_flat);
  abstol = N_VNew_Cuda(size_flat);

  // Initialize y, abstol from flattened state
  int nzone = 0;
  for (int i=lo[0]; i<=hi[0]; i++) {
    for (int j=lo[1]; j<=hi[1]; j++) {
      for (int k=lo[2]; k<=hi[2]; k++) {
        for (int n=1; n<=neqs; n++) {
          get_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n, &state_y[nzone*neqs + n - 1]);
	  if (n==1) abstol_values[nzone*neqs + n - 1] = 1.e-8;
	  else if (n==2) abstol_values[nzone*neqs + n - 1] = 1.e-14;
	  else abstol_values[nzone*neqs + n - 1] = 1.e-6;
        }
	nzone++;
      }
    }
  }
  set_nvector_cuda(y, state_y, size_flat);
  set_nvector_cuda(abstol, abstol_values, size_flat);

  // Initialize CVODE
  cvode_mem = CVodeCreate(CV_BDF);
  flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));
  flag = CVodeInit(cvode_mem, fun_rhs, time, y);
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  flag = CVodeSetMaxNumSteps(cvode_mem, 150000);

  // Initialize Linear Solver
  Linsol = SUNSPGMR(y, PREC_NONE, 0);
  flag = CVSpilsSetLinearSolver(cvode_mem, Linsol);
  flag = CVSpilsSetJacTimes(cvode_mem, NULL, fun_jac_times_vec);

  // Do Integration
  time = time + static_cast<realtype>(dt);
  flag = CVode(cvode_mem, time, yout, &tout, CV_NORMAL);
  if (flag != CV_SUCCESS) amrex::Abort("Failed integration");

  // Save Final State
  get_nvector_cuda(yout, state_y, size_flat);
  nzone = 0;
  for (int i=lo[0]; i<=hi[0]; i++) {
    for (int j=lo[1]; j<=hi[1]; j++) {
      for (int k=lo[2]; k<=hi[2]; k++) {
        for (int n=1; n<=neqs; n++) {
          set_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n, &state_y[nzone*neqs + n - 1]);
        }
	nzone++;
      }
    }
  }

  // Free Memory
  N_VDestroy(y);
  N_VDestroy(yout);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(Linsol);

}


static void set_nvector_cuda(N_Vector vec, realtype* data, sunindextype size)
{
  realtype* vec_host_ptr = N_VGetHostArrayPointer_Cuda(vec);
  for (sunindextype i = 0; i < size; i++) {
    vec_host_ptr[i] = data[i];
  }
  N_VCopyToDevice_Cuda(vec);
}


static void get_nvector_cuda(N_Vector vec, realtype* data, sunindextype size)
{
  N_VCopyFromDevice_Cuda(vec);  
  realtype* vec_host_ptr = N_VGetHostArrayPointer_Cuda(vec);
  for (sunindextype i = 0; i < size; i++) {
    data[i] = vec_host_ptr[i];
  }
}


static int fun_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype* ydot_d = N_VGetDeviceArrayPointer_Cuda(ydot);
  realtype* y_d = N_VGetDeviceArrayPointer_Cuda(y);
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  fun_rhs_kernel<<<numBlocks, numThreads>>>(t, y_d, ydot_d,
					    user_data);
  return 0;
}


static int fun_jac_times_vec(N_Vector v, N_Vector Jv, realtype t,
			     N_Vector y, N_Vector fy,
			     void *user_data, N_Vector tmp)
{
  realtype* v_d   = N_VGetDeviceArrayPointer_Cuda(v);
  realtype* Jv_d  = N_VGetDeviceArrayPointer_Cuda(Jv);
  realtype* y_d   = N_VGetDeviceArrayPointer_Cuda(y);
  realtype* fy_d  = N_VGetDeviceArrayPointer_Cuda(fy);
  realtype* tmp_d = N_VGetDeviceArrayPointer_Cuda(tmp);
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  fun_jtv_kernel<<<numBlocks, numThreads>>>(v_d, Jv_d, t,
					    y_d, fy_d,
					    user_data, tmp_d);
  return 0;
}

__global__ static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
				      void *user_data)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
#ifdef CPP_RHS
    ydot[offset] = -.04e0*y[offset] + 1.e4*y[offset+1]*y[offset+2];
    ydot[offset+2] = 3.e7*y[offset+1]*y[offset+1];
    ydot[offset+1] = -ydot[offset]-ydot[offset+2];
#else
    cv_f_rhs_device(&y[offset], &ydot[offset]);
#endif
  }
}


__global__ static void fun_jtv_kernel(realtype* v, realtype* Jv, realtype t,
				      realtype* y, realtype* fy,
				      void* user_data, realtype* tmp)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
#ifdef CPP_RHS
    Jv[offset] = -0.04e0*v[offset] + 1.e4*y[offset+2]*v[offset+1] + 1.e4*y[offset+1]*v[offset+2];
    Jv[offset+2] = 6.0e7*y[offset+1]*v[offset+1];
    Jv[offset+1] = 0.04e0*v[offset] + (-1.e4*y[offset+2]-6.0e7*y[offset+1])*v[offset+1] + (-1.e4*y[offset+1])*v[offset+2];
#else
    cv_f_jtv_device(&v[offset], &Jv[offset], &y[offset], &fy[offset]);
#endif
  }
}
