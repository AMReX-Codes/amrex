#include <AMReX_MultiFab.H>
#include "test_react.H"
#include "test_react_F.H"
#include <iostream>
#include <assert.h>

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

  realtype reltol=1.0e-4, start_time=0.0e0, end_time=0.0e0, tout=0.0e0;

  realtype abstol_values[size_flat];
  realtype state_y[size_flat];

  cuSolver_method LinearSolverMethod = QR;
  const int jac_number_nonzero = 8;
  int csr_row_count[neqs+1] = {1, 4, 7, 9};
  int csr_col_index[jac_number_nonzero] = {1, 2, 3, 1, 2, 3, 2, 3};

  N_Vector y = NULL, yout=NULL;
  N_Vector abstol = NULL;
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
  flag = CVodeInit(cvode_mem, fun_rhs, start_time, y);
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  flag = CVodeSetMaxNumSteps(cvode_mem, 150000);

  // Initialize cuSolver Linear Solver
  flag = cv_cuSolver_SetLinearSolver(cvode_mem, LinearSolverMethod, 1==1, 0);
  flag = cv_cuSolver_CSR_SetSizes(cvode_mem, neqs, jac_number_nonzero, size_state);
  flag = cv_cuSolver_SetJacFun(cvode_mem, &fun_csr_jac);
  flag = cv_cuSolver_SystemInitialize(cvode_mem, &csr_row_count[0], &csr_col_index[0]);

  // Do Integration to dt/2
  end_time = start_time + static_cast<realtype>(dt*0.5);
  std::cout << "Reinitializing CVODE ..." << std::endl;
  CVodeReInit(cvode_mem, start_time, y);
  std::cout << "Integrating to dt/2 ..." << std::endl;
  flag = CVode(cvode_mem, end_time, yout, &tout, CV_NORMAL);
  if (flag != CV_SUCCESS) amrex::Abort("Failed integration");
  else std::cout << "Integrated to dt/2 ..." << std::endl;

  // Do Integration to dt
  start_time = tout;
  end_time = start_time + static_cast<realtype>(dt*0.5);
  std::cout << "Reinitializing CVODE ..." << std::endl;
  CVodeReInit(cvode_mem, start_time, yout);
  std::cout << "Integrating to dt ..." << std::endl;
  flag = CVode(cvode_mem, end_time, yout, &tout, CV_NORMAL);
  if (flag != CV_SUCCESS) amrex::Abort("Failed integration");
  else std::cout << "Integrated to dt ..." << std::endl;

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



int fun_csr_jac(realtype t, N_Vector y, N_Vector fy,
		CV_cuSolver_csr_sys csr_sys, void* user_data)
{
  cudaError_t cuda_status = cudaSuccess;
  cuda_status = cudaGetLastError();

#ifdef DEBUG
  std::cout << "Got CUDA Last Error of: ";
  std::cout << cudaGetErrorString(cuda_status) << std::endl;
#endif

  assert(cuda_status == cudaSuccess);
  realtype* y_d   = N_VGetDeviceArrayPointer_Cuda(y);
  realtype* fy_d  = N_VGetDeviceArrayPointer_Cuda(fy);
  UserData udata = static_cast<CVodeUserData*>(user_data);

#ifdef DEBUG
  std::cout << "num_cells = " << udata->num_cells << std::endl;
  std::cout << "size_per_subsystem = " << csr_sys->size_per_subsystem << std::endl;
  std::cout << "csr_number_nonzero = " << csr_sys->csr_number_nonzero << std::endl;
  std::cout << "number_subsystems = " << csr_sys->number_subsystems << std::endl;
#endif

  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  fun_csr_jac_kernel<<<numBlocks, numThreads>>>(t, y_d, fy_d, csr_sys->d_csr_values,
						user_data,
						csr_sys->size_per_subsystem,
						csr_sys->csr_number_nonzero,
						csr_sys->number_subsystems);
  cuda_status = cudaDeviceSynchronize();

#ifdef DEBUG
  std::cout << "Got CUDA Synchronize return message of: ";
  std::cout << cudaGetErrorString(cuda_status) << std::endl;
#endif

  assert(cuda_status == cudaSuccess);
  return 0;
}


__global__ static void fun_csr_jac_kernel(realtype t, realtype* y, realtype* fy,
					  realtype* csr_jac, void* user_data,
					  const int size, const int nnz, const int nbatched)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;

  if (tid < nbatched) {
    int jac_offset = tid * nnz;
    int y_offset = tid * size;

    realtype* csr_jac_cell = csr_jac + jac_offset;
    realtype* actual_y = y + y_offset;

    csr_jac_cell[0] = -0.04e0;
    csr_jac_cell[1] = 1.e4*actual_y[2];
    csr_jac_cell[2] = 1.e4*actual_y[1];
    csr_jac_cell[6] = 6.0e7*actual_y[1];
    csr_jac_cell[3] = 0.04e0;
    csr_jac_cell[4] = -csr_jac_cell[1]-csr_jac_cell[6];
    csr_jac_cell[5] = -csr_jac_cell[2];
    csr_jac_cell[7] = 0.0e0;
  }
}
