//#include <cuda_device_runtime_api.h>
//#include <thrust/device_vector.h>

#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Geometry.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_BaseFab.H>
#include <AMReX_Gpu.H>

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

extern "C"
{

  AMREX_GPU_DEVICE
  long fort_return_test (long* mallocd);

}
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

using namespace amrex;


int main (int argc, char* argv[])
{
    std::cout << "**********************************\n";
    amrex::Initialize(argc, argv);

    // ===================================
    // Simple cuda action to make sure all tests have cuda.
    // Allows nvprof to return data.
    int devices;
    cudaGetDeviceCount(&devices);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
    amrex::Print() << "**********************************\n"; 

    // Malloc
    {

      long n = 0;
      long *n_d;
      cudaMalloc(&n_d, size_t(sizeof(long)));
      cudaMemcpy(n_d, &n, sizeof(long), cudaMemcpyHostToDevice);

      std::cout << "n before = " << n << std::endl << std::endl;

      amrex::launch_global<<<1,1>>>(
      [=] AMREX_GPU_DEVICE () mutable
      {

         long returned = 0;
         printf("before: return = %li, mallocd = %li\n", returned, *n_d);

         returned = fort_return_test(n_d);

         printf(" after: return = %li, mallocd = %li\n\n", returned, *n_d);

      });

      cudaMemcpy(&n, n_d, sizeof(long), cudaMemcpyDeviceToHost);
      cudaFree(n_d);
      amrex::Gpu::Device::synchronize();

      std::cout << "n after = " << n << std::endl; 
    }
    amrex::Print() << std::endl << "=============================================" << std::endl << std::endl;

    amrex::Finalize();
}
