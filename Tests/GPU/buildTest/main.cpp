#include <cuda_device_runtime_api.h>

#include <iostream>
#include "MyKernel.H"
#include "MyKernelB.H"

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

int main (int argc, char* argv[])
{
    // ===================================
    // Simple cuda action to make sure all tests have cuda.
    // Allows nvprof to return data.
    int devices;
    cudaGetDeviceCount(&devices);
    std::cout << "Hello from CUDA make test. GPU devices: " << devices << "\n";
    std::cout << "**********************************\n"; 

    {
      double n = 0;
      double *n_d;
      std::cout << "Before plus: " << n << std::endl; 

      cudaMalloc(&n_d, sizeof(double));
      cudaMemcpy(n_d, &n, sizeof(double), cudaMemcpyHostToDevice);

      plusone<<<1, 1>>>(n_d);

      cudaMemcpy(&n, n_d, sizeof(double), cudaMemcpyDeviceToHost);
      cudaFree(n_d);

      std::cout << "After plus: " << n << std::endl; 
    }

      std::cout << "-----------------------------" << std::endl;

    {
      double n = 0;
      double *n_d;
      std::cout << "Before minus: " << n << std::endl; 

      cudaMalloc(&n_d, sizeof(double));
      cudaMemcpy(n_d, &n, sizeof(double), cudaMemcpyHostToDevice);

      minusone<<<1, 1>>>(n_d);

      cudaMemcpy(&n, n_d, sizeof(double), cudaMemcpyDeviceToHost);
      cudaFree(n_d);

      std::cout << "After minus: " << n << std::endl; 
    }
    std::cout << "=============================================" << std::endl;

}
