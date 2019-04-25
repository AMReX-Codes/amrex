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

using namespace amrex;

template <class T, typename std::enable_if<std::is_pod<T>::value,int>::type = 0 > 
struct PinnedData
{
   PinnedData ()
   {
     cudaMallocHost(&d_d, std::size_t(sizeof(T)));
   }

   PinnedData (T const & h_d)
   : PinnedData()
   {
     *d_d = h_d;
   }

   ~PinnedData ()
   {
     cudaFreeHost(d_d);
   }

   T* devicePtr() &
   {
     return d_d;
   }

   T const * devicePtr() const&
   {
     return d_d;
   }
/*
   T hostValue () const
   {
     T t;
//     cudaMemcpy(&t, d_d, sizeof(T), cudaMemDeviceToHost);
     return t; 
   }
*/
   void updateDevice(const T &t)
   {
//       cudaMemcpy(d_d, &t, sizeof(T), cudaMemHostToDevice);
   }

   T* data() && = delete; 
   PinnedData(PinnedData const &) = delete;
   PinnedData(PinnedData &&) = delete;
   void operator = (PinnedData const &) = delete;
   void operator = (PinnedData &&) = delete; 

   private:
   T* d_d = nullptr;
};


int main (int argc, char* argv[])
{
    std::cout << "**********************************\n";
//    Device::initialize_device();
    amrex::Initialize(argc, argv);
    amrex::Print() << "amrex::Initialize complete." << "\n";

    // ===================================
    // Simple cuda action to make sure all tests have cuda.
    // Allows nvprof to return data.
    int devices;
    cudaGetDeviceCount(&devices);
    amrex::Print() << "Hello world from AMReX version " << amrex::Version() << ". GPU devices: " << devices << "\n";
    amrex::Print() << "**********************************\n"; 

    amrex::Cuda::ExecutionConfig simple_config(1,1);

    // Malloc
    {

      amrex::Print() << "Malloc version" << std::endl;
      amrex::Print() << "=======================" << std::endl;

      int n = 101;
      int *n_d;
      cudaMalloc(&n_d, size_t(sizeof(int)));
      cudaMemcpy(n_d, &n, sizeof(int), cudaMemcpyHostToDevice);

      std::cout << "n before = " << n << std::endl;

      amrex::launch_global<<<1,1>>>(
      [=] AMREX_GPU_DEVICE () mutable
      {
         *n_d = *n_d / 10;
         printf("n during = %i\n", *n_d);
      });

      cudaMemcpy(&n, n_d, sizeof(int), cudaMemcpyDeviceToHost);
      cudaFree(n_d);
      amrex::Gpu::Device::synchronize();

      std::cout << "n after = " << n << std::endl << std::endl;
    }

    // Managed
    {
      amrex::Print() << "Managed version" << std::endl;
      amrex::Print() << "=======================" << std::endl;

      int *n;
      cudaMallocManaged(&n, size_t(sizeof(int)));
      *n = 101;

      std::cout << "n before = " << *n << std::endl;

      amrex::launch_global<<<1,1>>>(
      [=] AMREX_GPU_DEVICE () mutable
      {
         *n = *n / 10;
         printf("n during = %i\n", *n);
      });

      amrex::Gpu::Device::synchronize();

      std::cout << "n after = " << *n << std::endl << std::endl;

      cudaFree(n);
    }



    // Pinned 
    {
      amrex::Print() << "Pinned version" << std::endl;
      amrex::Print() << "=======================" << std::endl;

      int *n;
      cudaMallocHost(&n, size_t(sizeof(int)));
      *n = 101;

      std::cout << "n before = " << *n << std::endl;

      amrex::launch_global<<<1,1>>>(
      [=] AMREX_GPU_DEVICE () mutable
      {
         *n = *n / 10;
         printf("n during = %i\n", *n);
      });

      amrex::Gpu::Device::synchronize();

      std::cout << "n after = " << *n << std::endl << std::endl;

      cudaFreeHost(n);
    }

    // Pinned Data 
    {
      amrex::Print() << "Macro version" << std::endl;
      amrex::Print() << "=======================" << std::endl;

      PinnedData<int> n;
      PinnedData<int[5]> n5;
      int *p  = n.devicePtr();
      auto q  = n5.devicePtr();
      std::cout << "n before = " << *p << std::endl;

      amrex::launch_global<<<1,1>>>(
      [=] AMREX_GPU_DEVICE () mutable
      {
         (*q)[0] = 4;
         (*q)[1] = 3;
         (*q)[2] = 2;
         (*q)[3] = 1;
         (*q)[4] = 0;

         *p = *p / 10;
         printf("n during = %i\n", *p);
      });

      amrex::Gpu::Device::synchronize();

      for (int i=0;i<5;i++)
      {
        std::cout << "q[" << i << "] = " << (*q)[i] << std::endl;
      }

      std::cout << "n after = " << *p << std::endl << std::endl;
    }

    amrex::Print() << std::endl << "=============================================" << std::endl << std::endl;

    amrex::Finalize();
}
