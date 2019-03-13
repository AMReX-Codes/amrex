#include <cstdio>
#include <iostream>
#include <AMReX.H>
#include <AMReX_Print.H>

#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_Vector.H>
#include <AMReX_IntVect.H>
#include <AMReX_Gpu.H>

// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

using namespace amrex;

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

    // ===================================

    int* n_cpy;
    cudaMallocManaged(&n_cpy, sizeof(int));
    *n_cpy = 2;

    amrex::Print() << "n_cpy before = " << *n_cpy << std::endl << std::endl;

    amrex::Gpu::setLaunchRegion(false);


    amrex::launch_global<<<1,1>>>(
    [=] AMREX_GPU_HOST_DEVICE ()
    {
       *n_cpy = 1;
       printf("n_cpy during CPU = %i\n", *n_cpy);
    });

    int n;
    amrex::Gpu::Device::synchronize();
    cudaMemcpy(&n, n_cpy, sizeof(int), cudaMemcpyDeviceToHost);
    amrex::Print() << "n after CPU = " << n << std::endl << std::endl;

    amrex::Gpu::setLaunchRegion(true);

    amrex::launch_global<<<1,1>>>(
    [=] AMREX_GPU_HOST_DEVICE ()
    {
       *n_cpy = 3;
       printf("n_cpy during GPU = %i\n", *n_cpy);
    });

    amrex::Gpu::Device::synchronize();
    cudaMemcpy(&n, n_cpy, sizeof(int), cudaMemcpyDeviceToHost);
    amrex::Print() << "n after GPU = " << n << std::endl << std::endl;

    cudaFree(n_cpy);
    amrex::Print() << std::endl << "*************************************\n";

    amrex::Finalize();
}
