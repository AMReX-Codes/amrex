#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Array.H>
#include <AMReX_CudaContainers.H>
#include <curand_kernel.h>

using namespace amrex;

void RandomNumGen ();
AMREX_GPU_DEVICE void InitializePositions(double *d_xpos_ptr, double *d_ypos_ptr, double *d_zpos_ptr, int idx, curandState *d_RandStates_ptr);
AMREX_GPU_DEVICE void SetRandSeedOnDevice(curandState *d_RandStates_ptr, int idx);

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    RandomNumGen();
    amrex::Finalize();
}

void RandomNumGen ()
{
    int N = 1E5;

    Gpu::DeviceVector<curandState> dev_RandStates_Seed(N*2);
    curandState *d_RandStates_ptr = dev_RandStates_Seed.dataPtr();

    Gpu::DeviceVector<double> d_xpos(N); 
    Gpu::DeviceVector<double> d_ypos(N); 
    Gpu::DeviceVector<double> d_zpos(N); 

    double *d_xpos_ptr = d_xpos.dataPtr();
    double *d_ypos_ptr = d_ypos.dataPtr();
    double *d_zpos_ptr = d_zpos.dataPtr();

    double *h_xpos = new double[N];
    double *h_ypos = new double[N];
    double *h_zpos = new double[N];
 
    AMREX_PARALLEL_FOR_1D (N, idx,
    {
       SetRandSeedOnDevice(d_RandStates_ptr, idx);
    });

    for (int i=0; i<10; i++){
        AMREX_PARALLEL_FOR_1D (N, idx,
        {
           InitializePositions(d_xpos_ptr,d_ypos_ptr,d_zpos_ptr,idx,d_RandStates_ptr);
        });
   
        cudaMemcpy( h_xpos, d_xpos_ptr, sizeof(double)*N, cudaMemcpyDeviceToHost); 
        cudaMemcpy( h_ypos, d_ypos_ptr, sizeof(double)*N, cudaMemcpyDeviceToHost); 
        cudaMemcpy( h_zpos, d_zpos_ptr, sizeof(double)*N, cudaMemcpyDeviceToHost); 

        for ( int i = 0; i < N; i++ ){
            amrex::Print() <<  i << " " << h_xpos[i] << " " << h_ypos[i] << " " << h_zpos[i] << "\n" ;
        }
    }

    delete[] h_xpos;
    delete[] h_ypos;
    delete[] h_zpos;
    h_xpos = NULL;
    h_ypos = NULL;
    h_zpos = NULL;
}


AMREX_GPU_DEVICE void InitializePositions(double *d_xpos_ptr, double *d_ypos_ptr, double *d_zpos_ptr, int idx, curandState *d_RandStates_ptr)
{
    double loc_rand; 

    loc_rand = curand_uniform_double( &d_RandStates_ptr[idx] );
    d_xpos_ptr[idx] = loc_rand;

    loc_rand = curand_uniform_double( &d_RandStates_ptr[idx] );
    d_ypos_ptr[idx] = loc_rand;

    loc_rand = curand_uniform_double( &d_RandStates_ptr[idx] );
    d_zpos_ptr[idx] = loc_rand;
}

AMREX_GPU_DEVICE void SetRandSeedOnDevice(curandState *d_RandStates_ptr, int idx)
{
    unsigned long seed = idx + 10*idx;
    curand_init( seed, seed, 0, &d_RandStates_ptr[idx]);
}



