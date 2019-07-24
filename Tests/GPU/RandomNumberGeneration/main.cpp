#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Array.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_ParallelDescriptor.H>
#include <iostream>
#include <ctime>
#include <ratio>
#include <chrono>
using namespace amrex;


void RandomNumGen();

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);
    {
	using namespace std::chrono;

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	RandomNumGen();

	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);

	std::cout << "It took me " << time_span.count() << " seconds.";
	std::cout << std::endl;
     
    }
    amrex::Finalize();

}

void RandomNumGen ()
{
    int N = 1E5/ParallelDescriptor::NProcs();
    int num_rngs = N;
    
#ifdef AMREX_USE_CUDA    
    amrex::InitRandSeedOnDevice(num_rngs);
    amrex::Print() << amrex::Gpu::Device::deviceId() << "\n";
    amrex::Print() << amrex::ParallelDescriptor::MyProc() << "\n";
    
    Gpu::DeviceVector<double> d_xpos(N);
    double *dxpos = d_xpos.dataPtr();    

#else
    amrex::InitRandom(1024UL,1);
#endif

    amrex::Vector<double> hx(N);
    double *hxpos = hx.dataPtr();    



#ifdef AMREX_USE_CUDA
    if ( d_xpos.size() < N){
      CheckSeedArraySizeAndResize(N);
      d_xpos.resize(N);

      dxpos = d_xpos.dataPtr();

    }
#endif 
    if ( hx.size() < N){
      hx.resize(N);

      hxpos = hx.dataPtr();

    }

              
        AMREX_PARALLEL_FOR_1D (N, idx,
        {
#ifdef AMREX_USE_CUDA    
	  dxpos[idx] = amrex::gpusafe_Random();
#else
           hx[idx] = amrex::Random();
#endif
        });
    
#ifdef AMREX_USE_CUDA    
        cudaMemcpy(hxpos,dxpos,sizeof(double)*N,cudaMemcpyDeviceToHost);
#endif


	amrex::Print() << "size of hx = " << hx.size() << std::endl;

	for (int w; w <= hx.size(); w++)
	 {
	   	    amrex::Print() << "hx[" << w << "] = " << hx[w] << std::endl;
	 }
}





