#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Array.H>
#include <AMReX_CudaContainers.H>

using namespace amrex;


void RandomNumGen();

int main (int argc, char* argv[])
{

    amrex::Initialize(argc,argv);
    RandomNumGen();
    amrex::Finalize();

}

void RandomNumGen ()
{
    int N = 1E5;
    int num_rngs = N/20;
    
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
        
        AMREX_PARALLEL_FOR_1D (N, idx,
        {
#ifdef AMREX_USE_CUDA    
           dxpos[idx] = amrex::gpusafe_Random();
#else
           hx[idx] = amrex::gpusafe_Random();
#endif
        });
   
#ifdef AMREX_USE_CUDA    
        cudaMemcpy(hxpos,dxpos,sizeof(double)*N,cudaMemcpyDeviceToHost);
#endif

        for (int i = 0; i < N; i++ )
	  {
	    amrex::Print() << i << " " << hx[i]  << "\n";
        }
     	
	hx.resize(0);
	hx.shrink_to_fit();
}





