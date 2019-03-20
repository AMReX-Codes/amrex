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

#ifdef AMREX_USE_CUDA    
    amrex::InitRandSeedOnDevice(N);

    Gpu::DeviceVector<double> d_xpos(N);
    Gpu::DeviceVector<double> d_ypos(N);
    Gpu::DeviceVector<double> d_zpos(N);

    double *dxpos = d_xpos.dataPtr();    
    double *dypos = d_ypos.dataPtr();    
    double *dzpos = d_zpos.dataPtr();    
#else
    amrex::InitRandom(1024UL,1);
#endif
 
    amrex::Vector<double> hx(N);
    amrex::Vector<double> hy(N);
    amrex::Vector<double> hz(N);
    double *hxpos = hx.dataPtr();    
    double *hypos = hy.dataPtr();    
    double *hzpos = hz.dataPtr();    


    int timesteps = 10; // just for example
    for (int i=0; i<timesteps; i++)
    {
        
        AMREX_PARALLEL_FOR_1D (N, idx,
        {
#ifdef AMREX_USE_CUDA    
           dxpos[idx] = amrex::Random();
           dypos[idx] = amrex::Random();
           dzpos[idx] = amrex::Random();
#else
           hx[idx] = amrex::Random();
           hy[idx] = amrex::Random();
           hz[idx] = amrex::Random();
#endif
        });
   
#ifdef AMREX_USE_CUDA    
        cudaMemcpy(hxpos,dxpos,sizeof(double)*N,cudaMemcpyDeviceToHost);
        cudaMemcpy(hypos,dypos,sizeof(double)*N,cudaMemcpyDeviceToHost);
        cudaMemcpy(hzpos,dzpos,sizeof(double)*N,cudaMemcpyDeviceToHost);
#endif

        for (int i = 0; i < N; i++ )
        {
           amrex::Print() << i << " " << hx[i]  << " " << hy[i] << " " << hz[i]<< "\n";
        }

    }

#ifdef AMREX_USE_CUDA    
    amrex::DeallocateRandomSeedDevArray();
#endif

}





