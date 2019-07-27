#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Array.H>
#include <AMReX_CudaContainers.H>
#include <AMReX_ParmParse.H>

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

    ParmParse pp;

    int Nstates;
    int Ndraw;

    pp.get("num_states", Nstates);
    pp.get("num_draw", Ndraw);    

#ifdef AMREX_USE_CUDA    
    //    amrex::InitRandSeedOnDevice(Nstates);  // This will set the number of RNGs

    amrex::Print() << "Generating random numbers using GPU ";
    amrex::Print() << amrex::Gpu::Device::deviceId() << " on rank ";
    amrex::Print() << amrex::ParallelDescriptor::MyProc() << "\n";
#else
    amrex::InitRandom(1024UL,1);
#endif

    Gpu::DeviceVector<Real> x(Ndraw);
    Gpu::DeviceVector<Real> y(Ndraw);
    Gpu::DeviceVector<Real> z(Ndraw);

    {

    BL_PROFILE_REGION("Draw");

    auto x_ptr = x.dataPtr();
    auto y_ptr = y.dataPtr();
    auto z_ptr = z.dataPtr(); 
    AMREX_PARALLEL_FOR_1D (Ndraw, idx,
    {
        x_ptr[idx] = amrex::Random();
        y_ptr[idx] = amrex::Random();
        z_ptr[idx] = amrex::Random();
    });
   
    Gpu::Device::synchronize();

    }
 
    // for (int i = 0; i < Ndraw; i++ )
    // {
    //     amrex::Print() << i << " " << x[i]  << " " << y[i] << " " << z[i]<< "\n";
    // }
}





