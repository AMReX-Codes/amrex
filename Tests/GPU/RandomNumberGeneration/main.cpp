#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Gpu.H>
#include <AMReX_Utility.H>
#include <AMReX_Array.H>
#include <AMReX_GpuContainers.H>
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

#ifdef AMREX_USE_GPU
    amrex::Print() << "Generating random numbers using GPU ";
    amrex::Print() << amrex::Gpu::Device::deviceId() << " on rank ";
    amrex::Print() << amrex::ParallelDescriptor::MyProc() << "\n";

#else
    amrex::InitRandom(1024UL,1);
#endif

    Gpu::HostVector<Real> x_h(Ndraw);
    Gpu::HostVector<Real> y_h(Ndraw);
    Gpu::HostVector<Real> z_h(Ndraw);

    Gpu::DeviceVector<Real> x_d(Ndraw);
    Gpu::DeviceVector<Real> y_d(Ndraw);
    Gpu::DeviceVector<Real> z_d(Ndraw);

    // Test for random numbers. 
    {

        BL_PROFILE_REGION("Draw");

        auto x_d_ptr = x_d.dataPtr();
        auto y_d_ptr = y_d.dataPtr();
        auto z_d_ptr = z_d.dataPtr(); 
        AMREX_PARALLEL_FOR_1D (Ndraw, idx,
        {
            x_d_ptr[idx] = amrex::Random();
            y_d_ptr[idx] = amrex::Random();
            z_d_ptr[idx] = amrex::Random();
        });
   
        Gpu::Device::synchronize();

    }

    Gpu::dtoh_memcpy(x_h.dataPtr(), x_d.dataPtr(), sizeof(Real)*Ndraw);
    Gpu::dtoh_memcpy(y_h.dataPtr(), y_d.dataPtr(), sizeof(Real)*Ndraw);
    Gpu::dtoh_memcpy(z_h.dataPtr(), z_d.dataPtr(), sizeof(Real)*Ndraw);
    
    // Output to check for random-ness
    // for (int i = 0; i < Ndraw; i++ )
    // {
    //     amrex::Print() << i << " " << x_h[i]  << " " << y_h[i] << " " << z_h[i] << "\n";
    // }

    // Test for a subset of threads calling amrex::Random().
    // Testing for a possible hang.
    {
        BL_PROFILE_REGION("Draw2");

        auto x_d_ptr = x_d.dataPtr();
        auto y_d_ptr = y_d.dataPtr();
        auto z_d_ptr = z_d.dataPtr(); 
        AMREX_PARALLEL_FOR_1D (Ndraw, idx,
        {
            if (idx % 2 == 0)
            {
                x_d_ptr[idx] = amrex::Random();
                y_d_ptr[idx] = amrex::Random();
                z_d_ptr[idx] = amrex::Random();
            }
        });
   
        Gpu::Device::synchronize();

    }
}





