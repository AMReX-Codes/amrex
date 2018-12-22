
#include <AMReX.H>
#include <AMReX_MultiFab.H>

#include "MyKernel.H"
#include "MyKernel_F.H"

using namespace amrex;

void main_main();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();
}


void main_main ()
{
    BoxArray ba;
    {
        Box domain_box(IntVect(0), IntVect(127));
        ba.define(domain_box);
        ba.maxSize(64);
    }

    MultiFab mf(ba,DistributionMapping{ba},1,0);
    mf.setVal(0.0);

    // launch CUDA C++ kernel to add 1
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        // Tiling is off in case of gpu.
        // In that case, tilebox simply return validbox
        const Box& bx = mfi.tilebox();
        // Use fabPtr function to get a managed pointer to fab
        FArrayBox* fab = mf.fabPtr(mfi);
        // Enough threads are launched to cover bx, and tbx is a thread's work box
        AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
        {
            plusone_cudacpp(tbx, *fab);
        });
    }

    // launch C++ kernel on vector to add 1
    {
        int size = 100;
        amrex::Gpu::ManagedVector<int> ones(size, 0); 
        const auto data = ones.dataPtr();
        AMREX_LAUNCH_DEVICE_LAMBDA(size, iter,
        {
            data[iter] = data[iter] + 1;
        });

        Gpu::Device::synchronize();
    }

#ifdef AMREX_USE_CUDA_FORTRAN
    // launch CUDA Fortran kernel to add 1 if supported
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        FArrayBox* fab = mf.fabPtr(mfi);
        AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
        {
            plusone_cudafort(BL_TO_FORTRAN_BOX(tbx),
                             BL_TO_FORTRAN_ANYD(*fab));
        });
    }
#endif

#ifdef AMREX_USE_ACC
    // launch OpenACC kernel to add 1 if supported
    for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        // Use operator[] to get a reference to fab.
        // The reference points to host memory, but the data pointer inside is managed.
        FArrayBox& fab = mf[mfi];
        plusone_acc(BL_TO_FORTRAN_BOX(bx),
                    BL_TO_FORTRAN_ANYD(fab));
    }
#endif

}
