
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

    // Multiple ways of kernel launch

    // (1) C++, AMREX_FOR_3D
    {
        BL_PROFILE("1-amrex_for_3d");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> fab = mf.array(mfi);
            // loop over bx
            AMREX_FOR_3D ( bx, i, j, k,
            {
                fab(i,j,k) += 1.;
                // or fab(i,j,k,0) += 1.;
            });
        }
    }

    // (2) C++, AMREX_FOR_4D
    {
        BL_PROFILE("2-amrex_for_4d");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> fab = mf.array(mfi);
            int ncomp = mf.nComp();
            // loop over bx and component.
            AMREX_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                fab(i,j,k,n) += 1.;
            });
        }
    }

    // (3) C++, AMREX_FOR_1D
    {
        BL_PROFILE("3-amrex_for_1d");
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = mf[mfi];
            Real* p = fab.dataPtr();
            const long nitems = fab.box().numPts() * mf.nComp();
            // Enough threads are launched to work over nitems.
            // This only works on a contiguous chunk of memory.
            AMREX_FOR_1D ( nitems, idx,
            {
                p[idx] += 1.;
            });
        }
    }

    // (4) C++, AMREX_LAUNCH_DEVICE_LAMBDA, Capture Array4
    {
        BL_PROFILE("4-amrex_launch_array4");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> fab = mf.array(mfi);
            // Enough threads are launched to work over bx,
            // and tbx is a thread's work box
            AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
            {
                // Array4<Real> fab is captured
                const auto lo = amrex::lbound(tbx);
                const auto hi = amrex::ubound(tbx);
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) {
                    fab(i,j,k) += 1.;
                }}}
            });
        }
    }

    // (5) C++, AMREX_LAUNCH_DEVICE_LAMBDA, Capture FArrayBox
    {
        BL_PROFILE("5-amrex_launch_fab");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Use fabPtr function to get a managed pointer to fab
            FArrayBox* fab = mf.fabPtr(mfi);
            // Enough threads are launched to work over bx,
            // and tbx is a thread's work box
            AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
            {
                // FArrayBox* fab is captured
                plusone_cudacpp(tbx, *fab);
            });
        }
    }

#ifdef AMREX_USE_CUDA_FORTRAN
    // (6) launch CUDA Fortran kernel to add 1 if supported
    {
        BL_PROFILE("6-fortran");
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
    }
#endif

#ifdef AMREX_USE_ACC
    // (7) launch OpenACC kernel to add 1 if supported
    {
        BL_PROFILE("7-acc");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            // Use operator[] to get a reference to fab.
            // The reference points to host memory, but the data pointer inside is managed.
            FArrayBox& fab = mf[mfi];
            plusone_acc(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_ANYD(fab));
        }
    }
#endif

#ifdef AMREX_OMP_OFFLOAD
    // (8) launch OpenOMP kernel to add 1 if supported
    {
        BL_PROFILE("8-omp");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            // Use operator[] to get a reference to fab.
            // The reference points to host memory, but the data pointer inside is managed.
            FArrayBox& fab = mf[mfi];
            plusone_omp(BL_TO_FORTRAN_BOX(bx),
                        BL_TO_FORTRAN_ANYD(fab));
        }
    }
#endif
}
