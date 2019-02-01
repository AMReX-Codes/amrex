
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

    // (1) C++, 3D amrex::ParallelFor function
    {
        BL_PROFILE("1-amrex_for_3d_func");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> const& fab = mf.array(mfi);
            // loop over bx
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                fab(i,j,k) += 1.;
            });
        }
    }

    // (2) C++, AMREX_PARALLEL_FOR_3D macro
    {
        BL_PROFILE("2-amrex_for_3d-macro");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> const& fab = mf.array(mfi);
            // loop over bx
            AMREX_PARALLEL_FOR_3D ( bx, i, j, k,
            {
                fab(i,j,k) += 1.;
                // or fab(i,j,k,0) += 1.;
            });
        }
    }

    // (3) C++, 4D amrex::ParallelFor function
    {
        BL_PROFILE("3-amrex_for_4d-func");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> const& fab = mf.array(mfi);
            int ncomp = mf.nComp();
            // loop over bx and component.
            amrex::ParallelFor(bx, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n)
            {
                fab(i,j,k,n) += 1.;
            });
        }
    }

    // (4) C++, AMREX_PARALLEL_FOR_4D macro
    {
        BL_PROFILE("4-amrex_for_4d-macro");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> const& fab = mf.array(mfi);
            int ncomp = mf.nComp();
            // loop over bx and component.
            AMREX_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                fab(i,j,k,n) += 1.;
            });
        }
    }

    // (5) C++, 1D amrex::ParallelFor function
    {
        BL_PROFILE("5-amrex_for_1d-func");
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = mf[mfi];
            Real* AMREX_RESTRICT p = fab.dataPtr();
            const long nitems = fab.box().numPts() * mf.nComp();
            // Enough threads are launched to work over nitems.
            // This only works on a contiguous chunk of memory.
            amrex::ParallelFor(nitems,
            [=] AMREX_GPU_DEVICE (long idx)
            {
                p[idx] += 1.;
            });
        }
    }

    // (6) C++, AMREX_PARALLEL_FOR_1D macro
    {
        BL_PROFILE("6-amrex_for_1d-macro");
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = mf[mfi];
            Real* AMREX_RESTRICT p = fab.dataPtr();
            const long nitems = fab.box().numPts() * mf.nComp();
            // Enough threads are launched to work over nitems.
            // This only works on a contiguous chunk of memory.
            AMREX_PARALLEL_FOR_1D ( nitems, idx,
            {
                p[idx] += 1.;
            });
        }
    }

    // (7) C++, amrex::launch function, Capture Array4
    {
        BL_PROFILE("7-amrex_launch_array4-func");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> const& fab = mf.array(mfi);
            // Enough threads are launched to work over bx,
            // and tbx is a thread's work box
            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                // Array4<Real> fab is captured
                const auto lo = amrex::lbound(tbx);
                const auto hi = amrex::ubound(tbx);
                for (int k = lo.z; k <= hi.z; ++k) {
                for (int j = lo.y; j <= hi.y; ++j) {
                AMREX_PRAGMA_SIMD
                for (int i = lo.x; i <= hi.x; ++i) {
                    fab(i,j,k) += 1.;
                }}}  
            });
        }
    }

    // (8) C++, AMREX_LAUNCH_DEVICE_LAMBDA macro
    {
        BL_PROFILE("8-amrex_launch_array4-macro");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Get Array4 object
            Array4<Real> const& fab = mf.array(mfi);
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

    // (9) C++, amrex::launch function, Capture FArrayBox*
    {
        BL_PROFILE("9-amrex_launch_fab-func");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Tiling is off in case of gpu.
            // In that case, tilebox simply return validbox
            const Box& bx = mfi.tilebox();
            // Use fabPtr function to get a managed pointer to fab
            FArrayBox* fab = mf.fabPtr(mfi);
            // Enough threads are launched to work over bx,
            // and tbx is a thread's work box
            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                // FArrayBox* fab is captured
                plusone_cudacpp(tbx, *fab);
            });
        }
    }

    // (10) C++, AMREX_LAUNCH_DEVICE_LAMBDA macro, Capture FArrayBox*
    {
        BL_PROFILE("10-amrex_launch_fab-macro");
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

#if defined(AMREX_USE_CUDA_FORTRAN) || !defined(AMREX_USE_GPU)
    // (11) launch CUDA Fortran kernel to add 1
    {
        BL_PROFILE("11-fortran-func");
        for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            FArrayBox* fab = mf.fabPtr(mfi);
            amrex::launch(bx,
            [=] AMREX_GPU_DEVICE (Box const& tbx)
            {
                plusone_cudafort(BL_TO_FORTRAN_BOX(tbx),
                                 BL_TO_FORTRAN_ANYD(*fab));
            });
        }
    }

    // (12) launch CUDA Fortran kernel to add 1
    {
        BL_PROFILE("12-fortran-macro");
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
    // (13) launch OpenACC kernel to add 1 if supported
    {
        BL_PROFILE("13-acc");
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
    // (14) launch OpenOMP kernel to add 1 if supported
    {
        BL_PROFILE("14-omp");
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
