
#include <AMReX.H>
#include <AMReX_MultiFab.H>

using namespace amrex;

int  main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    {
        Box domain(IntVect(0), IntVect(255));
        BoxArray ba(domain);
        ba.maxSize(64);
        DistributionMapping dm(ba);
        const int ncomp = 1;
        MultiFab mf(ba,dm,ncomp,0);
        mf.setVal(1.0);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::NotInLaunchRegion())
#endif
        {
            FArrayBox tmpfab;
            for (MFIter mfi(mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                // what we used to do:  tmpfab.resize(bx,ncomp);
                Gpu::DeviceFab device_fab(tmpfab, bx, ncomp);
                FArrayBox* tfab = device_fab.fabPtr();
                FArrayBox* mffab = mf.fabPtr(mfi);  // get pointer so that it can be capture by the lambda below

                AMREX_LAUNCH_DEVICE_LAMBDA ( bx, tbx,
                {
                    // tfab, mffab and ncomp are captured by value.
                    // tbx is the box a gpu thread works on, NOT bx
                    tfab->copy(*mffab, tbx, 0, tbx, 0, ncomp);
                    tfab->saxpy(3.0, *mffab);
                    mffab->mult(*tfab);
                });
            }
        }
    }

    amrex::Finalize();
}
