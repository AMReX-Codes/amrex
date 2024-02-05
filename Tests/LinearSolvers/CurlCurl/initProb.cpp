#include "initProb_K.H"

#include "MyTest.H"

using namespace amrex;

void
MyTest::initProb ()
{
    const auto prob_lo = geom.ProbLoArray();
    const auto dx      = geom.CellSizeArray();
    const auto a = alpha;
    const auto b = beta;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhs[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& gbx = mfi.tilebox(IntVect(1),IntVect(1));
        GpuArray<Array4<Real>,3> rhsfab{rhs[0].array(mfi),
                                        rhs[1].array(mfi),
                                        rhs[2].array(mfi)};
        GpuArray<Array4<Real>,3> solfab{solution[0].array(mfi),
                                        solution[1].array(mfi),
                                        solution[2].array(mfi)};
        amrex::ParallelFor(gbx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            actual_init_prob(i,j,k,rhsfab,solfab,prob_lo,dx,a,b);
        });
    }
}
