#include "Poisson.H"

using namespace amrex;

void InitData (MultiFab& State)
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(State,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox();
        const Array4<Real>& q = State.array(mfi);
        amrex::ParallelFor(bx,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (i==70 and j==70 and k==0) {
                q(i,j,k) = 1.0;
            } else {
                q(i,j,k) = 0.0;
            }
        });
    }
}
