#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void plusone_cudacpp (amrex::Box const& bx, amrex::FArrayBox& fab)
{
    auto data = fab.array();   // array inside FArrayBox
    amrex::LoopConcurrent(bx,
    [=] (int i, int j, int k)
    {
        data(i,j,k) += 1.0;
    });
}

#endif
