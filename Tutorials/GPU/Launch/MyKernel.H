#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE
AMREX_INLINE
void plusone_cudacpp (amrex::Box const& bx, amrex::FArrayBox& fab)
{
    const auto len = amrex::length(bx);  // length of box
    const auto lo  = amrex::lbound(bx);  // lower bound of box
    const auto data = fab.view(lo);  // a view starting from lo

    for         (int k = 0; k < len.z; ++k) {
        for     (int j = 0; j < len.y; ++j) {
            // We know this is safe for simd on cpu.  So let's give compiler some help.
            AMREX_PRAGMA_SIMD
#pragma unroll(1)
            for (int i = 0; i < len.x; ++i) {
                data(i,j,k) += 1.0;
            }
        }
    }
}

#endif
