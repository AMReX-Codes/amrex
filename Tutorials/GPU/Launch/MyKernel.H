#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE
AMREX_INLINE
void plusone_cudacpp (amrex::Box const& bx, amrex::FArrayBox& fab)
{
    const auto lo  = amrex::lbound(bx);  // lower bound of box
    const auto hi  = amrex::ubound(bx);  // upper bound of box
    auto data = fab.array();             // array inside FArrayBox

    for (int k = lo.z; k <= hi.z; ++k) {
    for (int j = lo.y; j <= hi.y; ++j) {
    // We know this is safe for simd on cpu.  So let's give compiler some help.
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
        data(i,j,k) += 1.0;
    }}}
}

#endif
