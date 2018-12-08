#ifndef MY_KERNEL_H_
#define MY_KERNEL_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE void plusone_cudacpp (amrex::Box const& bx, amrex::FArrayBox& fab);

#endif
