#ifndef AMREX_CUDA_UTILITY_H_
#define AMREX_CUDA_UTILITY_H_

#include <AMReX_Managed.H>
#include <AMReX_Box.H>

namespace amrex {

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    Box getThreadBox(const Box& bx);

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    Box getThreadBox(const Box& bx, const IntVect& typ);

    AMREX_CUDA_HOST AMREX_CUDA_DEVICE
    void getThreadIndex(int &index, int &size, const int num_particles);

}

#endif
