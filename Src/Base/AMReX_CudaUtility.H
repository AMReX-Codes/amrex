#ifndef AMREX_CUDA_UTILITY_H_
#define AMREX_CUDA_UTILITY_H_

#include <AMReX_CudaQualifiers.H>

namespace amrex {

namespace CudaAtomic {

    template<class T>
    AMREX_CUDA_HOST_DEVICE
    void Add (T* sum, T value)
    {
#if defined(__CUDA_ARCH__)
        atomicAdd(sum, value);
#else
        *sum += value;
#endif
    }

    template<class T>
    AMREX_CUDA_HOST_DEVICE
    void Min (T* m, T value)
    {
#if defined(__CUDA_ARCH__)
        atomicMin(m, value);
#else
        *m = (*m) < value ? (*m) : value;
#endif
    }

    template<class T>
    AMREX_CUDA_HOST_DEVICE
    void Max (T* m, T value)
    {
#if defined(__CUDA_ARCH__)
        atomicMax(m, value);
#else
        *m = (*m) > value ? (*m) : value;
#endif
    }
}

}

#endif
