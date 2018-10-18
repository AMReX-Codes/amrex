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

AMREX_CUDA_DEVICE
inline float atomicMax(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_i, assumed,
            __float_as_int(fmaxf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

AMREX_CUDA_DEVICE
inline double atomicMax(double* address, double val)
{
    unsigned long long int* address_as_ull = 
                  (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(fmax(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}


AMREX_CUDA_DEVICE
inline float atomicMin(float* address, float val)
{
    int* address_as_i = (int*) address;
    int old = *address_as_i, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_i, assumed,
            __float_as_int(fminf(val, __int_as_float(assumed))));
    } while (assumed != old);
    return __int_as_float(old);
}

AMREX_CUDA_DEVICE
inline double atomicMin(double* address, double val)
{
    unsigned long long int* address_as_ull = 
                  (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
            __double_as_longlong(fmin(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
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
