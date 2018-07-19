
#ifndef AMREX_ARRAY_H_
#define AMREX_ARRAY_H_

#include <array>
#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>

namespace amrex {

    template <class T, std::size_t N>
    using Array = std::array<T,N>;

    using RealArray = Array<Real, AMREX_SPACEDIM>;
    using IntArray  = Array<int , AMREX_SPACEDIM>;

// If on GPUs, use a simple array wrapper that 
// designates needed functions as __host__ __device__
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
    template <class T, std::size_t N>
    using CudaArray = std::array<T,N>;
#else
#include <AMReX_Managed.H>
    template <class T, std::size_t N>
    struct CudaArray
    {
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE
        inline const T& operator [] (int i) const { return arr[i]; }
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE
        inline const T* data() const { return arr; };
        AMREX_CUDA_HOST AMREX_CUDA_DEVICE
        inline const std::size_t size() const { return N; };

        T arr[N];
    };
#endif

}

namespace amrex
{
    template <class T>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (std::array<T,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }
}

#endif

