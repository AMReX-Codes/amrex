
#ifndef AMREX_ARRAY_H_
#define AMREX_ARRAY_H_

#include <array>
#include <memory>
#include <utility>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>

namespace amrex {

    template <class T, std::size_t N>
    using Array = std::array<T,N>;

    using RealArray = Array<Real, AMREX_SPACEDIM>;
    using IntArray  = Array<int , AMREX_SPACEDIM>;

    struct Dim3 { int x; int y; int z; };

    template <typename T>
    struct Array4
    {
        T* AMREX_RESTRICT p;
        Dim3 begin;
        Dim3 end;
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T const& operator() (int i, int j, int k) const {
            AMREX_ASSERT(i>=begin.x && i<end.x && j>=begin.y && j<end.y && k>=begin.z && k<end.z);
            return p[i-begin.x+(j-begin.y)*(end.x-begin.x)+(k-begin.z)*(end.x-begin.x)*(end.y-begin.y)];
        }
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) {
            AMREX_ASSERT(i>=begin.x && i<end.x && j>=begin.y && j<end.y && k>=begin.z && k<end.z);
            return p[i-begin.x+(j-begin.y)*(end.x-begin.x)+(k-begin.z)*(end.x-begin.x)*(end.y-begin.y)];
        }
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T const& operator() (int i, int j, int k, int n) const {
            AMREX_ASSERT(i>=begin.x && i<end.x && j>=begin.y && j<end.y && k>=begin.z && k<end.z);
            return p[(i-begin.x+(j-begin.y)*(end.x-begin.x)+(k-begin.z)*(end.x-begin.x)*(end.y-begin.y))+static_cast<long>(n)*(((end.x-begin.x)*(end.y-begin.y))*(end.z-begin.z))];
        }
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k, int n) {
            AMREX_ASSERT(i>=begin.x && i<end.x && j>=begin.y && j<end.y && k>=begin.z && k<end.z);
            return p[(i-begin.x+(j-begin.y)*(end.x-begin.x)+(k-begin.z)*(end.x-begin.x)*(end.y-begin.y))+static_cast<long>(n)*(((end.x-begin.x)*(end.y-begin.y))*(end.z-begin.z))];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T const* ptr (int i, int j, int k) const {
            return p + (i-begin.x+(j-begin.y)*(end.x-begin.x)+(k-begin.z)*(end.x-begin.x)*(end.y-begin.y));
        }
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T const* ptr (int i, int j, int k, int n) const {
            return p + ((i-begin.x+(j-begin.y)*(end.x-begin.x)+(k-begin.z)*(end.x-begin.x)*(end.y-begin.y))+static_cast<long>(n)*(((end.x-begin.x)*(end.y-begin.y))*(end.z-begin.z)));
        }
    };

}

// If on GPUs, use a simple array wrapper that 
// designates needed functions as __host__ __device__
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

namespace amrex {
    template <class T, std::size_t N>
    struct GpuArray
    {
        AMREX_GPU_HOST_DEVICE
        inline const T& operator [] (int i) const { return arr[i]; }
        AMREX_GPU_HOST_DEVICE
        inline T& operator [] (int i) { return arr[i]; }
        AMREX_GPU_HOST_DEVICE
        inline const T* data() const { return arr; };
        AMREX_GPU_HOST_DEVICE
        inline std::size_t size() const { return N; };

        T arr[N];
    };
}

#else

namespace amrex {
    template <class T, std::size_t N>
    using GpuArray = std::array<T,N>;
}

#endif

namespace amrex
{
    template <class T, typename = typename T::FABType>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (std::array<T,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }

    template <class T>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (const std::array<std::unique_ptr<T>,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(a[0].get(), a[1].get(), a[2].get())};
    }
    
    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T*,AMREX_SPACEDIM>& a)
    {
       return {AMREX_D_DECL(a[0], a[1], a[2])};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<std::unique_ptr<T>,AMREX_SPACEDIM>& a)
    {
        return {AMREX_D_DECL(a[0].get(), a[1].get(), a[2].get())};
    }
}

#endif

