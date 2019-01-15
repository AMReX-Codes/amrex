
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
        T* p;
        int jstride;
        int kstride;
        long nstride;
        Dim3 lo;
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T const& operator() (int i, int j, int k, int n=0) const {
            return *(p + ((i-lo.x)+(j-lo.y)*jstride+(k-lo.z)*kstride+n*nstride));
        }
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k, int n=0) {
            return *(p + ((i-lo.x)+(j-lo.y)*jstride+(k-lo.z)*kstride+n*nstride));
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

