#ifndef AMREX_ARRAY_H_
#define AMREX_ARRAY_H_

#include <array>
#include <memory>
#include <utility>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <AMReX.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>

namespace amrex {

    template <class T, std::size_t N>
    using Array = std::array<T,N>;

    using RealArray = Array<Real, AMREX_SPACEDIM>;
    using IntArray  = Array<int , AMREX_SPACEDIM>;

}

// If on GPUs, use a simple array wrapper that 
// designates needed functions as __host__ __device__
#ifdef AMREX_USE_GPU

namespace amrex {
    template <class T, std::size_t N>
    struct GpuArray
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator [] (int i) const noexcept { return arr[i]; }
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator [] (int i) noexcept { return arr[i]; }
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* data() const noexcept { return arr; };
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        std::size_t size() const noexcept { return N; };
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin() const noexcept { return arr; };
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end() const noexcept { return arr + N; };
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin() noexcept { return arr; };
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end() noexcept { return arr + N; };

        T arr[N];
    };
}

#else

namespace amrex {
    template <class T, std::size_t N>
    using GpuArray = std::array<T,N>;
}

#endif

namespace amrex {
    template <class T, int XLO, int XHI>
    struct Array1D
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i) const noexcept {
            return arr[i-XLO];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i) noexcept {
            return arr[i-XLO];
        }

        T arr[(XHI-XLO+1)];
    };

    template <class T, int XLO, int XHI, int YLO, int YHI>
    struct Array2D
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j) const noexcept {
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j) noexcept {
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        T arr[(XHI-XLO+1)*(YHI-YLO+1)];
    };

    template <class T, int XLO, int XHI, int YLO, int YHI, int ZLO, int ZHI>
    struct Array3D
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j, int k) const noexcept {
            return arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                       -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) noexcept {
            return arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                       -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
        }

        T arr[(XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1)];
    };
}

namespace amrex
{
    template <class T, typename = typename T::FABType>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (std::array<T,AMREX_SPACEDIM>& a) noexcept
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }

    template <class T>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (const std::array<std::unique_ptr<T>,AMREX_SPACEDIM>& a) noexcept
    {
        return {AMREX_D_DECL(a[0].get(), a[1].get(), a[2].get())};
    }
    
    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T,AMREX_SPACEDIM>& a) noexcept
    {
        return {AMREX_D_DECL(&a[0], &a[1], &a[2])};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T*,AMREX_SPACEDIM>& a) noexcept
    {
       return {AMREX_D_DECL(a[0], a[1], a[2])};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<std::unique_ptr<T>,AMREX_SPACEDIM>& a) noexcept
    {
        return {AMREX_D_DECL(a[0].get(), a[1].get(), a[2].get())};
    }

}

#endif

