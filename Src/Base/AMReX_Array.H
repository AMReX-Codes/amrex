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
        AMREX_GPU_HOST_DEVICE
        inline const T& operator [] (int i) const noexcept { return arr[i]; }
        AMREX_GPU_HOST_DEVICE
        inline T& operator [] (int i) noexcept { return arr[i]; }
        AMREX_GPU_HOST_DEVICE
        inline const T* data() const noexcept { return arr; };
        AMREX_GPU_HOST_DEVICE
        inline std::size_t size() const noexcept { return N; };

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

