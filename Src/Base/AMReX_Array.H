
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
#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>

namespace amrex {

    template <class T, std::size_t N>
    using Array = std::array<T,N>;

    using RealArray = Array<Real, AMREX_SPACEDIM>;
    using IntArray  = Array<int , AMREX_SPACEDIM>;

    struct Dim3 { int x; int y; int z; };
    struct XDim3 { Real x; Real y; Real z; };

    template <typename T>
    struct Array4
    {
        T* AMREX_RESTRICT p;
        long jstride;
        long kstride;
        long nstride;
        Dim3 begin;
        Dim3 end;

        template <class U=T, class = typename std::enable_if<std::is_const<U>::value>::type >
        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (Array4<typename std::remove_const<T>::type> const& rhs) noexcept
            : p(rhs.p),
              jstride(rhs.jstride),
              kstride(rhs.kstride),
              nstride(rhs.nstride),
              begin(rhs.begin),
              end(rhs.end)
            {}

        AMREX_GPU_HOST_DEVICE
        constexpr Array4 (T* a_p, Dim3 const& a_begin, Dim3 const& a_end) noexcept
            : p(a_p),
              jstride(a_end.x-a_begin.x),
              kstride(jstride*(a_end.y-a_begin.y)),
              nstride(kstride*(a_end.z-a_begin.z)),
              begin(a_begin),
              end(a_end)
            {}

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k);
#endif
            return p[(i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k);
#endif
            return p[(i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride+n*nstride];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k);
#endif
            return p + ((i-begin.x)+(j-begin.y)*jstride+(k-begin.z)*kstride);
        }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        AMREX_GPU_HOST_DEVICE inline
        void index_assert (int i, int j, int k) const
        {
            if (i<begin.x || i>=end.x || j<begin.y || j>=end.y || k<begin.z || k>=end.z) {
#ifdef __CUDA_ARCH__
                std::printf(" (%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d)\n",
                            i, j, k, begin.x, end.x-1, begin.y, end.y-1, begin.z, end.z-1);
                amrex::Abort();
#else
                std::stringstream ss;
                ss << " (" << i << "," << j << "," << k << ") is out of bound ("
                   << begin.x << ":" << end.x-1 << ","
                   << begin.y << ":" << end.y-1 << ","
                   << begin.z << ":" << end.z-1 << ")";
                amrex::Abort(ss.str());
#endif
            }
        }
#endif
    };

    inline std::ostream& operator<< (std::ostream& os, const Dim3& d) {
        os << '(' << d.x << ',' << d.y << ',' << d.z << ')';
        return os;
    }

    template <typename T>
    std::ostream& operator<< (std::ostream& os, const Array4<T>& a) {
        os << '(' << a.begin << ',' << a.end << ')';
        return os;
    }
}

// If on GPUs, use a simple array wrapper that 
// designates needed functions as __host__ __device__
#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

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

    inline
    XDim3 makeXDim3 (const Array<Real,AMREX_SPACEDIM>& a) noexcept
    {
#if (AMREX_SPACEDIM == 1)
        return XDim3{a[0], 0., 0.};
#elif (AMREX_SPACEDIM == 2)
        return XDim3{a[0], a[1], 0.};
#else
        return XDim3{a[0], a[1], a[2]};
#endif
    }
}

#endif

