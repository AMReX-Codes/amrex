#ifndef AMREX_ARRAY_H_
#define AMREX_ARRAY_H_

#include <array>
#include <memory>
#include <utility>
#include <string>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <type_traits>
#include <AMReX.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Dim3.H>

namespace amrex {

    template <class T, std::size_t N>
    using Array = std::array<T,N>;

    using RealArray = Array<Real, AMREX_SPACEDIM>;
    using IntArray  = Array<int , AMREX_SPACEDIM>;

}

namespace amrex {
    template <class T, std::size_t N>
    struct GpuArray
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator [] (int i) const noexcept { return arr[i]; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator [] (int i) noexcept { return arr[i]; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* data () const noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* data () noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        std::size_t size() const noexcept { return N; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin() const noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end() const noexcept { return arr + N; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin() noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end() noexcept { return arr + N; }

        T arr[amrex::max(N,std::size_t{1})];
    };
}

namespace amrex {

    namespace Order {
        struct C {};
        struct F {};
    }

    template <class T, int XLO, int XHI>
    struct Array1D
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i) const noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI);
            return arr[i-XLO];
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i) noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI);
            return arr[i-XLO];
        }

        T arr[(XHI-XLO+1)];
    };

    template <class T, int XLO, int XHI, int YLO, int YHI,
              class ORDER=Order::F> // Fortran order by default
    struct Array2D
    {
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j) const noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j) noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j) const noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
        }

        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j) noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
        }

        T arr[(XHI-XLO+1)*(YHI-YLO+1)];
    };

    template <class T, int XLO, int XHI, int YLO, int YHI, int ZLO, int ZHI,
              class ORDER=Order::F> // Fortran order by default
    struct Array3D
    {
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j, int k) const noexcept {
            return arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                       -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
        }

        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) noexcept {
            return arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                       -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
        }

        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j, int k) const noexcept {
            return arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                       -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
        }

        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) noexcept {
            return arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                       -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
        }

        T arr[(XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1)];
    };
}

namespace amrex
{
    template <class T, typename = typename T::FABType>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (std::array<T,AMREX_SPACEDIM>& a) noexcept
    {
        return {{AMREX_D_DECL(&a[0], &a[1], &a[2])}};
    }

    template <class T>
    std::array<T*,AMREX_SPACEDIM> GetArrOfPtrs (const std::array<std::unique_ptr<T>,AMREX_SPACEDIM>& a) noexcept
    {
        return {{AMREX_D_DECL(a[0].get(), a[1].get(), a[2].get())}};
    }
    
    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T,AMREX_SPACEDIM>& a) noexcept
    {
        return {{AMREX_D_DECL(&a[0], &a[1], &a[2])}};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<T*,AMREX_SPACEDIM>& a) noexcept
    {
        return {{AMREX_D_DECL(a[0], a[1], a[2])}};
    }

    template <class T>
    std::array<T const*,AMREX_SPACEDIM> GetArrOfConstPtrs (const std::array<std::unique_ptr<T>,AMREX_SPACEDIM>& a) noexcept
    {
        return {{AMREX_D_DECL(a[0].get(), a[1].get(), a[2].get())}};
    }

}

namespace amrex
{
    inline XDim3 makeXDim3 (const Array<Real,AMREX_SPACEDIM>& a) noexcept
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

