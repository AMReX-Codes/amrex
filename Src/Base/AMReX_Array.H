#ifndef AMREX_ARRAY_H_
#define AMREX_ARRAY_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_BLassert.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Dim3.H>

#include <array>
#include <memory>
#include <utility>
#include <string>
#include <type_traits>

namespace amrex {

    template <class T, std::size_t N>
    using Array = std::array<T,N>;

    using RealArray = Array<Real, AMREX_SPACEDIM>;
    using IntArray  = Array<int , AMREX_SPACEDIM>;

}

namespace amrex {
    template <class T, unsigned int N>
    struct GpuArray
    {
        using value_type = T;
        using reference_type = T&;

        // Array access

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator [] (int i) const noexcept { return arr[i]; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator [] (int i) noexcept { return arr[i]; }

        // Array pointer access

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* data () const noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* data () noexcept { return arr; }

        // Size
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr unsigned int size() const noexcept { return N; }

        // Pointer access

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin() const noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end() const noexcept { return arr + N; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin() noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end() noexcept { return arr + N; }

        // Fill
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void fill( const T& value ) noexcept
        { for (unsigned int i = 0; i < N; ++i) arr[i] = value; }

        // sum() and product() for floating point types

        template <typename F=T,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0.0;
            for (unsigned int i = 0; i < N; ++i) s += arr[i];
            return s;
        }

        template <typename F=T,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for (unsigned int i = 0; i < N; ++i) p *= arr[i];
            return p;
        }

        // sum() and product() for integral types

        template <typename I=T,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0;
            for (unsigned int i = 0; i < N; ++i) s += arr[i];
            return s;
        }

        template <typename I=T,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1;
            for (unsigned int i = 0; i < N; ++i) p *= arr[i];
            return p;
        }

        // Member array
        T arr[amrex::max(N,1u)];
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
        // Size
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr unsigned int size() const noexcept { return (XHI-XLO+1); }

        // Length

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 length() const noexcept { return Dim3{(XHI-XLO+1),0,0}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr unsigned int length(const int axis) const noexcept {
            amrex::ignore_unused(axis);
            AMREX_ASSERT(axis == 0);
            return (XHI-XLO+1);
        }

        // Bounds

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 lbound() const noexcept { return Dim3{XLO,0,0}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 ubound() const noexcept { return Dim3{XHI,0,0}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr int lbound(const int axis) const noexcept {
            amrex::ignore_unused(axis);
            AMREX_ASSERT(axis == 0);
            return XLO;
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr int ubound(const int axis) const noexcept {
            amrex::ignore_unused(axis);
            AMREX_ASSERT(axis == 0);
            return XHI;
        }

        // Pointer access

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin() const noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end() const noexcept { return arr + XHI-XLO+1; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin() noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end() noexcept { return arr + XHI-XLO+1; }

        // Array access

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

        // sum() and product() for floating point types

        template <typename F=T,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0.0;
            for (int i = XLO; i <= XHI; ++i) s += arr[i-XLO];
            return s;
        }

        template <typename F=T,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for (int i = XLO; i <= XHI; ++i) p *= arr[i-XLO];
            return p;
        }

        // sum() and product() for integral types

        template <typename I=T,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0;
            for (int i = XLO; i <= XHI; ++i) s += arr[i-XLO];
            return s;
        }

        template <typename I=T,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for (int i = XLO; i <= XHI; ++i) p *= arr[i-XLO];
            return p;
        }

        // Member array
        T arr[(XHI-XLO+1)];
    };

    template <class T, int XLO, int XHI, int YLO, int YHI,
              class ORDER=Order::F> // Fortran order by default
    struct Array2D
    {
        // Size
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr unsigned int size() const noexcept { return (XHI-XLO+1)*(YHI-YLO+1); }

        // Length

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 length() const noexcept { return Dim3{(XHI-XLO+1),(YHI-YLO+1),0}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr unsigned int length(const int axis) const noexcept {
            AMREX_ASSERT(axis == 0 || axis == 1);
            if      (axis == 0)    { return (XHI-XLO+1); }
            else /* (axis == 1) */ { return (YHI-YLO+1); }
        }

        // Bounds

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 lbound() const noexcept { return Dim3{XLO,YLO,0}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 ubound() const noexcept { return Dim3{XHI,YHI,0}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr int lbound(const int axis) const noexcept {
            AMREX_ASSERT(axis == 0 || axis == 1);
            if      (axis == 0)    { return XLO; }
            else /* (axis == 1) */ { return YLO; }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr int ubound(const int axis) const noexcept {
            AMREX_ASSERT(axis == 0 || axis == 1);
            if      (axis == 0)    { return XHI; }
            else /* (axis == 1) */ { return YHI; }
        }

        // Pointer access

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin() const noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end() const noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1); }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin() noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end() noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1); }

        // Array access, Fortran order

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

        // sum() and product() for floating point types, Fortran order

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0.0;
            for     (int j = YLO; j <= YHI; ++j)
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc)
        {
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for     (int j = YLO; j <= XHI; ++j)
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc)
        {
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            }
            return p;
        }

        // sum() and product() for integral types, Fortran order

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0;
            for     (int j = YLO; j <= YHI; ++j)
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc)
        {
            T s = (T)0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1;
            for     (int j = YLO; j <= YHI; ++j)
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc)
        {
            T p = (T)1;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
            }
            return p;
        }

        // Array access, C order

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

        // sum() and product() for floating point types, C order

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0.0;
            for     (int i = XLO; i <= XHI; ++i)
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc)
        {
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for     (int i = XLO; i <= XHI; ++i)
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc)
        {
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            }
            return p;
        }

        // sum() and product() for integral types, C order

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0;
            for     (int i = XLO; i <= XHI; ++i)
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc)
        {
            T s = (T)0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1;
                for     (int i = XLO; i <= XHI; ++i)
                    for (int j = YLO; j <= YHI; ++j)
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc)
        {
            T p = (T)1;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
            }
            return p;
        }

        // Member array
        T arr[(XHI-XLO+1)*(YHI-YLO+1)];
    };

    template <class T, int XLO, int XHI, int YLO, int YHI, int ZLO, int ZHI,
              class ORDER=Order::F> // Fortran order by default
    struct Array3D
    {
        // Size
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr unsigned int size() const noexcept { return (XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1); }

        // Length

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 length() const noexcept { return Dim3{(XHI-XLO+1),(YHI-YLO+1),(ZHI-ZLO+1)}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr unsigned int length(const int axis) const noexcept {
            AMREX_ASSERT(axis >= 0 || axis <= 2);
            if      (axis == 0)    { return (XHI-XLO+1); }
            else if (axis == 1)    { return (YHI-YLO+1); }
            else /* (axis == 2) */ { return (ZHI-ZLO+1); }
        }

        // Bounds

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 lbound() const noexcept { return Dim3{XLO,YLO,ZLO}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Dim3 ubound() const noexcept { return Dim3{XHI,YHI,ZHI}; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr int lbound(const int axis) const noexcept {
            AMREX_ASSERT(axis >= 0 || axis <= 2);
            if      (axis == 0)    { return XLO; }
            else if (axis == 1)    { return YLO; }
            else /* (axis == 2) */ { return ZLO; }
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr int ubound(const int axis) const noexcept {
            AMREX_ASSERT(axis >= 0 || axis <= 2);
            if      (axis == 0)    { return XHI; }
            else if (axis == 1)    { return YHI; }
            else /* (axis == 2) */ { return ZHI; }
        }

        // Pointer access

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin() const noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end() const noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1); }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin() noexcept { return arr; }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end() noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1); }

        // Array access, Fortran order

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

        // sum() and product() for floating point types, Fortran order

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0.0;
            for         (int k = ZLO; k <= ZHI; ++k)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int i = XLO; i <= XHI; ++i)
                        s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                                 -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc0, const int loc1)
        {
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for         (int k = ZLO; k <= ZHI; ++k)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int i = XLO; i <= XHI; ++i)
                        p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                                 -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc0, const int loc1)
        {
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            }
            return p;
        }

        // sum() and product() for integral types, Fortran order

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0;
            for         (int k = ZLO; k <= ZHI; ++k)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int i = XLO; i <= XHI; ++i)
                        s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                                 -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc0, const int loc1)
        {
            T s = (T)0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1;
            for         (int k = ZLO; k <= ZHI; ++k)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int i = XLO; i <= XHI; ++i)
                        p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                                 -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc0, const int loc1)
        {
            T p = (T)1;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
            }
            return p;
        }

        // Array access, C order

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

        // sum() and product() for floating point types, C order

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0.0;
            for         (int i = XLO; i <= XHI; ++i)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int k = ZLO; k <= ZHI; ++k)
                        s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                                 -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc0, const int loc1)
        {
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for         (int i = XLO; i <= XHI; ++i)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int k = ZLO; k <= ZHI; ++k)
                        p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                                 -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename F=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc0, const int loc1)
        {
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            }
            return p;
        }

        // sum() and product() for integral types, C order

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum() noexcept
        {
            T s = (T)0;
            for         (int i = XLO; i <= XHI; ++i)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int k = ZLO; k <= ZHI; ++k)
                        s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                                 -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum(const int axis, const int loc0, const int loc1)
        {
            T s = (T)0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            }
            return s;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1;
            for         (int i = XLO; i <= XHI; ++i)
                for     (int j = YLO; j <= YHI; ++j)
                    for (int k = ZLO; k <= ZHI; ++k)
                        p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                                 -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            return p;
        }

        template <typename O=ORDER,
                  typename I=T,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0,
                  typename std::enable_if<std::is_integral<I>::value,int>::type FOO = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product(const int axis, const int loc0, const int loc1)
        {
            T p = (T)1;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                for (int i = XLO; i <= XHI; ++i)
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                for (int j = YLO; j <= YHI; ++j)
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                for (int k = ZLO; k <= ZHI; ++k)
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
            }
            return p;
        }

        // Member array
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
