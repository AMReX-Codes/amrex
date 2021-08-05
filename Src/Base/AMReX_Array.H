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

        /**
         * GpuArray elements are indexed using square brackets, as with any
         * other array.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator [] (int i) const noexcept { return arr[i]; }

        /**
         * GpuArray elements are indexed using square brackets, as with any
         * other array.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator [] (int i) noexcept { return arr[i]; }

        /**
         * Returns a const pointer to the underlying data of a GpuArray object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* data () const noexcept { return arr; }

        /**
         * Returns a pointer to the underlying data of a GpuArray object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* data () noexcept { return arr; }

        /**
         * Returns the number of elements in the GpuArray object as an
         * unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int size () noexcept { return N; }

        /**
         * Returns a const pointer address to the first element of the
         * GpuArray object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin () const noexcept { return arr; }

        /**
         * Returns a const pointer address right after the last element of the
         * GpuArray object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end () const noexcept { return arr + N; }

        /**
         * Returns a pointer address to the first element of the
         * GpuArray object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin () noexcept { return arr; }

        /**
         * Returns a pointer address right after the last element of the
         * GpuArray object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end () noexcept { return arr + N; }

        /**
         * Fills in all of the elements in the GpuArray object to the same
         * value.
         *
         * \param value The fill value
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void fill ( const T& value ) noexcept
        { for (unsigned int i = 0; i < N; ++i) arr[i] = value; }

        /**
         * Returns the sum of all elements in the GpuArray object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename F=T,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum () noexcept
        {
            T s = (T)0.0;
            for (unsigned int i = 0; i < N; ++i) { s += arr[i]; }
            return s;
        }

        /**
         * Returns the product of all elements in the GpuArray object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename F=T,
                  typename std::enable_if<std::is_floating_point<F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product () noexcept
        {
            T p = (T)1.0;
            for (unsigned int i = 0; i < N; ++i) { p *= arr[i]; }
            return p;
        }

        /**
         * Returns the sum of all elements in the GpuArray object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename I=T,
                  typename std::enable_if<std::is_integral<I>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum () noexcept
        {
            T s = (T)0;
            for (unsigned int i = 0; i < N; ++i) s += arr[i];
            return s;
        }

        /**
         * Returns the product of all elements in the GpuArray object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename I=T,
                  typename std::enable_if<std::is_integral<I>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product () noexcept
        {
            T p = (T)1;
            for (unsigned int i = 0; i < N; ++i) p *= arr[i];
            return p;
        }

        /**
         * GpuArray is implemented as a fixed-size array. Hence, no constructor or
         * destructor is given.
         */
        T arr[amrex::max(N,1u)];
    };
}

namespace amrex {

    /**
     * Array2D and Array3D objects can be indexed according to Fortran
     * column-major order (first index moving fastest) or C/C++ row-major
     * order (last index moving fastest). If not specified, Fortran order is
     * assumed.
     */
    namespace Order {
        struct C {};
        struct F {};
    }

    /**
     * A GPU-compatible one-dimensional array.
     *
     * \tparam XLO Index for lower bound. Can be other than 0.
     * \tparam XHI Index for upper bound.
     */
    template <class T, int XLO, int XHI>
    struct Array1D
    {
        /**
         * Returns the number of elements in the Array1D object as an unsigned
         * integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int size () noexcept { return (XHI-XLO+1); }

        /**
         * Returns the index of the lower bound of the Array1D object.
         * Can be other than 0.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int lo () noexcept { return XLO; }

        /**
         * Returns the index of the upper bound of the Array1D object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int hi () noexcept { return XHI; }

        /**
         * Returns the number of elements in the Array1D object as an unsigned
         * integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int len () noexcept { return (XHI-XLO+1); }

        /**
         * Returns a const pointer address to the first element of the
         * Array1D object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin () const noexcept { return arr; }

        /**
         * Returns a const pointer address right after the last element of the
         * Array1D object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end () const noexcept { return arr + XHI-XLO+1; }

        /**
         * Returns a pointer address to the first element of the
         * Array1D object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin () noexcept { return arr; }

        /**
         * Returns a pointer address right after the last element of the
         * Array1D object.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end () noexcept { return arr + XHI-XLO+1; }

        /**
         * The elements of an Array1D object are accessed using parentheses,
         * e.g. array(i), instead of using square brackets.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i) const noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI);
            return arr[i-XLO];
        }

        /**
         * The elements of an Array1D object are accessed using parentheses,
         * e.g. array(i), instead of using square brackets.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i) noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI);
            return arr[i-XLO];
        }

        /**
         * Returns the sum of all elements in the Array1D object.
         * Only makes sense for integers and floating-point types.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum () noexcept
        {
            T s = (T)0.0;
            for (int i = XLO; i <= XHI; ++i) { s += arr[i-XLO]; }
            return s;
        }

        /**
         * Returns the product of all elements in the Array1D object.
         * Only makes sense for integers and floating-point types.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product() noexcept
        {
            T p = (T)1.0;
            for (int i = XLO; i <= XHI; ++i) { p *= arr[i-XLO]; }
            return p;
        }

        /**
         * Array1D is implemented as a fixed-size array.
         * Hence, no constructor or destructor is given.
         */
        T arr[(XHI-XLO+1)];
    };

    /**
     * A GPU-compatible two-dimensional array.
     *
     * \tparam XLO Index for lower bound in x dimension. Can be other than 0.
     * \tparam XHI Index for upper bound in x dimension.
     * \tparam YLO Index for lower bound in y dimension. Can be other than 0.
     * \tparam YHI Index for upper bound in y dimension.
     * \tparam ORDER Either Order::C (C/C++ row-major order) or
     *               Order::F (Fortran column-major order, which is the
     *               default if not given)
     */
    template <class T, int XLO, int XHI, int YLO, int YHI,
              class ORDER=Order::F>
    struct Array2D
    {
        /**
         * Returns the total number of elements of the Array2D object as an
         * unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int size() noexcept { return (XHI-XLO+1)*(YHI-YLO+1); }

        /**
         * Returns the index of the lower bound of the Array2D object in the
         * x direction.
         * Can be other than 0.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int xlo () noexcept { return XLO; }

        /**
         * Returns the index of the upper bound of the Array2D object in the
         * x direction.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int xhi () noexcept { return XHI; }


        /**
         * Returns the number of elements of the Array2D object in the
         * x direction as an unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int xlen () noexcept { return (XHI-XLO+1); }

        /**
         * Returns the index of the lower bound of the Array2D object in the
         * y direction.
         * Can be other than 0.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int ylo () noexcept { return YLO; }

        /**
         * Returns the index of the upper bound of the Array2D object in the
         * y direction.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int yhi () noexcept { return YHI; }


        /**
         * Returns the number of elements of the Array2D object in the
         * y direction as an unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int ylen () noexcept { return (YHI-YLO+1); }

        /**
         * Returns a const pointer address to the first element of the
         * Array2D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin () const noexcept { return arr; }

        /**
         * Returns a const pointer address right after the last element of the
         * Array2D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end () const noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1); }

        /**
         * Returns a pointer address to the first element of the
         * Array2D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin () noexcept { return arr; }

        /**
         * Returns a pointer address right after the last element of the
         * Array2D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end () noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1); }

        /**
         * The elements of an Array2D object are accessed using parentheses,
         * e.g. array(i,j), instead of using square brackets.
         * If the order is not specified, Fortran column-major order is assumed
         * (the index i moves the fastest)
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j) const noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        /**
         * The elements of an Array2D object are accessed using parentheses,
         * e.g. array(i,j), instead of using square brackets.
         * If the order is not specified, Fortran column-major order is assumed
         * (the index i moves the fastest)
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j) noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
        }

        /**
         * When called without any arguments, returns the sum of all
         * elements in the Array2D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum () noexcept
        {
            T s = (T)0.0;
            for     (int j = YLO; j <= YHI; ++j) {
                for (int i = XLO; i <= XHI; ++i) {
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
                }
            }
            return s;
        }

        /**
         * When called with two arguments, performs a sum reduction over
         * the specified axis, for a particular location index loc.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension)
         * \param loc  The appropriate location index
         *
         * This can be used, for instance, to calculate the sum over the y
         * dimension of an Array2D object that was instantiated as
         * \verbatim Array2D<amrex::Real, 1, M, 1, N> array; \endverbatim
         *
         * One could instantiate an Array1D object to hold the results,
         * \verbatim Array1D<amrex::Real, 1, M> vec; \endverbatim
         * and then perform the summation for each element of the resulting
         * vector.
         * \verbatim
         for (int i = 1; i <= M; ++i) {
             vec(i) = array.sum(1,i)
         }
         \endverbatim
         * In this example, the axis is 1 and the location index is \c i.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum (const int axis, const int loc)
        {
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i) {
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
                }
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j) {
                    s += arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
                }
            }
            return s;
        }

        /**
         * When called without any arguments, returns the product of all
         * elements in the Array2D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product () noexcept
        {
            T p = (T)1.0;
            for     (int j = YLO; j <= XHI; ++j) {
                for (int i = XLO; i <= XHI; ++i) {
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
                }
            }
            return p;
        }


        /**
         * When called with two arguments, performs a product reduction over
         * the specified axis, for a particular location index loc.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension)
         * \param loc  The appropriate location index
         *
         * This can be used, for instance, to calculate the product over the x
         * dimension of an Array2D object that was instantiated as
         * \verbatim Array2D<amrex::Real, 1, M, 1, N> array; \endverbatim
         *
         * One could instantiate an Array1D object to hold the results with,
         * \verbatim Array1D<amrex::Real, 1, N> vec; \endverbatim
         * and then perform the product for each element of the resulting
         * vector.
         * 
         * \verbatim
         for (int j = 1; j <= N; ++j) {
             vec(j) = array.sum(0,j)
         }
         \endverbatim
         * In this example, the axis is 0 and the location index is \c j.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product (const int axis, const int loc)
        {
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc;
                for (int i = XLO; i <= XHI; ++i) {
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
                }
            } else if (axis == 1) {
                int i = loc;
                for (int j = YLO; j <= YHI; ++j) {
                    p *= arr[i+j*(XHI-XLO+1)-(YLO*(XHI-XLO+1)+XLO)];
                }
            }
            return p;
        }

        /**
         * The elements of an Array2D object are accessed using parentheses,
         * e.g. array(i,j), instead of using square brackets.
         * When the order is manually specified as Order::C, row-major order
         * is used (the index j moves the fastest).
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j) const noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
        }

        /**
         * The elements of an Array2D object are accessed using parentheses,
         * e.g. array(i,j), instead of using square brackets.
         * When the order is manually specified as Order::C, row-major order
         * is used (the index j moves the fastest).
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j) noexcept {
            AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
            return arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
        }

        /**
         * When called without any arguments, returns the sum of all
         * elements in the Array2D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum () noexcept
        {
            T s = (T)0.0;
            for     (int i = XLO; i <= XHI; ++i) {
                for (int j = YLO; j <= YHI; ++j) {
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
                }
            }
            return s;
        }

        /**
         * When called with two arguments, performs a sum reduction over
         * the specified axis, for a particular location index loc.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension)
         * \param loc  The appropriate location index
         *
         * This can be used, for instance, to calculate the sum over the y
         * dimension of an Array2D object that was instantiated as
         * \verbatim Array2D<amrex::Real, 1, M, 1, N, Order::C> array; \endverbatim
         *
         * One could instantiate an Array1D object to hold the results,
         * \verbatim Array1D<amrex::Real, 1, M> vec; \endverbatim
         * and then perform the summation for each element of the resulting
         * vector.
         * \verbatim
         for (int i = 1; i <= M; ++i) {
             vec(i) = array.sum(1,i)
         }
         \endverbatim
         * In this example, the axis is 1 and the location index is i.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum (const int axis, const int loc)
        {
            AMREX_ASSERT(axis == 0 || axis == 1);
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc;
                AMREX_ASSERT(j >= YLO && j <= YHI);
                for (int i = XLO; i <= XHI; ++i) {
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
                }
            } else if (axis == 1) {
                int i = loc;
                AMREX_ASSERT(i >= XLO && i <= XHI);
                for (int j = YLO; j <= YHI; ++j) {
                    s += arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
                }
            }
            return s;
        }

        /**
         * When called without any arguments, returns the product of all
         * elements in the Array2D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product () noexcept
        {
            T p = (T)1.0;
            for     (int i = XLO; i <= XHI; ++i) {
                for (int j = YLO; j <= YHI; ++j) {
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
                }
            }
            return p;
        }

        /**
         * When called with two arguments, performs a product reduction over
         * the specified axis, for a particular location index loc.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension)
         * \param loc  The appropriate location index
         *
         * This can be used, for instance, to calculate the product over the x
         * dimension of an Array2D object that was instantiated as
         * \verbatim Array2D<amrex::Real, 1, M, 1, N, Order::C> array; \endverbatim
         *
         * One could instantiate an Array1D object to hold the results,
         * \verbatim Array1D<amrex::Real, 1, N> vec; \endverbatim
         * and then perform the summation for each element of the resulting
         * vector.
         * \verbatim
         for (int j = 1; j <= N; ++j) {
             vec(j) = array.product(0,j)
         }
         \endverbatim
         * In this example, the axis is 0 and the location index is \c j.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product (const int axis, const int loc)
        {
            AMREX_ASSERT(axis == 0 || axis == 1);
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc;
                AMREX_ASSERT(j >= YLO && j <= YHI);
                for (int i = XLO; i <= XHI; ++i) {
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
                }
            } else if (axis == 1) {
                int i = loc;
                AMREX_ASSERT(i >= XLO && i <= XHI);
                for (int j = YLO; j <= YHI; ++j) {
                    p *= arr[j+i*(YHI-YLO+1)-(XLO*(YHI-YLO+1)+YLO)];
                }
            }
            return p;
        }

        /**
         * Array2D is implemented as a fixed-size array.
         * Hence, no constructor or destructor is given.
         */
        T arr[(XHI-XLO+1)*(YHI-YLO+1)];
    };


    /**
     * A GPU-compatible three-dimensional array.
     *
     * \tparam XLO Index for lower bound in x dimension. Can be other than 0.
     * \tparam XHI Index for upper bound in x dimension.
     * \tparam YLO Index for lower bound in y dimension. Can be other than 0.
     * \tparam YHI Index for upper bound in y dimension.
     * \tparam ZLO Index for lower bound in z dimension. Can be other than 0.
     * \tparam ZHI Index for upper bound in z dimension.
     * \tparam ORDER Either Order::C (C/C++ row-major order) or
     *               Order::F (Fortran column-major order, which is the
     *               default if not given)
     */
    template <class T, int XLO, int XHI, int YLO, int YHI, int ZLO, int ZHI,
              class ORDER=Order::F>
    struct Array3D
    {
        /**
         * Returns the total number of elements in the Array3D object as an
         * unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int size () noexcept { return (XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1); }

        /**
         * Returns the index of the lower bound of the Array3D object in the
         * x direction.
         * Can be other than 0.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int xlo () noexcept { return XLO; }

        /**
         * Returns the index of the upper bound of the Array3D object in the
         * x direction.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int xhi () noexcept { return XHI; }

        /**
         * Returns the number of elements of the Array3D object in the
         * x direction as an unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int xlen () noexcept { return (XHI-XLO+1); }

        /**
         * Returns the index of the lower bound of the Array3D object in the
         * y direction.
         * Can be other than 0.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int ylo () noexcept { return YLO; }

        /**
         * Returns the index of the upper bound of the Array3D object in the
         * y direction.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int yhi () noexcept { return YHI; }


        /**
         * Returns the number of elements of the Array3D object in the
         * y direction as an unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int ylen () noexcept { return (YHI-YLO+1); }

        /**
         * Returns the index of the lower bound of the Array3D object in the
         * z direction.
         * Can be other than 0.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int zlo () noexcept { return ZLO; }

        /**
         * Returns the index of the upper bound of the Array3D object in the
         * z direction.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr int zhi () noexcept { return ZHI; }

        /**
         * Number of elements of the Array3D object in the z direction
         * as an unsigned integer.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr unsigned int zlen () noexcept { return (ZHI-ZLO+1); }

        /**
         * Returns a const pointer address to the first element of the
         * Array3D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* begin () const noexcept { return arr; }

        /**
         * Returns a const pointer address right after the last element of the
         * Array3D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T* end () const noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1); }

        /**
         * Returns a pointer address to the first element of the
         * Array3D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* begin () noexcept { return arr; }

        /**
         * Returns a pointer address right after the last element of the
         * Array3D object, as if the object is treated as one-dimensional.
         */
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* end () noexcept { return arr + (XHI-XLO+1)*(YHI-YLO+1)*(ZHI-ZLO+1); }

        /**
         * The elements of an Array3D object are accessed using parentheses,
         * e.g. array(i,j,k), instead of using square brackets.
         * If the order is not specified, Fortran column-major order is assumed
         * (the index i moves the fastest)
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j, int k) const noexcept {
            return arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                       -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
        }

        /**
         * The elements of an Array3D object are accessed using parentheses,
         * e.g. array(i,j,k), instead of using square brackets.
         * If the order is not specified, Fortran column-major order is assumed
         * (the index i moves the fastest)
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) noexcept {
            return arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                       -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
        }

        /**
         * When called without any arguments, returns the sum of all
         * elements in the Array3D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum () noexcept
        {
            T s = (T)0.0;
            for         (int k = ZLO; k <= ZHI; ++k) {
                for     (int j = YLO; j <= YHI; ++j) {
                    for (int i = XLO; i <= XHI; ++i) {
                        s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                                 -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                    }
                }
            }
            return s;
        }


        /**
         * When called with three arguments, performs a sum reduction over
         * the specified axis, for a particular set of location indices loc0
         * and loc1.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension, 2 for z dimension)
         * \param loc0 The appropriate location index (either i or j)
         * \param loc1 The appropriate location index (either j or k)
         *
         * This can be used, for instance, to calculate the sum over the x
         * dimension of an Array3D object that was instantiated as
         * \verbatim Array3D<amrex::Real, 1, M, 1, N, 1, K> array; \endverbatim
         *
         * One could instantiate an Array2D object to hold the results,
         * \verbatim Array2D<amrex::Real, 1, N, 1, K> mat; \endverbatim
         * and then perform the summation for each element of the resulting
         * matrix.
         * \verbatim
         for     (int j = 1; j <= N; ++j) {
             for (int k = 1; k <= K; ++k) {
                 mat(j,k) = array.sum(0,j,k)
             }
         }
         \endverbatim
         * In this example, the axis is 0 and the location indices are loc0 = j
         * and loc1 = k. For axis = 1, the location indices are treated as
         * loc0 = i and loc1 = k; for axis = 2, loc0 = j and loc1 = k.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum (const int axis, const int loc0, const int loc1)
        {
            AMREX_ASSERT(axis >= 0 && axis <= 2);
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                AMREX_ASSERT(j >= YLO && j <= YHI && k >= ZLO && k <= ZHI);
                for (int i = XLO; i <= XHI; ++i) {
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                }
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && k >= ZLO && k <= ZHI);
                for (int j = YLO; j <= YHI; ++j) {
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                }
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
                for (int k = ZLO; k <= ZHI; ++k) {
                    s += arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                }
            }
            return s;
        }

        /**
         * When called without any arguments, returns the product of all
         * elements in the Array3D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product () noexcept
        {
            T p = (T)1.0;
            for         (int k = ZLO; k <= ZHI; ++k) {
                for     (int j = YLO; j <= YHI; ++j) {
                    for (int i = XLO; i <= XHI; ++i) {
                        p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                                 -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                    }
                }
            }
            return p;
        }

        /**
         * When called with three arguments, performs a product reduction over
         * the specified axis, for a particular set of location indices loc0
         * and loc1.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension, 2 for z dimension)
         * \param loc0 The appropriate location index (either i or j)
         * \param loc1 The appropriate location index (either j or k)
         * 
         * This can be used, for instance, to calculate the sum over the z
         * dimension of an Array3D object that was instantiated as
         * \verbatim Array3D<amrex::Real, 1, M, 1, N, 1, K> array; \endverbatim
         *
         * One could instantiate an Array2D object to hold the results,
         * \verbatim Array2D<amrex::Real, 1, M, 1, N> mat; \endverbatim
         * and then perform the summation for each element of the resulting
         * matrix.
         * \verbatim
         for     (int j = 1; j <= N; ++j) {
             for (int i = 1; i <= M; ++i) {
                 mat(i,j) = array.sum(2,i,j)
             }
         }
         \endverbatim
         * In this example, the axis is 2 and the location indices are loc0 = i
         * and loc1 = j. For axis = 0, the location indices are treated as
         * loc0 = j and loc1 = k; for axis = 1, loc0 = i and loc1 = k.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::F>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product (const int axis, const int loc0, const int loc1)
        {
            AMREX_ASSERT(axis >= 0 && axis <= 2);
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                AMREX_ASSERT(j >= YLO && j <= YHI && k >= ZLO && k <= ZHI);
                for (int i = XLO; i <= XHI; ++i) {
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                }
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && k >= ZLO && k <= ZHI);
                for (int j = YLO; j <= YHI; ++j) {
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                }
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
                for (int k = ZLO; k <= ZHI; ++k) {
                    p *= arr[i+j*(XHI-XLO+1)+k*((XHI-XLO+1)*(YHI-YLO+1))
                             -(ZLO*((XHI-XLO+1)*(YHI-YLO+1))+YLO*(XHI-XLO+1)+XLO)];
                }
            }
            return p;
        }

        /**
         * The elements of an Array3D object are accessed using parentheses,
         * e.g. array(i,j,k), instead of using square brackets.
         * When the order is manually specified as Order::C, row-major order
         * is used (the index k moves the fastest).
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        const T& operator() (int i, int j, int k) const noexcept {
            return arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                       -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
        }

        /**
         * The elements of an Array3D object are accessed using parentheses,
         * e.g. array(i,j,k), instead of using square brackets.
         * When the order is manually specified as Order::C, row-major order
         * is used (the index k moves the fastest).
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) noexcept {
            return arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                       -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
        }

        /**
         * When called without any arguments, returns the sum of all
         * elements in the Array3D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum () noexcept
        {
            T s = (T)0.0;
            for         (int i = XLO; i <= XHI; ++i) {
                for     (int j = YLO; j <= YHI; ++j) {
                    for (int k = ZLO; k <= ZHI; ++k) {
                        s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                                 -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                    }
                }
            }
            return s;
        }

        /**
         * When called with three arguments, performs a sum reduction over
         * the specified axis, for a particular set of location indices loc0
         * and loc1.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension, 2 for z dimension)
         * \param loc0 The appropriate location index (either i or j)
         * \param loc1 The appropriate location index (either j or k)
         *
         * This can be used, for instance, to calculate the sum over the y
         * dimension of an Array3D object that was instantiated as
         * \verbatim Array3D<amrex::Real, 1, M, 1, N, 1, K, Order::C> array; \endverbatim
         *
         * One could instantiate an Array2D object to hold the results,
         * \verbatim Array2D<amrex::Real, 1, M, 1, K, Order::C> mat; \endverbatim
         * and then perform the summation for each element of the resulting
         * matrix.
         * \verbatim
          for     (int i = 1; i <= M; ++i) {
              for (int k = 1; k <= K; ++k) {
                  mat(i,k) = array.sum(1,i,k)
              }
          }
          \endverbatim
         * In this example, the axis is 1 and the location indices are loc0 = i
         * and loc1 = k. For axis = 0, the location indices are treated as
         * loc0 = i and loc1 = j; for axis = 2, loc0 = j and loc1 = k.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T sum (const int axis, const int loc0, const int loc1)
        {
            AMREX_ASSERT(axis >= 0 && axis <= 2);
            T s = (T)0.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                AMREX_ASSERT(j >= YLO && j <= YHI && k >= ZLO && k <= ZHI);
                for (int i = XLO; i <= XHI; ++i) {
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                }
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && k >= ZLO && k <= ZHI);
                for (int j = YLO; j <= YHI; ++j) {
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                }
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
                for (int k = ZLO; k <= ZHI; ++k) {
                    s += arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                }
            }
            return s;
        }

        /**
         * When called without any arguments, returns the product of all
         * elements in the Array3D object.
         * Only makes sense for integers and floating-point types.
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product () noexcept
        {
            T p = (T)1.0;
            for         (int i = XLO; i <= XHI; ++i) {
                for     (int j = YLO; j <= YHI; ++j) {
                    for (int k = ZLO; k <= ZHI; ++k) {
                        p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                                 -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                    }
                }
            }
            return p;
        }

        /**
         * When called with three arguments, performs a product reduction over
         * the specified axis, for a particular set of location indices loc0
         * and loc1.
         *
         * \param axis The dimension to reduce (0 for x dimension,
         *             1 for y dimension, 2 for z dimension)
         * \param loc0 The appropriate location index (either i or j)
         * \param loc1 The appropriate location index (either j or k)
         *
         * This can be used, for instance, to calculate the product over the x
         * dimension of an Array3D object that was instantiated as
         * \verbatim Array3D<amrex::Real, 1, M, 1, N, 1, K, Order::C> array; \endverbatim
         *
         * One could instantiate an Array2D object to hold the results,
         * \verbatim Array2D<amrex::Real, 1, N, 1, K, Order::C> mat; \endverbatim
         * and then perform the summation for each element of the resulting
         * matrix.
         * \verbatim
         for     (int j = 1; j <= N; ++j) {
             for (int k = 1; k <= K; ++k) {
                 mat(j,k) = array.sum(0,j,k)
             }
         }
         \endverbatim
         * In this example, the axis is 0 and the location indices are loc0 = j
         * and loc1 = k. For axis = 1, the location indices are treated as
         * loc0 = i and loc1 = k; for axis = 2, loc0 = j and loc1 = k.
         *
         */
        template <typename O=ORDER,
                  typename std::enable_if<std::is_same<O,Order::C>::value,int>::type=0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T product (const int axis, const int loc0, const int loc1)
        {
            AMREX_ASSERT(axis >= 0 && axis <= 2);
            T p = (T)1.0;
            if        (axis == 0) {
                int j = loc0;
                int k = loc1;
                AMREX_ASSERT(j >= YLO && j <= YHI && k >= ZLO && k <= ZHI);
                for (int i = XLO; i <= XHI; ++i) {
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                }
            } else if (axis == 1) {
                int i = loc0;
                int k = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && k >= ZLO && k <= ZHI);
                for (int j = YLO; j <= YHI; ++j) {
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                }
            } else if (axis == 2) {
                int i = loc0;
                int j = loc1;
                AMREX_ASSERT(i >= XLO && i <= XHI && j >= YLO && j <= YHI);
                for (int k = ZLO; k <= ZHI; ++k) {
                    p *= arr[k+j*(ZHI-ZLO+1)+i*((ZHI-ZLO+1)*(YHI-YLO+1))
                             -(XLO*((ZHI-ZLO+1)*(YHI-YLO+1))+YLO*(ZHI-ZLO+1)+ZLO)];
                }
            }
            return p;
        }

        /**
         * Array3D is implemented as a fixed-size array.
         * Hence, no constructor or destructor is given.
         */
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
