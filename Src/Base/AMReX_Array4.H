#ifndef AMREX_ARRAY4_H_
#define AMREX_ARRAY4_H_

#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_ConstexprFor.H>

#include <iostream>
#include <sstream>

/// \cond DOXYGEN_IGNORE
#define AMREX_ARRAY4_INDEX_ASSERT(i,j,k,n)         \
    if ((i)<begin.vect[0] || (i)>=end.vect[0] ||   \
        (j)<begin.vect[1] || (j)>=end.vect[1] ||   \
        (k)<begin.vect[2] || (k)>=end.vect[2] ||   \
        (n)<0             || (n)>=end.vect[3])     \
    {                                              \
        index_assert_print_error_message(i,j,k,n); \
    }
/// \endcond

namespace amrex {

    template<int dim>
    class BoxND;

    /**
     * \ingroup amrex_fab_containers
     * \class CellData
     * \brief Lightweight accessor for data associated with a single cell.
     * \tparam T Data type. May be const-qualified.
     *
     * \c CellData provides a non-owning view of data associated with a single
     * cell in an \c Array4. It's typically created by the \c Array4::cellData
     * member function.
     */
    template <typename T>
    struct CellData
    {
        /// \cond DOXYGEN_IGNORE
        T* AMREX_RESTRICT p = nullptr;
        Long stride = 0;
        int ncomp = 0;
        /// \endcond

        /**
         * \brief Construct a CellData
         *
         * \param a_p data Pointer pointing to component 0 of the cell data
         * \param a_stride Component stride, i.e., the distance (in number of
         *                 elements, not bytes) between successive components
         *                 \c (i,j,k,n) and \c (i,j,k,n+1) of Array4.
         * \param a_ncomp number of components.
         */
        AMREX_GPU_HOST_DEVICE
        constexpr CellData (T* a_p, Long a_stride, int a_ncomp)
            : p(a_p), stride(a_stride), ncomp(a_ncomp)
            {}

        /**
         * \brief Construct a const CellData from a non-const CellData.
         *
         * This enables implicit conversion from \c CellData<T> to
         * \c CellData<const T>.
         *
         * \param rhs Source CellData object.
         */
        template <class U=T,
                  std::enable_if_t<std::is_const_v<U>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr CellData (CellData<std::remove_const_t<T>> const& rhs) noexcept
            : p(rhs.p), stride(rhs.stride), ncomp(rhs.ncomp)
            {}

        /**
         * \brief Check whether the accessor refers to valid data.
         *
         * \return \c true if the internal pointer is non-null, otherwise \c false.
         */
        AMREX_GPU_HOST_DEVICE
        explicit operator bool() const noexcept { return p != nullptr; }

        /**
         * \brief Return the number of components.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE
        int nComp() const noexcept { return ncomp; }

        /**
         * \brief Access the \p n-th component of the cell.
         *
         * \param n Component index.
         * \return Reference to component \p n.
         *
         * \warning No bounds checking is performed, unless \c either AMREX_DEBUG
         *          or AMREX_BOUND_CHECK is defined.  The caller is responsible
         *          for ensuring that \p n is in the range \c [0, ncomp).
         */
        template <class U=T,
                  std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator[] (int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            if (n < 0 || n >= ncomp) {
                AMREX_IF_ON_DEVICE((
                    AMREX_DEVICE_PRINTF(" %d is out of bound (0:%d)", n, ncomp-1);
                ))
                AMREX_IF_ON_HOST((
                    std::stringstream ss;
                    ss << " " << n << " is out of bound: (0:" << ncomp-1 << ")";
                    amrex::Abort(ss.str());
                ))
            }
#endif
            return p[n*stride];
        }
    };

    /// \cond DOXYGEN_IGNORE
    namespace detail {

        template <int N> struct Stride { Long a[N] = {}; };
        template <> struct Stride<0> {};

        /**
         * \brief Helper struct to check if a type is a valid index type.
         *
         * Valid index types are integral types that can be safely converted to int without narrowing.
         */
        template <typename T>
        struct IsValidIndexType {
            static constexpr bool value = std::is_integral_v<T> && IsNonNarrowingConversion_v<T, int>;
        };

        /**
         * \brief Helper struct to check if index types are valid.
         *
         * Valid index types are integral types that can be safely converted to int
         * without narrowing. The last index type can also be an enum.
         */
        template <int N, bool last_dim_component, typename... idx>
        struct ArrayNDIndexCheck_impl {
        private:
            static constexpr int num_indices = sizeof...(idx);

            template <typename... Ts>
            struct AllExceptLastEnum;

            template <typename first, typename... Ts>
            struct AllExceptLastEnum<first, Ts...> {
                static constexpr bool value =
                    IsValidIndexType<std::decay_t<first>>::value
                    && AllExceptLastEnum<Ts...>::value;
            };

            // Check if the last index type is an enum or valid index type.
            // An enum is only allowed if last_dim_component is true.
            template <typename last>
            struct AllExceptLastEnum<last> {
                static constexpr bool value = IsValidIndexType<std::decay_t<last>>::value
                        || (std::is_enum_v<std::decay_t<last>>);
            };

            static constexpr bool index_with_maybe_enum = AllExceptLastEnum<idx...>::value;

            static constexpr bool index_all_values = Conjunction<
                IsValidIndexType<std::decay_t<idx>>...>::value;

        public:
            static constexpr bool value = ((num_indices == N) && index_all_values)      // N indices are integer types
                                || ((num_indices == N - 1) && index_all_values && last_dim_component) // N-1 indices are integer types, last is component
                                || ((num_indices == N) && index_with_maybe_enum && last_dim_component); // N indices, last dim is component and can be enum
        };

        template <int N, bool last_dim_component, class... idx>
        inline constexpr bool ArrayNDIndexCheck_impl_v = ArrayNDIndexCheck_impl<N, last_dim_component, idx...>::value;

        template<std::size_t... idx>
        constexpr auto make_oob_message_impl (std::index_sequence<idx...>) {
            constexpr std::size_t N = sizeof...(idx);
            constexpr char prefix[] = " (";
            constexpr char middle[] = ") is out of bound (";
            constexpr char suffix[] = ")\n";

            constexpr std::size_t size =
                sizeof(prefix) - 1 +
                (2*N + (N-1)) +
                sizeof(middle) - 1 +
                (5*N + (N-1)) +
                sizeof(suffix); // include null terminator

            std::array<char, size> buf{};

            std::size_t pos = 0;

            // prefix
            for (char c : prefix) {
                if (c != '\0') {
                    buf[pos++] = c;
                }
            }
            // %d,%d,...
            ((buf[pos++] = '%', buf[pos++] = 'd', idx + 1 < N ? buf[pos++] = ',' : 0), ...);

            // ") is out of bound ("
            for (char c : middle) {
                if (c != '\0') {
                    buf[pos++] = c;
                }
            }

            // %d:%d,%d:%d,...
            ((buf[pos++] = '%', buf[pos++] = 'd', buf[pos++] = ':', buf[pos++] = '%', buf[pos++] = 'd',
                idx + 1 < N ? buf[pos++] = ',' : 0), ...);

            // suffix
            for (char c : suffix) {
                buf[pos++] = c;
            }
            return buf;
        }

        template <std::size_t... idx, std::size_t... idx2x>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void device_print_impl2 (std::index_sequence<idx...>,
                                std::index_sequence<idx2x...>,
                                IntVectND<sizeof...(idx)> const& iv,
                                IntVectND<sizeof...(idx)> const& begin,
                                IntVectND<sizeof...(idx)> const& end)
        {
            constexpr auto msg = make_oob_message_impl(std::index_sequence<idx...>{});

            AMREX_DEVICE_PRINTF(msg.data(),
                iv.vect[idx]...,
                ((idx2x % 2 == 0) ? begin.vect[idx2x / 2] : end.vect[idx2x / 2])...
            );
        }

        template <int N, std::enable_if_t<(N>=1),int> = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void device_printf_impl (const IntVectND<N>& iv,
                                const IntVectND<N>& begin,
                                const IntVectND<N>& end)
        {
            device_print_impl2(
                std::make_index_sequence<N>{},
                std::make_index_sequence<N*2>{},
                iv, begin, end
            );
        }
    }
    /// \endcond

    /**
     * \ingroup amrex_fab_containers
     * \class ArrayND
     * \brief A multidimensional array accessor.
     * \tparam T The type of data (e.g., Real, int). May be const-qualified.
     * \tparam N The dimension of the array (rank).
     * \tparam last_dim_component If true, the last dimension represents
     *  components rather than a spatial index. In this case, the component
     *  index always starts at zero.
     *
     * ArrayND is a lightweight, non-owning wrapper around a pointer that provides
     * multidimensional indexing into a contiguous memory block. It is trivially
     * copyable and designed to be passed by value into GPU kernels. Note, ArrayND
     * follows fortran-style column-major ordering.
     *
     * It stores the pointer, stride information, and the bounds (begin and end)
     * of the array.
     *
     * \note \c Array4 is an alias of ArrayND with \p N=4 and
     *          \p last_dim_component=true. It represents three spatial dimensions
     *          plus one component dimension. When \c AMREX_SPACEDIM is less than 3,
     *          the extra spatial dimensions typically have bounds \c [0,1), and
     *          the correspding indices are zero.
     */
    template<typename T, int N, bool last_dim_component = false>
    struct ArrayND
    {
        static_assert(N >= 1, "ArrayND must have at least one dimension");
        static_assert(N > 1 || !last_dim_component, "ArrayND with N=1 cannot have last_dim_component=true");

        /// True if the last dimension is treated as components.
        static constexpr bool IsLastDimComponent_v = last_dim_component;
        /// True if this is an Array4 (\c N==4 and last dim is component).
        static constexpr bool IsArray4_v = (N==4 && last_dim_component);

        T* AMREX_RESTRICT p = nullptr; //!< Base pointer to the first element
        /// \cond DOXYGEN_IGNORE
        AMREX_NO_UNIQUE_ADDRESS detail::Stride<N-1> stride{};
        /// \endcond
        IntVectND<N> begin{1}; //!< Inclusive lower bounds
        IntVectND<N> end{0}; //!< Exclusive upper bounds

        /**
         * \brief Default-construct an empty accessor.
         *
         * The resulting accessor is invalid (\c operator \c bool() returns false).
         */
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND () noexcept : p(nullptr) {}

        /**
         * \brief Construct a const accessor from a non-const accessor.
         *
         * This enables implicit conversion from \c ArrayND<T,...> to
         * \c ArrayND<const T,...>.
         *
         * \param rhs Source accessor.
         */
        template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (ArrayND<std::remove_const_t<T>, N, last_dim_component> const& rhs) noexcept
            : p(rhs.p), stride(rhs.stride), begin(rhs.begin), end(rhs.end) {}

        /**
         * \brief Constructor using a BoxND.
         * \param a_p Pointer to data.
         * \param box The domain covered by this array.
         */
        template <bool C = last_dim_component, std::enable_if_t<!C,int> = 0>
        // TODO: Make BoxND functions constexpr to allow this constructor to be constexpr.
        AMREX_GPU_HOST_DEVICE
        ArrayND (T* a_p, BoxND<N> const& box) noexcept
            : ArrayND(a_p, box.smallEnd(), box.bigEnd() + 1)
        {}

        /**
         * \brief Constructor using a BoxND and the number of components.
         * \param a_p Pointer to data.
         * \param box BoxND<M> with M <= N-1 defining the spatial dimensions.
         * \param ncomp Number of components (size of last dimension).
         */
        template <int M, std::enable_if_t<((M+1==N) || (N == 4 && M == AMREX_SPACEDIM))
            && last_dim_component, int> = 0>
        // TODO: Make BoxND functions constexpr to allow this constructor to be constexpr.
        AMREX_GPU_HOST_DEVICE
        ArrayND (T* a_p, BoxND<M> const& box, int ncomp) noexcept
            : ArrayND(a_p, box.smallEnd(), box.bigEnd() + 1, ncomp)
        {}

        /**
         * \brief IntVectND<N> constructor.
         *
         * This constructor specifies begin and end indices using IntVectND<N> for all N dimensions.
         *
         * \param a_p Pointer to data.
         * \param a_begin Start index (inclusive).
         * \param a_end End index (exclusive).
         */
        template <bool C = last_dim_component, std::enable_if_t<!C,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (T* a_p, IntVectND<N> const& a_begin, IntVectND<N> const& a_end) noexcept
            : p(a_p), begin(a_begin), end(a_end)
        {
            set_stride();
        }

        /**
         * \brief Constructor for N=4 using Dim3.
         *
         * Typically used for standard (i, j, k, n) access in 3D simulations.
         *
         * \param a_p Pointer to data.
         * \param a_begin Start index.
         * \param a_end End index exclusive.
         * \param a_ncomp Number of components (4th dimension).
         */
        template <bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (T* a_p, Dim3 const& a_begin, Dim3 const& a_end, int a_ncomp) noexcept
            : p(a_p), begin(a_begin.x, a_begin.y, a_begin.z, 0), end(a_end.x, a_end.y, a_end.z, a_ncomp)
        {
            set_stride();
        }

        /**
         * \brief Reduced dimension constructor with component count.
         *
         * Pads intermediate dimensions with [0,1) and sets the last dimension to [0, ncomp).
         *
         * \param a_p Pointer to data.
         * \param a_begin Start index (inclusive).
         * \param a_end End index (exclusive).
         * \param ncomp Number of components (size of last dimension).
         */
        template <int M, std::enable_if_t<((M+1 == N) || (N == 4 && M == AMREX_SPACEDIM))
            && last_dim_component, int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (T* a_p, IntVectND<M> const& a_begin, IntVectND<M> const& a_end, int ncomp) noexcept
            : p(a_p)
        {
            constexpr_for<0, M>([&](int d) {
                begin.vect[d] = a_begin.vect[d];
                end.vect[d] = a_end.vect[d];
            });
            constexpr_for<M, N>([&](int d) {
                begin.vect[d] = 0;
                end.vect[d] = 1;
            });
            end.vect[N-1] = ncomp;

            set_stride();
        }

        /**
         * \brief Slicing constructor (Component subset).
         *
         * Creates a new ArrayND that is a slice of the input ArrayND along the last dimension (components).
         *
         * \param rhs The source ArrayND to slice from.
         * \param start_comp The starting component index for the slice.
         */
        template <class U,
            std::enable_if_t
                  <std::is_same_v<std::remove_const_t<T>, std::remove_const_t<U>>
                    && (N >= 2) && last_dim_component,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (ArrayND<U, N, last_dim_component> const& rhs, int start_comp) noexcept
            : p((T*)(rhs.p + start_comp*rhs.stride.a[N-2])),
              stride(rhs.stride),
              begin(rhs.begin),
              end(rhs.end)
        {
            begin.vect[N-1] = 0;
            end.vect[N-1] = rhs.end.vect[N-1] - start_comp;
        }

        /**
         * \brief Slicing constructor (Component subset with count).
         *
         * Creates a new ArrayND that is a slice of the input ArrayND along the last dimension (components).
         *
         * \param rhs The source ArrayND to slice from.
         * \param start_comp The starting component index for the slice.
         * \param num_comp The number of components to include in the slice.
         */
        template <class U,
            std::enable_if_t
                  <std::is_same_v<std::remove_const_t<T>, std::remove_const_t<U>>
                    && (N >= 2) && last_dim_component,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (ArrayND<U, N, last_dim_component> const& rhs, int start_comp, int num_comp) noexcept
            : p((T*)(rhs.p + start_comp*rhs.stride.a[N-2])),
              stride(rhs.stride),
              begin(rhs.begin),
              end(rhs.end)
        {
            begin.vect[N-1] = 0;
            end.vect[N-1] = num_comp;
        }

        /**
         * \brief Check if the ArrayND pointer is valid.
         * \return true if pointer is valid, false otherwise.
         */
        AMREX_GPU_HOST_DEVICE
        constexpr explicit operator bool() const noexcept { return p != nullptr; }

        /**
         * \brief Check if the ArrayND pointer is valid and bounds are valid.
         *
         * \return true if pointer and bounds are valid, false otherwise.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE
        constexpr bool ok () const noexcept { return p != nullptr && end.allGT(begin); }

        /**
         * \brief Multi-index operator() for accessing elements.
         *
         * This operator is only enabled when the number of indices provided matches the dimension N or
         * N-1 if last_dim_component is true, where in the latter case the last dimension is assumed to be 0.
         * The passed indices must be of integral types that can be converted to int
         * without narrowing. The last index can also be an enum type.
         *
         * \note This overload is disabled for Array4, which has its own overload of
         * \c operator() that takes three integer indices for spatial dimensions and
         * an optional integer component index.
         *
         * \tparam idx Variadic template for indices.
         * \param i Variadic list of indices.
         * \return Reference to the element.
         */
        template <typename... idx, class U=T,
            std::enable_if_t<!std::is_void_v<U> && !IsArray4_v
                && detail::ArrayNDIndexCheck_impl_v<N, last_dim_component, idx...>,int> = 0>
        // TODO: Use concepts when we move to C++20.
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (idx... i) const noexcept {
            constexpr auto nidx = sizeof...(i);
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(IntVectND<nidx>{i...});
#endif
            return p[get_offset(IntVectND<nidx>{i...})];
        }

        /**
         * \brief Access element by IntVectND.
         *
         * This function handles access via:
         * 1. Full dimension vector (M == N).
         * 2. Spatial dimension vector (M == N-1) && !IsArray4_v
         *    (only enabled if last_dim_component is true).
         *
         * \tparam M The dimension of the index vector.
         * \param iv The index vector.
         * \return Reference to the element.
         */
        template <int M, class U=T, std::enable_if_t<
            !std::is_void_v<U> &&
            ((M == N) || (!IsArray4_v && last_dim_component && (M + 1 == N))), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVectND<M> const& iv) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv);
#endif
            return p[get_offset(iv)];
        }

        /**
         * \brief Access element by spatial IntVectND and component index.
         *
         * This overload is enabled only when last_dim_component is true and M represents
         * the spatial dimensions (either N-1 or AMREX_SPACEDIM).
         *
         * \tparam M The dimension of the spatial index vector (typically N-1 or AMREX_SPACEDIM).
         * \param iv The spatial index vector.
         * \param n Component index (last dimension).
         * \return Reference to the element.
         */
        template <int M, class U=T, std::enable_if_t<
            !std::is_void_v<U> && last_dim_component && !IsArray4_v
            && (M + 1 == N), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVectND<M> const& iv, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv, n);
#endif
            return p[get_offset(iv, n)];
        }

        /**
         * \brief Multi-index ptr() for accessing pointer to element.
         * \tparam idx Variadic template for indices.
         * \param i Variadic list of indices.
         * \return Pointer to the element.
         */
        template <typename... idx, class U=T,
            std::enable_if_t<!std::is_void_v<U> && !IsArray4_v
            && detail::ArrayNDIndexCheck_impl_v<N, last_dim_component, idx...>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (idx... i) const noexcept {
            constexpr auto nidx = sizeof...(i);
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(IntVectND<nidx>{i...});
#endif
            return p + get_offset(IntVectND<nidx>{i...});
        }

        /**
         * \brief Access pointer by IntVectND.
         *
         * This function handles access via:
         * 1. Full dimension vector (M == N).
         * 2. Spatial dimension vector (M == N-1) && !IsArray4_v
         *    (only enabled if last_dim_component is true).
         *
         * \tparam M The dimension of the index vector.
         * \param iv The index vector.
         * \return Pointer to the element.
         */
        template <int M, std::enable_if_t<(M == N)
            || (!IsArray4_v && last_dim_component && (M + 1 == N)), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVectND<M> const& iv) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv);
#endif
            return p + get_offset(iv);
        }

        /**
         * \brief Access pointer by spatial IntVectND and component index.
         *
         * This overload is enabled only when last_dim_component is true and M represents
         * the spatial dimensions (either N-1 or AMREX_SPACEDIM).
         *
         * \tparam M The dimension of the spatial index vector (typically N-1 or AMREX_SPACEDIM).
         * \param iv The spatial index vector.
         * \param n Component index (last dimension).
         * \return Pointer to the element.
         */
        template <int M, std::enable_if_t<last_dim_component && !IsArray4_v
            && (M + 1 == N), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVectND<M> const& iv, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv, n);
#endif
            return p + get_offset(iv, n);
        }

        /**
         * \brief Get raw data pointer.
         * \return Raw data pointer.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr T* dataPtr () const noexcept {
            return this->p;
        }

        /**
         * \brief Get number of components
         * \return Number of components.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr int nComp () const noexcept {
            if constexpr (last_dim_component) {
                return end.vect[N-1];
            } else {
                return 1;
            }
        }

        /**
         * \brief Total number of elements in the ArrayND's index region.
         *
         * Returns the product of extents in each dimension. If the object
         * is not valid (see \c ok()), this returns 0.
         *
         * \return Number of elements in the ArrayND, or 0 if invalid.
         */
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr std::size_t size () const noexcept {
            if (ok()) {
                std::size_t s = 1;
                constexpr_for<0, N>([&](int d) {
                    s *= (end.vect[d] - begin.vect[d]);
                });
                return s;
            } else {
                return 0;
            }
        }

        /**
         * \brief Return the stride (in elements) for dimension \p d.
         *
         * This returns the multiplier used to advance the linear address when
         * incrementing index \p d by one. For \c d==0, the stride is 1.
         *
         * \tparam d Dimension index in the range \c [0, N).
         *
         * \return Stride for dimension \p d, in number of elements.
         */
        template <int d, std::enable_if_t<(d < N) && (d >= 0), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Long get_stride () const noexcept {
            if constexpr (N > 1 && d > 0) {
                return stride.a[d-1];
            } else {
                return 1;
            }
        }

        /**
         * \brief Test whether an index tuple lies inside the ArrayND bounds.
         *
         * This overload accepts a pack of integer indices and it is
         * disabled for Array4, which has its own overload of \c contains
         * that takes three integer indices for spatial dimensions.
         *
         * \tparam idx Index argument types.
         * \param i Indices to test.
         * \return \c true if all indices are within bounds.
         */
        template <typename... idx, std::enable_if_t<!IsArray4_v &&
            detail::ArrayNDIndexCheck_impl_v<N, last_dim_component, idx...>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr bool contains (idx... i) const noexcept {
            constexpr auto nidx = sizeof...(i);
            return this->contains(IntVectND<nidx>{i...});
        }

        /**
         * \brief Test whether an IntVectND lies inside the ArrayND bounds.
         *
         * Checks whether each component of \p iv is within the half-open bounds
         * \c [begin, end) for the corresponding dimensions.
         *
         * This overload supports two cases:
         * - \c M==N : the full index space is checked; if the last dimension is a
         *   component dimension, it is checked as well. This can be used \c Array4.
         * - \c (M+1==N) with \c last_dim_component=true : only the spatial dimensions
         *   are checked (component handled separately).
         *
         * \tparam M Dimension of the input IntVectND.
         * \param[in] iv IntVectND to test.
         * \return \c true if \p iv is inside the bounds, otherwise \c false.
         */
        template <int M, std::enable_if_t<(M == N)
            || (!IsArray4_v && last_dim_component && (M + 1 == N)), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr bool contains (IntVectND<M> const& iv) const noexcept {
            bool inside = true;
            constexpr_for<0, M>([&](int d) {
                inside = inside && (iv.vect[d] >= begin.vect[d]) && (iv.vect[d] < end.vect[d]);
            });
            return inside;
        }

        /**
         * \brief Test whether a spatial index and component lie inside the bounds.
         *
         * Checks that the spatial index \p iv is within \c [begin, end) for the
         * spatial dimensions and that the component index \p n satisfies
         * \c 0 <= n < end[N-1].
         *
         * This overload is enabled only when the last dimension is a component
         * dimension.
         *
         * \tparam M Spatial dimension of the input IntVectND.
         * \param[in] iv Spatial IntVectND to test.
         * \param[in] n  Component index to test.
         * \return \c true if both the spatial index and component are inside,
         *         otherwise \c false.
         */
        template <int M, std::enable_if_t<last_dim_component && !IsArray4_v
            && (M + 1 == N), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr bool contains (IntVectND<M> const& iv, int n) const noexcept {
            bool inside = true;
            constexpr_for<0, M>([&](int d) {
                inside = inside && (iv.vect[d] >= begin.vect[d]) && (iv.vect[d] < end.vect[d]);
            });
            inside = inside && (n >= 0) && (n < end.vect[N-1]);
            return inside;
        }

        /**
         * \brief Create a single-cell component accessor.
         *
         * Returns a \c CellData<T> view for the cell \f$(i,j,k)\f$ in an
         * Array4. The returned object provides access to all components
         * without recomputing full \c (i,j,k,n) indexing.
         *
         * \param i Cell index in x.
         * \param j Cell index in y.
         * \param k Cell index in z.
         * \return A \c CellData<T> accessor for cell \f$(i,j,k)\f$.
         *
         * \note This overload is enabled only for Array4.
         *
         * \sa CellData
         */
        template <int M=N, bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CellData<T> cellData (int i, int j, int k) const noexcept {
            int ncomp = end.vect[M-1];
            return CellData<T>(this->ptr(i,j,k), stride.a[M-2], ncomp);
        }

        /**
         * \brief Compute the linear offset (in elements) for an IntVectND.
         *
         * \tparam M Dimension of IntVectND.
         * \param[in] iv IntVectND.
         * \return Linear offset in number of elements from base pointer \c p.
         */
        template <int M, std::enable_if_t<(M == N)
            || (last_dim_component && (M + 1 == N || M == AMREX_SPACEDIM)), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Long get_offset (IntVectND<M> const& iv) const noexcept
        {
            Long offset = iv.vect[0] - begin.vect[0];
            // If M == N and we have a component at the end, we only loop up to N-1 for spatial.
            // Otherwise, we loop up to M.
            constexpr int idx = (last_dim_component && M == N) ? N - 1 : M;
            if constexpr (N > 1) {
                constexpr_for<1, idx>([&](int d) {
                    offset += (iv.vect[d] - begin.vect[d]) * stride.a[d-1];
                });
            }

            // Handle the component offset if the input is the full N dimensions
            // and the last dimension is a component.
            if constexpr (last_dim_component && M == N) {
                offset += iv.vect[N-1] * stride.a[N-2];
            }
            return offset;
        }

        /**
         * \brief Compute the linear offset (in elements) for an IntVectND and
         *        component index.
         *
         * \tparam M Dimension of IntVectND.
         * \param iv IntVectND for spatial part.
         * \param n  Component index.
         * \return Linear offset in number of elements from base pointer \c p.
         */
        template <int M, std::enable_if_t<last_dim_component
            && ((M + 1 == N) || (N == 4 && M == AMREX_SPACEDIM)), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Long get_offset (IntVectND<M> const& iv, int n) const noexcept
        {
            Long offset = iv.vect[0] - begin.vect[0];
            constexpr_for<1, M>([&](int d) {
                offset += (iv.vect[d] - begin.vect[d]) * stride.a[d-1];
            });
            offset += n * stride.a[N-2];
            return offset;
        }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE
#endif
        void index_assert (IntVectND<N> const& iv) const
        {
            bool out_of_bounds = false;
            for (int d = 0; d < N; ++d) {
                if (iv.vect[d] < begin.vect[d] || iv.vect[d] >= end.vect[d]) {
                    out_of_bounds = true;
                }
            }
            if (out_of_bounds) {
                AMREX_IF_ON_DEVICE((
                    detail::device_printf_impl(iv, begin, end);
                    amrex::Abort();
                ))
                AMREX_IF_ON_HOST((
                    std::stringstream ss;
                    ss << " (";
                    for (int d = 0; d < N; ++d) {
                        ss << iv.vect[d];
                        if (d + 1 < N) ss << ",";
                    }
                    ss << ") is out of bound (";
                    for (int d = 0; d < N; ++d) {
                        ss << begin.vect[d] << ":" << end.vect[d]-1;
                        if (d + 1 < N) ss << ",";
                    };
                    ss << ")";
                    amrex::Abort(ss.str());
                ))
            }
        }

        // index_assert overload for M == N-1 last index assumed 0
        // Only valid when last_dim_component == true
        template <int M, std::enable_if_t<((M+1 == N) || (M == AMREX_SPACEDIM))
            && last_dim_component,int> = 0>
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE
#endif
        void index_assert (IntVectND<M> const& iv) const
        {
            IntVectND<N> iv_full = iv.template expand<N>(0);
            for (int d = M; d < N; ++d) {
                iv_full.vect[d] = begin.vect[d];
            }
            index_assert(iv_full);
        }

        // index_assert overload for M == N-1 and last index n
        // Only valid when last_dim_component == true
        template <int M, std::enable_if_t<((M+1 == N) || (M == AMREX_SPACEDIM))
            && last_dim_component,int> = 0>
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE
#endif
        void index_assert (IntVectND<M> const& iv, int n) const
        {
            IntVectND<N> iv_full = iv.template expand<N>(0);
            for (int d = M; d < N-1; ++d) {
                iv_full.vect[d] = begin.vect[d];
            }
            iv_full.vect[N-1] = n;
            index_assert(iv_full);
        }
#endif

        //
        // Specialization for Array4
        //

        /**
         * \brief Access Array4 element using three integer spatial indices.
         *
         * The component index is implicitly assumed to be zero. When
         * \c AMREX_SPACEDIM is less than 3, the bounds for extra spatial dimensions
         * are typically \c [0,1). In particular, when launching a kernel with
         * \c amrex::ParallelFor(Box, f), where \c f is a lambda function with the
         * signature of \c void(int i,int j,int k), AMReX supplies zero for the
         * extra dimensions. Thus this operator is safe under the usual AMReX
         * execution model.
         *
         * \param i Cell index in x.
         * \param j Cell index in y.
         * \param k Cell index in z.
         *
         * \note This overload is enabled only for Array4.
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            AMREX_ARRAY4_INDEX_ASSERT(i,j,k,0);
#endif
#if defined(AMREX_USE_GPU) || defined(AMREX_DEBUG)
            return p[(i-begin.vect[0]) +
                     (j-begin.vect[1]) * stride.a[0] +
                     (k-begin.vect[2]) * stride.a[1]];
#else
            Long idx1 = i + j*stride.a[0] + k*stride.a[1];
            Long idx0 = begin.vect[0]
                +       begin.vect[1] * stride.a[0]
                +       begin.vect[2] * stride.a[1];
            return p[idx1-idx0];
#endif
        }

        /**
         * \brief Access Array4 element using three integer spatial indices and
         *        component index.
         *
         * When \c AMREX_SPACEDIM is less than 3, the bounds for extra spatial
         * dimensions are typically \c [0,1). For example, when AMREX_SPADIM is 2,
         * \p k should typically be zero.
         *
         * \param i Cell index in x.
         * \param j Cell index in y.
         * \param k Cell index in z.
         * \param n Component index.
         *
         * \note This overload is enabled only for Array4.
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            AMREX_ARRAY4_INDEX_ASSERT(i,j,k,n);
#endif
#if defined(AMREX_USE_GPU) || defined(AMREX_DEBUG)
            return p[(i-begin.vect[0]) +
                     (j-begin.vect[1]) * stride.a[0] +
                     (k-begin.vect[2]) * stride.a[1] +
                                     n * stride.a[2]];
#else
            Long idx1 = i + j*stride.a[0] + k*stride.a[1] + n*stride.a[2];
            Long idx0 = begin.vect[0]
                +       begin.vect[1] * stride.a[0]
                +       begin.vect[2] * stride.a[1];
            return p[idx1-idx0];
#endif
        }

        /**
         * \brief Access Array4 element by spatial indices in IntVectND.
         *
         * The component index is implicitly assumed to be zero. The dimension of
         * \p iv is either 3 or \c AMREX_SPACEDIM. If it's less than 3, zero is used
         * in the extra dimensions.
         *
         * \param iv Spatial indices.
         */
        template <int M, class U=T,
                  std::enable_if_t<!std::is_void_v<U> && IsArray4_v && (M == 3 || M == AMREX_SPACEDIM), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVectND<M> const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            if constexpr (M == 1) {
                return this->operator()(iv.vect[0],0,0);
            } else
#elif (AMREX_SPACEDIM == 2)
            if constexpr (M == 2) {
                return this->operator()(iv.vect[0],iv.vect[1],0);
            } else
#endif
            {
                return this->operator()(iv.vect[0],iv.vect[1],iv.vect[2]);
            }
        }

        /**
         * \brief Access Array4 element by spatial indices in IntVectND and
         *        component index.
         *
         * The dimension of \p iv is either 3 or \c AMREX_SPACEDIM. If it's less
         * than 3, zero is used in the extra dimensions.
         *
         * \param iv Spatial indices.
         * \param n  Component index.
         */
        template <int M, class U=T,
                  std::enable_if_t<!std::is_void_v<U> && IsArray4_v && (M == 3 || M == AMREX_SPACEDIM), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVectND<M> const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            if constexpr (M == 1) {
                return this->operator()(iv.vect[0],0,0,n);
            } else
#elif (AMREX_SPACEDIM == 2)
            if constexpr (M == 2) {
                return this->operator()(iv.vect[0],iv.vect[1],0,n);
            } else
#endif
            {
                return this->operator()(iv.vect[0],iv.vect[1],iv.vect[2],n);
            }
        }

        /**
         * \brief Access \c Array4 element by spatial indices in \c Dim3.
         *
         * The component index is implicitly assumed to be zero.
         *
         * \param cell Spatial indices.
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z);
        }

        /**
         * \brief Access \c Array4 element by spatial indices in \c Dim3 and
         *        component index.
         *
         * \param cell Spatial indices.
         * \param n    Component index.
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell, int n) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z,n);
        }

        /**
         * \brief Return the pointer to an element
         *
         * \param i Cell index in x.
         * \param j Cell index in y.
         * \param k Cell index in z.
         * \return Pointer to element at \c (i,j,k).
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            AMREX_ARRAY4_INDEX_ASSERT(i,j,k,0);
#endif
#if defined(AMREX_USE_GPU) || defined(AMREX_DEBUG)
            return p + ((i-begin.vect[0]) +
                        (j-begin.vect[1]) * stride.a[0] +
                        (k-begin.vect[2]) * stride.a[1]);
#else
            Long idx1 = i + j*stride.a[0] + k*stride.a[1];
            Long idx0 = begin.vect[0]
                +       begin.vect[1] * stride.a[0]
                +       begin.vect[2] * stride.a[1];
            return p + (idx1-idx0);
#endif
        }

        /**
         * \brief Return the pointer to an element
         *
         * \param i Cell index in x.
         * \param j Cell index in y.
         * \param k Cell index in z.
         * \param n Component index.
         * \return Pointer to element at \c (i,j,k,n).
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            AMREX_ARRAY4_INDEX_ASSERT(i,j,k,n);
#endif
#if defined(AMREX_USE_GPU) || defined(AMREX_DEBUG)
            return p + ((i-begin.vect[0]) +
                        (j-begin.vect[1]) * stride.a[0] +
                        (k-begin.vect[2]) * stride.a[1] +
                                        n * stride.a[2]);
#else
            Long idx1 = i + j*stride.a[0] + k*stride.a[1] + n*stride.a[2];
            Long idx0 = begin.vect[0]
                +       begin.vect[1] * stride.a[0]
                +       begin.vect[2] * stride.a[1];
            return p + (idx1-idx0);
#endif
        }

        /**
         * \brief Return the pointer to an element
         *
         * \param iv Spatial indices
         * \return Pointer to element at \c (iv).
         */
        template <int M, class U=T,
                  std::enable_if_t<!std::is_void_v<U> && IsArray4_v && (M == 3 || M == AMREX_SPACEDIM), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (IntVectND<M> const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            if constexpr (M == 1) {
                return this->ptr(iv.vect[0],0,0);
            } else
#elif (AMREX_SPACEDIM == 2)
            if constexpr (M == 2) {
                return this->ptr(iv.vect[0],iv.vect[1],0);
            } else
#endif
            {
                return this->ptr(iv.vect[0],iv.vect[1],iv.vect[2]);
            }
        }

        /**
         * \brief Return the pointer to an element
         *
         * \param iv Spatial indices.
         * \param n  Component index.
         * \return Pointer to element at \p (iv,n).
         */
        template <int M, class U=T,
                  std::enable_if_t<!std::is_void_v<U> && IsArray4_v && (M == 3 || M == AMREX_SPACEDIM), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (IntVectND<M> const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            if constexpr (M == 1) {
                return this->ptr(iv.vect[0],0,0,n);
            } else
#elif (AMREX_SPACEDIM == 2)
            if constexpr (M == 2) {
                return this->ptr(iv.vect[0],iv.vect[1],0,n);
            } else
#endif
            {
                return this->ptr(iv.vect[0],iv.vect[1],iv.vect[2],n);
            }
        }

        /**
         * \brief Return the pointer to an element
         *
         * \param cell Spatial indices
         * \return Pointer to element at \c (cell).
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (Dim3 const& cell) const noexcept {
            return this->ptr(cell.x,cell.y,cell.z);
        }

        /**
         * \brief Return the pointer to an element
         *
         * \param cell Spatial indices.
         * \param n    Component index.
         * \return Pointer to element at \c (cell,n).
         */
        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (Dim3 const& cell, int n) const noexcept {
            return this->ptr(cell.x,cell.y,cell.z,n);
        }

        /**
         * \brief Test whether the spatial indices are inside the Array4 bounds.
         *
         * \param i Cell index in x.
         * \param j Cell index in y.
         * \param k Cell index in z.
         * \return \c true if all indices are within bounds.
         */
        template <bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (int i, int j, int k) const noexcept {
            return (i>=begin.vect[0] && i<end.vect[0] &&
                    j>=begin.vect[1] && j<end.vect[1] &&
                    k>=begin.vect[2] && k<end.vect[2]);
        }

        /**
         * \brief Test whether the spatial indices are inside the Array4 bounds.
         *
         * \param iv Spatial indices.
         * \return \c true if all indices are within bounds.
         */
        template <int M, std::enable_if_t<IsArray4_v && (M == 3 || M == AMREX_SPACEDIM), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (IntVectND<M> const& iv) const noexcept {
#if (AMREX_SPACEDIM < 3)
            if constexpr (M == AMREX_SPACEDIM) {
                return AMREX_D_TERM((iv.vect[0]>=begin.vect[0]) && (iv.vect[0]<end.vect[0]),
                                 && (iv.vect[1]>=begin.vect[1]) && (iv.vect[1]<end.vect[1]),
                                 && (iv.vect[2]>=begin.vect[2]) && (iv.vect[2]<end.vect[2]));
            } else
#endif
            {
                return (iv.vect[0]>=begin.vect[0]) && (iv.vect[0]<end.vect[0])
                    && (iv.vect[1]>=begin.vect[1]) && (iv.vect[1]<end.vect[1])
                    && (iv.vect[2]>=begin.vect[2]) && (iv.vect[2]<end.vect[2]);
            }
        }

        /**
         * \brief Test whether the spatial indices are inside the Array4 bounds.
         *
         * \param cell Spatial indices.
         * \return \c true if all indices are within bounds.
         */
        template <bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (Dim3 const& cell) const noexcept {
            return (cell.x>=begin.vect[0] && cell.x<end.vect[0] &&
                    cell.y>=begin.vect[1] && cell.y<end.vect[1] &&
                    cell.z>=begin.vect[2] && cell.z<end.vect[2]);
        }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        template <bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE
#endif
        void index_assert_print_error_message (int i, int j, int k, int n) const
        {
            AMREX_IF_ON_DEVICE((
                    AMREX_DEVICE_PRINTF(" (%d,%d,%d,%d) is out of bound (%d:%d,%d:%d,%d:%d,0:%d)\n",
                                        i, j, k, n,
                                        begin.vect[0], end.vect[0]-1,
                                        begin.vect[1], end.vect[1]-1,
                                        begin.vect[2], end.vect[2]-1,
                                                       end.vect[3]-1);
                    amrex::Abort();
            ))
            AMREX_IF_ON_HOST((
                    std::stringstream ss;
                    ss << " (" << i << "," << j << "," << k << "," <<  n
                    << ") is out of bound ("
                    << begin.vect[0] << ":" << end.vect[0]-1 << ","
                    << begin.vect[1] << ":" << end.vect[1]-1 << ","
                    << begin.vect[2] << ":" << end.vect[2]-1 << ","
                    << "0:" << end.vect[3]-1 << ")";
                    amrex::Abort(ss.str());
            ))
        }
#endif

    private:
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr void set_stride () noexcept {
            if constexpr (N > 1) {
                Long current_stride = 1;
                constexpr_for<0, N-1>([&](int d) {
                    Long len = end.vect[d] - begin.vect[d];
                    current_stride *= len;
                    stride.a[d] = current_stride;
                });
            }
        }
    };

    // Deduction guides for ArrayND
    // 1: Matches ArrayND (T*, BoxND<N>) -> ArrayND<T,N, last_dim_component=false>
    template <typename T, int N>
    ArrayND (T*, BoxND<N> const&) -> ArrayND<T, N, false>;

    // 2: Matches ArrayND (T*, BoxND<N>, ncomp) -> ArrayND<T,N+1, last_dim_component=true>
    // This supports the "N+1" logic (Spatial + Component)
    template <typename T, int N>
    ArrayND (T*, BoxND<N> const&, int) -> ArrayND<T, N+1, true>;

    // 3: Matches ArrayND (T*, IntVectND<N>, IntVectND<N>) -> ArrayND<T,N, last_dim_component=false>
    template <typename T, int N>
    ArrayND (T*, IntVectND<N> const&, IntVectND<N> const&) -> ArrayND<T, N, false>;

    // 4: Matches ArrayND (T*, IntVectND<N>, IntVectND<N>, int) -> ArrayND<T,N+1, last_dim_component=true>
    template <typename T, int N>
    ArrayND (T*, IntVectND<N> const&, IntVectND<N> const&, int) -> ArrayND<T, N+1, true>;

    // 5: Matches ArrayND (T*, Dim3, Dim3) -> ArrayND<T,4, last_dim_component=true>
    template <typename T>
    ArrayND (T*, Dim3 const&, Dim3 const&, int) -> ArrayND<T, 4, true>;

    /**
     * \ingroup amrex_fab_containers
     * \brief Convenient alias for a 4D \c ArrayND with three spatial dimensions
     *        plus a component dimension.
     */
    template<typename T>
    using Array4 = ArrayND<T, 4, true>;

    /**
     * \ingroup amrex_fab_containers
     * \brief Return the inclusive lower bounds of an Array4 in \c Dim3 form.
     *
     * \param a Array4
     * \return Inclusive lower bounds.
     */
    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 lbound (Array4<T> const& a) noexcept
    {
        return Dim3{a.begin.vect[0],a.begin.vect[1],a.begin.vect[2]};
    }

    /**
     * \ingroup amrex_fab_containers
     * \brief Return the inclusive upper bounds of an Array4 in \c Dim3 form.
     *
     * \param a Array4
     * \return Inclusive upper bounds.
     */
    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 ubound (Array4<T> const& a) noexcept
    {
        return Dim3{a.end.vect[0]-1,a.end.vect[1]-1,a.end.vect[2]-1};
    }

    /**
     * \ingroup amrex_fab_containers
     * \brief Return the spatial extents of an Array4 in \c Dim3 form.
     *
     * \param a Array4
     * \return Lengths
     */
    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 length (Array4<T> const& a) noexcept
    {
        return Dim3{a.end.vect[0]-a.begin.vect[0],a.end.vect[1]-a.begin.vect[1],a.end.vect[2]-a.begin.vect[2]};
    }

    /**
     * \ingroup amrex_fab_containers
     * \brief Stream an ArrayND as its index bounds.
     *
     * Prints a compact representation of the index bounds, e.g.
     * \c ((lo...), (hi...)) with \c hi reported as inclusive (i.e., \c end-1).
     *
     * \param os Output stream.
     * \param a  ArrayND to print.
     * \return \p os.
     */
    template <typename T, int N, bool C>
    std::ostream& operator<< (std::ostream& os, const ArrayND<T,N,C>& a) {
        os << "(" << a.begin << ',' << a.end-1 << ")";
        return os;
    }

    //
    // Type traits for detecting if a class has a size() constexpr function.
    //
    template <class A, class Enable = void> struct HasMultiComp : std::false_type {};
    //
    template <class B>
    struct HasMultiComp<B, std::enable_if_t<B().size() >= 1>>
        : std::true_type {};

    //
    // PolymorphicArray4 can be used to access both AoS and SoA with
    // (i,j,k,n).  Here SoA refers multi-component BaseFab, and AoS refers
    // to single-component BaseFab of multi-component GpuArray.
    //
    template <typename T>
    struct PolymorphicArray4
        : public Array4<T>
    {
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        PolymorphicArray4 (Array4<T> const& a)
            : Array4<T>{a} {}

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T& operator() (int i, int j, int k) const noexcept {
            return this->Array4<T>::operator()(i,j,k);
        }

        template <class U=T, std::enable_if_t< amrex::HasMultiComp<U>::value,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        typename U::reference_type
        operator() (int i, int j, int k, int n) const noexcept {
            return this->Array4<T>::operator()(i,j,k,0)[n];
        }

        template <class U=T, std::enable_if_t<!amrex::HasMultiComp<U>::value,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k, int n) const noexcept {
            return this->Array4<T>::operator()(i,j,k,n);
        }
    };

    template <typename T>
    [[nodiscard]] PolymorphicArray4<T>
    makePolymorphic (Array4<T> const& a)
    {
        return PolymorphicArray4<T>(a);
    }
}

#endif
