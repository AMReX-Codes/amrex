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

    template <typename T>
    struct CellData // Data in a single cell
    {
        T* AMREX_RESTRICT p = nullptr;
        Long stride = 0;
        int ncomp = 0;

        AMREX_GPU_HOST_DEVICE
        constexpr CellData (T* a_p, Long a_stride, int a_ncomp)
            : p(a_p), stride(a_stride), ncomp(a_ncomp)
            {}

        template <class U=T,
                  std::enable_if_t<std::is_const_v<U>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr CellData (CellData<std::remove_const_t<T>> const& rhs) noexcept
            : p(rhs.p), stride(rhs.stride), ncomp(rhs.ncomp)
            {}

        AMREX_GPU_HOST_DEVICE
        explicit operator bool() const noexcept { return p != nullptr; }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE
        int nComp() const noexcept { return ncomp; }

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
                iv[idx]...,
                ((idx2x % 2 == 0) ? begin[idx2x / 2] : end[idx2x / 2])...
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

    /**
     * \brief A multidimensional array accessor.
     * \tparam T The type of data (e.g., Real, int).
     * \tparam N The dimension of the array (rank).
     * \tparam last_dim_component If true, the last dimension is treated as components.
     *
     * ArrayND is a lightweight, non-owning wrapper around a pointer that provides
     * multidimensional indexing into a contiguous memory block. It is trivially
     * copyable and designed to be passed by value into GPU kernels. Note, ArrayND
     * follows fortran-style column-major ordering.
     *
     * It stores the pointer, stride information, and the bounds (begin and end)
     * of the array.
     */
    template<typename T, int N, bool last_dim_component = false>
    struct ArrayND
    {
        static_assert(N >= 1, "ArrayND must have at least one dimension");
        static_assert(N > 1 || !last_dim_component, "ArrayND with N=1 cannot have last_dim_component=true");

        static constexpr bool IsLastDimComponent_v = last_dim_component;
        static constexpr bool IsArray4_v = (N==4 && last_dim_component);

        T* AMREX_RESTRICT p = nullptr;
        AMREX_NO_UNIQUE_ADDRESS detail::Stride<N-1> stride{};
        IntVectND<N> begin{1};
        IntVectND<N> end{0}; // end is hi+1

        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND () noexcept : p(nullptr) {}

        template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (ArrayND<std::remove_const_t<T>, N, last_dim_component> const& rhs) noexcept
            : p(rhs.p), stride(rhs.stride), begin(rhs.begin), end(rhs.end) {}

        /**
         * \brief Constructor using a BoxND.
         * \param a_p Pointer to data.
         * \param box The domain covered by this array.
         *
         * TODO: Make BoxND functions constexpr to allow this constructor to be constexpr.
         */
        AMREX_GPU_HOST_DEVICE
        ArrayND (T* a_p, BoxND<N> const& box) noexcept
            : ArrayND(a_p, box.smallEnd(), box.bigEnd() + 1)
        {}

        /**
         * \brief Constructor using a BoxND.
         * \param a_p Pointer to data.
         * \param box BoxND<M> with M <= N-1 defining the spatial dimensions.
         * \param ncomp Number of components (size of last dimension).
         *
         * TODO: Make BoxND functions constexpr to allow this constructor to be constexpr.
         */
        template <int M, std::enable_if_t<(M+1==N) && last_dim_component, int> = 0>
        AMREX_GPU_HOST_DEVICE
        ArrayND (T* a_p, BoxND<M> const& box, int ncomp) noexcept
            : ArrayND(a_p, box.smallEnd(), box.bigEnd() + 1, ncomp)
        {}


        /**
         * \brief IntVectND<N> constructor.
         * \param a_p Pointer to data.
         * \param a_begin Start index (inclusive).
         * \param a_end End index (exclusive).
         *
         * This constructor specifies begin and end indices using IntVectND<N> for all N dimensions.
         */
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (T* a_p, IntVectND<N> const& a_begin, IntVectND<N> const& a_end) noexcept
            : p(a_p), begin(a_begin), end(a_end)
        {
            if constexpr (N > 1) {
                Long current_stride = 1;
                constexpr_for<0, N-1>([&](int d) {
                    Long len = a_end[d] - a_begin[d];
                    current_stride *= len;
                    stride.a[d] = current_stride;
                });
            }
        }

        /**
         * \brief Constructor for N=4 using Dim3.
         * \param a_p Pointer to data.
         * \param a_begin Start index.
         * \param a_end End index exclusive.
         * \param a_ncomp Number of components (4th dimension).
         *
         * Typically used for standard (i, j, k, n) access in 3D simulations.
         */
        template <bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (T* a_p, Dim3 const& a_begin, Dim3 const& a_end, int a_ncomp) noexcept
            : p(a_p), begin(a_begin.x, a_begin.y, a_begin.z, 0), end(a_end.x, a_end.y, a_end.z, a_ncomp)
        {
            Long current_stride = 1;
            constexpr_for<0, N-1>([&](int d) {
                Long len = end[d] - begin[d];
                current_stride *= len;
                stride.a[d] = current_stride;
            });
        }

        /**
         * \brief Reduced dimension constructor with component count.
         * \param a_p Pointer to data.
         * \param a_begin Start index (inclusive).
         * \param a_end End index (exclusive).
         * \param ncomp Number of components (size of last dimension).
         *
         * Pads intermediate dimensions with [0,1) and sets the last dimension to [0, ncomp).
         */
        template <int M, std::enable_if_t<((M+1 == N) || (N == 4 && M == AMREX_SPACEDIM))
            && last_dim_component, int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (T* a_p, IntVectND<M> const& a_begin, IntVectND<M> const& a_end, int ncomp) noexcept
            : p(a_p)
        {
            constexpr_for<0, M>([&](int d) {
                begin[d] = a_begin[d];
                end[d] = a_end[d];
            });
            constexpr_for<M, N>([&](int d) {
                begin[d] = 0;
                end[d] = 1;
            });
            end[N-1] = ncomp;

            if constexpr (N > 1) {
                Long current_stride = 1;
                constexpr_for<0, N-1>([&](int d) {
                    Long len = end[d] - begin[d];
                    current_stride *= len;
                    stride.a[d] = current_stride;
                });
            }
        }

        /**
         * \brief Slicing constructor (Component subset).
         * \param rhs The source ArrayND to slice from.
         * \param start_comp The starting component index for the slice.
         *
         * Creates a new ArrayND that is a slice of the input ArrayND along the last dimension (components).
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
            begin[N-1] = 0;
            end[N-1] = rhs.end[N-1] - start_comp;
        }

        /**
         * \brief Slicing constructor (Component subset with count).
         * \param rhs The source ArrayND to slice from.
         * \param start_comp The starting component index for the slice.
         * \param num_comp The number of components to include in the slice.
         *
         * Creates a new ArrayND that is a slice of the input ArrayND along the last dimension (components).
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
            begin[N-1] = 0;
            end[N-1] = num_comp;
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
         * \tparam idx Variadic template for indices.
         * \param i Variadic list of indices.
         * \return Reference to the element.
         *
         * This operator is only enabled when the number of indices provided matches the dimension N or
         * N-1 if last_dim_component is true, where in the latter case the last dimension is assumed to be 0.
         * The passed indices must be of integral types that can be converted to int
         * without narrowing. The last index can also be an enum type.
         *
         * TODO: Use concepts when we move to C++20.
         */
        template <typename... idx, class U=T,
            std::enable_if_t<!std::is_void_v<U> && !IsArray4_v
                && detail::ArrayNDIndexCheck_impl_v<N, last_dim_component, idx...>,int> = 0>
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
         * \tparam M The dimension of the index vector.
         * \param iv The index vector.
         * \return Reference to the element.
         *
         * This function handles access via:
         * 1. Full dimension vector (M == N).
         * 2. Spatial dimension vector (M == N-1) && !IsArray4_v
         *    (only enabled if last_dim_component is true).
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
         * \tparam M The dimension of the spatial index vector (typically N-1 or AMREX_SPACEDIM).
         * \param iv The spatial index vector.
         * \param n Component index (last dimension).
         * \return Reference to the element.
         *
         * This overload is enabled only when last_dim_component is true and M represents
         * the spatial dimensions (either N-1 or AMREX_SPACEDIM).
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
         * \tparam M The dimension of the index vector.
         * \param iv The index vector.
         * \return Pointer to the element.
         *
         * This function handles access via:
         * 1. Full dimension vector (M == N).
         * 2. Spatial dimension vector (M == N-1) && !IsArray4_v
         *    (only enabled if last_dim_component is true).
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
         * \tparam M The dimension of the spatial index vector (typically N-1 or AMREX_SPACEDIM).
         * \param iv The spatial index vector.
         * \param n Component index (last dimension).
         * \return Pointer to the element.
         *
         * This overload is enabled only when last_dim_component is true and M represents
         * the spatial dimensions (either N-1 or AMREX_SPACEDIM).
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
                return end[N-1];
            } else {
                return 1;
            }
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr std::size_t size () const noexcept {
            if (ok()) {
                std::size_t s = 1;
                constexpr_for<0, N>([&](int d) {
                    s *= (end[d] - begin[d]);
                });
                return s;
            } else {
                return 0;
            }
        }

        template <typename... idx, std::enable_if_t<!IsArray4_v &&
            detail::ArrayNDIndexCheck_impl_v<N, last_dim_component, idx...>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr bool contains (idx... i) const noexcept {
            constexpr auto nidx = sizeof...(i);
            return this->contains(IntVectND<nidx>{i...});
        }

        template <int M, std::enable_if_t<(M == N)
            || (!IsArray4_v && last_dim_component && (M + 1 == N)), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr bool contains (IntVectND<M> const& iv) const noexcept {
            bool inside = true;
            constexpr_for<0, M>([&](int d) {
                inside = inside && (iv[d] >= begin[d]) && (iv[d] < end[d]);
            });
            return inside;
        }

        template <int M, std::enable_if_t<last_dim_component && !IsArray4_v
            && (M + 1 == N), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr bool contains (IntVectND<M> const& iv, int n) const noexcept {
            bool inside = true;
            constexpr_for<0, M>([&](int d) {
                inside = inside && (iv[d] >= begin[d]) && (iv[d] < end[d]);
            });
            inside = inside && (n >= 0) && (n < end[N-1]);
            return inside;
        }

        template <int M=N, bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CellData<T> cellData (int i, int j, int k) const noexcept {
            int ncomp = end[M-1];
            return CellData<T>(this->ptr(i,j,k), stride.a[M-2], ncomp);
        }

        template <int M, std::enable_if_t<(M == N)
            || (last_dim_component && (M + 1 == N || M == AMREX_SPACEDIM)), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Long get_offset (IntVectND<M> const& iv) const noexcept
        {
            Long offset = iv[0] - begin[0];
            // If M == N and we have a component at the end, we only loop up to N-1 for spatial.
            // Otherwise, we loop up to M.
            constexpr int idx = (last_dim_component && M == N) ? N - 1 : M;
            constexpr_for<1, idx>([&](int d) {
                offset += (iv[d] - begin[d]) * stride.a[d-1];
            });

            // Handle the component offset if the input is the full N dimensions
            // and the last dimension is a component.
            if constexpr (last_dim_component && M == N) {
                offset += iv[N-1] * stride.a[N-2];
            }
            return offset;
        }

        template <int M, std::enable_if_t<last_dim_component
            && ((M + 1 == N) || (N == 4 && M == AMREX_SPACEDIM)), int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        constexpr Long get_offset (IntVectND<M> const& iv, int n) const noexcept
        {
            Long offset = iv[0] - begin[0];
            constexpr_for<1, M>([&](int d) {
                offset += (iv[d] - begin[d]) * stride.a[d-1];
            });
            offset += n * stride.a[N-2];
            return offset;
        }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE inline
#endif
        void index_assert (IntVectND<N> const& iv) const
        {
            bool out_of_bounds = false;
            for (int d = 0; d < N; ++d) {
                if (iv[d] < begin[d] || iv[d] >= end[d]) {
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
                        ss << iv[d];
                        if (d + 1 < N) ss << ",";
                    }
                    ss << ") is out of bound (";
                    for (int d = 0; d < N; ++d) {
                        ss << begin[d] << ":" << end[d]-1;
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
        AMREX_GPU_HOST_DEVICE inline
#endif
        void index_assert (IntVectND<M> const& iv) const
        {
            IntVectND<N> iv_full = iv.template expand<N>(0);
            for (int d = M; d < N; ++d) {
                iv_full[d] = begin[d];
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
        AMREX_GPU_HOST_DEVICE inline
#endif
        void index_assert (IntVectND<M> const& iv, int n) const
        {
            IntVectND<N> iv_full = iv.template expand<N>(0);
            for (int d = M; d < N-1; ++d) {
                iv_full[d] = begin[d];
            }
            iv_full[N-1] = n;
            index_assert(iv_full);
        }
#endif

        //
        // Specialization for Array4
        //

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

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv.vect[0],0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv.vect[0],iv.vect[1],0);
#else
            return this->operator()(iv.vect[0],iv.vect[1],iv.vect[2]);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv.vect[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv.vect[0],iv.vect[1],0,n);
#else
            return this->operator()(iv.vect[0],iv.vect[1],iv.vect[2],n);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z);
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell, int n) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z,n);
        }

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

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv.vect[0],0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv.vect[0],iv.vect[1],0);
#else
            return this->ptr(iv.vect[0],iv.vect[1],iv.vect[2]);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv.vect[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv.vect[0],iv.vect[1],0,n);
#else
            return this->ptr(iv.vect[0],iv.vect[1],iv.vect[2],n);
#endif
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (Dim3 const& cell) const noexcept {
            return this->ptr(cell.x,cell.y,cell.z);
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U> && IsArray4_v, int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U* ptr (Dim3 const& cell, int n) const noexcept {
            return this->ptr(cell.x,cell.y,cell.z,n);
        }

        template <bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (int i, int j, int k) const noexcept {
            return (i>=begin.vect[0] && i<end.vect[0] &&
                    j>=begin.vect[1] && j<end.vect[1] &&
                    k>=begin.vect[2] && k<end.vect[2]);
        }

        template <bool B = IsArray4_v, std::enable_if_t<B,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (IntVect const& iv) const noexcept {
            return AMREX_D_TERM((iv.vect[0]>=begin.vect[0]) && (iv.vect[0]<end.vect[0]),
                             && (iv.vect[1]>=begin.vect[1]) && (iv.vect[1]<end.vect[1]),
                             && (iv.vect[2]>=begin.vect[2]) && (iv.vect[2]<end.vect[2]));
        }

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
        AMREX_GPU_HOST_DEVICE inline
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

    template<typename T>
    using Array4 = ArrayND<T, 4, true>;

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 lbound (Array4<T> const& a) noexcept
    {
        return Dim3{a.begin[0],a.begin[1],a.begin[2]};
    }

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 ubound (Array4<T> const& a) noexcept
    {
        return Dim3{a.end[0]-1,a.end[1]-1,a.end[2]-1};
    }

    template <class T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Dim3 length (Array4<T> const& a) noexcept
    {
        return Dim3{a.end[0]-a.begin[0],a.end[1]-a.begin[1],a.end[2]-a.begin[2]};
    }

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
