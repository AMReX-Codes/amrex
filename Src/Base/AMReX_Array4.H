#ifndef AMREX_ARRAY4_H_
#define AMREX_ARRAY4_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_GpuPrint.H>
#include <AMReX_ConstexprFor.H>

#include <iostream>
#include <sstream>

namespace amrex {

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
        // Helper to get N-th argument from a parameter pack
        template <int idx, typename... Ts>
        inline constexpr auto get_impl (Ts&&... ts) {                    // NOLINT (cppcoreguidelines-missing-std-forward)
            return std::get<idx>(std::forward_as_tuple(ts...));
        }

        /**
         * /brief Indices checker for ArrayND
         *
         * 1. All indices are int
         * 2. Number of indices matches dimension or is one less than dimension
         */
        template <int dim, typename... idx>
        struct ArrayNDIndexCheck_impl {
            constexpr static bool value =
                (sizeof...(idx) == dim || sizeof...(idx) == dim - 1)
                && (Conjunction<std::is_same<int, std::decay_t<idx>>...>::value);
        };

        template <int dim, class... idx>
        inline constexpr bool ArrayNDIndexCheck_impl_v = ArrayNDIndexCheck_impl<dim, idx...>::value;

        template<std::size_t... idx>
        constexpr auto make_oob_message_impl (std::index_sequence<idx...>) {
            constexpr std::size_t Dim = sizeof...(idx);
            constexpr char prefix[] = " (";
            constexpr char middle[] = ") is out of bound (";
            constexpr char suffix[] = ")\n";

            constexpr std::size_t size =
                sizeof(prefix) - 1 +
                (2*Dim + (Dim-1)) +
                sizeof(middle) - 1 +
                (5*Dim + (Dim-1)) +
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
            ((buf[pos++] = '%', buf[pos++] = 'd', idx + 1 < Dim ? buf[pos++] = ',' : 0), ...);

            // ") is out of bound ("
            for (char c : middle) {
                if (c != '\0') {
                    buf[pos++] = c;
                }
            }

            // %d:%d,%d:%d,...
            ((buf[pos++] = '%', buf[pos++] = 'd', buf[pos++] = ':', buf[pos++] = '%', buf[pos++] = 'd',
                idx + 1 < Dim ? buf[pos++] = ',' : 0), ...);

            // suffix
            for (char c : suffix) {
                buf[pos++] = c;
            }
            return buf;
        }

        template <std::size_t... idx, std::enable_if_t<sizeof...(idx)>=1,int> = 0>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        void device_printf_impl (std::index_sequence<idx...>&& index_seq,
                                   const IntVectND<sizeof...(idx)>& iv,
                                   const IntVectND<sizeof...(idx)>& begin,
                                   const IntVectND<sizeof...(idx)>& end) {
            constexpr auto msg = make_oob_message_impl(std::move(index_seq));
            AMREX_DEVICE_PRINTF(msg.data(),
                iv[idx]...,
                ((begin[idx]), (end[idx]))...);
        }
    }

    template<typename T, int dim>
    struct ArrayND
    {
        static_assert(dim >= 1, "ArrayND must have at least one dimension");

        T* AMREX_RESTRICT p;
        GpuArray<Long, dim-1> nstride = {0}; // GpuArray safeguards against zero-size array
        IntVectND<dim> begin{1};
        IntVectND<dim> end{0}; // end is hi+1

        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND () noexcept : p(nullptr) {}

        template <class U=T, std::enable_if_t<std::is_const_v<U>,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (ArrayND<std::remove_const_t<T>, dim> const& rhs) noexcept
            : p(rhs.p), nstride(rhs.nstride), begin(rhs.begin), end(rhs.end) {}

        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (T* a_p, IntVectND<dim> const& a_begin, IntVectND<dim> const& a_end) noexcept
            : p(a_p), begin(a_begin), end(a_end)
        {
            if constexpr (dim > 1) {
                Long current_stride = 1;
                constexpr_for<0, dim - 1>([&] (int d) {
                    Long len = a_end[d] - a_begin[d];
                    current_stride *= len;
                    nstride[d] = current_stride;
                });
            }
        }

        // Dim3 overload for dim=4
        template <int Dim=dim, std::enable_if_t<Dim==4,int> = 0>
        constexpr ArrayND (T* a_p, Dim3 const& a_begin, Dim3 const& a_end, int a_ncomp) noexcept
            : p(a_p), begin(a_begin.x, a_begin.y, a_begin.z, 0), end(a_end.x, a_end.y, a_end.z, a_ncomp)
        {
            Long current_stride = 1;
            constexpr_for<0, dim-1>([&] (int d) {
                Long len = end[d] - begin[d];
                current_stride *= len;
                nstride[d] = current_stride;
            });
        }

        template <class U,
            std::enable_if_t
                  <std::is_same_v<std::remove_const_t<T>,
                                std::remove_const_t<U>> && dim == 4,int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (ArrayND<U, dim> const& rhs, int start_comp) noexcept
            : p((T*)(rhs.p + start_comp*rhs.nstride[dim-2])),
              nstride(rhs.nstride),
              begin(rhs.begin),
              end(rhs.end)
        {
            begin[3] = 0;
            end[3] = rhs.end[3] - start_comp;
        }

        template <class U,
            std::enable_if_t
                  <std::is_same_v<std::remove_const_t<T>,
                                std::remove_const_t<U>> && (dim == 4),int> = 0>
        AMREX_GPU_HOST_DEVICE
        constexpr ArrayND (ArrayND<U, dim> const& rhs, int start_comp, int num_comp) noexcept
            : p((T*)(rhs.p + start_comp*rhs.nstride[dim-2])),
              nstride(rhs.nstride),
              begin(rhs.begin),
              end(rhs.end)
        {
            begin[3] = 0;
            end[3] = num_comp;
        }

        AMREX_GPU_HOST_DEVICE
        explicit operator bool() const noexcept { return p != nullptr; }

        template <typename... idx, class U=T,
            std::enable_if_t<!std::is_void_v<U> && detail::ArrayNDIndexCheck_impl_v<dim, idx...>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (idx... i) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            if constexpr (sizeof...(i) == dim) {
                index_assert(i...);
            } else if constexpr (sizeof...(i) == dim - 1) {
                // When number of indices is one less than dimension, we assume the last index is zero
                index_assert(i...,0);
            }
#endif
            if constexpr (sizeof...(i) == dim) {
                return p[get_offset(i...)];
            } else {
                return p[get_offset(i...,0)];
            }
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVectND<dim> const& iv) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(iv);
#endif
            return p[get_offset_iv(iv)];
        }

        template <class U=T, int Dim=dim,
            std::enable_if_t<!std::is_void_v<U> && (Dim==4),int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv[0],0,0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv[0],iv[1],0,0);
#elif (AMREX_SPACEDIM == 3)
            return this->operator()(iv[0],iv[1],iv[2],0);
#endif
        }

        template <class U=T, int Dim=dim,
            std::enable_if_t<!std::is_void_v<U> && (Dim==4),int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->operator()(iv[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->operator()(iv[0],iv[1],0,n);
#elif (AMREX_SPACEDIM == 3)
            return this->operator()(iv[0],iv[1],iv[2],n);
#endif
        }

        // Explicit overload for dim=4
        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0,
                    std::enable_if_t<dim==4,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,0);
#endif
            return p[get_offset(i,j,k,0)];
        }

        template <class U=T, std::enable_if_t<!std::is_void_v<U>,int> = 0,
                    std::enable_if_t<dim==4,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (int i, int j, int k, int n) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            index_assert(i,j,k,n);
#endif
            return p[get_offset(i,j,k,n)];
        }

        template <typename... idx, class U=T,
            std::enable_if_t<!std::is_void_v<U> && detail::ArrayNDIndexCheck_impl_v<dim, idx...>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (idx... i) const noexcept {
#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
            if constexpr (sizeof...(i) == dim) {
                index_assert(i...);
            } else if constexpr (sizeof...(i) == dim - 1) {
                // When number of indices is one less than dimension, we assume the last index is zero
                index_assert(i...,0);
            }
#endif
            if constexpr (sizeof...(i) == dim) {
                return p + get_offset(i...);
            } else {
                return p + get_offset(i...,0);
            }
        }

        template <class U=T, int Dim=dim, std::enable_if_t<!std::is_void_v<U>,int> = 0,
                    std::enable_if_t<dim==4,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv[0],0,0,0);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv[0],iv[1],0,0);
#elif (AMREX_SPACEDIM == 3)
            return this->ptr(iv[0],iv[1],iv[2],0);
#endif
        }

        template <class U=T, int Dim=dim, std::enable_if_t<!std::is_void_v<U>,int> = 0,
                    std::enable_if_t<dim==4,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* ptr (IntVect const& iv, int n) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->ptr(iv[0],0,0,n);
#elif (AMREX_SPACEDIM == 2)
            return this->ptr(iv[0],iv[1],0,n);
#elif (AMREX_SPACEDIM == 3)
            return this->ptr(iv[0],iv[1],iv[2],n);
#endif
        }

        template <class U=T, int Dim=dim,
            std::enable_if_t<!std::is_void_v<U> && (Dim==4),int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z,0);
        }

        template <class U=T, int Dim=dim,
            std::enable_if_t<!std::is_void_v<U> && (Dim==4),int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        U& operator() (Dim3 const& cell, int n) const noexcept {
            return this->operator()(cell.x,cell.y,cell.z,n);
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        T* dataPtr () const noexcept {
            return this->p;
        }

        template <int Dim=dim, std::enable_if_t<Dim==4,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        int nComp () const noexcept {
            return end[3] - begin[3];
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        std::size_t size () const noexcept {
            std::size_t s = 1;
            constexpr_for<0, dim>([&] (auto d) {
                s *= (end[d] - begin[d]);
            });
            return s;
        }

        template <typename... idx, std::enable_if_t<detail::ArrayNDIndexCheck_impl_v<dim, idx...>,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (idx... i) const noexcept {
            if constexpr (sizeof...(i) == dim) {
                return this->contains(IntVectND<dim>{i...});
            } else {
                // When number of indices is one less than dimension, we assume the last index is zero
                return this->contains(IntVectND<dim>{i...,0});
            }
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (IntVect const& iv) const noexcept {
#if (AMREX_SPACEDIM == 1)
            return this->contains(IntVectND<dim>{iv[0],0,0,0});
#elif (AMREX_SPACEDIM == 2)
            return this->contains(IntVectND<dim>{iv[0],iv[1],0,0});
#elif (AMREX_SPACEDIM == 3)
            return this->contains(IntVectND<dim>{iv[0],iv[1],iv[2],0});
#endif
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (IntVectND<dim> const& iv) const noexcept {
            bool inside = true;
            constexpr_for<0, dim>([&] (auto d) {
                inside = inside && (iv[d] >= begin[d]) && (iv[d] < end[d]);
            });
            return inside;
        }

        template <int Dim=dim, std::enable_if_t<Dim==4,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        bool contains (Dim3 const& cell) const noexcept {
            return this->contains(IntVectND<dim>{cell.x, cell.y, cell.z, 0});
        }

        template <int Dim=dim, std::enable_if_t<Dim==4,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        CellData<T> cellData (int i, int j, int k) const noexcept {
            int ncomp = end[3] - begin[3];
            return CellData<T>(this->ptr(i,j,k), nstride[dim-2], ncomp);
        }


        template <typename... idx, std::enable_if_t<sizeof...(idx) == dim,int> = 0>
        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Long get_offset (idx... indices) const noexcept
        {
            Long offset = detail::get_impl<0>(indices...) - begin[0];
            constexpr_for<1, dim>([&] (auto i) {
                offset += (detail::get_impl<i>(indices...) - begin[i]) * nstride[i - 1];
            });
            return offset;
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        Long get_offset_iv (IntVectND<dim> const& iv) const noexcept
        {
            Long offset = iv[0] - begin[0];
            constexpr_for<1, dim>([&] (auto d) {
                offset += (iv[d] - begin[d]) * nstride[d-1];
            });
            return offset;
        }

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE inline
#endif
        void index_assert (IntVectND<dim> const& iv) const
        {
            bool out_of_bounds = false;
            constexpr_for<0, dim>([&] (auto d) {
                if (iv[d] < begin[d] || iv[d] >= end[d]) {
                    out_of_bounds = true;
                }
            });
            if (out_of_bounds) {
                AMREX_IF_ON_DEVICE((
                    device_printf_impl(std::make_index_sequence<dim>{}, iv, begin, end);
                    amrex::Abort();
                ))
                AMREX_IF_ON_HOST((
                    std::stringstream ss;
                    ss << " (";
                    constexpr_for<0, dim>([&] (auto d) {
                        ss << iv[d];
                        if (d + 1 < dim) ss << ",";
                    });
                    ss << ") is out of bound (";
                    constexpr_for<0, dim>([&] (auto d) {
                        ss << begin[d] << ":" << end[d]-1;
                        if (d + 1 < dim) ss << ",";
                    });
                    ss << ")";
                    amrex::Abort(ss.str());
                ))
            }
        }
#endif

#if defined(AMREX_DEBUG) || defined(AMREX_BOUND_CHECK)
        template<typename... idx, std::enable_if_t<sizeof...(idx) == dim,int> = 0>
#if defined(AMREX_USE_HIP)
        AMREX_GPU_HOST_DEVICE AMREX_NO_INLINE
#else
        AMREX_GPU_HOST_DEVICE inline
#endif
        void index_assert (idx... indices) const
        {
            this->index_assert(IntVectND<dim>{indices...});
        }
#endif
    };

    template<typename T>
    using Array4 = ArrayND<T,4>;

    template <class Tto, class Tfrom>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Array4<Tto> ToArray4 (Array4<Tfrom> const& a_in) noexcept
    {
        return Array4<Tto>((Tto*)(a_in.p), a_in.begin, a_in.end);
    }

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

    template <typename T>
    std::ostream& operator<< (std::ostream& os, const Array4<T>& a) {
        os << "((" << lbound(a) << ',' << ubound(a) << ")," << a.nComp() << ')';
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
