
#ifndef AMREX_TUPLE_H_
#define AMREX_TUPLE_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_Tuple.H
 * \brief GPU-compatible tuple type and associated utilities.
 */

#include <AMReX_Array.H>
#include <AMReX_Functional.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_TypeList.H>
#include <AMReX_TypeTraits.H>

#include <array>
#include <functional>
#include <tuple>

namespace amrex {
    /// \ingroup amrex_linear_containers
    //! Alias for std::tuple.
    template <class... Ts>
    using Tuple = std::tuple<Ts...>;
}

namespace amrex {

/// \cond DOXYGEN_IGNORE
namespace detail {

template <std::size_t I, typename T>
struct gpu_tuple_element
{
    template <typename U=T, std::enable_if_t<std::is_default_constructible_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr gpu_tuple_element () {} // NOLINT

    explicit constexpr gpu_tuple_element (T const& a_value)
        : m_value(a_value)
        {}

    template <typename U, std::enable_if_t<std::is_convertible_v<U&&,T>,int> = 0>
    explicit constexpr gpu_tuple_element (U && a_value) // NOLINT(bugprone-forwarding-reference-overload)
        : m_value(std::forward<U>(a_value))
        {}

    T m_value{};
};

template <std::size_t I, typename... Ts> struct gpu_tuple_impl;

template <std::size_t I, typename Head, typename... Tail>
struct gpu_tuple_impl<I, Head, Tail...>
    : public gpu_tuple_impl<I+1, Tail...>,
      public gpu_tuple_element<I, Head>
{
    template<typename U=Head, std::enable_if_t<std::is_default_constructible_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr gpu_tuple_impl () {} // NOLINT

    constexpr gpu_tuple_impl (Head const& a_head, Tail const&... a_tail)
        : gpu_tuple_impl<I+1, Tail...>(a_tail...),
          gpu_tuple_element<I, Head>(a_head)
        {}

    template <typename UH, typename... UT, std::enable_if_t<std::is_convertible_v<UH&&,Head>,int> = 0>
    constexpr gpu_tuple_impl (UH&& a_head, UT &&... a_tail)
        : gpu_tuple_impl<I+1, Tail...>(std::forward<UT>(a_tail)...),
          gpu_tuple_element<I, Head>(std::forward<UH>(a_head))
        {}
};

template <std::size_t I, typename Head>
struct gpu_tuple_impl<I, Head>
    : public gpu_tuple_element<I, Head>
{

    template<typename U=Head, std::enable_if_t<std::is_default_constructible_v<U>,int> = 0>
    AMREX_GPU_HOST_DEVICE
    constexpr gpu_tuple_impl () {} // NOLINT

    explicit constexpr gpu_tuple_impl (Head const& a_head)
        : gpu_tuple_element<I, Head>(a_head)
        {}

    template <typename U, std::enable_if_t<std::is_convertible_v<U&&,Head>,int> = 0>
    explicit constexpr gpu_tuple_impl (U&& a_head) // NOLINT(bugprone-forwarding-reference-overload)
        : gpu_tuple_element<I, Head>(std::forward<U>(a_head))
        {}
};

} // detail
/// \endcond

/**
 * \class amrex::GpuTuple
 * \ingroup amrex_linear_containers
 * \brief GPU-compatible tuple
 */
template <typename... Ts>
class GpuTuple
    : public detail::gpu_tuple_impl<0, Ts...>
{
public:
    //! Default-construct all elements.
    AMREX_GPU_HOST_DEVICE // Some versions of nvcc require this in debug build
    constexpr GpuTuple () = default;

    //! Construct from a pack of const references.
    constexpr GpuTuple (Ts const&... args)
        : detail::gpu_tuple_impl<0, Ts...>(args...)
        {}

    //! Construct from a pack of forwarded arguments.
    template <typename... Us, std::enable_if_t<sizeof...(Us) == sizeof...(Ts),int> = 0>
    constexpr GpuTuple (Us&&... args)
        : detail::gpu_tuple_impl<0, Ts...>(std::forward<Us>(args)...)
        {}

    //! Copy-assign from a GpuTuple of compatible types.
    template <typename... Us, std::enable_if_t<sizeof...(Us) == sizeof...(Ts),int> = 0>
    AMREX_GPU_HOST_DEVICE
    inline GpuTuple<Ts...>&
    operator= (GpuTuple<Us...> const& rhs);

    //! Move-assign from a GpuTuple of compatible types.
    template <typename... Us, std::enable_if_t<sizeof...(Us) == sizeof...(Ts),int> = 0>
    AMREX_GPU_HOST_DEVICE
    inline GpuTuple<Ts...>&
    operator= (GpuTuple<Us...> && rhs);
};

/**
 * \class amrex::GpuTupleSize
 * \ingroup amrex_linear_containers
 * \brief Metafunction returning the number of elements in a GpuTuple
 *        (analogous to std::tuple_size).
 */
template <typename T> struct GpuTupleSize;

/// \cond DOXYGEN_IGNORE
template <typename... Ts>
struct GpuTupleSize<GpuTuple<Ts...> >
    : public std::integral_constant<std::size_t, sizeof...(Ts)> {};
/// \endcond

/**
 * \class amrex::GpuTupleElement
 * \ingroup amrex_linear_containers
 * \brief Metafunction returning the type of the I-th element of a GpuTuple
 *        (analogous to std::tuple_element).
 */
template <std::size_t I, typename T> struct GpuTupleElement;

/// \cond DOXYGEN_IGNORE
template <std::size_t I, typename Head, typename... Tail>
struct GpuTupleElement<I, GpuTuple<Head, Tail...> >
    : GpuTupleElement<I-1, GpuTuple<Tail...> > {};

template <typename Head, typename... Tail>
struct GpuTupleElement<0, GpuTuple<Head, Tail...> > {
    using type = Head;
};
/// \endcond

// get

/// \cond DOXYGEN_IGNORE
namespace detail {

template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type&
get_impl (detail::gpu_tuple_element
          <I, typename GpuTupleElement<I, GpuTuple<Ts...> >::type>& te) noexcept
{
    return te.m_value;
}

template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type const&
get_impl (detail::gpu_tuple_element
          <I, typename GpuTupleElement<I, GpuTuple<Ts...> >::type> const& te) noexcept
{
    return te.m_value;
}

template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type &&
get_impl (detail::gpu_tuple_element
          <I, typename GpuTupleElement<I, GpuTuple<Ts...> >::type> && te) noexcept
{
    return std::move(te).m_value;
}

} // detail
/// \endcond

/// \ingroup amrex_linear_containers
//! Return a mutable reference to the I-th element of \p tup.
template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type&
get (GpuTuple<Ts...>& tup) noexcept
{
    return detail::get_impl<I,Ts...>(tup);
}

/// \ingroup amrex_linear_containers
//! Return a const reference to the I-th element of \p tup.
template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type const&
get (GpuTuple<Ts...> const& tup) noexcept
{
    return detail::get_impl<I,Ts...>(tup);
}

/// \ingroup amrex_linear_containers
//! Return an rvalue reference to the I-th element of \p tup.
template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type &&
get (GpuTuple<Ts...> && tup) noexcept
{
    return detail::get_impl<I,Ts...>(std::move(tup));
}

/// \cond DOXYGEN_IGNORE
namespace detail {
    template <std::size_t I, std::size_t N, typename TP1, typename TP2>
    AMREX_GPU_HOST_DEVICE
    std::enable_if_t<(I<N-1),void>
    tuple_copy (TP1 & a, TP2 && b)
    {
        (amrex::get<I>(a) = amrex::get<I>(std::forward<TP2>(b)),
         tuple_copy<I+1,N>(a,std::forward<TP2>(b)));
    }

    template <std::size_t I, std::size_t N, typename TP1, typename TP2>
    AMREX_GPU_HOST_DEVICE
    std::enable_if_t<I==N-1,void>
    tuple_copy (TP1 & a, TP2 && b)
    {
        amrex::get<I>(a) = amrex::get<I>(std::forward<TP2>(b));
    }
}
/// \endcond

template <typename... Ts>
template <typename... Us, std::enable_if_t<sizeof...(Us) == sizeof...(Ts),int> >
AMREX_GPU_HOST_DEVICE inline GpuTuple<Ts...>&
GpuTuple<Ts...>::operator= (GpuTuple<Us...> const& rhs)
{
    detail::tuple_copy<0,sizeof...(Ts)>(*this, rhs);
    return *this;
}

template <typename... Ts>
template <typename... Us, std::enable_if_t<sizeof...(Us) == sizeof...(Ts),int> >
AMREX_GPU_HOST_DEVICE inline GpuTuple<Ts...>&
GpuTuple<Ts...>::operator= (GpuTuple<Us...> && rhs)
{
    detail::tuple_copy<0,sizeof...(Ts)>(*this, std::move(rhs));
    return *this;
}

// makeTuple

/// \cond DOXYGEN_IGNORE
namespace detail {
    template <typename T> struct unwrap { using type = T; };
    template <typename T> struct unwrap<std::reference_wrapper<T> > { using type = T&; };
    template <typename T>
    using tuple_decay_t = typename unwrap<std::decay_t<T>>::type;
}
/// \endcond

/// \ingroup amrex_linear_containers
//! Construct a GpuTuple deducing element types from the arguments.
template <typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
GpuTuple<detail::tuple_decay_t<Ts>...>
makeTuple (Ts &&... args)
{
    return GpuTuple<detail::tuple_decay_t<Ts>...>(std::forward<Ts>(args)...);
}

/// \cond DOXYGEN_IGNORE
namespace detail {
    template <typename...> struct tuple_cat_result {};

    template <typename... Ts>
    struct tuple_cat_result<GpuTuple<Ts...> >
    {
        using type = GpuTuple<Ts...>;
    };

    template <typename... T1s, typename... T2s, typename... TPs>
    struct tuple_cat_result<GpuTuple<T1s...>,GpuTuple<T2s...>,TPs...>
    {
        using type = typename tuple_cat_result<GpuTuple<T1s..., T2s...>, TPs...>::type;
    };

    template <typename R, typename TP1, typename TP2, std::size_t... N1, std::size_t... N2>
    AMREX_GPU_HOST_DEVICE constexpr R
    make_tuple (TP1 const& a, TP2 const& b,
                std::index_sequence<N1...> const& /*n1*/, std::index_sequence<N2...> const& /*n2*/)
    {
        return R(amrex::get<N1>(a)..., amrex::get<N2>(b)...);
    }
}
/// \endcond

// TupleCat

/// \ingroup amrex_linear_containers
//! Return the input GpuTuple unchanged.
template <typename TP>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP && a) -> typename detail::tuple_cat_result<detail::tuple_decay_t<TP> >::type
{
    using ReturnType = typename detail::tuple_cat_result<detail::tuple_decay_t<TP> >::type;
    return ReturnType(std::forward<TP>(a));
}

/// \ingroup amrex_linear_containers
//! Concatenate two GpuTuples into a single GpuTuple.
template <typename TP1, typename TP2>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP1 && a, TP2 && b) -> typename detail::tuple_cat_result<detail::tuple_decay_t<TP1>,
                                                                   detail::tuple_decay_t<TP2> >::type
{
    using ReturnType =  typename detail::tuple_cat_result<detail::tuple_decay_t<TP1>,
                                                          detail::tuple_decay_t<TP2> >::type;
    return detail::make_tuple<ReturnType>
        (std::forward<TP1>(a), std::forward<TP2>(b),
         std::make_index_sequence<GpuTupleSize<std::decay_t<TP1>>::value>(),
         std::make_index_sequence<GpuTupleSize<std::decay_t<TP2>>::value>());
}

/// \ingroup amrex_linear_containers
//! Concatenate three or more GpuTuples into a single GpuTuple.
template <typename TP1, typename TP2, typename... TPs>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP1&& a, TP2&& b, TPs&&... args)
    -> typename detail::tuple_cat_result<detail::tuple_decay_t<TP1>,
                                         detail::tuple_decay_t<TP2>,
                                         detail::tuple_decay_t<TPs>...>::type
{
    return TupleCat(TupleCat(std::forward<TP1>(a),std::forward<TP2>(b)),
                    std::forward<TPs>(args)...);
}

// TupleSplit

/// \cond DOXYGEN_IGNORE
namespace detail {

    template<std::size_t...Is>
    struct SplitIndexList {
        template<std::size_t J>
        AMREX_GPU_HOST_DEVICE
        static constexpr std::size_t get () noexcept {
            std::size_t arr[sizeof...(Is)] = {Is...};
            return arr[J];
        }

        template<std::size_t J>
        AMREX_GPU_HOST_DEVICE
        static constexpr std::size_t get_exclusive_sum () noexcept {
            std::size_t arr[sizeof...(Is)] = {Is...};
            std::size_t sum = 0;
            for (std::size_t k=0; k<J; ++k) {
                sum += arr[k];
            }
            return sum;
        }
    };

    template <std::size_t start, typename... Args, std::size_t... Is>
    AMREX_GPU_HOST_DEVICE
    constexpr auto
    GetSubTuple (const GpuTuple<Args...>& tup, std::index_sequence<Is...>) noexcept
    {
        return makeTuple(amrex::get<start+Is>(tup)...);
    }

    template <typename... Args, std::size_t... Is, typename SIL>
    AMREX_GPU_HOST_DEVICE
    constexpr auto
    TupleSplitImp (const GpuTuple<Args...>& tup, std::index_sequence<Is...>, SIL) noexcept
    {
        return makeTuple(
            GetSubTuple<(SIL::template get_exclusive_sum<Is>())>(
                tup,
                std::make_index_sequence<SIL::template get<Is>()>()
            )...
        );
    }
}
/// \endcond

/**
 * \ingroup amrex_linear_containers
 * \brief Returns a GpuTuple of GpuTuples obtained by splitting the input GpuTuple
 * according to the sizes specified by the template arguments.
 */
template <std::size_t... Is, typename... Args>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleSplit (const GpuTuple<Args...>& tup) noexcept
{
    static_assert((0 + ... + Is) == sizeof...(Args), "Incorrect total size in TupleSplit");
    return detail::TupleSplitImp(
        tup,
        std::make_index_sequence<sizeof...(Is)>(),
        detail::SplitIndexList<Is...>()
    );
}

// Apply

/// \cond DOXYGEN_IGNORE
namespace detail {

    template <typename F, typename... Args>
    AMREX_GPU_HOST_DEVICE
    auto INVOKE (F&& f, Args&&... args) -> decltype(f(std::forward<Args>(args)...));

    template <typename V, typename F, typename... Args> struct invoke_result {};

    template <typename F, typename... Args>
    struct invoke_result<decltype(void(INVOKE(std::declval<F>(), std::declval<Args>()...))),
                         F, Args...>
    {
        using type = decltype(INVOKE(std::declval<F>(), std::declval<Args>()...));
    };

    template <typename F, typename...> struct apply_result {};

    template <typename F, typename... Ts>
    struct apply_result<F, GpuTuple<Ts...> >
    {
        using type = typename invoke_result<void, F, Ts...>::type;
    };

    template <typename F, typename TP, std::size_t... N>
    AMREX_GPU_HOST_DEVICE
    constexpr auto
    apply_impl (F&& f, TP&& t, std::index_sequence<N...> /*is*/)
        -> typename detail::apply_result<F,detail::tuple_decay_t<TP> >::type
    {
        return std::forward<F>(f)(amrex::get<N>(std::forward<TP>(t))...);
    }
}
/// \endcond

/// \ingroup amrex_linear_containers
//! Invoke callable \p f with the elements of GpuTuple \p t as arguments.
template <typename F, typename TP>
AMREX_GPU_HOST_DEVICE
constexpr auto
Apply (F&& f, TP&& t) -> typename detail::apply_result<F,detail::tuple_decay_t<TP> >::type
{
    return detail::apply_impl(std::forward<F>(f), std::forward<TP>(t),
                              std::make_index_sequence<GpuTupleSize<std::decay_t<TP>>::value>());
}

// Tie

/// \ingroup amrex_linear_containers
//! Return a GpuTuple of lvalue references to the given arguments.
template <typename... Args>
AMREX_GPU_HOST_DEVICE
constexpr GpuTuple<Args&...>
Tie (Args&... args) noexcept
{
    return GpuTuple<Args&...>(args...);
}

// ForwardAsTuple

/// \ingroup amrex_linear_containers
//! Return a GpuTuple of forwarding references to the given arguments.
template <typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr GpuTuple<Ts&&...>
ForwardAsTuple (Ts&&... args) noexcept
{
    return GpuTuple<Ts&&...>(std::forward<Ts>(args)...);
}

// MakeZeroTuple

/**
 * \ingroup amrex_linear_containers
 * \brief Return a GpuTuple containing all zeros.
 * Note that a default-constructed GpuTuple can have uninitialized values.
 */
template <typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr GpuTuple<Ts...>
MakeZeroTuple (GpuTuple<Ts...>) noexcept
{
    return GpuTuple<Ts...>(static_cast<Ts>(0)...);
}

/// \cond DOXYGEN_IGNORE
namespace detail {
    template <typename T, std::size_t... I>
    AMREX_GPU_HOST_DEVICE constexpr
    auto tuple_to_array_helper (T const& tup, std::index_sequence<I...>) {
        return GpuArray<typename GpuTupleElement<0,T>::type, sizeof...(I)>{amrex::get<I>(tup)...};
    }
}
/// \endcond

/// \ingroup amrex_linear_containers
//! Convert a single-element GpuTuple<T> to a GpuArray<T,1>.
template <typename T>
AMREX_GPU_HOST_DEVICE constexpr
auto tupleToArray (GpuTuple<T> const& tup)
{
    return GpuArray<T,1>{amrex::get<0>(tup)};
}

/// \ingroup amrex_linear_containers
//! Convert GpuTuple<T,T2,Ts...> to GpuArray when all element types are the same.
template <typename T, typename T2, typename... Ts, std::enable_if_t<Same<T,T2,Ts...>::value, int> = 0>
AMREX_GPU_HOST_DEVICE constexpr
auto tupleToArray (GpuTuple<T,T2,Ts...> const& tup)
{
    return detail::tuple_to_array_helper(tup, std::index_sequence_for<T,T2,Ts...>{});
}

} // namespace amrex

// Specialize std::tuple_size for GpuTuple. Used by structured bindings.
template<typename... Ts>
struct std::tuple_size<amrex::GpuTuple<Ts...>> {
    static constexpr std::size_t value = sizeof...(Ts);
};

// Specialize std::tuple_element for GpuTuple. Used by structured bindings.
template<typename T, typename... Ts>
struct std::tuple_element<std::size_t{0}, amrex::GpuTuple<T, Ts...>> {
    using type = T;
};

template<std::size_t s, typename T, typename... Ts>
struct std::tuple_element<s, amrex::GpuTuple<T, Ts...>> {
    using type = typename std::tuple_element<s-1, amrex::GpuTuple<Ts...>>::type;
};

#endif /*AMREX_TUPLE_H_*/
