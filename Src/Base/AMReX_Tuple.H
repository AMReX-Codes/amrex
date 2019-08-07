
#ifndef AMREX_TUPLE_H_
#define AMREX_TUPLE_H_

#include <array>

#include <tuple>
#include <functional>
#include <type_traits>
#include <AMReX_IndexSequence.H>

#include <AMReX_GpuQualifiers.H>

namespace amrex {
    template <class... Ts>
    using Tuple = std::tuple<Ts...>;
}

namespace amrex {

namespace detail {

template <std::size_t I, typename T>
struct gpu_tuple_element
{
    constexpr gpu_tuple_element () = default;

    explicit constexpr gpu_tuple_element (T const& a_value)
        : m_value(a_value)
        {}

    explicit constexpr gpu_tuple_element (T && a_value)
        : m_value(std::move(a_value))
        {}

    T m_value;
};

template <std::size_t I, typename... Ts> struct gpu_tuple_impl;

template <std::size_t I, typename Head, typename... Tail>
struct gpu_tuple_impl<I, Head, Tail...>
    : public gpu_tuple_impl<I+1, Tail...>,
      public gpu_tuple_element<I, Head>
{
    constexpr gpu_tuple_impl () = default;

    constexpr gpu_tuple_impl (Head const& a_head, Tail const&... a_tail)
        : gpu_tuple_impl<I+1, Tail...>(a_tail...),
          gpu_tuple_element<I, Head>(a_head)
        {}

    constexpr gpu_tuple_impl (Head&& a_head, Tail &&... a_tail)
        : gpu_tuple_impl<I+1, Tail...>(std::forward<Tail>(a_tail)...),
          gpu_tuple_element<I, Head>(std::forward<Head>(a_head))
        {}
};

template <std::size_t I, typename Head>
struct gpu_tuple_impl<I, Head>
    : public gpu_tuple_element<I, Head>
{
    constexpr gpu_tuple_impl () = default;

    explicit constexpr gpu_tuple_impl (Head const& a_head)
        : gpu_tuple_element<I, Head>(a_head)
        {}

    explicit constexpr gpu_tuple_impl (Head&& a_head)
        : gpu_tuple_element<I, Head>(std::forward<Head>(a_head))
        {}
};

} // detail

// GpuTuple

template <typename... Ts>
class GpuTuple
    : public detail::gpu_tuple_impl<0, Ts...>
{
public:
    constexpr GpuTuple () = default;

    constexpr GpuTuple (Ts const&... args)
        : detail::gpu_tuple_impl<0, Ts...>(args...)
        {}

    constexpr GpuTuple (Ts&&... args)
        : detail::gpu_tuple_impl<0, Ts...>(std::forward<Ts>(args)...)
        {}
};

// GpuTupleSize

template <typename T> struct GpuTupleSize;

template <typename... Ts>
struct GpuTupleSize<GpuTuple<Ts...> >
    : public std::integral_constant<std::size_t, sizeof...(Ts)> {};

// GpuTupleElement

template <std::size_t I, typename T> struct GpuTupleElement;

template <std::size_t I, typename Head, typename... Tail>
struct GpuTupleElement<I, GpuTuple<Head, Tail...> >
    : GpuTupleElement<I-1, GpuTuple<Tail...> > {};

template <typename Head, typename... Tail>
struct GpuTupleElement<0, GpuTuple<Head, Tail...> > {
    using type = Head;
};

// get

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

} // detail

template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type&
get (GpuTuple<Ts...>& tup) noexcept
{
    return detail::get_impl<I,Ts...>(tup);
}

template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type const&
get (GpuTuple<Ts...> const& tup) noexcept
{
    return detail::get_impl<I,Ts...>(tup);
}

// makeTuple

namespace detail {
    template <typename T> struct unwrap { using type = T; };
    template <typename T> struct unwrap<std::reference_wrapper<T> > { using type = T&; };
    template <typename T>
    using tuple_decay_t = typename unwrap<typename std::decay<T>::type>::type;
}

template <typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
GpuTuple<Ts...>
makeTuple (Ts &&... args)
{
    return GpuTuple<detail::tuple_decay_t<Ts>...>(std::forward<Ts>(args)...);
}

namespace detail {
    template <typename...> struct tuple_cat_t {};

    template <typename... Ts>
    struct tuple_cat_t<GpuTuple<Ts...> >
    {
        typedef GpuTuple<Ts...> type;
    };

    template <typename... T1s, typename... T2s, typename... TPs>
    struct tuple_cat_t<GpuTuple<T1s...>,GpuTuple<T2s...>,TPs...>
    {
        typedef typename tuple_cat_t<GpuTuple<T1s..., T2s...>, TPs...>::type type;
    };

    template <typename R, typename TP1, typename TP2, std::size_t... N1, std::size_t... N2>
    AMREX_GPU_HOST_DEVICE constexpr R
    make_tuple (TP1 && a, TP2 && b,
                IndexSequence<N1...> const& n1, IndexSequence<N2...> const& n2)
    {
        return R(amrex::get<N1>(a)..., amrex::get<N2>(b)...);
    }
}

template <typename TP>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP && a) -> typename detail::tuple_cat_t<detail::tuple_decay_t<TP> >::type
{
    typedef typename detail::tuple_cat_t<detail::tuple_decay_t<TP> >::type ReturnType;
    return ReturnType(std::forward<TP>(a));
}

template <typename TP1, typename TP2>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP1 && a, TP2 && b) -> typename detail::tuple_cat_t<detail::tuple_decay_t<TP1>,
                                                              detail::tuple_decay_t<TP2> >::type
{
    typedef typename detail::tuple_cat_t<detail::tuple_decay_t<TP1>,
                                         detail::tuple_decay_t<TP2> >::type ReturnType;
    return detail::make_tuple<ReturnType>
        (a, b,
         makeIndexSequence<GpuTupleSize<typename std::decay<TP1>::type>::value>(),
         makeIndexSequence<GpuTupleSize<typename std::decay<TP2>::type>::value>());
}

template <typename TP1, typename TP2, typename... TPs>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP1&& a, TP2&& b, TPs&&... args)
    -> typename detail::tuple_cat_t<detail::tuple_decay_t<TP1>,
                                    detail::tuple_decay_t<TP2>,
                                    detail::tuple_decay_t<TPs>...>::type
{
    return TupleCat(TupleCat(std::forward<TP1>(a),std::forward<TP2>(b)),
                    std::forward<TPs>(args)...);
}

}

#endif /*AMREX_TUPLE_H_*/
