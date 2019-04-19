
#ifndef AMREX_TUPLE_H_
#define AMREX_TUPLE_H_

#include <array>

#include <tuple>
#include <functional>
#include <type_traits>

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

    explicit constexpr GpuTuple (Ts const&... args)
        : detail::gpu_tuple_impl<0, Ts...>(args...)
        {}

    explicit constexpr GpuTuple (Ts&&... args)
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
GpuTuple<Ts...>
makeTuple (Ts &&... args)
{
    return GpuTuple<detail::tuple_decay_t<Ts>...>(std::forward<Ts>(args)...);
}

}

#endif /*AMREX_TUPLE_H_*/
