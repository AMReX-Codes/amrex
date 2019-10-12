
#ifndef AMREX_TUPLE_H_
#define AMREX_TUPLE_H_

#include <array>

#include <tuple>
#include <functional>
#include <type_traits>
#include <utility>
#include <AMReX_TypeTraits.H>
#include <AMReX_IndexSequence.H>
#include <AMReX_Functional.H>
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
    template <typename U=T, amrex::EnableIf_t<std::is_default_constructible<U>::value,int> = 0>
    constexpr gpu_tuple_element () {}

    explicit constexpr gpu_tuple_element (T const& a_value)
        : m_value(a_value)
        {}

    template <typename U, amrex::EnableIf_t<std::is_convertible<U&&,T>::value,int> = 0>
    explicit constexpr gpu_tuple_element (U && a_value)
        : m_value(std::forward<U>(a_value))
        {}

    T m_value;
};

template <std::size_t I, typename... Ts> struct gpu_tuple_impl;

template <std::size_t I, typename Head, typename... Tail>
struct gpu_tuple_impl<I, Head, Tail...>
    : public gpu_tuple_impl<I+1, Tail...>,
      public gpu_tuple_element<I, Head>
{
    template<typename U=Head, amrex::EnableIf_t<std::is_default_constructible<U>::value,int> = 0>
    constexpr gpu_tuple_impl () {}

    constexpr gpu_tuple_impl (Head const& a_head, Tail const&... a_tail)
        : gpu_tuple_impl<I+1, Tail...>(a_tail...),
          gpu_tuple_element<I, Head>(a_head)
        {}

    template <typename UH, typename... UT, amrex::EnableIf_t<std::is_convertible<UH&&,Head>::value,int> = 0>
    constexpr gpu_tuple_impl (UH&& a_head, UT &&... a_tail)
        : gpu_tuple_impl<I+1, Tail...>(std::forward<UT>(a_tail)...),
          gpu_tuple_element<I, Head>(std::forward<UH>(a_head))
        {}
};

template <std::size_t I, typename Head>
struct gpu_tuple_impl<I, Head>
    : public gpu_tuple_element<I, Head>
{

    template<typename U=Head, amrex::EnableIf_t<std::is_default_constructible<U>::value,int> = 0>
    constexpr gpu_tuple_impl () {}

    explicit constexpr gpu_tuple_impl (Head const& a_head)
        : gpu_tuple_element<I, Head>(a_head)
        {}

    template <typename U, amrex::EnableIf_t<std::is_convertible<U&&,Head>::value,int> = 0>
    explicit constexpr gpu_tuple_impl (U&& a_head)
        : gpu_tuple_element<I, Head>(std::forward<U>(a_head))
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

    template <typename... Us, amrex::EnableIf_t<sizeof...(Us) == sizeof...(Ts),int> = 0>
    constexpr GpuTuple (Us&&... args)
        : detail::gpu_tuple_impl<0, Ts...>(std::forward<Us>(args)...)
        {}

    template <typename... Us, amrex::EnableIf_t<sizeof...(Us) == sizeof...(Ts),int> = 0>
    AMREX_GPU_HOST_DEVICE
    inline GpuTuple<Ts...>&
    operator= (GpuTuple<Us...> const& rhs);

    template <typename... Us, amrex::EnableIf_t<sizeof...(Us) == sizeof...(Ts),int> = 0>
    AMREX_GPU_HOST_DEVICE
    inline GpuTuple<Ts...>&
    operator= (GpuTuple<Us...> && rhs);
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

template <std::size_t I, typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr
typename GpuTupleElement<I, GpuTuple<Ts...> >::type &&
get (GpuTuple<Ts...> && tup) noexcept
{
    return detail::get_impl<I,Ts...>(std::move(tup));
}

namespace detail {
    template <std::size_t I, std::size_t N, typename TP1, typename TP2>
    AMREX_GPU_HOST_DEVICE
    amrex::EnableIf_t<(I<N-1),void>
    tuple_copy (TP1 & a, TP2 && b)
    {
        (amrex::get<I>(a) = amrex::get<I>(std::forward<TP2>(b)),
         tuple_copy<I+1,N>(a,std::forward<TP2>(b)));
    }

    template <std::size_t I, std::size_t N, typename TP1, typename TP2>
    AMREX_GPU_HOST_DEVICE
    amrex::EnableIf_t<I==N-1,void>
    tuple_copy (TP1 & a, TP2 && b)
    {
        amrex::get<I>(a) = amrex::get<I>(std::forward<TP2>(b));
    }
}

template <typename... Ts>
template <typename... Us, amrex::EnableIf_t<sizeof...(Us) == sizeof...(Ts),int> >
AMREX_GPU_HOST_DEVICE inline GpuTuple<Ts...>&
GpuTuple<Ts...>::operator= (GpuTuple<Us...> const& rhs)
{
    detail::tuple_copy<0,sizeof...(Ts)>(*this, rhs);
    return *this;
}

template <typename... Ts>
template <typename... Us, amrex::EnableIf_t<sizeof...(Us) == sizeof...(Ts),int> >
AMREX_GPU_HOST_DEVICE inline GpuTuple<Ts...>&
GpuTuple<Ts...>::operator= (GpuTuple<Us...> && rhs)
{
    detail::tuple_copy<0,sizeof...(Ts)>(*this, std::move(rhs));
    return *this;
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
GpuTuple<detail::tuple_decay_t<Ts>...>
makeTuple (Ts &&... args)
{
    return GpuTuple<detail::tuple_decay_t<Ts>...>(std::forward<Ts>(args)...);
}

namespace detail {
    template <typename...> struct tuple_cat_result {};

    template <typename... Ts>
    struct tuple_cat_result<GpuTuple<Ts...> >
    {
        typedef GpuTuple<Ts...> type;
    };

    template <typename... T1s, typename... T2s, typename... TPs>
    struct tuple_cat_result<GpuTuple<T1s...>,GpuTuple<T2s...>,TPs...>
    {
        typedef typename tuple_cat_result<GpuTuple<T1s..., T2s...>, TPs...>::type type;
    };

    template <typename R, typename TP1, typename TP2, std::size_t... N1, std::size_t... N2>
    AMREX_GPU_HOST_DEVICE constexpr R
    make_tuple (TP1 && a, TP2 && b,
                IndexSequence<N1...> const& n1, IndexSequence<N2...> const& n2)
    {
        return R(amrex::get<N1>(a)..., amrex::get<N2>(b)...);
    }
}

// TupleCat

template <typename TP>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP && a) -> typename detail::tuple_cat_result<detail::tuple_decay_t<TP> >::type
{
    typedef typename detail::tuple_cat_result<detail::tuple_decay_t<TP> >::type ReturnType;
    return ReturnType(std::forward<TP>(a));
}

template <typename TP1, typename TP2>
AMREX_GPU_HOST_DEVICE
constexpr auto
TupleCat (TP1 && a, TP2 && b) -> typename detail::tuple_cat_result<detail::tuple_decay_t<TP1>,
                                                                   detail::tuple_decay_t<TP2> >::type
{
    typedef typename detail::tuple_cat_result<detail::tuple_decay_t<TP1>,
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
    -> typename detail::tuple_cat_result<detail::tuple_decay_t<TP1>,
                                         detail::tuple_decay_t<TP2>,
                                         detail::tuple_decay_t<TPs>...>::type
{
    return TupleCat(TupleCat(std::forward<TP1>(a),std::forward<TP2>(b)),
                    std::forward<TPs>(args)...);
}

// Apply

namespace detail {

    template <typename F, typename... Args>
    AMREX_GPU_HOST_DEVICE
    auto INVOKE (F&& f, Args&&... args) -> decltype(f(std::forward<Args>(args)...));

    template <typename V, typename F, typename... Args> struct invoke_result {};

    template <typename F, typename... Args>
    struct invoke_result<decltype(void(INVOKE(std::declval<F>(), std::declval<Args>()...))),
                         F, Args...>
    {
        typedef decltype(INVOKE(std::declval<F>(), std::declval<Args>()...)) type;
    };

    template <typename F, typename...> struct apply_result {};

    template <typename F, typename... Ts>
    struct apply_result<F, GpuTuple<Ts...> >
    {
        typedef typename invoke_result<void, F, Ts...>::type type;
    };

    template <typename F, typename TP, std::size_t... N>
    AMREX_GPU_HOST_DEVICE
    constexpr auto
    apply_impl (F&& f, TP&& t, IndexSequence<N...> is)
        -> typename detail::apply_result<F,detail::tuple_decay_t<TP> >::type
    {
        return f(amrex::get<N>(std::forward<TP>(t))...);
    }
}

template <typename F, typename TP>
AMREX_GPU_HOST_DEVICE
constexpr auto
Apply (F&& f, TP&& t) -> typename detail::apply_result<F,detail::tuple_decay_t<TP> >::type
{
    return detail::apply_impl(std::forward<F>(f), std::forward<TP>(t),
                              makeIndexSequence<GpuTupleSize<typename std::decay<TP>::type>::value>());
}

// Tie

template <typename... Args>
AMREX_GPU_HOST_DEVICE
constexpr GpuTuple<Args&...>
Tie (Args&... args) noexcept
{
    return GpuTuple<Args&...>(args...);
}

// ForwardAsTuple

template <typename... Ts>
AMREX_GPU_HOST_DEVICE
constexpr GpuTuple<Ts&&...>
ForwardAsTuple (Ts&&... args) noexcept
{
    return GpuTuple<Ts&&...>(std::forward<Ts>(args)...);
}

}

#endif /*AMREX_TUPLE_H_*/
