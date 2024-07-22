#ifndef AMREX_TYPELIST_H_
#define AMREX_TYPELIST_H_
#include <AMReX_Config.H>

#include <utility>

namespace amrex {

//! Struct for holding types
template <class... Ts>
struct TypeList
{
    //! Number of types in the TypeList
    static constexpr std::size_t size () noexcept { return sizeof...(Ts); }
};

namespace detail {
    template <std::size_t I, typename T> struct TypeListGet;

    template <std::size_t I, typename Head, typename... Tail>
    struct TypeListGet<I, TypeList<Head, Tail...> >
        : TypeListGet<I-1, TypeList<Tail...> > {};

    template <typename Head, typename... Tail>
    struct TypeListGet<0, TypeList<Head, Tail...> > {
        using type = Head;
    };
}

//! Type at position I of a TypeList
template <std::size_t I, typename T>
using TypeAt = typename detail::TypeListGet<I,T>::type;

namespace detail
{
    template <typename TL, typename F, std::size_t...N>
    constexpr void for_each_impl (F const&f, std::index_sequence<N...>)
    {
        (f(TypeAt<N,TL>{}), ...);
    }

    template <typename TL, typename F, std::size_t...N>
    constexpr bool for_each_until_impl (F const&f, std::index_sequence<N...>)
    {
        return (f(TypeAt<N,TL>{}) || ...);
    }
}

/**
 * \brief For each type t in TypeList, call f(t)
 *
 * For example, instead of
 \verbatim
     int order = ...;
     if (order == 1) {
         interp<1>(...);
     } else if (order == 2) {
         interp<2>(...);
     } else if (order == 4) {
         interp<4>(...);
     }
 \endverbatim
 * we could have
 \verbatim
     int order = ...;
     ForEach(TypeList<std::integral_constant<int,1>,
                      std::integral_constant<int,2>,
                      std::integral_constant<int,4>>{},
             [&] (auto order_const) {
                 if (order_const() == order) {
                     interp<order_const()>(...);
                 }
             });
 \endverbatim
 */
template <typename... Ts, typename F>
constexpr void
ForEach (TypeList<Ts...>, F&& f)
{
    detail::for_each_impl<TypeList<Ts...>>
        (std::forward<F>(f), std::make_index_sequence<sizeof...(Ts)>());
}

/**
 * \brief For each type t in TypeList, call f(t) until true is returned.
 *
 * This behaves like return (f(t0) || f(t1) || f(t2) || ...).  Note that
 * shor-circuting occurs for the || operators.
 *
 * An example,
 \verbatim
     void AnyF (Any& dst, Any const& src) {
         // dst and src are either MultiFab or fMultiFab
         auto tt = CartesianProduct(TypeList<MultiFab,fMultiFab>{},
                                    TypeList<MultiFab,fMultiFab>{});
         bool r = ForEachUntil(tt, [&] (auto t) -> bool
         {
             using MF0 = TypeAt<0,decltype(t)>;
             using MF1 = TypeAt<1,decltype(t)>;
             if (dst.is<MF0>() && src.is<MF1>()) {
                 MF0      & dmf = dst.get<MF0>();
                 MF1 const& smf = src.get<MF1>();
                 f(dmf, smf);
                 return true;
             } else {
                 return false;
             }
         });
         if (!r) { amrex::Abort("Unsupported types"); }
     }
 \endverbatim
 */
template <typename... Ts, typename F>
constexpr bool
ForEachUntil (TypeList<Ts...>, F&& f)
{
    return detail::for_each_until_impl<TypeList<Ts...>>
        (std::forward<F>(f), std::make_index_sequence<sizeof...(Ts)>());
}

//! Concatenate two TypeLists
template <typename... As, typename... Bs>
constexpr auto operator+ (TypeList<As...>, TypeList<Bs...>) {
    return TypeList<As..., Bs...>{};
}

template <typename... Ls, typename A>
constexpr auto single_product (TypeList<Ls...>, A) {
    return TypeList<decltype(Ls{} + TypeList<A>{})...>{};
}

template <typename LLs, typename... As>
constexpr auto operator* (LLs, TypeList<As...>) {
    return (TypeList<>{} + ... + single_product(LLs{}, As{}));
}

/**
 * \brief Cartesian Product of TypeLists.
 *
 * For example,
 \verbatim
     CartesianProduct(TypeList<std::integral_constant<int,0>,
                               std::integral_constant<int,1>>{},
                      TypeList<std::integral_constant<int,2>,
                               std::integral_constant<int,3>>{});
 \endverbatim
 * returns TypeList of TypeList of integral_constants {{0,2},{1,2},{0,3},{1,3}}.
 */
template <typename... Ls>
constexpr auto CartesianProduct (Ls...) {
    return (TypeList<TypeList<>>{} * ... * Ls{});
}

namespace detail {
    // return TypeList<T, T, T, T, ... (N times)> by using the fast power algorithm
    template <class T, std::size_t N>
    constexpr auto SingleTypeMultiplier_impl () {
        if constexpr (N == 0) {
            return TypeList<>{};
        } else if constexpr (N == 1) {
            return TypeList<T>{};
        } else if constexpr (N % 2 == 0) {
            return SingleTypeMultiplier_impl<T, N / 2>() + SingleTypeMultiplier_impl<T, N / 2>();
        } else {
            return SingleTypeMultiplier_impl<T, N - 1>() + TypeList<T>{};
        }
    }

    // overload of SingleTypeMultiplier for multiple types:
    // convert T[N] to  T, T, T, T, ... (N times with N >= 1)
    template <class T, std::size_t N>
    constexpr auto SingleTypeMultiplier (const T (&)[N]) {
        return SingleTypeMultiplier_impl<T, N>();
    }

    // overload of SingleTypeMultiplier for one regular type
    template <class T>
    constexpr auto SingleTypeMultiplier (T) {
        return TypeList<T>{};
    }

    // apply the types of the input TypeList as template arguments to TParam
    template <template <class...> class TParam, class... Args>
    constexpr auto TApply (TypeList<Args...>) {
        return TypeList<TParam<Args...>>{};
    }
}

/**
 * \brief Return the first template argument with the later arguments applied to it.
 * Types of the form T[N] are expanded to T, T, T, T, ... (N times with N >= 1).
 *
 * For example, TypeMultiplier<ReduceData, Real[4], int[2], Long>
 * is an alias to the type ReduceData<Real, Real, Real, Real, int, int, Long>.
 */
template <template <class...> class TParam, class... Types>
using TypeMultiplier = TypeAt<0, decltype(detail::TApply<TParam>(
    (TypeList<>{} + ... + detail::SingleTypeMultiplier(Types{}))
))>;

}

#endif
