#ifndef AMREX_TYPELIST_H_
#define AMREX_TYPELIST_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_TypeList.H
 * \brief Template metaprogramming utilities: TypeList, ForEach, CartesianProduct, etc.
 */

#include <type_traits>
#include <utility>

namespace amrex {

/**
 * \class amrex::TypeList
 * \ingroup amrex_utilities
 * \brief A compile-time list of types used for template metaprogramming.
 *
 * \tparam Ts The types held in the list.
 */
template <class... Ts>
struct TypeList
{
    //! Number of types in the TypeList.
    static constexpr std::size_t size () noexcept { return sizeof...(Ts); }
};

/// \cond DOXYGEN_IGNORE
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
/// \endcond

/// \ingroup amrex_utilities
//! Type at position \p I of a TypeList.
template <std::size_t I, typename T>
using TypeAt = typename detail::TypeListGet<I,T>::type;

/// \cond DOXYGEN_IGNORE
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
/// \endcond

/**
 * \ingroup amrex_utilities
 * \brief For each type t in TypeList, call f(t)
 *
 * For example, instead of
 * \code{.cpp}
 *     int order = ...;
 *     if (order == 1) {
 *         interp<1>(...);
 *     } else if (order == 2) {
 *         interp<2>(...);
 *     } else if (order == 4) {
 *         interp<4>(...);
 *     }
 * \endcode
 * we could have
 * \code{.cpp}
 *     int order = ...;
 *     ForEach(TypeList<std::integral_constant<int,1>,
 *                      std::integral_constant<int,2>,
 *                      std::integral_constant<int,4>>{},
 *             [&] (auto order_const) {
 *                 if (order_const() == order) {
 *                     interp<order_const()>(...);
 *                 }
 *             });
 * \endcode
 */
template <typename... Ts, typename F>
constexpr void
ForEach (TypeList<Ts...>, F&& f)
{
    detail::for_each_impl<TypeList<Ts...>>
        (std::forward<F>(f), std::make_index_sequence<sizeof...(Ts)>());
}

/**
 * \ingroup amrex_utilities
 * \brief For each type t in TypeList, call f(t) until true is returned.
 *
 * This behaves like return (f(t0) || f(t1) || f(t2) || ...).  Note that
 * short-circuiting occurs for the || operators.
 *
 * An example,
 * \code{.cpp}
 *     void AnyF (Any& dst, Any const& src) {
 *         // dst and src are either MultiFab or fMultiFab
 *         auto tt = CartesianProduct(TypeList<MultiFab,fMultiFab>{},
 *                                    TypeList<MultiFab,fMultiFab>{});
 *         bool r = ForEachUntil(tt, [&] (auto t) -> bool
 *         {
 *             using MF0 = TypeAt<0,decltype(t)>;
 *             using MF1 = TypeAt<1,decltype(t)>;
 *             if (dst.is<MF0>() && src.is<MF1>()) {
 *                 MF0      & dmf = dst.get<MF0>();
 *                 MF1 const& smf = src.get<MF1>();
 *                 f(dmf, smf);
 *                 return true;
 *             } else {
 *                 return false;
 *             }
 *         });
 *         if (!r) { amrex::Abort("Unsupported types"); }
 *     }
 * \endcode
 */
template <typename... Ts, typename F>
constexpr bool
ForEachUntil (TypeList<Ts...>, F&& f)
{
    return detail::for_each_until_impl<TypeList<Ts...>>
        (std::forward<F>(f), std::make_index_sequence<sizeof...(Ts)>());
}

/// \ingroup amrex_utilities
//! Concatenate two TypeLists.
template <typename... As, typename... Bs>
constexpr auto operator+ (TypeList<As...>, TypeList<Bs...>) {
    return TypeList<As..., Bs...>{};
}

/// \cond DOXYGEN_IGNORE
template <typename... Ls, typename A>
constexpr auto single_product (TypeList<Ls...>, A) {
    return TypeList<decltype(Ls{} + TypeList<A>{})...>{};
}

template <typename LLs, typename... As>
constexpr auto operator* (LLs, TypeList<As...>) {
    return (TypeList<>{} + ... + single_product(LLs{}, As{}));
}
/// \endcond

/**
 * \ingroup amrex_utilities
 * \brief Cartesian Product of TypeLists.
 *
 * For example,
 * \code{.cpp}
 *     CartesianProduct(TypeList<std::integral_constant<int,0>,
 *                               std::integral_constant<int,1>>{},
 *                      TypeList<std::integral_constant<int,2>,
 *                               std::integral_constant<int,3>>{});
 * \endcode
 * returns TypeList of TypeList of integral_constants {{0,2},{1,2},{0,3},{1,3}}.
 */
template <typename... Ls>
constexpr auto CartesianProduct (Ls...) {
    return (TypeList<TypeList<>>{} * ... * Ls{});
}

/**
 * \class amrex::TypeArray
 * \ingroup amrex_utilities
 * \brief A compile-time descriptor for N copies of type T, used with TypeMultiplier.
 *
 * \tparam T The element type.
 * \tparam N The number of repetitions.
 */
template <class T, std::size_t N>
struct TypeArray {
    using Type = T; ///< The element type.
    //! Return the number of repetitions.
    static constexpr std::size_t size () noexcept { return N; }
};

/// \cond DOXYGEN_IGNORE
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

    // overload of SingleTypeMultiplier for multiple types:
    // convert TypeArray<T,N> to  T, T, T, T, ... (N times with N >= 0)
    // Note that only this overload works with N = 0
    template <class T, std::size_t N>
    constexpr auto SingleTypeMultiplier (TypeArray<T, N>) {
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
/// \endcond

/**
 * \ingroup amrex_utilities
 * \brief Return the first template argument with the later arguments applied to it.
 * Types of the form T[N] are expanded to T, T, T, T, ... (N times with N >= 1).
 * Types of the form TypeArray<T,N> are expanded to T, T, T, T, ... (N times with N >= 0).
 *
 * For example, TypeMultiplier<ReduceData, Real[4], int[2], Long>
 * is an alias to the type ReduceData<Real, Real, Real, Real, int, int, Long>.
 */
template <template <class...> class TParam, class... Types>
using TypeMultiplier = TypeAt<0, decltype(detail::TApply<TParam>(
    (TypeList<>{} + ... + detail::SingleTypeMultiplier(std::declval<Types>()))
))>;

/**
 * \class amrex::ToTypeList
 * \ingroup amrex_utilities
 * \brief Wrap a non-TypeList type into a TypeList; leave a TypeList unchanged.
 *
 * \tparam T The type to wrap.
 */
template <typename T>
struct ToTypeList { using type = TypeList<T>; };

template <typename... Ts>
struct ToTypeList<TypeList<Ts...>> { using type = TypeList<Ts...>; };

/// \ingroup amrex_utilities
//! Alias for the wrapped TypeList type (see amrex::ToTypeList).
template <typename T>
using ToTypeList_t = typename ToTypeList<T>::type;

/**
 * \class amrex::IsTypeList
 * \ingroup amrex_utilities
 * \brief Type trait that is true when T is a TypeList specialization.
 *
 * \tparam T The type to test.
 */
template <typename T> struct IsTypeList : std::false_type {};
template <typename... Ts> struct IsTypeList<TypeList<Ts...>> : std::true_type {};
/// \ingroup amrex_utilities
//! True when T is a TypeList specialization.
template <typename T> inline constexpr bool IsTypeList_v = IsTypeList<T>::value;

}

#endif
