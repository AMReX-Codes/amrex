#ifndef AMREX_EB2_IF_UNION_H_
#define AMREX_EB2_IF_UNION_H_

#include <AMReX_Array.H>
#include <AMReX_IndexSequence.H>

#include <type_traits>
#include <algorithm>
#include <tuple>
#include <utility>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

// Union of bodies

namespace UIF_detail {
    template <typename F>
    Real do_max (const RealArray& p, F&& f)
    {
        return f(p);
    }

    template <typename F, typename... Fs>
    Real do_max (const RealArray& p, F&& f, Fs&... fs)
    {
        return std::max(f(p), do_max(p, std::forward<Fs>(fs)...));
    }
}

template <class... Fs>
class UnionIF
    : public std::tuple<Fs...>
{
public:
    using std::tuple<Fs...>::tuple;

    Real operator() (const RealArray& p) const
    {
        constexpr std::size_t n = std::tuple_size<std::tuple<Fs...> >::value;
        return op_impl(p, makeIndexSequence<n>());
    }

protected:

    template <std::size_t... Is>
    Real op_impl (const RealArray& p, IndexSequence<Is...>) const
    {
        return UIF_detail::do_max(p, std::get<Is>(*this)...);
    }
};

template <class... Fs>
constexpr UnionIF<typename std::decay<Fs>::type ...>
makeUnion (Fs&&... fs)
{
    return UnionIF<typename std::decay<Fs>::type ...>(std::forward<Fs>(fs)...);
}

}}

#endif
