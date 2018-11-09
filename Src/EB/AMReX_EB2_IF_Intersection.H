#ifndef AMREX_EB2_IF_INTERSECTION_H_
#define AMREX_EB2_IF_INTERSECTION_H_

#include <AMReX_Array.H>
#include <AMReX_IndexSequence.H>

#include <type_traits>
#include <algorithm>
#include <tuple>
#include <utility>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

// Intersection of bodies

namespace IIF_detail {
    template <typename F>
    Real do_min (const RealArray& p, F&& f)
    {
        return f(p);
    }

    template <typename F, typename... Fs>
    Real do_min (const RealArray& p, F&& f, Fs&... fs)
    {
        return std::min(f(p), do_min(p, std::forward<Fs>(fs)...));
    }
}

template <class... Fs>
class IntersectionIF
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
        return IIF_detail::do_min(p, std::get<Is>(*this)...);
    }
};

template <class... Fs>
constexpr IntersectionIF<typename std::decay<Fs>::type ...>
makeIntersection (Fs&&... fs)
{
    return IntersectionIF<typename std::decay<Fs>::type ...>(std::forward<Fs>(fs)...);
}

}}

#endif
