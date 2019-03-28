#ifndef AMREX_EB2_IF_INTERSECTION_H_
#define AMREX_EB2_IF_INTERSECTION_H_

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <AMReX_IndexSequence.H>
#include <AMReX_Tuple.H>

#include <algorithm>
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

    template <typename F>
    Real do_min (AMREX_D_DECL(Real x, Real y, Real z), F&& f)
    {
        return f(AMREX_D_DECL(x,y,z));
    }

    template <typename F, typename... Fs>
    Real do_min (AMREX_D_DECL(Real x, Real y, Real z), F&& f, Fs&... fs)
    {
        return amrex::min(f(AMREX_D_DECL(x,y,z)), do_min(AMREX_D_DECL(x,y,z), std::forward<Fs>(fs)...));
    }
}

template <class... Fs>
class IntersectionIF
    : public GpuTuple<Fs...>
{
public:
    using GpuTuple<Fs...>::GpuTuple;

    Real operator() (const RealArray& p) const
    {
        return op_impl(p, makeIndexSequence<sizeof...(Fs)>());
    }

protected:

    template <std::size_t... Is>
    Real op_impl (const RealArray& p, IndexSequence<Is...>) const
    {
        return IIF_detail::do_min(p, amrex::get<Is>(*this)...);
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
