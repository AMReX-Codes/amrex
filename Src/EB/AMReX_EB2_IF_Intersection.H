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
    inline Real do_min (const RealArray& p, F&& f) noexcept
    {
        return f(p);
    }

    template <typename F, typename... Fs>
    inline Real do_min (const RealArray& p, F&& f, Fs&... fs) noexcept
    {
        return std::min(f(p), do_min(p, std::forward<Fs>(fs)...));
    }

    template <typename F>
    AMREX_GPU_HOST_DEVICE inline
    Real do_min (AMREX_D_DECL(Real x, Real y, Real z), F&& f) noexcept
    {
        return f(AMREX_D_DECL(x,y,z));
    }

    template <typename F, typename... Fs>
    AMREX_GPU_HOST_DEVICE inline
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

    inline Real operator() (const RealArray& p) const noexcept
    {
        return op_impl(p, makeIndexSequence<sizeof...(Fs)>());
    }

    template <class U=IntersectionIF<Fs...>, class = typename std::enable_if<IsGPUable<U>::value>::type>
    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return op_impl(AMREX_D_DECL(x,y,z), makeIndexSequence<sizeof...(Fs)>());
    }

protected:

    template <std::size_t... Is>
    inline Real op_impl (const RealArray& p, IndexSequence<Is...>) const noexcept
    {
        return IIF_detail::do_min(p, amrex::get<Is>(*this)...);
    }

    template <std::size_t... Is>
    AMREX_GPU_HOST_DEVICE inline
    Real op_impl (AMREX_D_DECL(Real x, Real y, Real z), IndexSequence<Is...>) const noexcept
    {
        return IIF_detail::do_min(AMREX_D_DECL(x,y,z), amrex::get<Is>(*this)...);
    }
};

template <class Head, class... Tail>
struct IsGPUable<IntersectionIF<Head, Tail...>, typename std::enable_if<IsGPUable<Head>::value>::type>
    : IsGPUable<IntersectionIF<Tail...> > {};

template <class F>
struct IsGPUable<IntersectionIF<F>, typename std::enable_if<IsGPUable<F>::value>::type>
    : std::true_type {};

template <class... Fs>
constexpr IntersectionIF<typename std::decay<Fs>::type ...>
makeIntersection (Fs&&... fs)
{
    return IntersectionIF<typename std::decay<Fs>::type ...>(std::forward<Fs>(fs)...);
}

}}

#endif
