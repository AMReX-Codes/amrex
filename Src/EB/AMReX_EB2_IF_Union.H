#ifndef AMREX_EB2_IF_UNION_H_
#define AMREX_EB2_IF_UNION_H_

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <AMReX_IndexSequence.H>
#include <AMReX_Tuple.H>

#include <algorithm>
#include <utility>

namespace amrex { namespace EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

// Union of bodies

namespace UIF_detail {
    template <typename F>
    inline Real do_max (const RealArray& p, F&& f) noexcept
    {
        return f(p);
    }

    template <typename F, typename... Fs>
    inline Real do_max (const RealArray& p, F&& f, Fs&... fs) noexcept
    {
        return std::max(f(p), do_max(p, std::forward<Fs>(fs)...));
    }

    template <typename F>
    AMREX_GPU_HOST_DEVICE inline
    Real do_max (AMREX_D_DECL(Real x, Real y, Real z), F&& f) noexcept
    {
        return f(AMREX_D_DECL(x,y,z));
    }

    template <typename F, typename... Fs>
    AMREX_GPU_HOST_DEVICE inline
    Real do_max (AMREX_D_DECL(Real x, Real y, Real z), F&& f, Fs&... fs) noexcept
    {
        return amrex::max(f(AMREX_D_DECL(x,y,z)), do_max(AMREX_D_DECL(x,y,z), std::forward<Fs>(fs)...));
    }
}

template <class... Fs>
class UnionIF
    : public GpuTuple<Fs...>
{
public:
    using GpuTuple<Fs...>::GpuTuple;

    inline Real operator() (const RealArray& p) const noexcept
    {
        return op_impl(p, makeIndexSequence<sizeof...(Fs)>());
    }

    template <class U=UnionIF<Fs...>, class = typename std::enable_if<IsGPUable<U>::value>::type>
    AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return op_impl(AMREX_D_DECL(x,y,z), makeIndexSequence<sizeof...(Fs)>());
    }

protected:

    template <std::size_t... Is>
    inline Real op_impl (const RealArray& p, IndexSequence<Is...>) const noexcept
    {
        return UIF_detail::do_max(p, amrex::get<Is>(*this)...);
    }

    template <std::size_t... Is>
    AMREX_GPU_HOST_DEVICE inline
    Real op_impl (AMREX_D_DECL(Real x, Real y, Real z), IndexSequence<Is...>) const noexcept
    {
        return UIF_detail::do_max(AMREX_D_DECL(x,y,z), amrex::get<Is>(*this)...);
    }
};

template <class Head, class... Tail>
struct IsGPUable<UnionIF<Head, Tail...>, typename std::enable_if<IsGPUable<Head>::value>::type>
    : IsGPUable<UnionIF<Tail...> > {};

template <class F>
struct IsGPUable<UnionIF<F>, typename std::enable_if<IsGPUable<F>::value>::type>
    : std::true_type {};

template <class... Fs>
constexpr UnionIF<typename std::decay<Fs>::type ...>
makeUnion (Fs&&... fs)
{
    return UnionIF<typename std::decay<Fs>::type ...>(std::forward<Fs>(fs)...);
}

}}

#endif
