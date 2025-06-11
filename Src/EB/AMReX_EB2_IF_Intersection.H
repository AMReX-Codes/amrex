#ifndef AMREX_EB2_IF_INTERSECTION_H_
#define AMREX_EB2_IF_INTERSECTION_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <AMReX_Tuple.H>

#include <algorithm>
#include <utility>

namespace amrex::EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

// Intersection of bodies

namespace IIF_detail {
    template <typename F>
    [[nodiscard]] inline Real do_min (const RealArray& p, F&& f) noexcept
    {
        return std::forward<F>(f)(p);
    }

    template <typename F, typename... Fs>
    [[nodiscard]] inline Real do_min (const RealArray& p, F&& f, Fs&... fs) noexcept
    {
        return amrex::min(std::forward<F>(f)(p), do_min(p, std::forward<Fs>(fs)...));
    }

    template <typename F>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real do_min (AMREX_D_DECL(Real x, Real y, Real z), F&& f) noexcept
    {
        return std::forward<F>(f)(AMREX_D_DECL(x,y,z));
    }

    template <typename F, typename... Fs>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real do_min (AMREX_D_DECL(Real x, Real y, Real z), F&& f, Fs&... fs)
    {
        return amrex::min(std::forward<F>(f)(AMREX_D_DECL(x,y,z)), do_min(AMREX_D_DECL(x,y,z), std::forward<Fs>(fs)...));
    }
}

template <class... Fs>
class IntersectionIF
    : public GpuTuple<Fs...>
{
public:
    using GpuTuple<Fs...>::GpuTuple;

    [[nodiscard]] inline Real operator() (const RealArray& p) const noexcept
    {
        return op_impl(p, std::make_index_sequence<sizeof...(Fs)>());
    }

    template <class U=IntersectionIF<Fs...>, std::enable_if_t<IsGPUable<U>::value,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return op_impl(AMREX_D_DECL(x,y,z), std::make_index_sequence<sizeof...(Fs)>());
    }

protected:

    template <std::size_t... Is>
    [[nodiscard]] inline Real op_impl (const RealArray& p, std::index_sequence<Is...>) const noexcept
    {
        return IIF_detail::do_min(p, amrex::get<Is>(*this)...);
    }

    template <std::size_t... Is>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real op_impl (AMREX_D_DECL(Real x, Real y, Real z), std::index_sequence<Is...>) const noexcept
    {
        return IIF_detail::do_min(AMREX_D_DECL(x,y,z), amrex::get<Is>(*this)...);
    }
};

template <class Head, class... Tail>
struct IsGPUable<IntersectionIF<Head, Tail...>, std::enable_if_t<IsGPUable<Head>::value>>
    : IsGPUable<IntersectionIF<Tail...> > {};

template <class F>
struct IsGPUable<IntersectionIF<F>, std::enable_if_t<IsGPUable<F>::value>>
    : std::true_type {};

template <class... Fs>
constexpr IntersectionIF<std::decay_t<Fs> ...>
makeIntersection (Fs&&... fs)
{
    return IntersectionIF<std::decay_t<Fs> ...>(std::forward<Fs>(fs)...);
}

}

#endif
