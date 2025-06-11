#ifndef AMREX_EB2_IF_UNION_H_
#define AMREX_EB2_IF_UNION_H_
#include <AMReX_Config.H>

#include <AMReX_EB2_IF_Base.H>
#include <AMReX_Array.H>
#include <AMReX_Tuple.H>

#include <algorithm>
#include <utility>

namespace amrex::EB2 {

// For all implicit functions, >0: body; =0: boundary; <0: fluid

// Union of bodies

namespace UIF_detail {
    template <typename F>
    [[nodiscard]] inline Real do_max (const RealArray& p, F&& f) noexcept
    {
        return std::forward<F>(f)(p);
    }

    template <typename F, typename... Fs>
    [[nodiscard]] inline Real do_max (const RealArray& p, F&& f, Fs&... fs) noexcept
    {
        return amrex::max(std::forward<F>(f)(p), do_max(p, std::forward<Fs>(fs)...));
    }

    template <typename F>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real do_max (AMREX_D_DECL(Real x, Real y, Real z), F&& f) noexcept
    {
        return std::forward<F>(f)(AMREX_D_DECL(x,y,z));
    }

    template <typename F, typename... Fs>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real do_max (AMREX_D_DECL(Real x, Real y, Real z), F&& f, Fs&... fs) noexcept
    {
        return amrex::max(std::forward<F>(f)(AMREX_D_DECL(x,y,z)), do_max(AMREX_D_DECL(x,y,z), std::forward<Fs>(fs)...));
    }
}

template <class... Fs>
class UnionIF
    : public GpuTuple<Fs...>
{
public:
    using GpuTuple<Fs...>::GpuTuple;

    [[nodiscard]] inline Real operator() (const RealArray& p) const noexcept
    {
        return op_impl(p, std::make_index_sequence<sizeof...(Fs)>());
    }

    template <class U=UnionIF<Fs...>, std::enable_if_t<IsGPUable<U>::value,int> = 0>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        return op_impl(AMREX_D_DECL(x,y,z), std::make_index_sequence<sizeof...(Fs)>());
    }

protected:

    template <std::size_t... Is>
    [[nodiscard]] inline Real op_impl (const RealArray& p, std::index_sequence<Is...>) const noexcept
    {
        return UIF_detail::do_max(p, amrex::get<Is>(*this)...);
    }

    template <std::size_t... Is>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE inline
    Real op_impl (AMREX_D_DECL(Real x, Real y, Real z), std::index_sequence<Is...>) const noexcept
    {
        return UIF_detail::do_max(AMREX_D_DECL(x,y,z), amrex::get<Is>(*this)...);
    }
};

template <class Head, class... Tail>
struct IsGPUable<UnionIF<Head, Tail...>, std::enable_if_t<IsGPUable<Head>::value>>
    : IsGPUable<UnionIF<Tail...> > {};

template <class F>
struct IsGPUable<UnionIF<F>, std::enable_if_t<IsGPUable<F>::value>>
    : std::true_type {};

template <class... Fs>
constexpr UnionIF<std::decay_t<Fs> ...>
makeUnion (Fs&&... fs)
{
    return UnionIF<std::decay_t<Fs> ...>(std::forward<Fs>(fs)...);
}

}

#endif
