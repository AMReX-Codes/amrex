#ifndef AMREX_GPU_RANGE_H_
#define AMREX_GPU_RANGE_H_
#include <AMReX_Config.H>

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuError.H>
#include <AMReX_Box.H>
#include <AMReX_TypeTraits.H>

namespace amrex {

template <typename T, std::enable_if_t<std::is_integral_v<T>,int> = 0>
bool isEmpty (T n) noexcept { return n <= 0; }

template<int dim>
AMREX_FORCE_INLINE bool isEmpty (BoxND<dim> const& b) noexcept { return b.isEmpty(); }

namespace Gpu {

namespace range_detail {

//! integer version
template <typename T, std::enable_if_t<std::is_integral_v<T>,int> = 0>
AMREX_GPU_HOST_DEVICE
Long size (T const& b) noexcept { return static_cast<Long>(b); }

template <typename T, std::enable_if_t<std::is_integral_v<T>,int> = 0>
AMREX_GPU_HOST_DEVICE
Long at (T const& /*b*/, Long offset) noexcept { return offset; }

//! BoxND version
template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE Long size (BoxND<dim> const& b) noexcept
{
    AMREX_IF_ON_DEVICE((return b.numPts();))
    AMREX_IF_ON_HOST((
        amrex::ignore_unused(b);
        return 1;
    ))
}

template<int dim>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE BoxND<dim> at (BoxND<dim> const& b, Long offset) noexcept
{
    AMREX_IF_ON_DEVICE((
        auto iv = b.atOffset(offset);
        return BoxND<dim>(iv,iv,b.ixType());
    ))
    AMREX_IF_ON_HOST((
        amrex::ignore_unused(offset);
        return b;
    ))
}

template <typename T>
struct range_impl
{
    AMREX_GPU_HOST_DEVICE
    explicit range_impl (T const& b) noexcept : m_b(b), m_n(range_detail::size(b)) {}

#ifdef AMREX_USE_SYCL
    range_impl (T const& b, Long gid, Long grange) noexcept
        : m_b(b), m_n(range_detail::size(b)), m_gid(gid), m_grange(grange) {}
#endif

    struct iterator
    {
        AMREX_GPU_HOST_DEVICE
        iterator (T const& b, Long i, Long s) noexcept : mi_b(&b), mi_i(i), mi_s(s)  {}

        AMREX_GPU_HOST_DEVICE
        void operator++ () noexcept { mi_i += mi_s; }

        AMREX_GPU_HOST_DEVICE
        bool operator!= (iterator const& rhs) const noexcept { return mi_i < rhs.mi_i; }

        AMREX_GPU_HOST_DEVICE
        T operator* () const noexcept { return range_detail::at(*mi_b,mi_i); }

    private:
        T const* mi_b;
        Long mi_i;
        Long mi_s;
    };

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    iterator begin () const noexcept {
#if defined (__SYCL_DEVICE_ONLY__)
        return iterator(m_b, m_gid, m_grange);
#else
        AMREX_IF_ON_DEVICE((
            return iterator(m_b, blockDim.x*blockIdx.x+threadIdx.x, blockDim.x*gridDim.x);
        ))
        AMREX_IF_ON_HOST((
            return iterator(m_b,0,1);
        ))
#endif
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    iterator end () const noexcept { return iterator(m_b,m_n,0); }

private:
    T m_b;
    Long m_n;
#ifdef AMREX_USE_SYCL
    Long m_gid;
    Long m_grange;
#endif
};
}

#ifdef AMREX_USE_SYCL
template <typename T>
range_detail::range_impl<T> Range (T const& b, Long gid, Long grange) noexcept {
    return range_detail::range_impl<T>(b,gid,grange);
}
#endif

template <typename T>
AMREX_GPU_HOST_DEVICE
range_detail::range_impl<T> Range (T const& b) noexcept { return range_detail::range_impl<T>(b); }

}}

#endif
