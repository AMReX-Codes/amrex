#ifndef AMREX_GPU_ASYNC_ARRAY_H_
#define AMREX_GPU_ASYNC_ARRAY_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuDevice.H>

#include <cstddef>
#include <cstring>
#include <cstdlib>
#include <memory>

namespace amrex {
namespace Gpu {

template <typename T>
requires (std::is_trivially_copyable_v<T>)
class AsyncArray
{
public:

    AsyncArray (T const* h_p, const std::size_t n)
    {
        if (n == 0) { return; }
        h_data = static_cast<T*>(The_Pinned_Arena()->alloc(n*sizeof(T)));
        std::memcpy(h_data, h_p, n*sizeof(T));
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Arena()->alloc(n*sizeof(T)));
            Gpu::htod_memcpy_async(d_data, h_data, n*sizeof(T));
        }
#endif
    }

    explicit AsyncArray (const std::size_t n)
        requires (std::is_standard_layout_v<T> && std::is_trivial_v<T>)

    {
        if (n == 0) { return; }
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Arena()->alloc(n*sizeof(T)));
        }
        else
#endif
        {
            h_data = static_cast<T*>(The_Pinned_Arena()->alloc(n*sizeof(T)));
        }
    }

    ~AsyncArray () { clear(); }

    AsyncArray (AsyncArray const&) = delete;
    AsyncArray (AsyncArray &&) = delete;
    void operator= (AsyncArray const&) = delete;
    void operator= (AsyncArray &&) = delete;

    [[nodiscard]] T const* data () const noexcept { return (d_data != nullptr) ? d_data : h_data; }
    [[nodiscard]] T* data () noexcept { return (d_data != nullptr) ? d_data : h_data; }
    void clear ()
    {
        if (d_data != nullptr) { Gpu::freeAsync(The_Arena(), d_data); }
        if (h_data != nullptr) { Gpu::freeAsync(The_Pinned_Arena(), h_data); }
        d_data = nullptr;
        h_data = nullptr;
    }

    void copyToHost (T* h_p, std::size_t n) const
    {
        if (n == 0) { return; }
#ifdef AMREX_USE_GPU
        if (d_data)
        {
            Gpu::dtoh_memcpy(h_p, d_data, n*sizeof(T));
        }
        else
#endif
        if (h_data)
        {
            std::memcpy(h_p, h_data, n*sizeof(T));
        }
    }

private:
    T* d_data = nullptr;
    T* h_data = nullptr;
};

}

using Gpu::AsyncArray;
}

#endif
