#ifndef AMREX_GPU_DEVICE_BUFFER_H_
#define AMREX_GPU_DEVICE_BUFFER_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuDevice.H>
#include <cstring>
#include <cstdlib>
#include <initializer_list>
#include <memory>

namespace amrex {
namespace Gpu {

template <typename T, std::enable_if_t<std::is_trivially_copyable<T>::value,int> = 0>
class Buffer
{
public:

    Buffer (std::initializer_list<T> init)
        : m_size(init.size())
    {
        if (m_size == 0) return;
#ifdef AMREX_USE_GPU
        h_data = static_cast<T*>(The_Pinned_Arena()->alloc(m_size*sizeof(T)));
#else
        h_data = static_cast<T*>(std::malloc(m_size*sizeof(T)));
#endif
        std::memcpy(h_data, init.begin(), m_size*sizeof(T));
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Arena()->alloc(m_size*sizeof(T)));
            Gpu::htod_memcpy_async(d_data, h_data, m_size*sizeof(T));
#ifdef AMREX_USE_DPCPP
            if (Gpu::onNullStream()) Gpu::synchronize();
#endif
        }
#endif
    }

    Buffer (T const* h_p, const std::size_t n)
        : m_size(n)
    {
        if (m_size == 0) return;
#ifdef AMREX_USE_GPU
        h_data = static_cast<T*>(The_Pinned_Arena()->alloc(m_size*sizeof(T)));
#else
        h_data = static_cast<T*>(std::malloc(m_size*sizeof(T)));
#endif
        std::memcpy(h_data, h_p, m_size*sizeof(T));
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Arena()->alloc(m_size*sizeof(T)));
            Gpu::htod_memcpy_async(d_data, h_data, m_size*sizeof(T));
#ifdef AMREX_USE_DPCPP
            if (Gpu::onNullStream()) Gpu::synchronize();
#endif
        }
#endif
    }

    ~Buffer () { clear(); }

    Buffer (Buffer const&) = delete;
    Buffer (Buffer &&) = delete;
    void operator= (Buffer const&) = delete;
    void operator= (Buffer &&) = delete;

    T const* data () const noexcept { return (d_data != nullptr) ? d_data : h_data; }
    T* data () noexcept { return (d_data != nullptr) ? d_data : h_data; }

    T const* hostData () const noexcept { return h_data; }
    T* hostDatat () noexcept { return h_data; }

    std::size_t size () const noexcept { return m_size; }

    void clear ()
    {
#ifdef AMREX_USE_GPU
        if (d_data) The_Arena()->free(d_data);
        if (h_data) The_Pinned_Arena()->free(h_data);
#else
        std::free(h_data);
#endif
        d_data = nullptr;
        h_data = nullptr;
    }

    T* copyToHost ()
    {
#ifdef AMREX_USE_GPU
        if (d_data)
        {
            Gpu::dtoh_memcpy_async(h_data, d_data, m_size*sizeof(T));
            Gpu::streamSynchronize();
        }
#endif
        return h_data;
    }

private:
    std::size_t m_size;
    T* d_data = nullptr;
    T* h_data = nullptr;
};

}
}

#endif
