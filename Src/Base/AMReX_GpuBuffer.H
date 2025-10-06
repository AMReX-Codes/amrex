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

namespace amrex::Gpu {

template <typename T, std::enable_if_t<std::is_trivially_copyable_v<T>,int> = 0>
class Buffer
{
public:

    Buffer (std::initializer_list<T> init) {
        if (init.size() == 0) { return; }

        define(init.size());
        std::memcpy(h_data, init.begin(), m_size*sizeof(T));
        copyToDeviceAsync();
    }

    Buffer (T const* h_p, const std::size_t n) {
        if (n == 0) { return; }

        define(n);
        std::memcpy(h_data, h_p, m_size*sizeof(T));
        copyToDeviceAsync();
    }

    Buffer (const std::size_t n) {
        define(n);
    }

    Buffer () = default;

    ~Buffer () { clear(); }

    Buffer (Buffer const&) = default; // FIX
    Buffer (Buffer &&) = default; // FIX
    Buffer& operator= (Buffer const&) = default; // FIX
    Buffer& operator= (Buffer &&) = default; // FIX

    [[nodiscard]] T const* data () const noexcept { return (d_data != nullptr) ? d_data : h_data; }
    [[nodiscard]] T* data () noexcept { return (d_data != nullptr) ? d_data : h_data; }

    [[nodiscard]] T const* hostData () const noexcept { return h_data; }
    [[nodiscard]] T* hostData () noexcept { return h_data; }

    [[nodiscard]] T& operator[] (const std::size_t i) noexcept {
        AMREX_ASSERT(i < m_size);
        return h_data[i];
    }

    [[nodiscard]] const T& operator[] (const std::size_t i) const noexcept {
        AMREX_ASSERT(i < m_size);
        return h_data[i];
    }

    [[nodiscard]] std::size_t size () const noexcept { return m_size; }

    void define (const std::size_t n) noexcept {
        if (m_size != 0) {
            clear();
        }
        m_size = n;
        if (n == 0) {
            return;
        }
#ifdef AMREX_USE_GPU
        h_data = static_cast<T*>(The_Pinned_Arena()->alloc(m_size*sizeof(T)));
#else
        h_data = static_cast<T*>(std::malloc(m_size*sizeof(T)));
#endif
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            d_data = static_cast<T*>(The_Arena()->alloc(m_size*sizeof(T)));
        }
#endif
    }

    void clear ()
    {
#ifdef AMREX_USE_GPU
        if (d_data) { The_Arena()->free(d_data); }
        if (h_data) { The_Pinned_Arena()->free(h_data); }
#else
        std::free(h_data);
#endif
        d_data = nullptr;
        h_data = nullptr;
    }

    T* copyToDeviceAsync () noexcept {
#ifdef AMREX_USE_GPU
        if (d_data && Gpu::inLaunchRegion())
        {
            Gpu::htod_memcpy_async(d_data, h_data, m_size*sizeof(T));
            return d_data;
        }
#endif
        return h_data;
    }

    T* copyToHost () noexcept {
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
    std::size_t m_size = 0;
    T* d_data = nullptr;
    T* h_data = nullptr;
};

}

#endif
