#ifndef AMREX_GPU_DEVICE_BUFFER_H_
#define AMREX_GPU_DEVICE_BUFFER_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuContainers.H>
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
        resize(init.size());

        if (init.size() > 0) {
            std::memcpy(h_vect.data(), init.begin(), init.size()*sizeof(T));
            copyToDeviceAsync();
        }
    }

    Buffer (T const* h_p, const std::size_t n) {
        resize(n);

        if (n > 0 && h_p != nullptr) {
            std::memcpy(h_vect.data(), h_p, n*sizeof(T));
            copyToDeviceAsync();
        }
    }

    Buffer (const std::size_t n) {
        resize(n);
    }

    Buffer () = default;

    [[nodiscard]] T const* data () const noexcept {
        return (useDVect() && !d_vect.empty()) ? d_vect.data() : h_vect.data();
    }
    [[nodiscard]] T* data () noexcept {
        return (useDVect() && !d_vect.empty()) ? d_vect.data() : h_vect.data();
    }

    [[nodiscard]] T const* hostData () const noexcept { return h_vect.data(); }
    [[nodiscard]] T* hostData () noexcept { return h_vect.data(); }

    /**
    * \brief Changes the value of an element of the host (CPU) vector.
    * Does not update the device (GPU) vector, so copyToDeviceAsync()
    * needs to be called before accessing the data on the GPU.
    * \code{.cpp}
    *   amrex::Gpu::Buffer<int> buf;
    *   buf.resize(n);
    *   for (int i=0; i<n; ++i) {
    *       buf[i] = i*i;
    *   }
    *   buf.copyToDeviceAsync();
    *   int * ptr = buf.data();
    *   // Use ptr inside ParallelFor
    *   // optional:
    *   // Change values of ptr inside ParallelFor
    *   buf.copyToHost();
    *   // Use buf.hostData() or buf[] on the CPU
    * \endcode
    */
    [[nodiscard]] T& operator[] (const std::size_t i) noexcept {
        return h_vect[i];
    }

    [[nodiscard]] const T& operator[] (const std::size_t i) const noexcept {
        return h_vect[i];
    }

    [[nodiscard]] std::size_t size () const noexcept { return h_vect.size(); }

    [[nodiscard]] bool empty () const noexcept { return h_vect.size() == 0; }

    void resize (const std::size_t n) noexcept {
        h_vect.resize(n);
        if (useDVect()) {
            d_vect.resize(n);
        }
    }

    void clear () noexcept {
        h_vect.clear();
        d_vect.clear();
    }

    void shrink_to_fit () noexcept {
        h_vect.shrink_to_fit();
        d_vect.shrink_to_fit();
    }

    void reserve (const std::size_t n) noexcept {
        h_vect.reserve(n);
        if (useDVect()) {
            d_vect.reserve(n);
        }
    }

    /**
    * \brief Adds an element to the back of the host (CPU) vector.
    * Does not update the device (GPU) vector, so copyToDeviceAsync()
    * needs to be called before accessing the data on the GPU.
    * \code{.cpp}
    *   amrex::Gpu::Buffer<int> buf;
    *   buf.reserve(n);
    *   for (int i=0; i<n; ++i) {
    *       buf.push_back(i*i);
    *   }
    *   buf.copyToDeviceAsync();
    *   int * ptr = buf.data();
    *   // Use ptr inside ParallelFor
    *   // optional:
    *   // Change values of ptr inside ParallelFor
    *   buf.copyToHost();
    *   // Use buf.hostData() or buf[] on the CPU
    * \endcode
    */
    void push_back (const T& value) noexcept {
        h_vect.push_back(value);
    }

    T* copyToDeviceAsync () noexcept {
#ifdef AMREX_USE_GPU
        if (useDVect() && !h_vect.empty())
        {
            d_vect.resize(h_vect.size());
            Gpu::htod_memcpy_async(d_vect.data(), h_vect.data(), h_vect.size()*sizeof(T));
            return d_vect.data();
        }
#endif
        return h_vect.data();
    }

    T* copyToHost () noexcept {
#ifdef AMREX_USE_GPU
        if (useDVect() && !d_vect.empty())
        {
            h_vect.resize(d_vect.size());
            Gpu::dtoh_memcpy_async(h_vect.data(), d_vect.data(), d_vect.size()*sizeof(T));
            Gpu::streamSynchronize();
        }
#endif
        return h_vect.data();
    }

private:

    [[nodiscard]] bool useDVect () const noexcept {
        return Gpu::inLaunchRegion() /* && !use_unified_gpu_memory */;
    }

    DeviceVector<T> d_vect;
    PinnedVector<T> h_vect;
};

}

#endif
