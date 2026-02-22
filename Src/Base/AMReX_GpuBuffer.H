#ifndef AMREX_GPU_DEVICE_BUFFER_H_
#define AMREX_GPU_DEVICE_BUFFER_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuBuffer.H
 * \brief A mirrored host+device buffer for trivially-copyable data.
 */

#include <AMReX_Arena.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuContainers.H>
#include <cstring>
#include <cstdlib>
#include <initializer_list>
#include <memory>

namespace amrex::Gpu {

/**
 * \class amrex::Gpu::Buffer
 * \ingroup amrex_linear_containers
 * \brief A mirrored host+device buffer for trivially-copyable data.
 *
 * Keeps a pinned host copy and a device copy in sync.  The device copy is
 * only allocated when inside a GPU launch region.  \c data() returns the
 * device pointer inside a launch region and the host pointer otherwise.
 *
 * \tparam T Element type (must be trivially copyable).
 */
template <typename T, std::enable_if_t<std::is_trivially_copyable_v<T>,int> = 0>
class Buffer
{
public:

    //! Construct from an initializer list; data are copied to both host and device.
    Buffer (std::initializer_list<T> init) {
        resize(init.size());

        if (init.size() > 0) {
            std::memcpy(h_vect.data(), init.begin(), init.size()*sizeof(T));
            copyToDeviceAsync();
        }
    }

    //! Construct from a host pointer \p h_p and size \p n; data are copied to both host and device.
    Buffer (T const* h_p, const std::size_t n) {
        resize(n);

        if (n > 0 && h_p != nullptr) {
            std::memcpy(h_vect.data(), h_p, n*sizeof(T));
            copyToDeviceAsync();
        }
    }

    //! Construct with \p n default-initialized elements (no values written).
    Buffer (const std::size_t n) {
        resize(n);
    }

    //! Default-construct an empty buffer.
    Buffer () = default;

    //! Return a const pointer to the data (device pointer inside a launch region, host otherwise).
    [[nodiscard]] T const* data () const noexcept {
        return (useDVect() && !d_vect.empty()) ? d_vect.data() : h_vect.data();
    }
    //! Return a pointer to the data (device pointer inside a launch region, host otherwise).
    [[nodiscard]] T* data () noexcept {
        return (useDVect() && !d_vect.empty()) ? d_vect.data() : h_vect.data();
    }

    //! Return a const pointer to the host (CPU) data.
    [[nodiscard]] T const* hostData () const noexcept { return h_vect.data(); }
    //! Return a pointer to the host (CPU) data.
    [[nodiscard]] T* hostData () noexcept { return h_vect.data(); }

    /**
    * \brief Return a mutable reference to element \p i of the host (CPU) vector.
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

    //! Return a const reference to element \p i of the host (CPU) vector.
    [[nodiscard]] const T& operator[] (const std::size_t i) const noexcept {
        return h_vect[i];
    }

    //! Return the number of elements.
    [[nodiscard]] std::size_t size () const noexcept { return h_vect.size(); }

    //! Return true if the buffer contains no elements.
    [[nodiscard]] bool empty () const noexcept { return h_vect.size() == 0; }

    //! Resize both host and device vectors to \p n elements.
    void resize (const std::size_t n) noexcept {
        h_vect.resize(n);
        if (useDVect()) {
            d_vect.resize(n);
        }
    }

    //! Clear both host and device vectors.
    void clear () noexcept {
        h_vect.clear();
        d_vect.clear();
    }

    //! Release excess capacity in both host and device vectors.
    void shrink_to_fit () noexcept {
        h_vect.shrink_to_fit();
        d_vect.shrink_to_fit();
    }

    //! Reserve capacity for at least \p n elements in both host and device vectors.
    void reserve (const std::size_t n) noexcept {
        h_vect.reserve(n);
        if (useDVect()) {
            d_vect.reserve(n);
        }
    }

    /**
    * \brief Append \p value to the back of the host (CPU) vector.
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

    //! Copy host data to the device asynchronously; returns a pointer to the device data.
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

    //! Copy device data back to the host synchronously; returns a pointer to the host data.
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

    //! Return true if the buffer should access the device mirror (inside GPU launch regions).
    [[nodiscard]] bool useDVect () const noexcept {
        return Gpu::inLaunchRegion() /* && !use_unified_gpu_memory */;
    }

    DeviceVector<T> d_vect;
    PinnedVector<T> h_vect;
};

}

#endif
