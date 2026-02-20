#ifndef AMREX_GPU_UTILITY_H_
#define AMREX_GPU_UTILITY_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuUtility.H
 * \brief Miscellaneous GPU utility functions and the StreamIter class.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuTypes.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_Extension.H>
#include <AMReX_REAL.H>
#include <AMReX_INT.H>
#include <AMReX_Array.H>
#include <AMReX_Array4.H>
#include <iosfwd>
#include <cmath>
#include <cstring>

#ifdef AMREX_USE_CUDA
#include <cuda.h>
#include <curand_kernel.h>   // Is this needed here?
#endif

namespace amrex {
namespace Gpu {

    /// \ingroup amrex_execution
    //! Load element (i,j,k) of \p a using a read-only cache (L1 texture on CUDA).
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T LDG (Array4<T> const& a, int i, int j, int k) noexcept {
#if defined(AMREX_USE_CUDA)
        AMREX_IF_ON_DEVICE((return __ldg(a.ptr(i,j,k));))
        AMREX_IF_ON_HOST((return a(i,j,k);))
#else
        return a(i,j,k);
#endif
    }

    /// \ingroup amrex_execution
    //! Load component \p n at (i,j,k) of \p a using a read-only cache (L1 texture on CUDA).
    template <typename T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T LDG (Array4<T> const& a, int i, int j, int k, int n) noexcept {
#if defined(AMREX_USE_CUDA)
        AMREX_IF_ON_DEVICE((return __ldg(a.ptr(i,j,k,n));))
        AMREX_IF_ON_HOST((return a(i,j,k,n);))
#else
        return a(i,j,k,n);
#endif
    }

    /// \ingroup amrex_execution
    //! Return true if \p p points to CUDA/HIP/SYCL managed (unified) memory.
    inline bool isManaged (void const* p) noexcept {
#ifdef AMREX_USE_CUDA
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_IS_MANAGED;
        unsigned int is_managed = 0;
        void* data[] = { (void*)(&is_managed) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS && is_managed;
#elif defined(AMREX_USE_SYCL)
        auto type = sycl::get_pointer_type(p,Device::syclContext());
        return type == sycl::usm::alloc::shared;
#else
        amrex::ignore_unused(p);
        return false;
#endif
    }

    /// \ingroup amrex_execution
    //! Return true if \p p points to device memory.
    inline bool isDevicePtr (void const* p) noexcept {
#if defined(AMREX_USE_HIP)
        hipPointerAttribute_t attrib;
        hipError_t r = hipPointerGetAttributes(&attrib, p);
#if defined(HIP_VERSION_MAJOR) && (HIP_VERSION_MAJOR < 6)
        return r == hipSuccess && attrib.memoryType == hipMemoryTypeDevice;
#else
        return r == hipSuccess && attrib.type == hipMemoryTypeDevice;
#endif // (HIP_VERSION_MAJOR) && (HIP_VERSION_MAJOR < 6)
#elif defined(AMREX_USE_CUDA)
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_MEMORY_TYPE;
        CUmemorytype mem_type = static_cast<CUmemorytype>(0);
        void* data[] = { (void*)(&mem_type) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS && mem_type == CU_MEMORYTYPE_DEVICE;
#elif defined(AMREX_USE_SYCL)
        auto type = sycl::get_pointer_type(p,Device::syclContext());
        return type == sycl::usm::alloc::device;
#else
        amrex::ignore_unused(p);
        return false;
#endif
    }

    /// \ingroup amrex_execution
    //! Return true if \p p points to pinned (page-locked) host memory.
    inline bool isPinnedPtr (void const* p) noexcept {
#if defined(AMREX_USE_HIP)
        hipPointerAttribute_t attrib;
        hipError_t r = hipPointerGetAttributes(&attrib, p);
#if defined(HIP_VERSION_MAJOR) && (HIP_VERSION_MAJOR < 6)
        return r == hipSuccess && attrib.memoryType == hipMemoryTypeHost;
#else
        return r == hipSuccess && attrib.type == hipMemoryTypeHost;
#endif // (HIP_VERSION_MAJOR) && (HIP_VERSION_MAJOR < 6)
#elif defined(AMREX_USE_CUDA)
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_MEMORY_TYPE;
        CUmemorytype mem_type = static_cast<CUmemorytype>(0);
        void* data[] = { (void*)(&mem_type) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS && mem_type == CU_MEMORYTYPE_HOST;
#elif defined(AMREX_USE_SYCL)
        auto type = sycl::get_pointer_type(p,Device::syclContext());
        return type == sycl::usm::alloc::host;
#else
        amrex::ignore_unused(p);
        return false;
#endif
    }

    /// \ingroup amrex_execution
    //! Return true if \p p is any GPU-accessible pointer (device, host, unified, or array).
    inline bool isGpuPtr (void const* p) noexcept {
#if defined(AMREX_USE_HIP)
        if (isManaged(p)) { // We might be using CUDA/NVCC
            return true;
        } else {
            hipPointerAttribute_t attrib;
            hipError_t r = hipPointerGetAttributes(&attrib, p);
#if defined(HIP_VERSION_MAJOR) && (HIP_VERSION_MAJOR < 6)
            return r == hipSuccess &&
                (attrib.memoryType == hipMemoryTypeHost   ||
                 attrib.memoryType == hipMemoryTypeDevice);
#else
            return r == hipSuccess &&
                (attrib.type == hipMemoryTypeHost   ||
                 attrib.type == hipMemoryTypeDevice);
#endif // (HIP_VERSION_MAJOR) && (HIP_VERSION_MAJOR < 6)
        }
#elif defined(AMREX_USE_CUDA)
        CUpointer_attribute attrib = CU_POINTER_ATTRIBUTE_MEMORY_TYPE;
        CUmemorytype mem_type = static_cast<CUmemorytype>(0);
        void* data[] = { (void*)(&mem_type) };
        CUresult r = cuPointerGetAttributes(1, &attrib, data, (CUdeviceptr)p);
        return r == CUDA_SUCCESS &&
            (mem_type == CU_MEMORYTYPE_HOST   ||
             mem_type == CU_MEMORYTYPE_DEVICE ||
             mem_type == CU_MEMORYTYPE_ARRAY  ||
             mem_type == CU_MEMORYTYPE_UNIFIED);
#elif defined(AMREX_USE_SYCL)
        auto type = sycl::get_pointer_type(p,Device::syclContext());
        return type != sycl::usm::alloc::unknown;
#else
        amrex::ignore_unused(p);
        return false;
#endif
    }

    /// \ingroup amrex_execution
    //! Return true if \p m is a NaN; works on both host and device.
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isnan (T m) noexcept
    {
#if defined(_WIN32)
        AMREX_IF_ON_DEVICE((return m != m;))
        AMREX_IF_ON_HOST((return std::isnan(m);))
#elif defined(AMREX_USE_SYCL)
        return sycl::isnan(m);
#else
        return std::isnan(m);
#endif
    }

    /// \ingroup amrex_execution
    //! Return true if \p m is infinite; works on both host and device.
    template <class T>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool isinf (T m) noexcept
    {
#if defined(_WIN32)
        AMREX_IF_ON_DEVICE((return (2*m == m) && (m != 0);))
        AMREX_IF_ON_HOST((return std::isinf(m);))
#elif defined(AMREX_USE_SYCL)
        return sycl::isinf(m);
#else
        return std::isinf(m);
#endif
    }

    /**
     * \class amrex::Gpu::StreamItInfo
     * \ingroup amrex_execution
     * \brief Configuration options for a StreamIter.
     */
    struct StreamItInfo
    {
        bool device_sync; ///< Whether to synchronize the device after each stream iteration.
        //! Default-construct with device_sync set from the current no-sync region state.
        StreamItInfo () noexcept
            : device_sync(!Gpu::inNoSyncRegion()) {}
        //! Disable device synchronization for this StreamIter.
        StreamItInfo& DisableDeviceSync () noexcept {
            device_sync = false;
            return *this;
        }
    };

    /**
     * \class amrex::Gpu::StreamIter
     * \ingroup amrex_execution
     * \brief Iterator over GPU streams; advances to the next stream on each increment.
     *
     * Used to distribute work across multiple GPU streams.  The iterator
     * optionally synchronizes the device between streams.
     */
    class StreamIter
    {
    public:
        //! Construct iterating over \p n streams.
        StreamIter (int n, bool is_thread_safe=true) noexcept;
        //! Construct iterating over \p n streams with the given \p info options.
        StreamIter (int n, const StreamItInfo& info, bool is_thread_safe=true) noexcept;

        //! Destructor: synchronizes the device if requested.
        ~StreamIter ();

        //! StreamIter is not copyable.
        StreamIter (StreamIter const&) = delete;
        //! StreamIter is not movable.
        StreamIter (StreamIter &&) = delete;
        //! StreamIter is not copy-assignable.
        void operator= (StreamIter const&) = delete;
        //! StreamIter is not move-assignable.
        void operator= (StreamIter &&) = delete;

        //! Return the current stream index.
        int operator() () const noexcept { return m_i; }

        //! Return true if there are more streams to iterate over.
        [[nodiscard]] bool isValid () const noexcept { return m_i < m_n; }

#if !defined(AMREX_USE_GPU)
        //! Advance to the next stream.
        void operator++ () noexcept { ++m_i; }
#else
        //! Advance to the next stream.
        void operator++ () noexcept;
#endif

    private:
        void init () noexcept; // NOLINT

        int m_n;
        int m_i;
        bool m_threadsafe;
        bool m_sync;
    };

/// \ingroup amrex_execution
//! GPU-compatible memcpy; falls back to std::memcpy on the host.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void* memcpy (void* dest, const void* src, std::size_t count)
{
#if defined(__HIP_DEVICE_COMPILE__) && defined(AMREX_USE_HIP)
    return ::memcpy(dest, src, count);
#else
    return std::memcpy(dest, src, count);
#endif
}

} // namespace Gpu

#ifdef AMREX_USE_GPU
/// \ingroup amrex_execution
//! Print a dim3 value to an output stream.
std::ostream& operator<< (std::ostream& os, const dim3& d);
#endif

using Gpu::isnan;
using Gpu::isinf;

} // namespace amrex

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
/* Somehow I cannot get cudaMemcpyFromSymbol to work. */
#define AMREX_GET_DEVICE_FUNC_PTR(FP_t, sym) \
    [] () { \
        auto* amrex_i_hp = (FP_t*) amrex::The_Pinned_Arena()->alloc(sizeof(FP_t)); \
        amrex::single_task([=] AMREX_GPU_DEVICE () { *amrex_i_hp = sym; }); \
        amrex::Gpu::streamSynchronize(); \
        auto amrex_i_r = *amrex_i_hp; \
        amrex::The_Pinned_Arena()->free(amrex_i_hp); \
        return amrex_i_r; \
    } ()
#elif !defined(AMREX_USE_SYCL)
#define AMREX_GET_DEVICE_FUNC_PTR(FP_t, sym) sym
#endif

#endif
