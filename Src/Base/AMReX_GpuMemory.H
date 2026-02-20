#ifndef AMREX_GPU_MEMORY_H_
#define AMREX_GPU_MEMORY_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_GpuMemory.H
 * \brief GPU memory helpers: Managed, Pinned, Deleter, DeviceScalar, and SharedMemory.
 */

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Arena.H>
#include <cstdlib>

namespace amrex::Gpu {

/**
 * \class amrex::Gpu::Managed
 * \ingroup amrex_execution
 * \brief Base class that overrides operator new/delete to use managed (unified) GPU memory.
 *
 * Deriving from this class causes heap allocation to go through
 * The_Managed_Arena() when GPU support is enabled.
 */
struct Managed {

#ifdef AMREX_USE_GPU

    //! Allocate \p len bytes from the managed arena.
    [[nodiscard]] void *operator new (std::size_t len)
    {
        return The_Managed_Arena()->alloc(len);
    }

    //! Free \p ptr back to the managed arena.
    void operator delete (void *ptr)
    {
        if (ptr) { The_Managed_Arena()->free(ptr); }
    }

#endif
};

/**
 * \class amrex::Gpu::Pinned
 * \ingroup amrex_execution
 * \brief Base class that overrides operator new/delete to use pinned host memory.
 *
 * Deriving from this class causes heap allocation to go through
 * The_Pinned_Arena() when GPU support is enabled.
 */
struct Pinned {

#ifdef AMREX_USE_GPU

    //! Allocate \p len bytes from the pinned arena.
    [[nodiscard]] void *operator new (std::size_t len)
    {
        return The_Pinned_Arena()->alloc(len);
    }

    //! Free \p ptr back to the pinned arena.
    void operator delete (void *ptr)
    {
        if (ptr) { The_Pinned_Arena()->free(ptr); }
    }

#endif
};

/**
 * \class amrex::Gpu::Deleter
 * \ingroup amrex_execution
 * \brief A callable deleter that frees memory through a given Arena.
 *
 * Suitable for use as a custom deleter with std::unique_ptr.
 */
struct Deleter {
    Arena* m_arena = nullptr; ///< The arena used for deallocation.
    //! Construct with a pointer to the arena to use for deallocation.
    Deleter (Arena* ar) noexcept : m_arena(ar) {}
    //! Free \p pt through the stored arena.
    void operator() (void* pt) const noexcept { m_arena->free(pt); }
};

/**
 * \class amrex::Gpu::DeviceScalar
 * \ingroup amrex_execution
 * \brief A single scalar value stored in GPU device memory (or on the host when GPU is off).
 *
 * \tparam T Value type (must be trivially copyable).
 */
template <class T, std::enable_if_t<std::is_trivially_copyable_v<T>,int> = 0>
struct DeviceScalar
{
    //! DeviceScalar is not copyable.
    DeviceScalar (DeviceScalar const&) = delete;
    //! DeviceScalar is not movable.
    DeviceScalar (DeviceScalar &&) = delete;
    //! DeviceScalar is not copy-assignable.
    void operator= (DeviceScalar const&) = delete;
    //! DeviceScalar is not move-assignable.
    void operator= (DeviceScalar &&) = delete;

#ifdef AMREX_USE_GPU

    //! Default-construct: allocate device storage without initializing the value.
    DeviceScalar () {
        if (Gpu::inLaunchRegion()) {
            dp = (T*)(The_Arena()->alloc(sizeof(T)));
        } else {
            dp = (T*)(std::malloc(sizeof(T)));
        }
    }

    //! Construct and initialize the scalar to \p init_val.
    explicit DeviceScalar (T init_val) {
        if (Gpu::inLaunchRegion()) {
            dp = (T*)(The_Arena()->alloc(sizeof(T)));
            Gpu::htod_memcpy(dp, &init_val, sizeof(T));
        } else {
            dp = (T*)(std::malloc(sizeof(T)));
            *dp = init_val;
        }
    }

    //! Destructor: free device storage.
    ~DeviceScalar () {
        if (Gpu::inLaunchRegion()) {
            The_Arena()->free(dp);
        } else {
            std::free(dp);
        }
    }

    //! Return a pointer to the device scalar (usable in GPU kernels).
    [[nodiscard]] T* dataPtr () { return dp; }
    //! Return a const pointer to the device scalar.
    [[nodiscard]] T const* dataPtr () const { return dp; }
    //! Copy the scalar value back to the host and return it.
    [[nodiscard]] T dataValue () const {
        if (Gpu::inLaunchRegion()) {
            T r;
            Gpu::dtoh_memcpy(&r, dp, sizeof(T));
            return r;
        } else {
            return *dp;
        }
    }

private:
    T* dp;

#else

    //! Construct and initialize the scalar to \p init_val (CPU build).
    DeviceScalar (T const& init_val) : d(init_val) {}
    //! Default-construct (CPU build).
    DeviceScalar () = default;
    ~DeviceScalar () = default;

    //! Return a pointer to the scalar.
    [[nodiscard]] T* dataPtr () { return &d; }
    //! Return a const pointer to the scalar.
    [[nodiscard]] T const* dataPtr () const { return &d; }
    //! Return the scalar value.
    [[nodiscard]] T dataValue () const { return d; }

private:
    T d;

#endif
};

#ifdef AMREX_USE_GPU

/**
 * \class amrex::Gpu::SharedMemory
 * \ingroup amrex_execution
 * \brief Provides a pointer to GPU shared (local) memory of type T inside a kernel.
 *
 * Specializations are provided for all common arithmetic types.
 * Using an unsupported type will trigger a compile-time error.
 *
 * \tparam T Element type.
 */
template <class T>
struct SharedMemory
{
    //! Return a pointer to the dynamically-sized shared memory array.
    [[nodiscard]] AMREX_GPU_DEVICE T* dataPtr () noexcept {
        static_assert(sizeof(T) < 0, "We must specialize struct SharedMemory");
        return nullptr;
    }
};

/// \cond DOXYGEN_IGNORE
#ifndef AMREX_USE_SYCL
// xxxxx SYCL todo: extern __shared__

template <>
struct SharedMemory<double>
{
    [[nodiscard]] AMREX_GPU_DEVICE double* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(double,amrex_sm_double);,
                          extern __shared__  double amrex_sm_double[];)
        return amrex_sm_double;
    }
};

template <>
struct SharedMemory<float>
{
    [[nodiscard]] AMREX_GPU_DEVICE float* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(float,amrex_sm_float);,
                          extern __shared__  float amrex_sm_float[];)
        return amrex_sm_float;
    }
};

template <>
struct SharedMemory<long>
{
    [[nodiscard]] AMREX_GPU_DEVICE long* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(long,amrex_sm_long);,
                          extern __shared__  long amrex_sm_long[];)
        return amrex_sm_long;
    }
};

template <>
struct SharedMemory<long long>
{
    [[nodiscard]] AMREX_GPU_DEVICE long long* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(long long,amrex_sm_long_long);,
                          extern __shared__  long long amrex_sm_long_long[];)
        return amrex_sm_long_long;
    }
};

template <>
struct SharedMemory<int>
{
    [[nodiscard]] AMREX_GPU_DEVICE int* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(int,amrex_sm_int);,
                          extern __shared__  int amrex_sm_int[];)
        return amrex_sm_int;
    }
};

template <>
struct SharedMemory<short>
{
    [[nodiscard]] AMREX_GPU_DEVICE short* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(short,amrex_sm_short);,
                          extern __shared__  short amrex_sm_short[];)
        return amrex_sm_short;
    }
};

template <>
struct SharedMemory<char>
{
    [[nodiscard]] AMREX_GPU_DEVICE char* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(char,amrex_sm_char);,
                          extern __shared__  char amrex_sm_char[];)
        return amrex_sm_char;
    }
};

template <>
struct SharedMemory<unsigned long>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned long* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned long,amrex_sm_ulong);,
                          extern __shared__  unsigned long amrex_sm_ulong[];)
        return amrex_sm_ulong;
    }
};

template <>
struct SharedMemory<unsigned long long>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned long long* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned long long,amrex_sm_ulonglong);,
                          extern __shared__  unsigned long long amrex_sm_ulonglong[];)
        return amrex_sm_ulonglong;
    }
};

template <>
struct SharedMemory<unsigned int>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned int* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned int,amrex_sm_uint);,
                          extern __shared__  unsigned int amrex_sm_uint[];)
        return amrex_sm_uint;
    }
};

template <>
struct SharedMemory<unsigned short>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned short* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned short,amrex_sm_ushort);,
                          extern __shared__  unsigned short amrex_sm_ushort[];)
        return amrex_sm_ushort;
    }
};

template <>
struct SharedMemory<unsigned char>
{
    [[nodiscard]] AMREX_GPU_DEVICE unsigned char* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned char,amrex_sm_uchar);,
                          extern __shared__  unsigned char amrex_sm_uchar[];)
        return amrex_sm_uchar;
    }
};

template <>
struct SharedMemory<bool>
{
    [[nodiscard]] AMREX_GPU_DEVICE bool* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(bool,amrex_sm_bool);,
                          extern __shared__  bool amrex_sm_bool[];)
        return amrex_sm_bool;
    }
};

#endif
/// \endcond

#endif

}  //namespace


#endif
