#ifndef AMREX_GPU_MEMORY_H_
#define AMREX_GPU_MEMORY_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Arena.H>
#include <cstdlib>

namespace amrex {
namespace Gpu {

struct Managed {

#ifdef AMREX_USE_GPU

    void *operator new (std::size_t len)
    {
        return The_Managed_Arena()->alloc(len);
    }
    
    void operator delete (void *ptr)
    {
        The_Managed_Arena()->free(ptr);
    }

#endif
};

struct Pinned {
    
#ifdef AMREX_USE_GPU

    void *operator new (std::size_t len)
    {
        return The_Pinned_Arena()->alloc(len);
    }
    
    void operator delete (void *ptr)
    {
        The_Pinned_Arena()->free(ptr);
    }
    
#endif
};

struct Deleter {
    Arena* m_arena = nullptr;
    Deleter (Arena* ar) noexcept : m_arena(ar) {}
    void operator() (void* pt) const noexcept { m_arena->free(pt); }
};

template <class T, typename = amrex::EnableIf_t<AMREX_IS_TRIVIALLY_COPYABLE(T)> >
struct DeviceScalar
{
    DeviceScalar (DeviceScalar const&) = delete;
    DeviceScalar (DeviceScalar &&) = delete;
    void operator= (DeviceScalar const&) = delete;
    void operator= (DeviceScalar &&) = delete;

#ifdef AMREX_USE_GPU

    DeviceScalar () {
        if (Gpu::inLaunchRegion()) {
            dp = (T*)(The_Device_Arena()->alloc(sizeof(T)));
        } else {
            dp = (T*)(std::malloc(sizeof(T)));
        }
    }

    explicit DeviceScalar (T init_val) {
        if (Gpu::inLaunchRegion()) {
            dp = (T*)(The_Device_Arena()->alloc(sizeof(T)));
            Gpu::htod_memcpy(dp, &init_val, sizeof(T));
        } else {
            dp = (T*)(std::malloc(sizeof(T)));
            *dp = init_val;
        }
    }

    ~DeviceScalar () {
        if (Gpu::inLaunchRegion()) {
            The_Device_Arena()->free(dp);
        } else {
            std::free(dp);
        }
    }

    T* dataPtr () { return dp; }
    T const* dataPtr () const { return dp; }
    T dataValue () const {
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

    DeviceScalar (T init_val) : d(init_val) {}
    DeviceScalar () {}
    ~DeviceScalar () {}

    T* dataPtr () { return &d; }
    T const* dataPtr () const { return &d; }
    T dataValue () const { return d; }

private:
    T d;

#endif
};

#ifdef AMREX_USE_GPU

template <class T>
struct SharedMemory
{
    AMREX_GPU_DEVICE T* dataPtr () noexcept {
        static_assert(sizeof(T) < 0, "We must specialize struct SharedMemory");
        return nullptr;
    }
};

template <>
struct SharedMemory<double>
{
    AMREX_GPU_DEVICE double* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(double,amrex_sm_double);,
                          extern __shared__  double amrex_sm_double[];)
        return amrex_sm_double;
    }
};

template <>
struct SharedMemory<float>
{
    AMREX_GPU_DEVICE float* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(float,amrex_sm_float);,
                          extern __shared__  float amrex_sm_float[];)
        return amrex_sm_float;
    }
};

template <>
struct SharedMemory<long>
{
    AMREX_GPU_DEVICE long* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(long,amrex_sm_long);,
                          extern __shared__  long amrex_sm_long[];)
        return amrex_sm_long;
    }
};

template <>
struct SharedMemory<int>
{
    AMREX_GPU_DEVICE int* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(int,amrex_sm_int);,
                          extern __shared__  int amrex_sm_int[];)
        return amrex_sm_int;
    }
};

template <>
struct SharedMemory<short>
{
    AMREX_GPU_DEVICE short* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(short,amrex_sm_short);,
                          extern __shared__  short amrex_sm_short[];)
        return amrex_sm_short;
    }
};

template <>
struct SharedMemory<char>
{
    AMREX_GPU_DEVICE char* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(char,amrex_sm_char);,
                          extern __shared__  char amrex_sm_char[];)
        return amrex_sm_char;
    }
};

template <>
struct SharedMemory<unsigned long>
{
    AMREX_GPU_DEVICE unsigned long* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned long,amrex_sm_ulong);,
                          extern __shared__  unsigned long amrex_sm_ulong[];)
        return amrex_sm_ulong;
    }
};

template <>
struct SharedMemory<unsigned int>
{
    AMREX_GPU_DEVICE unsigned int* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned int,amrex_sm_uint);,
                          extern __shared__  unsigned int amrex_sm_uint[];)
        return amrex_sm_uint;
    }
};

template <>
struct SharedMemory<unsigned short>
{
    AMREX_GPU_DEVICE unsigned short* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned short,amrex_sm_ushort);,
                          extern __shared__  unsigned short amrex_sm_ushort[];)
        return amrex_sm_ushort;
    }
};

template <>
struct SharedMemory<unsigned char>
{
    AMREX_GPU_DEVICE unsigned char* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(unsigned char,amrex_sm_uchar);,
                          extern __shared__  unsigned char amrex_sm_uchar[];)
        return amrex_sm_uchar;
    }
};

template <>
struct SharedMemory<bool>
{
    AMREX_GPU_DEVICE bool* dataPtr () noexcept {
        AMREX_HIP_OR_CUDA(HIP_DYNAMIC_SHARED(bool,amrex_sm_bool);,
                          extern __shared__  bool amrex_sm_bool[];)
        return amrex_sm_bool;
    }
};

#endif

}}  //namespace


#endif
