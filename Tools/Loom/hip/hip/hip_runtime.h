// HIP source dialect for AMReX_GPU_BACKEND=LOOM.
//
// AMReX's HIP code paths are reused unchanged. This header replaces ROCm's
// <hip/hip_runtime.h>:
//  * host runtime API (the subset AMReX calls) implemented on native libhrx in
//    Src/Base/AMReX_Loom.cpp;
//  * device vocabulary (thread ids, barriers, atomics, shuffles) implemented
//    on Loom C++ operations when the Loom device pass reads the file
//    (AMREX_LOOM_DEVICE_PASS), and as inert host definitions otherwise;
// Kernel launches live in Src/Base/AMReX_LoomLaunch.H.
#ifndef AMREX_LOOM_HIP_RUNTIME_H_
#define AMREX_LOOM_HIP_RUNTIME_H_

#include <cstddef>
#include <cstdint>
#include <cstring>
#include <type_traits>

// AMREX_LOOM_DEVICE_PASS is defined exactly when Loom's C++ importer reads
// the translation unit. This is the only test of the importer's predefined
// __loom__: the rest of the LOOM backend uses AMREX_LOOM_DEVICE_PASS.
// AMReX sources under Src/ test neither; they see the device pass only
// through the HIP vocabulary below (__HIP_DEVICE_COMPILE__). Replace this
// test when Loom no longer predefines __loom__.
#if defined(__loom__)
#define AMREX_LOOM_DEVICE_PASS 1
#endif

#if defined(AMREX_LOOM_DEVICE_PASS)
#include <loomcxx/atomic.h>
#include <loomcxx/buffer.h>
#include <loomcxx/kernel.h>
#define __HIP_DEVICE_COMPILE__ 1
#define AMREX_LOOM_DEVICE_INLINE [[loom::inline]] inline
#else
#include <cstdio>
#include <cstdlib>
#include <string>
#include <typeinfo>
#define AMREX_LOOM_DEVICE_INLINE inline
#endif

#ifndef __HIP_PLATFORM_AMD__
#define __HIP_PLATFORM_AMD__ 1
#endif

#define HIP_VERSION_MAJOR 7
#define HIP_VERSION_MINOR 0
#define HIP_VERSION_PATCH 0
#define HIP_VERSION (HIP_VERSION_MAJOR * 10000000 + HIP_VERSION_MINOR * 100000)

// ---------------------------------------------------------------------------
// Qualifiers
// ---------------------------------------------------------------------------
#define __host__
#define __device__
#define __global__
#define __constant__
// Documented exclusion: the LOOM backend has no managed memory.
#define __managed__ _Pragma("GCC error \"AMReX LOOM backend: __managed__ variables are not supported\"")
#define __forceinline__ inline __attribute__((always_inline))
#define __launch_bounds__(...)
#if defined(AMREX_LOOM_DEVICE_PASS)
#define __shared__ [[loom::workgroup]]
#else
#define __shared__ static
#endif
#define HIP_DYNAMIC_SHARED(type, var) \
  __shared__ type var[AMREX_LOOM_DYNAMIC_SHARED_BYTES / sizeof(type)];
#define AMREX_LOOM_DYNAMIC_SHARED_BYTES 4096

// ---------------------------------------------------------------------------
// Runtime types
// ---------------------------------------------------------------------------
// A stream is a native HRX stream (hrx_stream_t), as AMReX's gpuStream_t is
// under AMREX_USE_LOOM. The other handle tag names match ROCm's.
typedef struct hrx_stream_s* hipStream_t;
typedef struct ihipEvent_t* hipEvent_t;
typedef int hipDevice_t;
typedef struct ihipMemPoolHandle_t* hipMemPool_t;

enum hipError_t {
  hipSuccess = 0,
  hipErrorInvalidValue = 1,
  hipErrorOutOfMemory = 2,
  hipErrorNotInitialized = 3,
  hipErrorInvalidDevice = 101,
  hipErrorNotReady = 600,
  hipErrorNotSupported = 801,
  hipErrorUnknown = 999,
};
#define hipCpuDeviceId (-1)

enum hipMemcpyKind {
  hipMemcpyHostToHost = 0,
  hipMemcpyHostToDevice = 1,
  hipMemcpyDeviceToHost = 2,
  hipMemcpyDeviceToDevice = 3,
  hipMemcpyDefault = 4,
};

enum hipMemoryType {
  hipMemoryTypeUnregistered = 0,
  hipMemoryTypeHost = 1,
  hipMemoryTypeDevice = 2,
  hipMemoryTypeManaged = 3,
  hipMemoryTypeArray = 10,
  hipMemoryTypeUnified = 11,
};

struct hipPointerAttribute_t {
  hipMemoryType type;
  int device;
  void* devicePointer;
  void* hostPointer;
  int isManaged;
  unsigned allocationFlags;
};

struct hipUUID {
  char bytes[16];
};

struct hipDeviceProp_t {
  char name[256];
  hipUUID uuid;
  std::size_t totalGlobalMem;
  std::size_t sharedMemPerBlock;
  int regsPerBlock;
  int warpSize;
  int maxThreadsPerBlock;
  int maxThreadsDim[3];
  int maxGridSize[3];
  int clockRate;
  std::size_t totalConstMem;
  int major;
  int minor;
  int multiProcessorCount;
  int maxThreadsPerMultiProcessor;
  int managedMemory;
  int concurrentManagedAccess;
  int integrated;
  int pciBusID;
  int pciDeviceID;
  int pciDomainID;
  char gcnArchName[256];
};

enum hipDeviceAttribute_t {
  hipDeviceAttributeMemoryPoolsSupported = 1,
  hipDeviceAttributeConcurrentManagedAccess = 2,
  hipDeviceAttributeManagedMemory = 3,
};

enum hipMemPoolAttr {
  hipMemPoolAttrReleaseThreshold = 1,
};

enum hipMemoryAdvise {
  hipMemAdviseSetReadMostly = 1,
  hipMemAdviseUnsetReadMostly = 2,
  hipMemAdviseSetPreferredLocation = 3,
  hipMemAdviseUnsetPreferredLocation = 4,
  hipMemAdviseSetAccessedBy = 5,
  hipMemAdviseUnsetAccessedBy = 6,
  hipMemAdviseSetCoarseGrain = 100,
  hipMemAdviseUnsetCoarseGrain = 101,
};

#define hipHostMallocDefault 0x0u
#define hipHostMallocPortable 0x1u
#define hipHostMallocMapped 0x2u
#define hipHostMallocWriteCombined 0x4u
#define hipHostMallocNonCoherent 0x40000000u
#define hipStreamDefault 0x0u
#define hipStreamNonBlocking 0x1u
#define hipMemAttachGlobal 0x1u

typedef void (*hipStreamCallback_t)(hipStream_t stream, hipError_t status,
                                    void* user_data);


// ---------------------------------------------------------------------------
// Runtime API (host; declarations are also visible to the device pass)
// ---------------------------------------------------------------------------
hipError_t hipGetLastError();
hipError_t hipPeekAtLastError();
const char* hipGetErrorString(hipError_t error);
const char* hipGetErrorName(hipError_t error);
hipError_t hipGetDeviceCount(int* count);
hipError_t hipSetDevice(int device);
hipError_t hipGetDevice(int* device);
hipError_t hipGetDeviceProperties(hipDeviceProp_t* prop, int device);
hipError_t hipDeviceGetAttribute(int* value, hipDeviceAttribute_t attribute,
                                 int device);
hipError_t hipDeviceGetUuid(hipUUID* uuid, hipDevice_t device);
hipError_t hipDeviceSynchronize();
hipError_t hipDeviceReset();
hipError_t hipMemGetInfo(std::size_t* free_bytes, std::size_t* total_bytes);
hipError_t hipMalloc(void** pointer, std::size_t bytes);
hipError_t hipMallocAsync(void** pointer, std::size_t bytes,
                          hipStream_t stream);
hipError_t hipMallocAsync(void** pointer, std::size_t bytes, hipMemPool_t pool,
                          hipStream_t stream);
hipError_t hipMallocManaged(void** pointer, std::size_t bytes,
                            unsigned flags = hipMemAttachGlobal);
hipError_t hipHostMalloc(void** pointer, std::size_t bytes,
                         unsigned flags = hipHostMallocDefault);
hipError_t hipFree(void* pointer);
hipError_t hipFreeAsync(void* pointer, hipStream_t stream);
hipError_t hipHostFree(void* pointer);
hipError_t hipMemcpy(void* dst, const void* src, std::size_t bytes,
                     hipMemcpyKind kind);
hipError_t hipMemcpyAsync(void* dst, const void* src, std::size_t bytes,
                          hipMemcpyKind kind, hipStream_t stream = nullptr);
hipError_t hipMemset(void* dst, int value, std::size_t bytes);
hipError_t hipMemsetAsync(void* dst, int value, std::size_t bytes,
                          hipStream_t stream = nullptr);
hipError_t hipMemcpyToSymbolAsync(const void* symbol, const void* src,
                                  std::size_t bytes, std::size_t offset,
                                  hipMemcpyKind kind, hipStream_t stream);
hipError_t hipMemcpyFromSymbolAsync(void* dst, const void* symbol,
                                    std::size_t bytes, std::size_t offset,
                                    hipMemcpyKind kind, hipStream_t stream);
hipError_t hipMemAdvise(const void* pointer, std::size_t bytes,
                        hipMemoryAdvise advice, int device);
hipError_t hipMemPrefetchAsync(const void* pointer, std::size_t bytes,
                               int device, hipStream_t stream = nullptr);
hipError_t hipPointerGetAttributes(hipPointerAttribute_t* attributes,
                                   const void* pointer);
hipError_t hipDeviceGetMemPool(hipMemPool_t* pool, int device);
hipError_t hipMemPoolSetAttribute(hipMemPool_t pool, hipMemPoolAttr attribute,
                                  void* value);
hipError_t hipMemPoolGetAttribute(hipMemPool_t pool, hipMemPoolAttr attribute,
                                  void* value);
hipError_t hipStreamCreate(hipStream_t* stream);
hipError_t hipStreamCreateWithFlags(hipStream_t* stream, unsigned flags);
hipError_t hipStreamDestroy(hipStream_t stream);
hipError_t hipStreamSynchronize(hipStream_t stream);
hipError_t hipStreamQuery(hipStream_t stream);
hipError_t hipStreamGetDevice(hipStream_t stream, hipDevice_t* device);
hipError_t hipStreamAddCallback(hipStream_t stream,
                                hipStreamCallback_t callback, void* user_data,
                                unsigned flags);
hipError_t hipStreamWaitEvent(hipStream_t stream, hipEvent_t event,
                              unsigned flags);
hipError_t hipEventCreate(hipEvent_t* event);
hipError_t hipEventDestroy(hipEvent_t event);
hipError_t hipEventRecord(hipEvent_t event, hipStream_t stream = nullptr);
hipError_t hipEventSynchronize(hipEvent_t event);
hipError_t hipEventElapsedTime(float* ms, hipEvent_t start, hipEvent_t stop);


// ---------------------------------------------------------------------------
// Device vocabulary
// ---------------------------------------------------------------------------
struct amrex_loom_uint3 {
  unsigned int x, y, z;
};

#if defined(AMREX_LOOM_DEVICE_PASS)
AMREX_LOOM_DEVICE_INLINE amrex_loom_uint3 amrex_loom_thread_idx() {
  return {static_cast<unsigned>(loom::kernel::workitem::id.x),
          static_cast<unsigned>(loom::kernel::workitem::id.y),
          static_cast<unsigned>(loom::kernel::workitem::id.z)};
}
AMREX_LOOM_DEVICE_INLINE amrex_loom_uint3 amrex_loom_block_idx() {
  return {static_cast<unsigned>(loom::kernel::workgroup::id.x),
          static_cast<unsigned>(loom::kernel::workgroup::id.y),
          static_cast<unsigned>(loom::kernel::workgroup::id.z)};
}
AMREX_LOOM_DEVICE_INLINE amrex_loom_uint3 amrex_loom_block_dim() {
  return {static_cast<unsigned>(loom::kernel::workgroup::size.x),
          static_cast<unsigned>(loom::kernel::workgroup::size.y),
          static_cast<unsigned>(loom::kernel::workgroup::size.z)};
}
AMREX_LOOM_DEVICE_INLINE amrex_loom_uint3 amrex_loom_grid_dim() {
  return {static_cast<unsigned>(loom::kernel::workgroup::count.x),
          static_cast<unsigned>(loom::kernel::workgroup::count.y),
          static_cast<unsigned>(loom::kernel::workgroup::count.z)};
}
AMREX_LOOM_DEVICE_INLINE void __syncthreads() {
  loom::kernel::workgroup::barrier();
}
AMREX_LOOM_DEVICE_INLINE void __threadfence() {
  loom::buffer::fence<loom::atomic::ordering::seq_cst,
                      loom::atomic::scope::device>();
}
AMREX_LOOM_DEVICE_INLINE void __threadfence_block() {
  loom::buffer::fence<loom::atomic::ordering::acq_rel,
                      loom::atomic::scope::workgroup>();
}
#else
inline amrex_loom_uint3 amrex_loom_thread_idx() { return {0, 0, 0}; }
inline amrex_loom_uint3 amrex_loom_block_idx() { return {0, 0, 0}; }
inline amrex_loom_uint3 amrex_loom_block_dim() { return {1, 1, 1}; }
inline amrex_loom_uint3 amrex_loom_grid_dim() { return {1, 1, 1}; }
inline void __syncthreads() {}
inline void __threadfence() {}
inline void __threadfence_block() {}
#endif
#define threadIdx (::amrex_loom_thread_idx())
#define blockIdx (::amrex_loom_block_idx())
#define blockDim (::amrex_loom_block_dim())
#define gridDim (::amrex_loom_grid_dim())
inline constexpr int warpSize = 64;
#define __builtin_amdgcn_s_sleep(n) ((void)0)

// HIP device math beyond <cmath>. Only device code (AMREX_IF_ON_DEVICE)
// references these.
#if defined(AMREX_LOOM_DEVICE_PASS)
AMREX_LOOM_DEVICE_INLINE double sinpi(double x) { return __builtin_sin(3.141592653589793238462643 * x); }
AMREX_LOOM_DEVICE_INLINE float sinpif(float x) { return __builtin_sinf(3.14159265358979323846f * x); }
AMREX_LOOM_DEVICE_INLINE double cospi(double x) { return __builtin_cos(3.141592653589793238462643 * x); }
AMREX_LOOM_DEVICE_INLINE float cospif(float x) { return __builtin_cosf(3.14159265358979323846f * x); }
AMREX_LOOM_DEVICE_INLINE void sincospi(double x, double* s, double* c) { *s = sinpi(x); *c = cospi(x); }
AMREX_LOOM_DEVICE_INLINE void sincospif(float x, float* s, float* c) { *s = sinpif(x); *c = cospif(x); }
AMREX_LOOM_DEVICE_INLINE double rsqrt(double x) { return 1.0 / __builtin_sqrt(x); }
AMREX_LOOM_DEVICE_INLINE float rsqrtf(float x) { return 1.0f / __builtin_sqrtf(x); }
AMREX_LOOM_DEVICE_INLINE double powi(double x, int n) {
  double r = 1.0;
  unsigned m = n < 0 ? 0u - static_cast<unsigned>(n) : static_cast<unsigned>(n);
  double b = x;
  while (m) { if (m & 1u) r *= b; b *= b; m >>= 1; }
  return n < 0 ? 1.0 / r : r;
}
AMREX_LOOM_DEVICE_INLINE float powif(float x, int n) {
  return static_cast<float>(powi(static_cast<double>(x), n));
}
AMREX_LOOM_DEVICE_INLINE unsigned long long __umul64hi(unsigned long long a,
                                                       unsigned long long b) {
  unsigned long long a_lo = a & 0xffffffffull, a_hi = a >> 32;
  unsigned long long b_lo = b & 0xffffffffull, b_hi = b >> 32;
  unsigned long long lo_lo = a_lo * b_lo, hi_lo = a_hi * b_lo;
  unsigned long long lo_hi = a_lo * b_hi, hi_hi = a_hi * b_hi;
  unsigned long long cross = (lo_lo >> 32) + (hi_lo & 0xffffffffull) + lo_hi;
  return hi_hi + (hi_lo >> 32) + (cross >> 32);
}
#endif

namespace amrex_loom_detail {
template <class T>
inline constexpr bool atomic_float = std::is_floating_point_v<T>;
}  // namespace amrex_loom_detail

#if defined(AMREX_LOOM_DEVICE_PASS)
#define AMREX_LOOM_RMW(KIND, value, pointer)                                  \
  loom::view::atomic::rmw<loom::atomic::kind::KIND,                           \
                          loom::atomic::ordering::relaxed,                    \
                          loom::atomic::scope::device>(value, pointer)
namespace amrex_loom_detail {
// Float RMW kinds are not admitted by the importer's atomic binding; emulate
// them with a compare-exchange loop on the integer bit pattern.
template <class T>
using atomic_bits = std::conditional_t<sizeof(T) == 8, unsigned long long, unsigned int>;
enum class float_op { add, min, max, exch };
template <float_op Op, class T>
AMREX_LOOM_DEVICE_INLINE T atomic_float_cas(T* address, T value) {
  using U = atomic_bits<T>;
  U* bits = reinterpret_cast<U*>(address);
  U old = *bits;
  U assumed;
  do {
    assumed = old;
    T current = __builtin_bit_cast(T, assumed);
    T desired = Op == float_op::add   ? current + value
                : Op == float_op::min ? (value < current ? value : current)
                : Op == float_op::max ? (current < value ? value : current)
                                      : value;
    old = loom::view::atomic::cmpxchg<loom::atomic::ordering::relaxed,
                                      loom::atomic::ordering::relaxed,
                                      loom::atomic::scope::device>(
        assumed, __builtin_bit_cast(U, desired), bits);
  } while (old != assumed);
  return __builtin_bit_cast(T, old);
}
}  // namespace amrex_loom_detail
template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicAdd(T* address, T value) {
  if constexpr (amrex_loom_detail::atomic_float<T>) {
    return amrex_loom_detail::atomic_float_cas<amrex_loom_detail::float_op::add>(address, value);
  } else {
    return AMREX_LOOM_RMW(addi, value, address);
  }
}
template <class T>
AMREX_LOOM_DEVICE_INLINE void atomicAddNoRet(T* address, T value) {
  (void)atomicAdd(address, value);
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicMin(T* address, T value) {
  if constexpr (amrex_loom_detail::atomic_float<T>) {
    return amrex_loom_detail::atomic_float_cas<amrex_loom_detail::float_op::min>(address, value);
  } else if constexpr (std::is_signed_v<T>) {
    return AMREX_LOOM_RMW(minsi, value, address);
  } else {
    return AMREX_LOOM_RMW(minui, value, address);
  }
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicMax(T* address, T value) {
  if constexpr (amrex_loom_detail::atomic_float<T>) {
    return amrex_loom_detail::atomic_float_cas<amrex_loom_detail::float_op::max>(address, value);
  } else if constexpr (std::is_signed_v<T>) {
    return AMREX_LOOM_RMW(maxsi, value, address);
  } else {
    return AMREX_LOOM_RMW(maxui, value, address);
  }
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicExch(T* address, T value) {
  if constexpr (amrex_loom_detail::atomic_float<T>) {
    return amrex_loom_detail::atomic_float_cas<amrex_loom_detail::float_op::exch>(address, value);
  } else {
    return AMREX_LOOM_RMW(xchgi, value, address);
  }
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicAnd(T* address, T value) {
  return AMREX_LOOM_RMW(andi, value, address);
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicOr(T* address, T value) {
  return AMREX_LOOM_RMW(ori, value, address);
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicCAS(T* address, T compare, T value) {
  return loom::view::atomic::cmpxchg<loom::atomic::ordering::relaxed,
                                     loom::atomic::ordering::relaxed,
                                     loom::atomic::scope::device>(
      compare, value, address);
}
#undef AMREX_LOOM_RMW
template <class T>
AMREX_LOOM_DEVICE_INLINE T __shfl_down(T value, unsigned int delta,
                                       int width = warpSize) {
  auto result = loom::kernel::subgroup::shuffle_down(
      value, delta, static_cast<unsigned>(width));
  return result.valid ? result.value : value;
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T __shfl_xor(T value, int mask, int width = warpSize) {
  auto result = loom::kernel::subgroup::shuffle_xor(
      value, static_cast<unsigned>(mask), static_cast<unsigned>(width));
  return result.valid ? result.value : value;
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T __shfl(T value, int lane, int width = warpSize) {
  auto result = loom::kernel::subgroup::shuffle_index(
      value, static_cast<unsigned>(lane), static_cast<unsigned>(width));
  return result.valid ? result.value : value;
}
AMREX_LOOM_DEVICE_INLINE unsigned long long __ballot(int predicate) {
  return loom::kernel::subgroup::vote::ballot<unsigned long long>(predicate !=
                                                                  0);
}
#else
template <class T> inline T atomicAdd(T* a, T v) { T o = *a; *a = o + v; return o; }
template <class T> inline void atomicAddNoRet(T* a, T v) { *a += v; }
template <class T> inline T atomicMin(T* a, T v) { T o = *a; if (v < o) *a = v; return o; }
template <class T> inline T atomicMax(T* a, T v) { T o = *a; if (v > o) *a = v; return o; }
template <class T> inline T atomicExch(T* a, T v) { T o = *a; *a = v; return o; }
template <class T> inline T atomicAnd(T* a, T v) { T o = *a; *a = o & v; return o; }
template <class T> inline T atomicOr(T* a, T v) { T o = *a; *a = o | v; return o; }
template <class T> inline T atomicCAS(T* a, T c, T v) { T o = *a; if (o == c) *a = v; return o; }
template <class T> inline T __shfl_down(T v, unsigned int, int = warpSize) { return v; }
template <class T> inline T __shfl_xor(T v, int, int = warpSize) { return v; }
template <class T> inline T __shfl(T v, int, int = warpSize) { return v; }
inline unsigned long long __ballot(int p) { return p ? 1ull : 0ull; }
#endif
// HIP's typed atomic builtins as used by AMReX (relaxed, agent scope), routed
// to the relaxed device-scope atomics above. Macros, because clang reserves
// the builtin names; the ordering and scope arguments are not evaluated.
#define __hip_atomic_fetch_add(pointer, value, order, scope) atomicAdd((pointer), (value))
#define __hip_atomic_fetch_min(pointer, value, order, scope) atomicMin((pointer), (value))
#define __hip_atomic_fetch_max(pointer, value, order, scope) atomicMax((pointer), (value))
// HIP's bit reinterpretations between floating-point and integer values.
AMREX_LOOM_DEVICE_INLINE int __float_as_int(float x) { return __builtin_bit_cast(int, x); }
AMREX_LOOM_DEVICE_INLINE float __int_as_float(int x) { return __builtin_bit_cast(float, x); }
AMREX_LOOM_DEVICE_INLINE long long __double_as_longlong(double x) {
  return __builtin_bit_cast(long long, x);
}
AMREX_LOOM_DEVICE_INLINE double __longlong_as_double(long long x) {
  return __builtin_bit_cast(double, x);
}
template <class T>
AMREX_LOOM_DEVICE_INLINE T __shfl_down_sync(unsigned long long, T value,
                                            unsigned int delta,
                                            int width = warpSize) {
  return __shfl_down(value, delta, width);
}

#endif  // AMREX_LOOM_HIP_RUNTIME_H_
