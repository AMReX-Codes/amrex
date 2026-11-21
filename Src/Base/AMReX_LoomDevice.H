#ifndef AMREX_LOOM_DEVICE_H_
#define AMREX_LOOM_DEVICE_H_
#include <AMReX_Config.H>

// Device vocabulary for AMReX_GPU_BACKEND=LOOM, written on Loom's loomcxx API.
//
// One parse of each translation unit feeds both the host compiler and the
// Loom importer, so everything here is valid C++ for both. The Loom
// operations are declarations only (loomcxx attributes them as Loom
// operations); the host never emits a call to one because it only appears in
// device code: kernel bodies, which host code never reaches, and the device
// arms of loom::target::is_device() branches (AMREX_IF_ON_DEVICE).
//
// Provided:
//  * amrex::threadIdx, blockIdx, blockDim, gridDim: the kernel topology
//    under the CUDA names that AMReX's shared CUDA/HIP/LOOM device code uses,
//    read from loom::kernel::workitem::id and loom::kernel::workgroup::{id,
//    size, count}. They live in namespace amrex only; nothing is a macro.
//  * amrex::LoomDevice: workgroup barrier and memory fences, relaxed
//    device-scope atomics (loom::view::atomic), subgroup shuffles and ballot.

#ifdef AMREX_USE_LOOM

#include <loomcxx/atomic.h>
#include <loomcxx/buffer.h>
#include <loomcxx/kernel.h>
#include <loomcxx/scalar.h>
#include <loomcxx/target.h>

#include <bit>
#include <type_traits>

// A helper that uses topology queries or collectives must be inlined into
// its kernel ([[loom::inline]]); host compilers ignore the attribute.
#define AMREX_LOOM_DEVICE_INLINE [[loom::inline]] inline

// Bytes of workgroup memory behind Gpu::SharedMemory, a fixed workgroup array:
// Loom has no dynamically sized workgroup memory. Launches may request at
// most this much shared memory.
#define AMREX_LOOM_DYNAMIC_SHARED_BYTES 4096

namespace amrex {

namespace LoomDevice {

// ---------------------------------------------------------------------------
// Kernel topology
// ---------------------------------------------------------------------------

enum class Topology { workitem_id, workgroup_id, workgroup_size, workgroup_count };

//! One axis of a topology query. It converts to unsigned, like the
//! components of CUDA's uint3 index variables.
template <Topology Q, int Axis>
struct TopologyAxis
{
    AMREX_LOOM_DEVICE_INLINE operator unsigned () const noexcept
    {
        if constexpr (Q == Topology::workitem_id) {
            if constexpr (Axis == 0) { return loom::kernel::workitem::id.x; }
            else if constexpr (Axis == 1) { return loom::kernel::workitem::id.y; }
            else { return loom::kernel::workitem::id.z; }
        } else if constexpr (Q == Topology::workgroup_id) {
            if constexpr (Axis == 0) { return loom::kernel::workgroup::id.x; }
            else if constexpr (Axis == 1) { return loom::kernel::workgroup::id.y; }
            else { return loom::kernel::workgroup::id.z; }
        } else if constexpr (Q == Topology::workgroup_size) {
            if constexpr (Axis == 0) { return loom::kernel::workgroup::size.x; }
            else if constexpr (Axis == 1) { return loom::kernel::workgroup::size.y; }
            else { return loom::kernel::workgroup::size.z; }
        } else {
            if constexpr (Axis == 0) { return loom::kernel::workgroup::count.x; }
            else if constexpr (Axis == 1) { return loom::kernel::workgroup::count.y; }
            else { return loom::kernel::workgroup::count.z; }
        }
    }
};

// Compound assignment of an axis into an arithmetic variable, spelled out
// because the host compiler does not apply the conversion to built-in
// compound assignment.
template <class T, Topology Q, int Axis>
requires std::is_arithmetic_v<T>
AMREX_LOOM_DEVICE_INLINE T& operator+= (T& lhs, TopologyAxis<Q,Axis> rhs) noexcept
{
    lhs += static_cast<unsigned>(rhs);
    return lhs;
}

template <class T, Topology Q, int Axis>
requires std::is_arithmetic_v<T>
AMREX_LOOM_DEVICE_INLINE T& operator-= (T& lhs, TopologyAxis<Q,Axis> rhs) noexcept
{
    lhs -= static_cast<unsigned>(rhs);
    return lhs;
}

template <class T, Topology Q, int Axis>
requires std::is_arithmetic_v<T>
AMREX_LOOM_DEVICE_INLINE T& operator*= (T& lhs, TopologyAxis<Q,Axis> rhs) noexcept
{
    lhs *= static_cast<unsigned>(rhs);
    return lhs;
}

template <class T, Topology Q, int Axis>
requires std::is_arithmetic_v<T>
AMREX_LOOM_DEVICE_INLINE T& operator/= (T& lhs, TopologyAxis<Q,Axis> rhs) noexcept
{
    lhs /= static_cast<unsigned>(rhs);
    return lhs;
}

template <class T, Topology Q, int Axis>
requires std::is_arithmetic_v<T>
AMREX_LOOM_DEVICE_INLINE T& operator%= (T& lhs, TopologyAxis<Q,Axis> rhs) noexcept
{
    lhs %= static_cast<unsigned>(rhs);
    return lhs;
}

template <Topology Q>
struct TopologyQuery
{
    TopologyAxis<Q,0> x;
    TopologyAxis<Q,1> y;
    TopologyAxis<Q,2> z;
};

// ---------------------------------------------------------------------------
// Synchronization
// ---------------------------------------------------------------------------

//! Workgroup barrier (CUDA's __syncthreads).
AMREX_LOOM_DEVICE_INLINE void syncThreads () noexcept
{
    loom::kernel::workgroup::barrier();
}

//! Device-scope sequentially consistent fence (CUDA's __threadfence).
AMREX_LOOM_DEVICE_INLINE void fenceDevice () noexcept
{
    loom::buffer::fence<loom::atomic::ordering::seq_cst, loom::atomic::scope::device>();
}

//! Workgroup-scope fence (CUDA's __threadfence_block).
AMREX_LOOM_DEVICE_INLINE void fenceBlock () noexcept
{
    loom::buffer::fence<loom::atomic::ordering::acq_rel, loom::atomic::scope::workgroup>();
}

// ---------------------------------------------------------------------------
// Atomics: relaxed, device scope, on pointers to global memory.
// ---------------------------------------------------------------------------

template <loom::atomic::kind K, class T>
AMREX_LOOM_DEVICE_INLINE T rmw (T* p, T v) noexcept
{
    return loom::view::atomic::rmw<K, loom::atomic::ordering::relaxed,
                                   loom::atomic::scope::device>(v, p);
}

//! Stores v at *p if *p equals expected; returns the old value.
template <class T>
AMREX_LOOM_DEVICE_INLINE T cas (T* p, T expected, T v) noexcept
{
    return loom::view::atomic::cmpxchg<loom::atomic::ordering::relaxed,
                                       loom::atomic::ordering::relaxed,
                                       loom::atomic::scope::device>(expected, v, p);
}

template <class T>
using AtomicBits = std::conditional_t<sizeof(T) == 8, unsigned long long, unsigned int>;

//! Atomically replaces *p by f(*p, v) with a compare-exchange loop on the bit
//! pattern; returns the old value. For types without a native RMW kind, such
//! as floating point (the importer admits integer RMW kinds only).
template <class T, class F>
AMREX_LOOM_DEVICE_INLINE T casLoop (T* p, T v, F const& f) noexcept
{
    static_assert(sizeof(T) == 4 || sizeof(T) == 8);
    using U = AtomicBits<T>;
    U* const bits = reinterpret_cast<U*>(p);
    U old = *bits;
    U assumed;
    do {
        assumed = old;
        T const desired = f(std::bit_cast<T>(assumed), v);
        old = cas(bits, assumed, std::bit_cast<U>(desired));
    } while (old != assumed);
    return std::bit_cast<T>(old);
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicAdd (T* p, T v) noexcept
{
    if constexpr (std::is_floating_point_v<T>) {
        return casLoop(p, v, [] (T a, T b) { return a + b; });
    } else {
        return rmw<loom::atomic::kind::addi>(p, v);
    }
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicMin (T* p, T v) noexcept
{
    if constexpr (std::is_floating_point_v<T>) {
        return casLoop(p, v, [] (T a, T b) { return b < a ? b : a; });
    } else if constexpr (std::is_signed_v<T>) {
        return rmw<loom::atomic::kind::minsi>(p, v);
    } else {
        return rmw<loom::atomic::kind::minui>(p, v);
    }
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicMax (T* p, T v) noexcept
{
    if constexpr (std::is_floating_point_v<T>) {
        return casLoop(p, v, [] (T a, T b) { return a < b ? b : a; });
    } else if constexpr (std::is_signed_v<T>) {
        return rmw<loom::atomic::kind::maxsi>(p, v);
    } else {
        return rmw<loom::atomic::kind::maxui>(p, v);
    }
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicExch (T* p, T v) noexcept
{
    if constexpr (std::is_integral_v<T>) {
        return rmw<loom::atomic::kind::xchgi>(p, v);
    } else {
        return casLoop(p, v, [] (T, T b) { return b; });
    }
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicAnd (T* p, T v) noexcept
{
    return rmw<loom::atomic::kind::andi>(p, v);
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicOr (T* p, T v) noexcept
{
    return rmw<loom::atomic::kind::ori>(p, v);
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T atomicCAS (T* p, T expected, T v) noexcept
{
    if constexpr (std::is_integral_v<T>) {
        return cas(p, expected, v);
    } else {
        using U = AtomicBits<T>;
        return std::bit_cast<T>(cas(reinterpret_cast<U*>(p), std::bit_cast<U>(expected),
                                    std::bit_cast<U>(v)));
    }
}

// ---------------------------------------------------------------------------
// Subgroup (wavefront) operations
// ---------------------------------------------------------------------------

//! Value from lane (this lane + delta), or this lane's value when that lane
//! is outside the logical width (CUDA's __shfl_down).
template <class T>
AMREX_LOOM_DEVICE_INLINE T shuffleDown (T v, unsigned delta, unsigned width) noexcept
{
    auto r = loom::kernel::subgroup::shuffle_down(v, delta, width);
    return r.valid ? r.value : v;
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T shuffleXor (T v, unsigned mask, unsigned width) noexcept
{
    auto r = loom::kernel::subgroup::shuffle_xor(v, mask, width);
    return r.valid ? r.value : v;
}

template <class T>
AMREX_LOOM_DEVICE_INLINE T shuffle (T v, unsigned lane, unsigned width) noexcept
{
    auto r = loom::kernel::subgroup::shuffle_index(v, lane, width);
    return r.valid ? r.value : v;
}

AMREX_LOOM_DEVICE_INLINE unsigned long long ballot (bool predicate) noexcept
{
    return loom::kernel::subgroup::vote::ballot<unsigned long long>(predicate);
}

} // namespace LoomDevice

// CUDA-named topology for the device code AMReX shares between CUDA, HIP and
// LOOM (namespace amrex only).
inline constexpr LoomDevice::TopologyQuery<LoomDevice::Topology::workitem_id> threadIdx{};
inline constexpr LoomDevice::TopologyQuery<LoomDevice::Topology::workgroup_id> blockIdx{};
inline constexpr LoomDevice::TopologyQuery<LoomDevice::Topology::workgroup_size> blockDim{};
inline constexpr LoomDevice::TopologyQuery<LoomDevice::Topology::workgroup_count> gridDim{};

} // namespace amrex

#endif // AMREX_USE_LOOM
#endif // AMREX_LOOM_DEVICE_H_
