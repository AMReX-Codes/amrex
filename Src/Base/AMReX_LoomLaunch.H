#ifndef AMREX_LOOM_LAUNCH_H_
#define AMREX_LOOM_LAUNCH_H_
#include <AMReX_Config.H>

// Kernel launch for AMReX_GPU_BACKEND=LOOM.
//
// AMREX_LAUNCH_KERNEL(MT, blocks, threads, shared, stream, f...) calls
// amrex_loom::launch_fixed<MT>. The callables f... are packed into a launch
// tag, amrex_loom::Tag<F...>, run by the kernel template instance
// kernel<MT, Tag<F...>>. launch_fixed odr-uses that instance's Loom launch
// descriptor, loom::kernel::kernel_ref<&kernel<MT, Tag<F...>>>, which
// loom-cxx sees from both sides of its one parse:
//  * the descriptors host code references are the Loom device pass's kernel
//    roots, so it compiles the kernel instance; loom-cxx defines the
//    descriptor inside the object, with local linkage (the kernel's export in
//    the object's own code object, and a function packing the kernel
//    arguments in the compiled kernel's layout);
//  * host code passes the descriptor and the tag to amrex_loom::dispatch
//    (AMReX_Loom.cpp), which packs the kernel arguments from the tag with the
//    descriptor's packing function and dispatches on the HRX stream.
// Launches need a compile-time block size MT.

#ifdef AMREX_USE_LOOM

#include <AMReX_GpuQualifiers.H>  // AMREX_LOOM_DEVICE_INLINE
#include <AMReX_GpuTypes.H>       // amrex::dim3, hrx_stream_t
#include <AMReX_Loom.H>           // amrex_loom::dispatch

#include <loomcxx/kernel.h>
#include <loomcxx/kernel_ref.h>

#include <cstddef>
#include <type_traits>

namespace amrex_loom {

// Launch tag: the callables run in order by each invocation, and the kernel's
// only parameter. Built as a recursive record: the first callable, then a tag
// for the rest.
template <class... F>
struct Tag;
template <class F>
struct Tag<F> {
  F f;
  AMREX_LOOM_DEVICE_INLINE void operator()() const { f(); }
};
template <class F, class... Rest>
struct Tag<F, Rest...> {
  F f;
  Tag<Rest...> rest;
  AMREX_LOOM_DEVICE_INLINE void operator()() const { f(); rest(); }
};

// Fixed workgroup sizes: blockDim/gridDim queries lower only when the
// workgroup size is static, so every launch names one. Both passes see the
// kernel definition. The host only names an instance, through kernel_ref
// (whose template argument is not a reference that emits code), so it never
// emits a kernel body or the device code the body reaches; the device pass
// compiles the instances whose descriptors host code references.
// kernel<MT, T>'s launch configuration function, kernel_config<MT>, takes its
// workgroup size from the template argument and its workgroup count from the
// launch; the importer requires it to have internal linkage.
template <int MT>
static loom::kernel::configuration kernel_config(unsigned gx, unsigned gy,
                                                 unsigned gz) {
  return {{gx, gy, gz}, {MT, 1, 1}};
}
template <int MT, class T>
[[loom::kernel(kernel_config<MT>)]]
void kernel(T tag) {
  tag();
}

template <int MT>
inline constexpr bool valid_block_size = MT >= 1 && MT <= 1024;

inline amrex::dim3 to_dim3(amrex::dim3 d) { return d; }
template <class I>
requires std::is_integral_v<I>
inline amrex::dim3 to_dim3(I n) { return amrex::dim3(static_cast<unsigned>(n)); }

template <class... A>
inline constexpr bool dependent_false = false;

// A launch with a run-time block size (AMREX_LAUNCH_KERNEL_NOBOUND, the
// amrex::launch(nblocks, nthreads_per_block, ...) overloads). LOOM kernels
// have a compile-time workgroup size.
template <class... A>
void launch(A const&...) {
  static_assert(dependent_false<A...>,
                "AMReX LOOM backend: kernel launches need a compile-time block "
                "size; use amrex::launch<MT>(...) or AMREX_LAUNCH_KERNEL(MT, ...) "
                "with 1 <= MT <= 1024");
}

// Launch with a compile-time workgroup size MT (threads must equal MT).
template <int MT, class B, class T, class... F>
void launch_fixed(B blocks, T threads, std::size_t shared_bytes,
                  hrx_stream_t stream, F const&... fs) {
  static_assert(valid_block_size<MT>,
                "AMReX LOOM backend: kernel block size MT must be in [1, 1024]");
  Tag<F...> tag{fs...};
  const amrex::dim3 g = to_dim3(blocks);
  const amrex::dim3 b = to_dim3(threads);
  const unsigned grid[3] = {g.x, g.y, g.z};
  const unsigned block[3] = {b.x, b.y, b.z};
  dispatch(loom::kernel::kernel_ref<&kernel<MT, Tag<F...>>>, &tag, grid, block, MT,
           shared_bytes, stream);
}

}  // namespace amrex_loom

#endif  // AMREX_USE_LOOM
#endif  // AMREX_LOOM_LAUNCH_H_
