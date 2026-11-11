#ifndef AMREX_LOOM_LAUNCH_H_
#define AMREX_LOOM_LAUNCH_H_
#include <AMReX_Config.H>

// Kernel launch for AMReX_GPU_BACKEND=LOOM.
//
// AMREX_LAUNCH_KERNEL(MT, blocks, threads, shared, stream, f...) calls
// amrex_loom::launch_fixed<MT>. The callables f... are packed into a launch
// tag, amrex_loom::Tag<F...>, run by a fixed-size kernel template instance
// kernel<MT><Tag<F...>>:
//  * in the Loom device pass, launch_fixed names that instance, so the
//    importer emits it with export name kernel<MT>__<Itanium name of the tag
//    type> (loom::kernel_type_symbol);
//  * in the host pass, launch_fixed passes the same export name and the tag
//    bytes to amrex_loom::launch_named (AMReX_Loom.cpp), which packs the
//    kernel arguments from the tag following the kernel's manifest and
//    dispatches on the HRX stream with hrx_stream_dispatch.
// Launches need a compile-time block size MT.

#ifdef AMREX_USE_LOOM

#include <AMReX_GpuQualifiers.H>  // AMREX_LOOM_DEVICE_PASS, AMREX_LOOM_DEVICE_INLINE
#include <AMReX_GpuTypes.H>       // amrex::dim3, hrx_stream_t

#include <cstddef>
#include <type_traits>

#if defined(AMREX_LOOM_DEVICE_PASS)
#include <loomcxx/kernel.h>
#else
#include <cstdio>
#include <cstdlib>
#include <typeinfo>
#endif

namespace amrex_loom {

// Launch tag: the callables run in order by each invocation. Its Itanium
// type name identifies the kernel export on both sides of the build. Built
// as a recursive record: the first callable, then a tag for the rest.
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

#if defined(AMREX_LOOM_DEVICE_PASS)
// Fixed workgroup sizes: blockDim/gridDim queries lower only when the
// workgroup size is static, so every launch names one.
#define AMREX_LOOM_FIXED_KERNEL(N)                                            \
  static loom::kernel::configuration kernel_config##N(                        \
      unsigned gx, unsigned gy, unsigned gz) {                                \
    return {{gx, gy, gz}, {N, 1, 1}};                                         \
  }                                                                           \
  template <class T>                                                          \
  [[loom::kernel(kernel_config##N), loom::kernel_type_symbol]]                \
  void kernel##N(T tag) {                                                     \
    tag();                                                                    \
  }
AMREX_LOOM_FIXED_KERNEL(1)
AMREX_LOOM_FIXED_KERNEL(64)
AMREX_LOOM_FIXED_KERNEL(128)
AMREX_LOOM_FIXED_KERNEL(256)
AMREX_LOOM_FIXED_KERNEL(512)
AMREX_LOOM_FIXED_KERNEL(1024)
#undef AMREX_LOOM_FIXED_KERNEL
#endif

template <int MT>
inline constexpr bool has_fixed_kernel =
    MT == 1 || MT == 64 || MT == 128 || MT == 256 || MT == 512 || MT == 1024;

// Export-name prefix of the fixed-size kernel template for MT.
template <int MT>
inline constexpr const char* fixed_kernel_name =
    MT == 1 ? "kernel1" : MT == 64 ? "kernel64" : MT == 128 ? "kernel128"
    : MT == 256 ? "kernel256" : MT == 512 ? "kernel512" : "kernel1024";

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
                "with MT in {1, 64, 128, 256, 512, 1024}");
}

// Launch with a compile-time workgroup size MT (threads must equal MT).
template <int MT, class B, class T, class... F>
void launch_fixed(B blocks, T threads, std::size_t shared_bytes,
                  hrx_stream_t stream, F const&... fs) {
  static_assert(has_fixed_kernel<MT>,
                "AMReX LOOM backend: kernel block size MT must be one of "
                "1, 64, 128, 256, 512, 1024");
  Tag<F...> tag{fs...};
  amrex::dim3 g = to_dim3(blocks);
#if defined(AMREX_LOOM_DEVICE_PASS)
  (void)threads; (void)shared_bytes; (void)stream;
  if constexpr (MT == 1) loom::kernel::launch<kernel1<Tag<F...>>>(loom::kernel::workload(g.x, g.y, g.z), tag);
  else if constexpr (MT == 64) loom::kernel::launch<kernel64<Tag<F...>>>(loom::kernel::workload(g.x, g.y, g.z), tag);
  else if constexpr (MT == 128) loom::kernel::launch<kernel128<Tag<F...>>>(loom::kernel::workload(g.x, g.y, g.z), tag);
  else if constexpr (MT == 256) loom::kernel::launch<kernel256<Tag<F...>>>(loom::kernel::workload(g.x, g.y, g.z), tag);
  else if constexpr (MT == 512) loom::kernel::launch<kernel512<Tag<F...>>>(loom::kernel::workload(g.x, g.y, g.z), tag);
  else loom::kernel::launch<kernel1024<Tag<F...>>>(loom::kernel::workload(g.x, g.y, g.z), tag);
#else
  amrex::dim3 b = to_dim3(threads);
  if (b.x * b.y * b.z != static_cast<unsigned>(MT) || b.y != 1 || b.z != 1) {
    std::fprintf(stderr, "amrex loom: fixed-size launch MT=%d with block %ux%ux%u\n",
                 MT, b.x, b.y, b.z);
    std::abort();
  }
  const unsigned grid[3] = {g.x, g.y, g.z};
  const unsigned block[3] = {b.x, b.y, b.z};
  launch_named(fixed_kernel_name<MT>, typeid(Tag<F...>).name(), &tag, sizeof(tag), grid,
               block, shared_bytes, stream);
#endif
}

}  // namespace amrex_loom

#endif  // AMREX_USE_LOOM
#endif  // AMREX_LOOM_LAUNCH_H_
