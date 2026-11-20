// Host runtime for AMReX_GPU_BACKEND=LOOM.
//
// loom-cxx, the C++ compiler of LOOM builds, embeds each translation unit's
// Loom code object. A launch names its kernel by a Loom kernel descriptor
// (loom::kernel::kernel_ref, defined locally in each object by loom-cxx),
// loads the code object the descriptor names the first time any of its
// kernels launches, dispatches the export the descriptor numbers there, and
// packs the kernel arguments from the host launch tag with the descriptor's
// packing function, which loom-cxx generates from the compiled kernel's
// argument layout. Nothing is loaded at startup. This file implements the
// amrex_loom runtime functions declared in AMReX_Loom.H and AMReX_LoomRandom.H
// on libhrx.
// Single device; streams are HRX streams, each bound to its own HRX queue
// where the device allows.
#include <AMReX_Config.H>

#if defined(AMREX_USE_LOOM)

#include <AMReX_GpuQualifiers.H>  // AMREX_LOOM_DYNAMIC_SHARED_BYTES
#include <AMReX_Loom.H>
#include <AMReX_LoomLaunch.H>
#include <AMReX_LoomRandom.H>
#include <hrx_runtime.h>
#include <loomcxx/kernel_ref.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#if defined(__linux__)
#include <dirent.h>
#include <sched.h>
#endif

namespace {

// --- launch checks -------------------------------------------------------------

struct Dim3 {
  unsigned x, y, z;
};

// Checks a launch against its kernel and returns false for an empty grid.
bool check_launch(const loom::kernel::descriptor& kernel, Dim3 grid, Dim3 block,
                  unsigned fixed_block, std::size_t shared_bytes) {
  // Every descriptor version keeps the version and the export name first.
  if (kernel.version != loom::kernel::descriptor_version) {
    std::fprintf(stderr, "amrex loom: kernel %s has a version %u descriptor; this runtime "
                 "reads version %u only: rebuild with the loom-cxx matching its "
                 "loomcxx/kernel_ref.h\n",
                 kernel.name, kernel.version, loom::kernel::descriptor_version);
    std::abort();
  }
  if (block.x != fixed_block || block.y != 1 || block.z != 1) {
    std::fprintf(stderr, "amrex loom: fixed-size launch MT=%u with block %ux%ux%u\n",
                 fixed_block, block.x, block.y, block.z);
    std::abort();
  }
  if (grid.x == 0 || grid.y == 0 || grid.z == 0) return false;
  if (!(kernel.flags & loom::kernel::descriptor_flag_available)) {
    std::fprintf(stderr, "amrex loom: no device kernel for %s\n", kernel.name);
    std::abort();
  }
  if (kernel.parameter_count != 1) {
    std::fprintf(stderr, "amrex loom: kernel %s has %u parameters, expected one\n",
                 kernel.name, kernel.parameter_count);
    std::abort();
  }
  if (shared_bytes > AMREX_LOOM_DYNAMIC_SHARED_BYTES) {
    std::fprintf(stderr, "amrex loom: %zu bytes of dynamic shared memory requested, "
                 "%d supported\n", shared_bytes, AMREX_LOOM_DYNAMIC_SHARED_BYTES);
    std::abort();
  }
  return true;
}

bool tracing() {
  static const bool trace = std::getenv("AMREX_LOOM_TRACE") != nullptr;
  return trace;
}

void trace_launch(const char* name, Dim3 grid, Dim3 block, std::size_t kernarg_size) {
  if (tracing()) {
    std::fprintf(stderr, "amrex loom: launch %s grid %ux%ux%u block %ux%ux%u kernarg %zu\n",
                 name, grid.x, grid.y, grid.z, block.x, block.y, block.z,
                 kernarg_size);
  }
}

bool sync_after_launch() {
  static const bool sync = std::getenv("AMREX_LOOM_SYNC_LAUNCH") != nullptr;
  return sync;
}

// How many operations a stream records before submitting them as one
// recording (AMREX_LOOM_FLUSH_EVERY, default 32): kernel launches in batch
// mode, and copies and fills in either mode. Completion points (stream and
// device synchronization, stream queries, frees and stream destruction) flush
// every stream first (flush_all), whatever the setting. 0 submits only at
// completion points; 1 submits each operation on its own.
unsigned flush_every() {
  static const unsigned every = [] {
    const char* value = std::getenv("AMREX_LOOM_FLUSH_EVERY");
    return value && *value ? static_cast<unsigned>(std::strtoul(value, nullptr, 10)) : 32u;
  }();
  return every;
}

// How kernel launches reach the GPU (AMREX_LOOM_SUBMIT):
//   immediate (default)  each launch is its own hrx_queue_dispatch on the
//                        stream's queue, chained on the stream's timeline. With
//                        no extra barriers it takes HRX's agent-scope argument
//                        acquisition, as HIP's dispatch packets do, and starts
//                        as soon as it is launched;
//   batch                launches are recorded on the stream and submitted as
//                        one recording per AMREX_LOOM_FLUSH_EVERY operations.
// Both pay a system-scope acquire and release at each completion point. On
// CNS (MI250X), immediate gave 1.22x/1.14x HIP's step time at 64^3/128^3 and
// batch 1.38x/1.20x.
enum class Submit { immediate, batch };
Submit submit_mode() {
  static const Submit mode = [] {
    const char* value = std::getenv("AMREX_LOOM_SUBMIT");
    if (!value || !*value || std::strcmp(value, "immediate") == 0) return Submit::immediate;
    if (std::strcmp(value, "batch") == 0) return Submit::batch;
    std::fprintf(stderr, "amrex loom: AMREX_LOOM_SUBMIT=%s is not immediate or batch\n", value);
    std::abort();
  }();
  return mode;
}

}  // namespace

// --- runtime state on libhrx ---------------------------------------------------

namespace {

struct Allocation {
  hrx_buffer_t buffer;
  std::size_t size;
};

std::recursive_mutex g_mutex;
bool g_initialized = false;
hrx_runtime_t g_runtime = nullptr;
hrx_device_group_t g_group = nullptr;
hrx_device_t g_device = nullptr;  // borrowed from g_group
// Memory scopes for address queries: kernels on g_queue_family and the host.
hrx_memory_scope_t g_program_scope{};
hrx_memory_scope_t g_host_scope{};
// Streams submit to exact HRX queues of one family that supports dispatch,
// transfer and host calls. The default stream uses provisioned queue 0. Each
// created stream gets a queue of its own where the device allows (an unused
// provisioned queue, else a newly acquired one), so kernels on different
// streams can run concurrently. HRX does not multiplex virtual streams onto
// queues, so when acquisition fails the stream shares the least-used
// provisioned queue: it keeps its own timeline, but its work serializes with
// the other streams on that queue. AMREX_LOOM_SHARED_QUEUE=1 forces sharing.
struct StreamQueue {
  hrx_queue_t queue = nullptr;
  bool acquired = false;           // owning reference from hrx_queue_acquire
  uint32_t provisioned = 0;        // provisioned ordinal when !acquired
};
uint32_t g_queue_family = 0;
bool g_queue_family_acquirable = false;
bool g_shared_queue = false;
std::vector<uint32_t> g_provisioned_users;  // streams bound to each provisioned queue
hrx_stream_t g_default_stream = nullptr;
std::map<hrx_stream_t, StreamQueue> g_streams;
std::map<uintptr_t, Allocation> g_allocations;
// The executable of each code object loaded so far, by its loom-cxx record.
std::unordered_map<const loom::kernel::embedded_code_object*, hrx_executable_t> g_executables;

[[noreturn]] void fatal(const char* what, hrx_status_t status) {
  char* message = nullptr;
  std::size_t length = 0;
  if (status) hrx_status_to_string(status, &message, &length);
  std::fprintf(stderr, "amrex loom: %s failed: %s\n", what,
               message ? message : "(no status)");
  std::abort();
}

void check(hrx_status_t status, const char* what) {
  if (!hrx_status_is_ok(status)) fatal(what, status);
}

void shutdown();

// NUMA node of the CPU the calling thread runs on, or -1 if unknown. HRX pins
// its device threads to node 0 unless the device URI names a host node, and
// fails device creation when the process may not run there (a batch job whose
// CPUs are all on another node).
int current_numa_node() {
#if defined(__linux__)
  const int cpu = sched_getcpu();
  if (cpu < 0) return -1;
  char path[64];
  std::snprintf(path, sizeof(path), "/sys/devices/system/cpu/cpu%d", cpu);
  DIR* dir = opendir(path);
  if (!dir) return -1;
  int node = -1;
  while (const dirent* entry = readdir(dir)) {
    char tail = 0;
    if (std::sscanf(entry->d_name, "node%d%c", &node, &tail) == 1) break;
    node = -1;
  }
  closedir(dir);
  return node;
#else
  return -1;
#endif
}

void initialize() {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (g_initialized) return;
  // One sealed group holding the first amdgpu device.
  const hrx_host_allocator_t allocator = hrx_host_allocator_system();
  hrx_runtime_options_t options{};
  check(hrx_runtime_create(&options, allocator, &g_runtime), "hrx_runtime_create");
  std::size_t device_count = 0;
  hrx_device_info_t* infos = nullptr;
  check(hrx_runtime_query_available_devices(g_runtime, "amdgpu", allocator, &device_count,
                                            &infos),
        "hrx_runtime_query_available_devices");
  if (device_count == 0) {
    std::fprintf(stderr, "amrex loom: no amdgpu devices available\n");
    std::abort();
  }
  // Place HRX's device threads on this thread's NUMA node, which the process
  // is allowed to use.
  std::string uri = infos[0].uri;
  if (const int node = current_numa_node(); node >= 0) {
    uri += (uri.find('?') == std::string::npos ? '?' : '&');
    uri += "host_numa_node=" + std::to_string(node);
  }
  const hrx_device_create_params_t params = {uri.c_str(), HRX_WAIT_FLAG_YIELD};
  hrx_status_t created = hrx_device_group_create(g_runtime, 1, &params, &g_group);
  hrx_host_allocator_free(allocator, infos);
  check(created, "hrx_device_group_create");
  g_device = hrx_device_group_get_device(g_group, 0);
  constexpr hrx_queue_role_flags_t roles =
      HRX_QUEUE_ROLE_DISPATCH | HRX_QUEUE_ROLE_TRANSFER | HRX_QUEUE_ROLE_HOST_CALL;
  bool found = false;
  const auto family_count = static_cast<uint32_t>(hrx_device_queue_family_count(g_device));
  for (uint32_t family = 0; family < family_count && !found; ++family) {
    hrx_queue_family_info_t info{};
    check(hrx_device_query_queue_family(g_device, family, &info),
          "hrx_device_query_queue_family");
    if ((info.roles & roles) != roles || info.provisioned_queue_count == 0) continue;
    found = true;
    g_queue_family = family;
    g_queue_family_acquirable = info.supports_acquisition;
    g_provisioned_users.assign(info.provisioned_queue_count, 0);
  }
  if (!found) {
    std::fprintf(stderr, "amrex loom: no HRX queue family supports dispatch, transfer and "
                         "host calls\n");
    std::abort();
  }
  check(hrx_device_query_queue_family_memory_scope(g_device, g_queue_family,
                                                   HRX_MEMORY_SITE_PROGRAM, &g_program_scope),
        "hrx_device_query_queue_family_memory_scope");
  check(hrx_device_query_host_memory_scope(g_device, &g_host_scope),
        "hrx_device_query_host_memory_scope");
  const char* shared = std::getenv("AMREX_LOOM_SHARED_QUEUE");
  g_shared_queue = shared && *shared && std::strcmp(shared, "0") != 0;
  check(hrx_stream_create(hrx_device_get_queue(g_device, g_queue_family, 0), 0,
                          &g_default_stream),
        "hrx_stream_create");
  g_provisioned_users[0] = 1;
  g_initialized = true;
  // Registered after the runtime's globals are constructed, so it runs before
  // they are destroyed.
  std::atexit(shutdown);
}

// Address kernels on the dispatch family use for |buffer|.
void* device_address(hrx_buffer_t buffer) {
  uint64_t address = 0;
  check(hrx_buffer_query_address(buffer, g_program_scope, HRX_BUFFER_INTERFACE_DEVICE_ADDRESS,
                                 &address),
        "hrx_buffer_query_address");
  return reinterpret_cast<void*>(static_cast<uintptr_t>(address));
}

// A null stream is the default stream.
hrx_stream_t resolve(hrx_stream_t stream) {
  initialize();
  return stream ? stream : g_default_stream;
}

// Operations recorded on each stream since its last flush when batching
// (AMREX_LOOM_FLUSH_EVERY != 1). Includes external streams (Gpu::setExternalGpuStream),
// which are not in g_streams.
std::unordered_map<hrx_stream_t, unsigned> g_unflushed;
// Streams with immediate dispatches since their last completion point.
std::unordered_map<hrx_stream_t, bool> g_unreleased;

// Submits every stream's pending recording. Completion points call this before
// waiting on or querying any one stream: HRX's query sees only submitted work,
// and AMReX synchronizes its streams one at a time (MFIter::Finalize), so
// submitting all of them first keeps them running concurrently. Flushing a
// stream with nothing pending does nothing.
void flush_all() {
  initialize();
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  // Immediate dispatches carry no system-scope barriers. A stream recording
  // begins with a system-scope acquire and ends with a system-scope release
  // (hrx_stream_flush), so an empty one after them publishes their writes,
  // including those cached in L2 for pinned host memory, before the stream's
  // timeline point, and discards stale L2 lines before the host's next access.
  for (auto& [s, pending] : g_unreleased) {
    if (pending) check(hrx_stream_barrier(s, nullptr), "hrx_stream_barrier");
  }
  g_unreleased.clear();
  check(hrx_stream_flush(g_default_stream), "hrx_stream_flush");
  for (auto& [s, q] : g_streams) check(hrx_stream_flush(s), "hrx_stream_flush");
  for (auto& [s, count] : g_unflushed) {
    if (count) check(hrx_stream_flush(s), "hrx_stream_flush");
  }
  g_unflushed.clear();
}

// Called after recording an operation on |s|: submits per AMREX_LOOM_FLUSH_EVERY.
void recorded(hrx_stream_t s) {
  const unsigned every = flush_every();
  if (every == 1) {
    check(hrx_stream_flush(s), "hrx_stream_flush");
    return;
  }
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  unsigned& count = g_unflushed[s];
  if (++count >= every && every != 0) {
    check(hrx_stream_flush(s), "hrx_stream_flush");
    count = 0;
  }
}

// Host copies of pageable sources of host-to-device transfers still recorded or
// in flight on each stream (memcpy_async), freed once the stream is known to
// have reached them.
struct Staging {
  std::vector<std::unique_ptr<uint8_t[]>> buffers;
  std::size_t bytes = 0;
};
std::unordered_map<hrx_stream_t, Staging> g_staging;

// Pinned staging for small pageable host-to-device copies under immediate
// submission (copy_immediate), one per stream, bump-allocated. A range is
// reused only after the stream is known complete: its completion point has
// then run the system-scope barrier that discards stale L2 lines.
struct PinnedRing {
  uint8_t* base = nullptr;
  std::size_t used = 0;
};
constexpr std::size_t kPinnedRingBytes = std::size_t(8) << 20;
// Largest pageable source staged through the ring; larger ones use
// hrx_stream_copy_h2d_linear.
constexpr std::size_t kPinnedRingCopyLimit = std::size_t(1) << 20;
std::unordered_map<hrx_stream_t, PinnedRing> g_rings;

// Frees |s|'s staging copies; the caller has waited for |s|.
void release_staging(hrx_stream_t s) {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  g_staging.erase(s);
  auto it = g_rings.find(s);
  if (it != g_rings.end()) it->second.used = 0;
}

// Submits one kernel directly to |s|'s queue, ordered after everything already
// on |s| and before everything later, by waiting on the stream's current
// timeline point and signaling the next one (AMREX_LOOM_SUBMIT=immediate).
void dispatch_immediate(hrx_stream_t s, hrx_executable_t executable, uint32_t ordinal,
                        const hrx_dispatch_config_t& config, const void* kernarg,
                        std::size_t kernarg_size) {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  // Copies or fills still recorded on |s| come first.
  check(hrx_stream_flush(s), "hrx_stream_flush");
  g_unflushed.erase(s);
  hrx_timeline_point_t after{};
  check(hrx_stream_get_timeline_position(s, &after), "hrx_stream_get_timeline_position");
  hrx_semaphore_t timeline = nullptr;
  check(hrx_stream_get_semaphore(s, &timeline), "hrx_stream_get_semaphore");
  uint64_t done = 0;
  check(hrx_stream_advance_timeline(s, &done), "hrx_stream_advance_timeline");
  hrx_semaphore_t wait_semaphore = after.semaphore;
  uint64_t wait_value = after.value;
  const hrx_semaphore_list_t waits{&wait_semaphore, &wait_value, after.value ? 1u : 0u};
  const hrx_semaphore_list_t signals{&timeline, &done, 1};

  // Empty barrier lists request no extra cache operations, which keeps the
  // dispatch eligible for HRX's agent-scope acquisition of its arguments.
  // Ordering on the queue comes from the timeline wait; host visibility comes
  // from the next completion point (flush_all).
  static const hrx_barrier_list_t no_barriers{0, nullptr};
  static const hrx_queue_barriers_t barriers{&no_barriers, &no_barriers};
  g_unreleased[s] = true;
  check(hrx_queue_dispatch(hrx_stream_get_queue(s), &waits, &signals, executable, ordinal,
                           &config, kernarg, kernarg_size, nullptr, 0, &barriers,
                           HRX_DISPATCH_FLAG_CUSTOM_DIRECT_ARGUMENTS),
        "hrx_queue_dispatch");
}

// Submits one copy directly to |s|'s queue, ordered like dispatch_immediate.
// hrx_stream_copy_h2d_linear/d2h_linear would make HRX wait on the host for earlier work
// before reading or writing host memory, and a stream recording carries
// system-scope barriers; either leaves the GPU idle for tens of microseconds.
// Both ends here are HRX buffers (device or pinned), so the copy needs neither:
// host visibility comes from the next completion point (flush_all).
void copy_immediate(hrx_stream_t s, hrx_buffer_t src, std::size_t src_offset, hrx_buffer_t dst,
                    std::size_t dst_offset, std::size_t bytes) {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  // Copies or fills still recorded on |s| come first.
  check(hrx_stream_flush(s), "hrx_stream_flush");
  g_unflushed.erase(s);
  hrx_timeline_point_t after{};
  check(hrx_stream_get_timeline_position(s, &after), "hrx_stream_get_timeline_position");
  hrx_semaphore_t timeline = nullptr;
  check(hrx_stream_get_semaphore(s, &timeline), "hrx_stream_get_semaphore");
  uint64_t done = 0;
  check(hrx_stream_advance_timeline(s, &done), "hrx_stream_advance_timeline");
  hrx_semaphore_t wait_semaphore = after.semaphore;
  uint64_t wait_value = after.value;
  const hrx_semaphore_list_t waits{&wait_semaphore, &wait_value, after.value ? 1u : 0u};
  const hrx_semaphore_list_t signals{&timeline, &done, 1};
  static const hrx_barrier_list_t no_barriers{0, nullptr};
  static const hrx_queue_barriers_t barriers{&no_barriers, &no_barriers};
  g_unreleased[s] = true;
  check(hrx_queue_copy_linear(hrx_stream_get_queue(s), &waits, &signals, src, src_offset, dst,
                       dst_offset, bytes, &barriers),
        "hrx_queue_copy_linear");
}

// How completion points wait for a stream (AMREX_LOOM_SYNC_WAIT):
//   poll (default)  spin on hrx_semaphore_query until the stream's timeline
//                   point is reached, then hrx_stream_synchronize, which then
//                   returns at once. HIP's stream synchronization spins too;
//   block           hrx_stream_synchronize alone. It sleeps in
//                   hrx_semaphore_wait, whose wake-up costs 15-25 us more than
//                   HIP's synchronization after short waits and about 120 us
//                   more after waits of a few milliseconds (MI250X).
// On CNS (MI250X), poll is 0.1-1.0 ms per step faster than block.
bool sync_poll() {
  static const bool poll = [] {
    const char* value = std::getenv("AMREX_LOOM_SYNC_WAIT");
    if (!value || !*value || std::strcmp(value, "poll") == 0) return true;
    if (std::strcmp(value, "block") == 0) return false;
    std::fprintf(stderr, "amrex loom: AMREX_LOOM_SYNC_WAIT=%s is not poll or block\n", value);
    std::abort();
  }();
  return poll;
}

void wait_stream(hrx_stream_t s) {
  if (sync_poll()) {
    hrx_timeline_point_t point{};
    check(hrx_stream_get_timeline_position(s, &point), "hrx_stream_get_timeline_position");
    if (point.semaphore && point.value) {
      uint64_t value = 0;
      do {
        check(hrx_semaphore_query(point.semaphore, &value), "hrx_semaphore_query");
      } while (value < point.value);
    }
  }
  check(hrx_stream_synchronize(s), "hrx_stream_synchronize");
}

void sync_all() {
  flush_all();
  wait_stream(g_default_stream);
  for (auto& [s, q] : g_streams) wait_stream(s);
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  // External streams' staging is released when they are synchronized.
  g_staging.erase(g_default_stream);
  for (auto& [s, q] : g_streams) g_staging.erase(s);
  for (auto& [s, ring] : g_rings) ring.used = 0;
}

// Returns the allocation containing |address|, or nullptr for host memory.
const std::pair<const uintptr_t, Allocation>* find_allocation(const void* address) {
  auto a = reinterpret_cast<uintptr_t>(address);
  auto it = g_allocations.upper_bound(a);
  if (it == g_allocations.begin()) return nullptr;
  --it;
  if (a >= it->first + it->second.size && !(a == it->first && it->second.size == 0)) {
    return nullptr;
  }
  return &*it;
}

}  // namespace

// --- launch (amrex_loom) ------------------------------------------------------

namespace {

// The executable of the code object |kernel| names, loaded the first time
// one of its kernels launches. Each translation unit's descriptors name the
// code object of that unit, which carries the kernels it launches: kernels of
// different units may share an export name (lambdas with internal linkage,
// for example) and only the unit's own code object has the right one. The
// descriptor names the kernel by its export ordinal there.
hrx_executable_t kernel_executable(const loom::kernel::descriptor& kernel) {
  const loom::kernel::embedded_code_object* code_object = kernel.code_object;
  if (!code_object) {
    std::fprintf(stderr, "amrex loom: kernel descriptor of %s names no code object\n",
                 kernel.name);
    std::abort();
  }
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto found = g_executables.find(code_object);
  if (found != g_executables.end()) return found->second;
  initialize();
  hrx_executable_t executable = nullptr;
  check(hrx_executable_load_data(g_device, code_object->data, code_object->size,
                                 code_object->target_family, code_object->target_key,
                                 &executable),
        "hrx_executable_load_data");
  if (tracing()) {
    std::fprintf(stderr, "amrex loom: loaded code object %p (%s:%s, %zu bytes) for %s\n",
                 static_cast<const void*>(code_object), code_object->target_family,
                 code_object->target_key, static_cast<std::size_t>(code_object->size),
                 kernel.name);
  }
  g_executables.emplace(code_object, executable);
  return executable;
}

// Checks, while tracing, that the export the descriptor numbers is its kernel.
void check_export(hrx_executable_t executable, const loom::kernel::descriptor& kernel) {
  hrx_executable_export_info_t info{};
  check(hrx_executable_export_info(executable, kernel.export_ordinal, &info),
        "hrx_executable_export_info");
  if (std::strcmp(info.name, kernel.name) != 0) {
    std::fprintf(stderr, "amrex loom: export %u of the code object of %s is %s\n",
                 kernel.export_ordinal, kernel.name, info.name);
    std::abort();
  }
}

}  // namespace

void amrex_loom::dispatch(const loom::kernel::descriptor& descriptor, const void* parameter,
                          const unsigned grid_size[3], const unsigned block_size[3],
                          unsigned fixed_block, std::size_t shared_bytes,
                          hrx_stream_t stream) {
  const Dim3 grid{grid_size[0], grid_size[1], grid_size[2]};
  const Dim3 block{block_size[0], block_size[1], block_size[2]};
  if (!check_launch(descriptor, grid, block, fixed_block, shared_bytes)) return;
  hrx_executable_t executable = kernel_executable(descriptor);
  if (tracing()) check_export(executable, descriptor);
  // The descriptor's packing function writes the whole kernarg segment from
  // the launch tag; HRX copies it at dispatch.
  alignas(16) uint8_t local[1024];
  std::vector<uint8_t> heap;
  uint8_t* kernarg = local;
  if (descriptor.kernarg_size > sizeof(local)) {
    heap.resize(descriptor.kernarg_size);
    kernarg = heap.data();
  }
  descriptor.pack(kernarg, &parameter);
  trace_launch(descriptor.name, grid, block, descriptor.kernarg_size);
  hrx_dispatch_config_t config{};
  config.workgroup_count[0] = grid.x;
  config.workgroup_count[1] = grid.y;
  config.workgroup_count[2] = grid.z;
  config.workgroup_size[0] = block.x;
  config.workgroup_size[1] = block.y;
  config.workgroup_size[2] = block.z;
  hrx_stream_t s = resolve(stream);
  if (submit_mode() == Submit::immediate) {
    dispatch_immediate(s, executable, descriptor.export_ordinal, config, kernarg,
                       descriptor.kernarg_size);
    if (sync_after_launch()) {
      flush_all();
      wait_stream(s);
    }
    return;
  }
  check(hrx_stream_dispatch(s, executable, descriptor.export_ordinal, &config,
                            kernarg, descriptor.kernarg_size, nullptr, 0,
                            HRX_DISPATCH_FLAG_CUSTOM_DIRECT_ARGUMENTS),
        "hrx_stream_dispatch");
  if (sync_after_launch()) {
    flush_all();
    wait_stream(s);
  } else {
    recorded(s);
  }
}

// --- device (amrex_loom) -------------------------------------------------------

namespace {
// Picks the queue for a new stream: an unused provisioned queue, else a newly
// acquired one, else the least-used provisioned queue (shared).
StreamQueue select_stream_queue() {
  auto least = std::min_element(g_provisioned_users.begin(), g_provisioned_users.end());
  auto ordinal = static_cast<uint32_t>(least - g_provisioned_users.begin());
  if (!g_shared_queue && *least != 0 && g_queue_family_acquirable) {
    hrx_queue_t queue = nullptr;
    hrx_status_t status = hrx_queue_acquire(g_device, g_queue_family, nullptr, &queue);
    if (hrx_status_is_ok(status)) return {queue, true, 0};
    hrx_status_ignore(status);
  }
  ++g_provisioned_users[ordinal];
  return {hrx_device_get_queue(g_device, g_queue_family, ordinal), false, ordinal};
}
// Drops this runtime's use of |q| once its stream has been released.
void release_stream_queue(const StreamQueue& q) {
  if (q.acquired) {
    hrx_queue_release(q.queue);
  } else {
    --g_provisioned_users[q.provisioned];
  }
}
}  // namespace

namespace amrex_loom {

int device_count() {
  initialize();
  return 1;
}

namespace {

// A 32-bit device property from HRX, or |fallback| when HRX does not report it.
int device_property(hrx_device_property_t property, int fallback) {
  uint32_t value = 0;
  hrx_status_t status = hrx_device_get_property(g_device, property, &value, sizeof(value));
  if (!hrx_status_is_ok(status)) {
    hrx_status_ignore(status);
    return fallback;
  }
  return value ? static_cast<int>(value) : fallback;
}

}  // namespace

// HRX reports memory, compute units and the maximum workgroup size. It does
// not report the wavefront size or shared memory, so those, and the limits
// with no HRX query, describe an MI200-class GCD.
void device_properties(DeviceProperties* prop) {
  initialize();
  std::memset(prop, 0, sizeof(*prop));
  std::snprintf(prop->name, sizeof(prop->name), "AMD GPU (HRX/Loom)");
  std::snprintf(prop->gcnArchName, sizeof(prop->gcnArchName), "%s",
                std::getenv("AMREX_LOOM_TARGET") ? std::getenv("AMREX_LOOM_TARGET") : "gfx90a");
  std::size_t free_bytes = 0, total_bytes = 0;
  memory_info(&free_bytes, &total_bytes);
  prop->totalGlobalMem = total_bytes;
  prop->sharedMemPerBlock = 65536;
  prop->warpSize = 64;
  prop->maxThreadsPerBlock = device_property(HRX_DEVICE_PROPERTY_MAX_WORKGROUP_SIZE, 1024);
  prop->maxThreadsDim[0] = prop->maxThreadsDim[1] = prop->maxThreadsDim[2] =
      prop->maxThreadsPerBlock;
  prop->maxGridSize[0] = 2147483647;
  prop->maxGridSize[1] = prop->maxGridSize[2] = 65535;
  prop->multiProcessorCount = device_property(HRX_DEVICE_PROPERTY_COMPUTE_UNITS, 104);
  prop->maxThreadsPerMultiProcessor = 2048;
  prop->managedMemory = 0;
  prop->concurrentManagedAccess = 0;
  prop->integrated = 0;
}

void device_uuid(char bytes[16]) {
  initialize();
  std::memset(bytes, 0, 16);
}

void memory_info(std::size_t* free_bytes, std::size_t* total_bytes) {
  initialize();
  std::size_t f = 0, t = 0;
  if (!hrx_status_is_ok(hrx_device_memory_info(g_device, &f, &t))) {
    f = t = std::size_t(64) << 30;
  }
  *free_bytes = f;
  *total_bytes = t;
}

void device_synchronize() { sync_all(); }

std::size_t loaded_code_object_count() {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  return g_executables.size();
}

// --- streams (amrex_loom) ------------------------------------------------------

hrx_stream_t stream_create() {
  initialize();
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  StreamQueue q = select_stream_queue();
  hrx_stream_t s = nullptr;
  hrx_status_t status = hrx_stream_create(q.queue, 0, &s);
  if (!hrx_status_is_ok(status)) {
    release_stream_queue(q);
    fatal("hrx_stream_create", status);
  }
  g_streams[s] = q;
  return s;
}

void stream_destroy(hrx_stream_t stream) {
  if (!stream) return;
  flush_all();
  wait_stream(stream);
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto it = g_streams.find(stream);
  if (it == g_streams.end()) {
    std::fprintf(stderr, "amrex loom: stream_destroy of a stream not from stream_create\n");
    std::abort();
  }
  StreamQueue q = it->second;
  g_streams.erase(it);
  g_unflushed.erase(stream);
  g_staging.erase(stream);
  g_unreleased.erase(stream);
  if (auto ring = g_rings.find(stream); ring != g_rings.end()) {
    pinned_free(ring->second.base);
    g_rings.erase(ring);
  }
  hrx_stream_release(stream);
  release_stream_queue(q);
}

// A cross-stream wait is a barrier on |waiter|'s queue that waits for both
// streams' timeline points and advances |waiter|'s timeline, so the GPU orders
// the streams and the host does not wait. Visibility (AMREX_LOOM_JOIN_VISIBILITY):
//   agent (default)  no extra cache operations: both queues run on one GPU,
//                    whose L2 holds the writes, as for consecutive dispatches;
//   system           HRX's conservative system acquire and release.
bool join_system_visibility() {
  static const bool system = [] {
    const char* v = std::getenv("AMREX_LOOM_JOIN_VISIBILITY");
    if (!v || !*v || std::strcmp(v, "agent") == 0) return false;
    if (std::strcmp(v, "system") == 0) return true;
    std::fprintf(stderr, "amrex loom: AMREX_LOOM_JOIN_VISIBILITY=%s is not agent or system\n", v);
    std::abort();
  }();
  return system;
}

void stream_wait_streams(hrx_stream_t waiter, const hrx_stream_t* waited, int n) {
  hrx_stream_t w = resolve(waiter);
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  std::vector<hrx_semaphore_t> wait_semaphores;
  std::vector<uint64_t> wait_values;
  bool unreleased = false;
  for (int i = 0; i < n; ++i) {
    hrx_stream_t v = resolve(waited[i]);
    if (v == w) continue;
    check(hrx_stream_flush(v), "hrx_stream_flush");
    g_unflushed.erase(v);
    hrx_timeline_point_t vp{};
    check(hrx_stream_get_timeline_position(v, &vp), "hrx_stream_get_timeline_position");
    if (!vp.semaphore || !vp.value) continue;  // nothing has been queued on it
    wait_semaphores.push_back(vp.semaphore);
    wait_values.push_back(vp.value);
    if (g_unreleased.count(v) && g_unreleased[v]) unreleased = true;
  }
  if (wait_semaphores.empty()) return;
  check(hrx_stream_flush(w), "hrx_stream_flush");
  g_unflushed.erase(w);
  hrx_timeline_point_t wp{};
  check(hrx_stream_get_timeline_position(w, &wp), "hrx_stream_get_timeline_position");
  if (wp.semaphore && wp.value) {
    wait_semaphores.push_back(wp.semaphore);
    wait_values.push_back(wp.value);
  }
  hrx_semaphore_t timeline = nullptr;
  check(hrx_stream_get_semaphore(w, &timeline), "hrx_stream_get_semaphore");
  uint64_t done = 0;
  check(hrx_stream_advance_timeline(w, &done), "hrx_stream_advance_timeline");
  const hrx_semaphore_list_t waits{wait_semaphores.data(), wait_values.data(),
                                   static_cast<uint32_t>(wait_semaphores.size())};
  const hrx_semaphore_list_t signals{&timeline, &done, 1};
  static const hrx_barrier_list_t no_barriers{0, nullptr};
  static const hrx_queue_barriers_t agent{&no_barriers, &no_barriers};
  // The waited streams' immediate work still needs its system-scope release at
  // the next completion point, which now also covers |waiter|.
  if (unreleased) g_unreleased[w] = true;
  check(hrx_queue_barrier(hrx_stream_get_queue(w), &waits, &signals,
                          join_system_visibility() ? nullptr : &agent),
        "hrx_queue_barrier");
}

void stream_wait_stream(hrx_stream_t waiter, hrx_stream_t waited) {
  stream_wait_streams(waiter, &waited, 1);
}

void stream_synchronize(hrx_stream_t stream) {
  flush_all();
  hrx_stream_t s = resolve(stream);
  wait_stream(s);
  release_staging(s);
}

// HRX's query compares only submitted work, so pending work is flushed first.
bool stream_query(hrx_stream_t stream) {
  flush_all();
  hrx_stream_t s = resolve(stream);
  bool complete = false;
  check(hrx_stream_query(s, &complete), "hrx_stream_query");
  if (complete) release_staging(s);
  return complete;
}

}  // namespace amrex_loom

// --- host random generator (amrex_loom) ----------------------------------------

struct amrex_loom::RandGenerator {
  std::mt19937_64 engine;
};

namespace {
// Draws n values on the host, copies them to device memory at p, and waits.
template <class T, class D>
void rand_generate(amrex_loom::RandGenerator* g, T* p, std::size_t n, D dist) {
  std::vector<T> host(n);
  for (auto& v : host) v = static_cast<T>(dist(g->engine));
  amrex_loom::memcpy_async(p, host.data(), n * sizeof(T), nullptr);
  amrex_loom::device_synchronize();
}
}  // namespace

namespace amrex_loom {
RandGenerator* rand_generator_create(unsigned long long seed) {
  return new RandGenerator{std::mt19937_64(seed)};
}
void rand_generator_destroy(RandGenerator* generator) { delete generator; }
void rand_generate_uniform(RandGenerator* g, float* p, std::size_t n) {
  rand_generate(g, p, n, std::uniform_real_distribution<float>(0.0f, 1.0f));
}
void rand_generate_uniform_double(RandGenerator* g, double* p, std::size_t n) {
  rand_generate(g, p, n, std::uniform_real_distribution<double>(0.0, 1.0));
}
void rand_generate_normal(RandGenerator* g, float* p, std::size_t n, float mean, float stddev) {
  rand_generate(g, p, n, std::normal_distribution<float>(mean, stddev));
}
void rand_generate_normal_double(RandGenerator* g, double* p, std::size_t n, double mean,
                                 double stddev) {
  rand_generate(g, p, n, std::normal_distribution<double>(mean, stddev));
}
}  // namespace amrex_loom

// --- memory (amrex_loom) -------------------------------------------------------

// Pinned host memory: host-local, device-visible HRX buffers. Device code
// dereferences these directly (e.g. ReduceOps::value writes the host tuple),
// so the mapped host address must equal the device address.
static std::map<uintptr_t, hrx_buffer_t> g_host_allocations;

namespace {

// Allocates on the default stream and waits for the backing; nullptr when
// HRX cannot allocate.
hrx_buffer_t allocate_buffer(std::size_t bytes, hrx_memory_type_t type, const char* what) {
  hrx_buffer_t buffer = nullptr;
  hrx_status_t status = hrx_buffer_allocate(g_default_stream, bytes, type,
                                            HRX_BUFFER_USAGE_DEFAULT, &buffer);
  if (!hrx_status_is_ok(status)) {
    hrx_status_ignore(status);
    return nullptr;
  }
  check(hrx_stream_synchronize(g_default_stream), what);
  return buffer;
}

}  // namespace

namespace amrex_loom {

void* device_malloc(std::size_t bytes) {
  initialize();
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (bytes == 0) bytes = 1;
  hrx_buffer_t buffer = allocate_buffer(bytes, HRX_MEMORY_TYPE_DEVICE_LOCAL, "allocation");
  if (!buffer) return nullptr;
  void* pointer = device_address(buffer);
  g_allocations[reinterpret_cast<uintptr_t>(pointer)] = {buffer, bytes};
  return pointer;
}

void device_free(void* pointer) {
  if (!pointer) return;
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto it = g_allocations.find(reinterpret_cast<uintptr_t>(pointer));
  if (it == g_allocations.end()) {
    std::fprintf(stderr, "amrex loom: device_free of %p, not from device_malloc\n", pointer);
    std::abort();
  }
  sync_all();
  hrx_buffer_release(it->second.buffer);
  g_allocations.erase(it);
}

void* pinned_malloc(std::size_t bytes) {
  initialize();
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (bytes == 0) bytes = 1;
  hrx_buffer_t buffer = allocate_buffer(
      bytes, HRX_MEMORY_TYPE_HOST_LOCAL | HRX_MEMORY_TYPE_DEVICE_VISIBLE, "host allocation");
  if (!buffer) return nullptr;
  void* device = device_address(buffer);
  void* host = nullptr;
  check(hrx_buffer_map(buffer, HRX_MAP_READ | HRX_MAP_WRITE | HRX_MAP_MAY_ALIAS, 0, bytes, &host),
        "hrx_buffer_map");
  if (host != device) {
    std::fprintf(stderr, "amrex loom: pinned host address %p differs from device address %p\n",
                 host, device);
    std::abort();
  }
  g_host_allocations[reinterpret_cast<uintptr_t>(host)] = buffer;
  return host;
}

void pinned_free(void* pointer) {
  if (!pointer) return;
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto it = g_host_allocations.find(reinterpret_cast<uintptr_t>(pointer));
  if (it == g_host_allocations.end()) {
    std::fprintf(stderr, "amrex loom: pinned_free of %p, not from pinned_malloc\n", pointer);
    std::abort();
  }
  sync_all();
  hrx_status_ignore(hrx_buffer_unmap(it->second));
  hrx_buffer_release(it->second);
  g_host_allocations.erase(it);
}

namespace {

// Returns the pinned allocation containing |address|, or nullptr.
const std::pair<const uintptr_t, hrx_buffer_t>* find_pinned(const void* address) {
  auto a = reinterpret_cast<uintptr_t>(address);
  auto it = g_host_allocations.upper_bound(a);
  if (it == g_host_allocations.begin()) return nullptr;
  --it;
  std::size_t size = 0;
  if (!hrx_status_is_ok(hrx_buffer_get_size(it->second, &size)) || a >= it->first + size) {
    return nullptr;
  }
  return &*it;
}

}  // namespace

// Host-to-device copies of at most AMREX_LOOM_COPY_KERNEL_MAX bytes (default
// 64 KiB; 0 disables) from pinned memory, including pageable sources staged in
// the stream's pinned ring, run as a copy kernel on the stream instead of an
// HRX queue copy, as HIP copies small transfers with a blit kernel. The
// kernel reads the pinned source through its device address.
std::size_t copy_kernel_max() {
  static const std::size_t limit = [] {
    const char* value = std::getenv("AMREX_LOOM_COPY_KERNEL_MAX");
    return value && *value ? static_cast<std::size_t>(std::strtoull(value, nullptr, 10))
                           : std::size_t(64) << 10;
  }();
  return limit;
}

void copy_kernel(hrx_stream_t s, void* dst, const void* src, std::size_t bytes) {
  constexpr int kThreads = 256;
  const bool words = ((reinterpret_cast<uintptr_t>(dst) | reinterpret_cast<uintptr_t>(src) |
                       bytes) & 7u) == 0;
  if (words) {
    auto* d = static_cast<uint64_t*>(dst);
    const auto* r = static_cast<const uint64_t*>(src);
    const unsigned n = static_cast<unsigned>(bytes / 8);
    amrex_loom::launch_fixed<kThreads>((n + kThreads - 1) / kThreads, kThreads, 0, s,
        [=] AMREX_GPU_DEVICE () {
          const unsigned i = loom::kernel::workgroup::id.x * kThreads + loom::kernel::workitem::id.x;
          if (i < n) d[i] = r[i];
        });
  } else {
    auto* d = static_cast<unsigned char*>(dst);
    const auto* r = static_cast<const unsigned char*>(src);
    const unsigned n = static_cast<unsigned>(bytes);
    amrex_loom::launch_fixed<kThreads>((n + kThreads - 1) / kThreads, kThreads, 0, s,
        [=] AMREX_GPU_DEVICE () {
          const unsigned i = loom::kernel::workgroup::id.x * kThreads + loom::kernel::workitem::id.x;
          if (i < n) d[i] = r[i];
        });
  }
}

void memcpy_async(void* dst, const void* src, std::size_t bytes, hrx_stream_t stream) {
  if (bytes == 0) return;
  hrx_stream_t s = resolve(stream);
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  const auto* d = find_allocation(dst);
  const auto* r = find_allocation(src);
  auto offset = [](const auto* a, const void* p) {
    return reinterpret_cast<uintptr_t>(p) - a->first;
  };
  if (d && r) {
    check(hrx_stream_copy_buffer_linear(s, r->second.buffer, offset(r, src), d->second.buffer,
                                 offset(d, dst), bytes),
          "hrx_stream_copy_buffer_linear");
    recorded(s);
  } else if (const auto* pin = find_pinned(src);
             d && pin && submit_mode() == Submit::immediate) {
    if (bytes <= copy_kernel_max()) {
      copy_kernel(s, dst, src, bytes);
    } else {
      copy_immediate(s, pin->second, offset(pin, src), d->second.buffer, offset(d, dst), bytes);
    }
  } else if (d && bytes <= kPinnedRingCopyLimit && submit_mode() == Submit::immediate) {
    // A small pageable source is copied into |s|'s pinned staging now, so the
    // caller may reuse it on return, as hipMemcpyAsync allows.
    PinnedRing& ring = g_rings[s];
    if (!ring.base) ring.base = static_cast<uint8_t*>(pinned_malloc(kPinnedRingBytes));
    std::size_t at = (ring.used + 255) & ~std::size_t(255);
    if (at + bytes > kPinnedRingBytes) {
      stream_synchronize(s);  // resets ring.used
      at = 0;
    }
    std::memcpy(ring.base + at, src, bytes);
    ring.used = at + bytes;
    if (bytes <= copy_kernel_max()) {
      copy_kernel(s, dst, ring.base + at, bytes);
    } else {
      const auto* staging = find_pinned(ring.base);
      copy_immediate(s, staging->second, offset(staging, ring.base + at), d->second.buffer,
                     offset(d, dst), bytes);
    }
  } else if (d) {
    // HRX reads the host range when the transfer executes. A pinned source
    // (pinned_malloc) stays asynchronous, as with HIP; the caller keeps it
    // unmodified until the stream reaches the transfer. A pageable source is
    // copied first, so the caller may reuse it on return, as hipMemcpyAsync
    // allows.
    const void* from = src;
    if (pointer_kind(src) != PointerKind::pinned) {
      constexpr std::size_t kStagingLimit = std::size_t(256) << 20;
      if (g_staging[s].bytes + bytes > kStagingLimit) stream_synchronize(s);
      Staging& staging = g_staging[s];
      auto copy = std::make_unique<uint8_t[]>(bytes);
      std::memcpy(copy.get(), src, bytes);
      from = copy.get();
      staging.buffers.push_back(std::move(copy));
      staging.bytes += bytes;
    }
    check(hrx_stream_copy_h2d_linear(s, from, d->second.buffer, offset(d, dst), bytes),
          "hrx_stream_copy_h2d_linear");
  } else if (const auto* pin = find_pinned(dst);
             r && pin && submit_mode() == Submit::immediate) {
    copy_immediate(s, r->second.buffer, offset(r, src), pin->second, offset(pin, dst), bytes);
  } else if (r) {
    // HRX writes the host range when the transfer executes. A pinned
    // destination stays asynchronous; a pageable one is complete on return, as
    // hipMemcpyAsync makes it.
    check(hrx_stream_copy_d2h_linear(s, r->second.buffer, offset(r, src), dst, bytes),
          "hrx_stream_copy_d2h_linear");
    if (pointer_kind(dst) != PointerKind::pinned) stream_synchronize(s);
  } else {
    std::memmove(dst, src, bytes);
  }
}

void memset_async(void* dst, int value, std::size_t bytes, hrx_stream_t stream) {
  if (bytes == 0) return;
  hrx_stream_t s = resolve(stream);
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  const auto* d = find_allocation(dst);
  if (!d) {
    std::memset(dst, value, bytes);
    return;
  }
  uint8_t pattern = static_cast<uint8_t>(value);
  check(hrx_stream_fill_buffer_linear(s, d->second.buffer,
                                      reinterpret_cast<uintptr_t>(dst) - d->first, bytes, &pattern, 1),
        "hrx_stream_fill_buffer_linear");
  recorded(s);
}

PointerKind pointer_kind(const void* pointer) {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (find_allocation(pointer)) return PointerKind::device;
  auto a = reinterpret_cast<uintptr_t>(pointer);
  auto it = g_host_allocations.upper_bound(a);
  if (it != g_host_allocations.begin()) {
    --it;
    std::size_t size = 0;
    if (hrx_status_is_ok(hrx_buffer_get_size(it->second, &size)) && a < it->first + size) {
      return PointerKind::pinned;
    }
  }
  return PointerKind::unregistered;
}

}  // namespace amrex_loom

namespace {

// Releases everything this runtime holds, then shuts HRX's GPU accelerator
// down, which also ends an HRX_PROFILE_FILE profiling session.
void shutdown() {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (!g_initialized) return;
  sync_all();
  for (auto& [host, buffer] : g_host_allocations) {
    hrx_status_ignore(hrx_buffer_unmap(buffer));
    hrx_buffer_release(buffer);
  }
  g_host_allocations.clear();
  g_rings.clear();
  for (auto& [address, allocation] : g_allocations) hrx_buffer_release(allocation.buffer);
  g_allocations.clear();
  for (auto& [code_object, executable] : g_executables) hrx_executable_release(executable);
  g_executables.clear();
  // Streams first (each final release joins its work), then the queues this
  // runtime acquired; HRX requires both before the group's final release.
  for (auto& [s, q] : g_streams) {
    hrx_stream_release(s);
    release_stream_queue(q);
  }
  g_streams.clear();
  hrx_stream_release(g_default_stream);
  g_default_stream = nullptr;
  g_provisioned_users.clear();
  g_device = nullptr;
  hrx_device_group_release(g_group);
  g_group = nullptr;
  hrx_runtime_release(g_runtime);
  g_runtime = nullptr;
  g_initialized = false;
}

}  // namespace

#endif  // AMREX_USE_LOOM
