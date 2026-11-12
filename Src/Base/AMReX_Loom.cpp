// Host runtime for AMReX_GPU_BACKEND=LOOM.
//
// loom-cxx, the C++ compiler of LOOM builds, embeds each translation unit's
// Loom code object and registers it before main(). A launch names its kernel by a Loom kernel
// descriptor (loom::kernel::kernel_ref, defined locally in each object by
// loom-cxx), dispatches the export the descriptor numbers in the code object
// it names, and packs the kernel arguments from the host launch tag following
// the descriptor's leaves. This file implements the
// amrex_loom runtime functions declared in AMReX_Loom.H and AMReX_LoomRandom.H
// on libhrx.
// Single device; streams are HRX streams, each bound to its own HRX queue
// where the device allows.
#include <AMReX_Config.H>

#if defined(AMREX_USE_LOOM)

#include <AMReX_GpuQualifiers.H>  // AMREX_LOOM_DYNAMIC_SHARED_BYTES
#include <AMReX_Loom.H>
#include <AMReX_LoomRandom.H>
#include <hrx_runtime.h>
#include <loomcxx/kernel_ref.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <random>
#include <unordered_map>
#include <vector>

namespace {

// --- launch registry ---------------------------------------------------------

struct Image {
  const void* data;
  std::size_t size;
};

std::vector<Image>& images() {
  static std::vector<Image> list;
  return list;
}

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

// Packs the kernel's explicit kernarg segment from the host launch tag at the
// offsets the descriptor gives.
std::vector<uint8_t> pack_kernargs(const loom::kernel::descriptor& kernel,
                                   const void* parameter) {
  const auto* bytes = static_cast<const uint8_t*>(parameter);
  std::vector<uint8_t> kernarg(kernel.kernarg_size, 0);
  for (unsigned i = 0; i < kernel.leaf_count; ++i) {
    const loom::kernel::descriptor_leaf& leaf = kernel.leaves[i];
    std::memcpy(kernarg.data() + leaf.kernarg_offset, bytes + leaf.source_offset, leaf.size);
  }
  return kernarg;
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

}  // namespace

extern "C" void amrex_loom_register_executable(const void* data, std::size_t size) {
  images().push_back({data, size});
}

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
// The loaded executable of each registered code object, by its address.
std::unordered_map<const void*, hrx_executable_t> g_executables;
bool g_images_loaded = false;

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
  const hrx_device_create_params_t params = {infos[0].uri, HRX_WAIT_FLAG_YIELD};
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

void sync_all() {
  initialize();
  check(hrx_stream_synchronize(g_default_stream), "hrx_stream_synchronize");
  for (auto& [s, q] : g_streams) check(hrx_stream_synchronize(s), "hrx_stream_synchronize");
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

void load_images() {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (g_images_loaded) return;
  initialize();
  const char* target = std::getenv("AMREX_LOOM_TARGET");
  if (!target || !*target) target = "gfx90a";
  for (const Image& image : images()) {
    hrx_executable_t executable = nullptr;
    check(hrx_executable_load_data(g_device, image.data, image.size, "amdgpu", target,
                                   &executable),
          "hrx_executable_load_data");
    g_executables.emplace(image.data, executable);
  }
  g_images_loaded = true;
}

}  // namespace

// --- launch (amrex_loom) ------------------------------------------------------

namespace {

// The executable of the code object |kernel| names. Each translation unit's
// descriptors name the code object of that unit, which carries the kernels it
// launches: kernels of different units may share an export name (lambdas with
// internal linkage, for example) and only the unit's own code object has the
// right one. The descriptor names the kernel by its export ordinal there.
hrx_executable_t kernel_executable(const loom::kernel::descriptor& kernel) {
  if (!kernel.code_object) {
    std::fprintf(stderr, "amrex loom: kernel descriptor of %s names no code object\n",
                 kernel.name);
    std::abort();
  }
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  load_images();
  auto found = g_executables.find(kernel.code_object);
  if (found == g_executables.end()) {
    std::fprintf(stderr, "amrex loom: the code object of %s is not registered\n", kernel.name);
    std::abort();
  }
  return found->second;
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
  std::vector<uint8_t> kernarg = pack_kernargs(descriptor, parameter);
  trace_launch(descriptor.name, grid, block, kernarg.size());
  hrx_dispatch_config_t config{};
  config.workgroup_count[0] = grid.x;
  config.workgroup_count[1] = grid.y;
  config.workgroup_count[2] = grid.z;
  config.workgroup_size[0] = block.x;
  config.workgroup_size[1] = block.y;
  config.workgroup_size[2] = block.z;
  hrx_stream_t s = resolve(stream);
  check(hrx_stream_dispatch(s, executable, descriptor.export_ordinal, &config,
                            kernarg.data(), kernarg.size(), nullptr, 0,
                            HRX_DISPATCH_FLAG_CUSTOM_DIRECT_ARGUMENTS),
        "hrx_stream_dispatch");
  check(hrx_stream_flush(s), "hrx_stream_flush");
  if (sync_after_launch()) {
    check(hrx_stream_synchronize(s), "hrx_stream_synchronize");
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

// HRX does not report these; values describe an MI200-class GCD.
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
  prop->maxThreadsPerBlock = 1024;
  prop->maxThreadsDim[0] = prop->maxThreadsDim[1] = prop->maxThreadsDim[2] = 1024;
  prop->maxGridSize[0] = 2147483647;
  prop->maxGridSize[1] = prop->maxGridSize[2] = 65535;
  prop->multiProcessorCount = 104;
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
  check(hrx_stream_synchronize(stream), "hrx_stream_synchronize");
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto it = g_streams.find(stream);
  if (it == g_streams.end()) {
    std::fprintf(stderr, "amrex loom: stream_destroy of a stream not from stream_create\n");
    std::abort();
  }
  StreamQueue q = it->second;
  g_streams.erase(it);
  hrx_stream_release(stream);
  release_stream_queue(q);
}

void stream_synchronize(hrx_stream_t stream) {
  check(hrx_stream_synchronize(resolve(stream)), "hrx_stream_synchronize");
}

bool stream_query(hrx_stream_t stream) {
  bool complete = false;
  check(hrx_stream_query(resolve(stream), &complete), "hrx_stream_query");
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
    check(hrx_stream_copy_buffer(s, r->second.buffer, offset(r, src), d->second.buffer,
                                 offset(d, dst), bytes),
          "hrx_stream_copy_buffer");
    check(hrx_stream_flush(s), "hrx_stream_flush");
  } else if (d) {
    // Host<->device copies are ordered on their stream: the host range stays
    // live until the stream reaches the transfer.
    check(hrx_stream_copy_h2d(s, src, d->second.buffer, offset(d, dst), bytes),
          "hrx_stream_copy_h2d");
  } else if (r) {
    check(hrx_stream_copy_d2h(s, r->second.buffer, offset(r, src), dst, bytes),
          "hrx_stream_copy_d2h");
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
  check(hrx_stream_fill_buffer(s, d->second.buffer,
                               reinterpret_cast<uintptr_t>(dst) - d->first, bytes, &pattern, 1),
        "hrx_stream_fill_buffer");
  check(hrx_stream_flush(s), "hrx_stream_flush");
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
  for (auto& [address, allocation] : g_allocations) hrx_buffer_release(allocation.buffer);
  g_allocations.clear();
  for (auto& [data, executable] : g_executables) hrx_executable_release(executable);
  g_executables.clear();
  g_images_loaded = false;
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
