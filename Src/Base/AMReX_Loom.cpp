// Host runtime for AMReX_GPU_BACKEND=LOOM.
//
// The launch registry is shared: the amrex-loom-cxx driver embeds each
// translation unit's Loom code object and kernel manifest, registers them
// before main(), and launches pack kernel arguments from the host launch tag
// following the manifest. This file also implements the HIP runtime subset
// declared by Tools/Loom/hip/hip/hip_runtime.h on libhrx. Single device;
// streams are HRX streams, each bound to its own HRX queue where the device
// allows. Needs no HIP headers.
#include <AMReX_Config.H>

// The Loom device pass never needs this host-only file.
#if defined(AMREX_USE_LOOM) && !defined(__loom__)

#include <hip/hip_runtime.h>
#include <hrx_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

// --- launch registry ---------------------------------------------------------

// One kernel ABI argument, copied from the host launch tag.
struct Argument {
  char kind;           // 'c' by-value bytes, 'b' device pointer (buffer)
  std::size_t offset;  // byte offset in the launch tag object
  std::size_t size;
};

// One kernel's manifest block.
struct ManifestKernel {
  std::string name;          // export name: <kernel>__<Itanium tag type name>
  std::size_t tag_size = 0;  // sizeof(Tag) seen by the device pass
  std::vector<Argument> arguments;
};

struct Image {
  const void* data;
  std::size_t size;
  const char* manifest;
};

std::vector<Image>& images() {
  static std::vector<Image> list;
  return list;
}

// Parses the importer's kernel manifest:
//   kernel <export> <parameter_count>
//   size <sizeof> <alignof>
//   arg <parameter> <c|b> <byte_offset> <byte_size>
std::vector<ManifestKernel> parse_manifest(const char* manifest) {
  std::vector<ManifestKernel> kernels;
  std::istringstream in(manifest);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream row(line);
    std::string word;
    row >> word;
    if (word == "kernel") {
      kernels.emplace_back();
      row >> kernels.back().name;
    } else if (kernels.empty()) {
      continue;
    } else if (word == "size") {
      row >> kernels.back().tag_size;
    } else if (word == "arg") {
      int parameter = 0;
      Argument a{};
      row >> parameter >> a.kind >> a.offset >> a.size;
      kernels.back().arguments.push_back(a);
    }
  }
  return kernels;
}

// Checks a launch against its kernel and returns false for an empty grid.
bool check_launch(const ManifestKernel* kernel, const char* tag_type_name,
                  std::size_t tag_size, dim3 grid, std::size_t shared_bytes) {
  if (grid.x == 0 || grid.y == 0 || grid.z == 0) return false;
  if (!kernel) {
    std::fprintf(stderr, "amrex loom: no device kernel for launch tag %s\n", tag_type_name);
    std::abort();
  }
  if (shared_bytes > AMREX_LOOM_DYNAMIC_SHARED_BYTES) {
    std::fprintf(stderr, "amrex loom: %zu bytes of dynamic shared memory requested, "
                 "%d supported\n", shared_bytes, AMREX_LOOM_DYNAMIC_SHARED_BYTES);
    std::abort();
  }
  // Both passes use one frontend, so the sizes agree unless the device pass
  // saw a different translation unit (for example, a differing macro).
  if (kernel->tag_size != 0 && kernel->tag_size != tag_size) {
    std::fprintf(stderr, "amrex loom: launch tag size mismatch for %s: host %zu, device pass %zu\n",
                 tag_type_name, tag_size, kernel->tag_size);
    std::abort();
  }
  return true;
}

// Packs the kernel's explicit kernarg segment from the host launch tag. Loom
// AMDGPU kernargs align buffers to 8 bytes and by-value scalars to 4 bytes
// regardless of size (the code object's .args .align), not C alignment.
std::vector<uint8_t> pack_kernargs(const ManifestKernel& kernel, const char* tag_type_name,
                                   const void* tag, std::size_t tag_size) {
  const auto* bytes = static_cast<const uint8_t*>(tag);
  std::vector<uint8_t> kernarg;
  kernarg.reserve(256);
  for (const Argument& a : kernel.arguments) {
    if (a.offset + a.size > tag_size) {
      std::fprintf(stderr, "amrex loom: manifest/host layout mismatch for %s\n",
                   tag_type_name);
      std::abort();
    }
    std::size_t align = a.kind == 'b' ? 8 : std::min<std::size_t>(a.size, 4);
    if (align == 0) align = 1;
    std::size_t offset = (kernarg.size() + align - 1) / align * align;
    kernarg.resize(offset + a.size, 0);
    std::memcpy(kernarg.data() + offset, bytes + a.offset, a.size);
  }
  kernarg.resize((kernarg.size() + 7) / 8 * 8, 0);
  return kernarg;
}

void trace_launch(const std::string& name, dim3 grid, dim3 block, std::size_t kernarg_size) {
  static const bool trace = std::getenv("AMREX_LOOM_TRACE") != nullptr;
  if (trace) {
    std::fprintf(stderr, "amrex loom: launch %s grid %ux%ux%u block %ux%ux%u kernarg %zu\n",
                 name.c_str(), grid.x, grid.y, grid.z, block.x, block.y, block.z,
                 kernarg_size);
  }
}

bool sync_after_launch() {
  static const bool sync = std::getenv("AMREX_LOOM_SYNC_LAUNCH") != nullptr;
  return sync;
}

}  // namespace

extern "C" void amrex_loom_register_executable(const void* data, std::size_t size,
                                               const char* manifest) {
  images().push_back({data, size, manifest});
}

void amrex_loom_launch(const char* tag_type_name, const void* tag, std::size_t tag_size,
                       dim3 grid, dim3 block, std::size_t shared_bytes,
                       hipStream_t stream) {
  amrex_loom_launch_named("kernel", tag_type_name, tag, tag_size, grid, block,
                          shared_bytes, stream);
}

// --- HIP runtime subset on libhrx ---------------------------------------------

namespace {

struct Kernel {
  ManifestKernel manifest;
  hrx_executable_t executable = nullptr;
  uint32_t ordinal = 0;
};

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
std::unordered_map<std::string, Kernel> g_kernels;
std::vector<hrx_executable_t> g_executables;
bool g_images_loaded = false;
hipError_t g_last_error = hipSuccess;

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

hrx_stream_t resolve(hipStream_t stream) {
  initialize();
  return stream ? reinterpret_cast<hrx_stream_t>(stream) : g_default_stream;
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
    g_executables.push_back(executable);
    for (auto& entry : parse_manifest(image.manifest)) {
      if (g_kernels.count(entry.name)) continue;
      Kernel k;
      k.executable = executable;
      check(hrx_executable_lookup_export_by_name(executable, entry.name.c_str(), &k.ordinal),
            "hrx_executable_lookup_export_by_name");
      std::string name = entry.name;
      k.manifest = std::move(entry);
      g_kernels.emplace(std::move(name), std::move(k));
    }
  }
  g_images_loaded = true;
}

}  // namespace

void amrex_loom_launch_named(const char* kernel_name, const char* tag_type_name,
                             const void* tag, std::size_t tag_size, dim3 grid,
                             dim3 block, std::size_t shared_bytes,
                             hipStream_t stream) {
  load_images();
  std::string name = std::string(kernel_name) + "__" + tag_type_name;
  auto found = g_kernels.find(name);
  const Kernel* kernel = found == g_kernels.end() ? nullptr : &found->second;
  if (!check_launch(kernel ? &kernel->manifest : nullptr, tag_type_name, tag_size, grid,
                    shared_bytes)) {
    return;
  }
  std::vector<uint8_t> kernarg = pack_kernargs(kernel->manifest, tag_type_name, tag, tag_size);
  trace_launch(name, grid, block, kernarg.size());
  hrx_dispatch_config_t config{};
  config.workgroup_count[0] = grid.x;
  config.workgroup_count[1] = grid.y;
  config.workgroup_count[2] = grid.z;
  config.workgroup_size[0] = block.x;
  config.workgroup_size[1] = block.y;
  config.workgroup_size[2] = block.z;
  hrx_stream_t s = resolve(stream);
  check(hrx_stream_dispatch(s, kernel->executable, kernel->ordinal, &config,
                            kernarg.data(), kernarg.size(), nullptr, 0,
                            HRX_DISPATCH_FLAG_CUSTOM_DIRECT_ARGUMENTS),
        "hrx_stream_dispatch");
  check(hrx_stream_flush(s), "hrx_stream_flush");
  if (sync_after_launch()) {
    check(hrx_stream_synchronize(s), "hrx_stream_synchronize");
  }
}

// --- errors -----------------------------------------------------------------
hipError_t hipGetLastError() { hipError_t e = g_last_error; g_last_error = hipSuccess; return e; }
hipError_t hipPeekAtLastError() { return g_last_error; }
const char* hipGetErrorString(hipError_t error) {
  switch (error) {
    case hipSuccess: return "hipSuccess";
    case hipErrorOutOfMemory: return "out of memory";
    case hipErrorNotReady: return "not ready";
    case hipErrorNotSupported: return "not supported by the LOOM backend";
    default: return "LOOM backend error";
  }
}
const char* hipGetErrorName(hipError_t error) { return hipGetErrorString(error); }

// --- device -----------------------------------------------------------------
hipError_t hipGetDeviceCount(int* count) {
  initialize();
  *count = 1;
  return hipSuccess;
}
hipError_t hipSetDevice(int device) { initialize(); return device == 0 ? hipSuccess : hipErrorInvalidDevice; }
hipError_t hipGetDevice(int* device) { initialize(); *device = 0; return hipSuccess; }

hipError_t hipGetDeviceProperties(hipDeviceProp_t* prop, int) {
  initialize();
  std::memset(prop, 0, sizeof(*prop));
  std::snprintf(prop->name, sizeof(prop->name), "AMD GPU (HRX/Loom)");
  std::snprintf(prop->gcnArchName, sizeof(prop->gcnArchName), "%s",
                std::getenv("AMREX_LOOM_TARGET") ? std::getenv("AMREX_LOOM_TARGET") : "gfx90a");
  std::size_t free_bytes = 0, total_bytes = 0;
  hipMemGetInfo(&free_bytes, &total_bytes);
  prop->totalGlobalMem = total_bytes;
  prop->sharedMemPerBlock = 65536;
  prop->regsPerBlock = 65536;
  prop->warpSize = 64;
  prop->maxThreadsPerBlock = 1024;
  prop->maxThreadsDim[0] = prop->maxThreadsDim[1] = prop->maxThreadsDim[2] = 1024;
  prop->maxGridSize[0] = 2147483647;
  prop->maxGridSize[1] = prop->maxGridSize[2] = 65535;
  prop->clockRate = 1700000;
  prop->totalConstMem = 65536;
  prop->major = 9;
  prop->minor = 0;
  prop->multiProcessorCount = 104;
  prop->maxThreadsPerMultiProcessor = 2048;
  prop->managedMemory = 0;
  prop->concurrentManagedAccess = 0;
  return hipSuccess;
}
hipError_t hipDeviceGetAttribute(int* value, hipDeviceAttribute_t, int) { *value = 0; return hipSuccess; }
hipError_t hipDeviceGetUuid(hipUUID* uuid, hipDevice_t) { std::memset(uuid, 0, sizeof(*uuid)); return hipSuccess; }
hipError_t hipDeviceSynchronize() { sync_all(); return hipSuccess; }
hipError_t hipDeviceReset() { sync_all(); return hipSuccess; }

hipError_t hipMemGetInfo(std::size_t* free_bytes, std::size_t* total_bytes) {
  initialize();
  std::size_t f = 0, t = 0;
  if (!hrx_status_is_ok(hrx_device_memory_info(g_device, &f, &t))) {
    f = t = std::size_t(64) << 30;
  }
  *free_bytes = f;
  *total_bytes = t;
  return hipSuccess;
}

// --- memory -----------------------------------------------------------------
hipError_t hipMalloc(void** pointer, std::size_t bytes) {
  initialize();
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (bytes == 0) bytes = 1;
  hrx_buffer_t buffer = nullptr;
  hrx_status_t status = hrx_buffer_allocate(g_default_stream, bytes, HRX_MEMORY_TYPE_DEVICE_LOCAL,
                                            HRX_BUFFER_USAGE_DEFAULT, &buffer);
  if (!hrx_status_is_ok(status)) {
    hrx_status_ignore(status);
    *pointer = nullptr;
    return g_last_error = hipErrorOutOfMemory;
  }
  check(hrx_stream_synchronize(g_default_stream), "allocation");
  *pointer = device_address(buffer);
  g_allocations[reinterpret_cast<uintptr_t>(*pointer)] = {buffer, bytes};
  return hipSuccess;
}
hipError_t hipMallocAsync(void** pointer, std::size_t bytes, hipStream_t) { return hipMalloc(pointer, bytes); }
hipError_t hipMallocAsync(void** pointer, std::size_t bytes, hipMemPool_t, hipStream_t) {
  return hipMalloc(pointer, bytes);
}
hipError_t hipMallocManaged(void** pointer, std::size_t, unsigned) {
  *pointer = nullptr;
  return g_last_error = hipErrorNotSupported;
}
// Pinned host memory: host-local, device-visible HRX buffers. Device code
// dereferences these directly (e.g. ReduceOps::value writes the host tuple),
// so the mapped host address must equal the device address.
static std::map<uintptr_t, hrx_buffer_t> g_host_allocations;
hipError_t hipHostMalloc(void** pointer, std::size_t bytes, unsigned) {
  initialize();
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (bytes == 0) bytes = 1;
  hrx_buffer_t buffer = nullptr;
  hrx_status_t status = hrx_buffer_allocate(
      g_default_stream, bytes, HRX_MEMORY_TYPE_HOST_LOCAL | HRX_MEMORY_TYPE_DEVICE_VISIBLE,
      HRX_BUFFER_USAGE_DEFAULT, &buffer);
  if (!hrx_status_is_ok(status)) {
    hrx_status_ignore(status);
    *pointer = nullptr;
    return g_last_error = hipErrorOutOfMemory;
  }
  check(hrx_stream_synchronize(g_default_stream), "host allocation");
  void* device = nullptr;
  void* host = nullptr;
  device = device_address(buffer);
  check(hrx_buffer_map(buffer, HRX_MAP_READ | HRX_MAP_WRITE | HRX_MAP_MAY_ALIAS, 0, bytes, &host),
        "hrx_buffer_map");
  if (host != device) {
    std::fprintf(stderr, "amrex loom: pinned host address %p differs from device address %p\n",
                 host, device);
    std::abort();
  }
  g_host_allocations[reinterpret_cast<uintptr_t>(host)] = buffer;
  *pointer = host;
  return hipSuccess;
}
hipError_t hipHostFree(void* pointer) {
  if (!pointer) return hipSuccess;
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto it = g_host_allocations.find(reinterpret_cast<uintptr_t>(pointer));
  if (it == g_host_allocations.end()) return g_last_error = hipErrorInvalidValue;
  check(hrx_stream_synchronize(g_default_stream), "host free");
  hrx_status_ignore(hrx_buffer_unmap(it->second));
  hrx_buffer_release(it->second);
  g_host_allocations.erase(it);
  return hipSuccess;
}
hipError_t hipFree(void* pointer) {
  if (!pointer) return hipSuccess;
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto it = g_allocations.find(reinterpret_cast<uintptr_t>(pointer));
  if (it == g_allocations.end()) return g_last_error = hipErrorInvalidValue;
  sync_all();
  hrx_buffer_release(it->second.buffer);
  g_allocations.erase(it);
  return hipSuccess;
}
hipError_t hipFreeAsync(void* pointer, hipStream_t) { return hipFree(pointer); }

static hipError_t copy(void* dst, const void* src, std::size_t bytes, hrx_stream_t s) {
  if (bytes == 0) return hipSuccess;
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  const auto* d = find_allocation(dst);
  const auto* r = find_allocation(src);
  auto doff = [&](const auto* a, const void* p) {
    return reinterpret_cast<uintptr_t>(p) - a->first;
  };
  if (d && r) {
    check(hrx_stream_copy_buffer(s, r->second.buffer, doff(r, src), d->second.buffer,
                                 doff(d, dst), bytes),
          "hrx_stream_copy_buffer");
    check(hrx_stream_flush(s), "hrx_stream_flush");
  } else if (d) {
    // Host<->device copies are ordered on their stream, as in HIP: the host
    // range stays live until the stream reaches the transfer.
    check(hrx_stream_copy_h2d(s, src, d->second.buffer, doff(d, dst), bytes),
          "hrx_stream_copy_h2d");
  } else if (r) {
    check(hrx_stream_copy_d2h(s, r->second.buffer, doff(r, src), dst, bytes),
          "hrx_stream_copy_d2h");
  } else {
    std::memmove(dst, src, bytes);
  }
  return hipSuccess;
}
hipError_t hipMemcpy(void* dst, const void* src, std::size_t bytes, hipMemcpyKind) {
  hipError_t e = copy(dst, src, bytes, resolve(nullptr));
  sync_all();
  return e;
}
hipError_t hipMemcpyAsync(void* dst, const void* src, std::size_t bytes, hipMemcpyKind,
                          hipStream_t stream) {
  return copy(dst, src, bytes, resolve(stream));
}
hipError_t hipMemsetAsync(void* dst, int value, std::size_t bytes, hipStream_t stream) {
  if (bytes == 0) return hipSuccess;
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  const auto* d = find_allocation(dst);
  if (!d) { std::memset(dst, value, bytes); return hipSuccess; }
  hrx_stream_t s = resolve(stream);
  uint8_t pattern = static_cast<uint8_t>(value);
  check(hrx_stream_fill_buffer(s, d->second.buffer,
                               reinterpret_cast<uintptr_t>(dst) - d->first, bytes, &pattern, 1),
        "hrx_stream_fill_buffer");
  check(hrx_stream_flush(s), "hrx_stream_flush");
  return hipSuccess;
}
hipError_t hipMemset(void* dst, int value, std::size_t bytes) {
  hipError_t e = hipMemsetAsync(dst, value, bytes, nullptr);
  sync_all();
  return e;
}
hipError_t hipMemcpyToSymbolAsync(const void*, const void*, std::size_t, std::size_t,
                                  hipMemcpyKind, hipStream_t) {
  return g_last_error = hipErrorNotSupported;
}
hipError_t hipMemcpyFromSymbolAsync(void*, const void*, std::size_t, std::size_t,
                                    hipMemcpyKind, hipStream_t) {
  return g_last_error = hipErrorNotSupported;
}
hipError_t hipMemAdvise(const void*, std::size_t, hipMemoryAdvise, int) { return hipSuccess; }
hipError_t hipMemPrefetchAsync(const void*, std::size_t, int, hipStream_t) { return hipSuccess; }
hipError_t hipPointerGetAttributes(hipPointerAttribute_t* attributes, const void* pointer) {
  std::memset(attributes, 0, sizeof(*attributes));
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (find_allocation(pointer)) {
    attributes->type = hipMemoryTypeDevice;
    attributes->devicePointer = const_cast<void*>(pointer);
  } else {
    attributes->type = hipMemoryTypeHost;
    attributes->hostPointer = const_cast<void*>(pointer);
    attributes->device = hipCpuDeviceId;
  }
  return hipSuccess;
}
hipError_t hipDeviceGetMemPool(hipMemPool_t* pool, int) { *pool = nullptr; return hipErrorNotSupported; }
hipError_t hipMemPoolSetAttribute(hipMemPool_t, hipMemPoolAttr, void*) { return hipErrorNotSupported; }
hipError_t hipMemPoolGetAttribute(hipMemPool_t, hipMemPoolAttr, void*) { return hipErrorNotSupported; }

// --- streams and events -----------------------------------------------------
hipError_t hipStreamCreate(hipStream_t* stream) { return hipStreamCreateWithFlags(stream, 0); }
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
hipError_t hipStreamCreateWithFlags(hipStream_t* stream, unsigned) {
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
  *stream = reinterpret_cast<hipStream_t>(s);
  return hipSuccess;
}
hipError_t hipStreamDestroy(hipStream_t stream) {
  if (!stream) return hipSuccess;
  auto s = reinterpret_cast<hrx_stream_t>(stream);
  check(hrx_stream_synchronize(s), "hrx_stream_synchronize");
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  auto it = g_streams.find(s);
  if (it == g_streams.end()) return g_last_error = hipErrorInvalidValue;
  StreamQueue q = it->second;
  g_streams.erase(it);
  hrx_stream_release(s);
  release_stream_queue(q);
  return hipSuccess;
}
hipError_t hipStreamSynchronize(hipStream_t stream) {
  check(hrx_stream_synchronize(resolve(stream)), "hrx_stream_synchronize");
  return hipSuccess;
}
hipError_t hipStreamQuery(hipStream_t stream) {
  bool complete = false;
  check(hrx_stream_query(resolve(stream), &complete), "hrx_stream_query");
  return complete ? hipSuccess : hipErrorNotReady;
}
hipError_t hipStreamGetDevice(hipStream_t, hipDevice_t* device) { *device = 0; return hipSuccess; }
namespace {
struct LoomCallbackThunk { hipStream_t stream; hipStreamCallback_t callback; void* user_data; };
hrx_status_t loom_callback_thunk(void* p) {
  auto* thunk = static_cast<LoomCallbackThunk*>(p);
  thunk->callback(thunk->stream, hipSuccess, thunk->user_data);
  delete thunk;
  return hrx_ok_status();
}
}  // namespace
// HIP semantics: the callback runs on a runtime thread once all prior work on
// the stream has completed, and later work on the stream waits until it
// returns. It is a host call that waits on the stream's current timeline point
// and signals the next one. As in HIP, the callback must not call back into
// the stream API.
hipError_t hipStreamAddCallback(hipStream_t stream, hipStreamCallback_t callback,
                                void* user_data, unsigned) {
  hrx_stream_t s = resolve(stream);
  check(hrx_stream_flush(s), "hrx_stream_flush");
  hrx_timeline_point_t before{};
  check(hrx_stream_get_timeline_position(s, &before), "hrx_stream_get_timeline_position");
  uint64_t after = 0;
  check(hrx_stream_advance_timeline(s, &after), "hrx_stream_advance_timeline");
  hrx_semaphore_t semaphore = before.semaphore;
  uint64_t wait_value = before.value, signal_value = after;
  hrx_semaphore_list_t wait{&semaphore, &wait_value, before.value > 0 ? std::size_t(1) : std::size_t(0)};
  hrx_semaphore_list_t signal{&semaphore, &signal_value, 1};
  // The host call goes to the stream's own queue, ordered by its timeline.
  check(hrx_queue_host_call(hrx_stream_get_queue(s), &wait, &signal, loom_callback_thunk,
                            new LoomCallbackThunk{stream, callback, user_data}),
        "hrx_queue_host_call");
  return hipSuccess;
}
hipError_t hipStreamWaitEvent(hipStream_t, hipEvent_t, unsigned) { sync_all(); return hipSuccess; }
hipError_t hipEventCreate(hipEvent_t* event) { *event = nullptr; return hipSuccess; }
hipError_t hipEventDestroy(hipEvent_t) { return hipSuccess; }
hipError_t hipEventRecord(hipEvent_t, hipStream_t stream) { return hipStreamSynchronize(stream); }
hipError_t hipEventSynchronize(hipEvent_t) { return hipSuccess; }
hipError_t hipEventElapsedTime(float* ms, hipEvent_t, hipEvent_t) { *ms = 0.0f; return hipSuccess; }

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
  g_kernels.clear();
  for (auto* executable : g_executables) hrx_executable_release(executable);
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
