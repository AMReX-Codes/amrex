// Native HRX runtime for AMReX_GPU_BACKEND=LOOM.
//
// Implements the HIP runtime subset declared by Tools/Loom/hip/hip/hip_runtime.h
// on libhrx, plus the registry and launcher for Loom-compiled kernels embedded
// by the amrex-loom-cxx driver. Single device; streams are HRX streams.
#include <AMReX_Config.H>

// The Loom device pass never needs this host-only file.
#if defined(AMREX_USE_LOOM) && !defined(__loom__) && defined(AMREX_LOOM_HRX_HIP)

// AMReX_LOOM_RUNTIME=HRX_HIP: HRX's libamdhip64.so implements the whole HIP
// runtime. This file only registers the Loom code objects embedded by the
// amrex-loom-cxx driver and launches their kernels with a packed kernarg
// buffer through the HIP module API.
#include <hip/hip_runtime.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Argument {
  char kind;           // 'c' scalar bytes, 'b' raw device pointer
  std::size_t offset;  // byte offset in the launch tag object
  std::size_t size;
};

struct Kernel {
  hipFunction_t function = nullptr;
  std::vector<Argument> arguments;
  std::size_t tag_size = 0;  // device pass's sizeof(Tag); 0 if not recorded
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

std::mutex g_mutex;
std::unordered_map<std::string, Kernel> g_kernels;
std::vector<hipModule_t> g_modules;
bool g_images_loaded = false;

void hip_check(hipError_t error, const char* what) {
  if (error == hipSuccess) return;
  std::fprintf(stderr, "amrex loom: %s failed: %s\n", what, hipGetErrorString(error));
  std::abort();
}

void load_images() {
  std::lock_guard<std::mutex> lock(g_mutex);
  if (g_images_loaded) return;
  for (const Image& image : images()) {
    hipModule_t module = nullptr;
    hip_check(hipModuleLoadData(&module, image.data), "hipModuleLoadData");
    g_modules.push_back(module);
    std::istringstream in(image.manifest);
    std::string line, current;
    std::vector<Argument> args;
    std::size_t tag_size = 0;
    auto flush = [&]() {
      if (current.empty() || g_kernels.count(current)) return;
      Kernel k;
      hip_check(hipModuleGetFunction(&k.function, module, current.c_str()),
                "hipModuleGetFunction");
      k.arguments = args;
      k.tag_size = tag_size;
      g_kernels.emplace(current, std::move(k));
    };
    while (std::getline(in, line)) {
      std::istringstream row(line);
      std::string word;
      row >> word;
      if (word == "kernel") {
        flush();
        row >> current;
        args.clear();
        tag_size = 0;
      } else if (word == "size") {
        row >> tag_size;
      } else if (word == "arg") {
        int parameter = 0;
        Argument a{};
        row >> parameter >> a.kind >> a.offset >> a.size;
        args.push_back(a);
      }
    }
    flush();
  }
  g_images_loaded = true;
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

void amrex_loom_launch_named(const char* kernel_name, const char* tag_type_name,
                             const void* tag, std::size_t tag_size, dim3 grid,
                             dim3 block, std::size_t shared_bytes,
                             hipStream_t stream) {
  load_images();
  if (grid.x == 0 || grid.y == 0 || grid.z == 0) return;
  std::string name = std::string(kernel_name) + "__" + tag_type_name;
  auto found = g_kernels.find(name);
  if (found == g_kernels.end()) {
    std::fprintf(stderr, "amrex loom: no device kernel for launch tag %s\n", tag_type_name);
    std::abort();
  }
  if (shared_bytes > AMREX_LOOM_DYNAMIC_SHARED_BYTES) {
    std::fprintf(stderr, "amrex loom: %zu bytes of dynamic shared memory requested, "
                 "%d supported\n", shared_bytes, AMREX_LOOM_DYNAMIC_SHARED_BYTES);
    std::abort();
  }
  if (found->second.tag_size != 0 && found->second.tag_size != tag_size) {
    std::fprintf(stderr, "amrex loom: launch tag size mismatch for %s: host %zu, device pass %zu\n",
                 tag_type_name, tag_size, found->second.tag_size);
    std::abort();
  }
  const auto* bytes = static_cast<const uint8_t*>(tag);
  std::vector<uint8_t> kernarg;
  kernarg.reserve(256);
  for (const Argument& a : found->second.arguments) {
    if (a.offset + a.size > tag_size) {
      std::fprintf(stderr, "amrex loom: manifest/host layout mismatch for %s\n",
                   tag_type_name);
      std::abort();
    }
    // Loom AMDGPU kernargs: buffers 8-byte aligned, by-value scalars 4-byte
    // aligned regardless of size (see the code object's .args .align).
    std::size_t align = a.kind == 'b' ? 8 : std::min<std::size_t>(a.size, 4);
    if (align == 0) align = 1;
    std::size_t offset = (kernarg.size() + align - 1) / align * align;
    kernarg.resize(offset + a.size, 0);
    std::memcpy(kernarg.data() + offset, bytes + a.offset, a.size);
  }
  kernarg.resize((kernarg.size() + 7) / 8 * 8, 0);
  std::size_t kernarg_size = kernarg.size();
  // The buffer is already in the kernel's native kernarg layout.
  void* extra[] = {HIP_LAUNCH_PARAM_BUFFER_POINTER, kernarg.data(),
                   HIP_LAUNCH_PARAM_BUFFER_SIZE, &kernarg_size, HIP_LAUNCH_PARAM_END};
  static const bool trace = std::getenv("AMREX_LOOM_TRACE") != nullptr;
  if (trace) {
    std::fprintf(stderr, "amrex loom: launch %s grid %ux%ux%u block %ux%ux%u kernarg %zu\n",
                 name.c_str(), grid.x, grid.y, grid.z, block.x, block.y, block.z,
                 kernarg_size);
  }
  // Dynamic shared memory is a fixed workgroup buffer inside the kernel.
  hip_check(hipModuleLaunchKernel(found->second.function, grid.x, grid.y, grid.z,
                                  block.x, block.y, block.z, 0, stream, nullptr, extra),
            "hipModuleLaunchKernel");
  if (std::getenv("AMREX_LOOM_SYNC_LAUNCH")) {
    hip_check(hipStreamSynchronize(stream), "hipStreamSynchronize");
  }
}

#elif defined(AMREX_USE_LOOM) && !defined(__loom__)

#include <hip/hip_runtime.h>
#include <hrx_runtime.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

namespace {

struct Argument {
  char kind;           // 'c' scalar bytes, 'b' raw device pointer
  std::size_t offset;  // byte offset in the launch tag object
  std::size_t size;
};

struct Kernel {
  hrx_executable_t executable = nullptr;
  uint32_t ordinal = 0;
  std::vector<Argument> arguments;
  std::size_t tag_size = 0;  // device pass's sizeof(Tag); 0 if not recorded
};

struct Allocation {
  hrx_buffer_t buffer;
  std::size_t size;
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

std::recursive_mutex g_mutex;
bool g_initialized = false;
hrx_device_t g_device = nullptr;
hrx_stream_t g_default_stream = nullptr;
std::set<hrx_stream_t> g_streams;
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

void initialize() {
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  if (g_initialized) return;
  check(hrx_gpu_initialize(0), "hrx_gpu_initialize");
  check(hrx_gpu_device_get(0, &g_device), "hrx_gpu_device_get");
  check(hrx_stream_create(g_device, 0, &g_default_stream), "hrx_stream_create");
  g_initialized = true;
}

hrx_stream_t resolve(hipStream_t stream) {
  initialize();
  return stream ? reinterpret_cast<hrx_stream_t>(stream) : g_default_stream;
}

void sync_all() {
  initialize();
  check(hrx_stream_synchronize(g_default_stream), "hrx_stream_synchronize");
  for (auto* s : g_streams) check(hrx_stream_synchronize(s), "hrx_stream_synchronize");
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
    std::istringstream in(image.manifest);
    std::string line, current;
    std::vector<Argument> args;
    std::size_t tag_size = 0;
    auto flush = [&]() {
      if (current.empty() || g_kernels.count(current)) return;
      Kernel k;
      k.executable = executable;
      check(hrx_executable_lookup_export_by_name(executable, current.c_str(), &k.ordinal),
            "hrx_executable_lookup_export_by_name");
      k.arguments = args;
      k.tag_size = tag_size;
      g_kernels.emplace(current, std::move(k));
    };
    while (std::getline(in, line)) {
      std::istringstream row(line);
      std::string word;
      row >> word;
      if (word == "kernel") {
        flush();
        row >> current;
        args.clear();
        tag_size = 0;
      } else if (word == "size") {
        row >> tag_size;
      } else if (word == "arg") {
        int parameter = 0;
        Argument a{};
        row >> parameter >> a.kind >> a.offset >> a.size;
        args.push_back(a);
      }
    }
    flush();
  }
  g_images_loaded = true;
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

void amrex_loom_launch_named(const char* kernel_name, const char* tag_type_name,
                             const void* tag, std::size_t tag_size, dim3 grid,
                             dim3 block, std::size_t shared_bytes,
                             hipStream_t stream) {
  load_images();
  if (grid.x == 0 || grid.y == 0 || grid.z == 0) return;
  std::string name = std::string(kernel_name) + "__" + tag_type_name;
  auto found = g_kernels.find(name);
  if (found == g_kernels.end()) {
    std::fprintf(stderr, "amrex loom: no device kernel for launch tag %s\n", tag_type_name);
    std::abort();
  }
  if (shared_bytes > AMREX_LOOM_DYNAMIC_SHARED_BYTES) {
    std::fprintf(stderr, "amrex loom: %zu bytes of dynamic shared memory requested, "
                 "%d supported\n", shared_bytes, AMREX_LOOM_DYNAMIC_SHARED_BYTES);
    std::abort();
  }
  // Raw explicit kernarg segment: arguments in manifest order at natural
  // alignment; pointers are raw device addresses.
  // The device pass computed the tag layout independently of the host
  // compiler; a size disagreement means the argument offsets cannot be trusted.
  if (found->second.tag_size != 0 && found->second.tag_size != tag_size) {
    std::fprintf(stderr, "amrex loom: launch tag size mismatch for %s: host %zu, device pass %zu\n",
                 tag_type_name, tag_size, found->second.tag_size);
    std::abort();
  }
  const auto* bytes = static_cast<const uint8_t*>(tag);
  std::vector<uint8_t> kernarg;
  kernarg.reserve(256);
  for (const Argument& a : found->second.arguments) {
    if (a.offset + a.size > tag_size) {
      std::fprintf(stderr, "amrex loom: manifest/host layout mismatch for %s\n",
                   tag_type_name);
      std::abort();
    }
    // Loom AMDGPU kernargs: buffers 8-byte aligned, by-value scalars 4-byte
    // aligned regardless of size (see the code object's .args .align).
    std::size_t align = a.kind == 'b' ? 8 : std::min<std::size_t>(a.size, 4);
    if (align == 0) align = 1;
    std::size_t offset = (kernarg.size() + align - 1) / align * align;
    kernarg.resize(offset + a.size, 0);
    std::memcpy(kernarg.data() + offset, bytes + a.offset, a.size);
  }
  kernarg.resize((kernarg.size() + 7) / 8 * 8, 0);
  hrx_dispatch_config_t config{};
  config.workgroup_count[0] = grid.x;
  config.workgroup_count[1] = grid.y;
  config.workgroup_count[2] = grid.z;
  config.workgroup_size[0] = block.x;
  config.workgroup_size[1] = block.y;
  config.workgroup_size[2] = block.z;
  hrx_stream_t s = resolve(stream);
  static const bool trace = std::getenv("AMREX_LOOM_TRACE") != nullptr;
  if (trace) {
    std::fprintf(stderr, "amrex loom: launch %s grid %ux%ux%u block %ux%ux%u kernarg %zu\n",
                 name.c_str(), grid.x, grid.y, grid.z, block.x, block.y, block.z,
                 kernarg.size());
  }
  check(hrx_stream_dispatch(s, found->second.executable, found->second.ordinal, &config,
                            kernarg.data(), kernarg.size(), nullptr, 0,
                            HRX_DISPATCH_FLAG_CUSTOM_DIRECT_ARGUMENTS),
        "hrx_stream_dispatch");
  check(hrx_stream_flush(s), "hrx_stream_flush");
  if (std::getenv("AMREX_LOOM_SYNC_LAUNCH")) {
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
  check(hrx_buffer_get_device_ptr(buffer, pointer), "hrx_buffer_get_device_ptr");
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
  check(hrx_buffer_get_device_ptr(buffer, &device), "hrx_buffer_get_device_ptr");
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
    sync_all();
    check(hrx_synchronous_h2d(g_device, src, d->second.buffer, doff(d, dst), bytes),
          "hrx_synchronous_h2d");
  } else if (r) {
    sync_all();
    check(hrx_synchronous_d2h(g_device, r->second.buffer, doff(r, src), dst, bytes),
          "hrx_synchronous_d2h");
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
hipError_t hipStreamCreateWithFlags(hipStream_t* stream, unsigned) {
  initialize();
  hrx_stream_t s = nullptr;
  check(hrx_stream_create(g_device, 0, &s), "hrx_stream_create");
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  g_streams.insert(s);
  *stream = reinterpret_cast<hipStream_t>(s);
  return hipSuccess;
}
hipError_t hipStreamDestroy(hipStream_t stream) {
  if (!stream) return hipSuccess;
  auto s = reinterpret_cast<hrx_stream_t>(stream);
  check(hrx_stream_synchronize(s), "hrx_stream_synchronize");
  std::lock_guard<std::recursive_mutex> lock(g_mutex);
  g_streams.erase(s);
  hrx_stream_release(s);
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
// Callbacks run once all prior work on the stream has completed.
hipError_t hipStreamAddCallback(hipStream_t stream, hipStreamCallback_t callback,
                                void* user_data, unsigned) {
  check(hrx_stream_synchronize(resolve(stream)), "hrx_stream_synchronize");
  callback(stream, hipSuccess, user_data);
  return hipSuccess;
}
hipError_t hipStreamWaitEvent(hipStream_t, hipEvent_t, unsigned) { sync_all(); return hipSuccess; }
hipError_t hipEventCreate(hipEvent_t* event) { *event = nullptr; return hipSuccess; }
hipError_t hipEventDestroy(hipEvent_t) { return hipSuccess; }
hipError_t hipEventRecord(hipEvent_t, hipStream_t stream) { return hipStreamSynchronize(stream); }
hipError_t hipEventSynchronize(hipEvent_t) { return hipSuccess; }
hipError_t hipEventElapsedTime(float* ms, hipEvent_t, hipEvent_t) { *ms = 0.0f; return hipSuccess; }

#endif  // AMREX_USE_LOOM
