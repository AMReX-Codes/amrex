#ifndef AMREX_LOOM_H_
#define AMREX_LOOM_H_
#include <AMReX_Config.H>

// Host runtime of AMReX_GPU_BACKEND=LOOM: thin functions over native libhrx,
// implemented in AMReX_Loom.cpp. AMReX's GPU layer calls these under
// AMREX_USE_LOOM. Streams are HRX streams; device and pinned memory are HRX
// buffers addressed by their device address. Every function initializes the
// runtime on first use and aborts on an HRX failure, except where noted.
//
// The Loom device pass sees these declarations but never calls them, so the
// header names the HRX handle types without including <hrx_runtime.h>.

#ifdef AMREX_USE_LOOM

#include <cstddef>

typedef struct hrx_stream_s* hrx_stream_t;

namespace loom::kernel { struct descriptor; }  // <loomcxx/kernel_ref.h>

namespace amrex_loom {

// Device properties AMReX reads (amrex::gpuDeviceProp_t under LOOM).
struct DeviceProperties {
    char name[256];
    char gcnArchName[256];
    std::size_t totalGlobalMem;
    std::size_t sharedMemPerBlock;
    int warpSize;
    int maxThreadsPerBlock;
    int maxThreadsDim[3];
    int maxGridSize[3];
    int multiProcessorCount;
    int maxThreadsPerMultiProcessor;
    int managedMemory;
    int concurrentManagedAccess;
    int integrated;
};

// --- device ------------------------------------------------------------------

//! Number of devices the runtime exposes (one).
int device_count ();
//! Properties of device 0.
void device_properties (DeviceProperties* properties);
//! 16-byte device identifier (all zero: HRX does not report one).
void device_uuid (char bytes[16]);
//! Free and total device memory in bytes.
void memory_info (std::size_t* free_bytes, std::size_t* total_bytes);
//! Waits for the work of every stream.
void device_synchronize ();

// --- streams -----------------------------------------------------------------

//! Creates a stream on an HRX queue of its own where the device allows.
hrx_stream_t stream_create ();
//! Waits for |stream|'s work, then releases it and its queue.
void stream_destroy (hrx_stream_t stream);
//! Flushes |stream| and waits for its work.
void stream_synchronize (hrx_stream_t stream);
//! True when all work submitted to |stream| has completed.
bool stream_query (hrx_stream_t stream);

// --- launch ------------------------------------------------------------------

//! Dispatches |kernel| (a loom::kernel::kernel_ref descriptor, defined by
//! the Loom device pass) on |stream| with a grid of grid[0..2] workgroups of
//! block[0..2] work-items. block must be fixed_block x 1 x 1. The kernel
//! arguments are packed from its one parameter object at |parameter|
//! following the descriptor's leaves. An empty grid launches nothing.
//! AMReX_LoomLaunch.H calls this; see there.
void dispatch (const loom::kernel::descriptor& kernel, const void* parameter,
               const unsigned grid[3], const unsigned block[3], unsigned fixed_block,
               std::size_t shared_bytes, hrx_stream_t stream);

// --- memory ------------------------------------------------------------------

//! Device-local buffer; returns its device address, or nullptr when out of
//! device memory.
void* device_malloc (std::size_t bytes);
//! Waits for all streams, then releases the device_malloc buffer at |pointer|.
void device_free (void* pointer);
//! Host-local, device-visible buffer whose host and device addresses agree;
//! returns nullptr when out of memory.
void* pinned_malloc (std::size_t bytes);
//! Waits for all streams, then releases the pinned_malloc buffer at |pointer|.
void pinned_free (void* pointer);
//! Copies |bytes| after prior work on |stream|. Either side may be device
//! memory (device_malloc) or host memory, including pinned memory; host
//! ranges stay live until the stream reaches the copy.
void memcpy_async (void* dst, const void* src, std::size_t bytes, hrx_stream_t stream);
//! Fills |bytes| with the byte |value| after prior work on |stream|.
void memset_async (void* dst, int value, std::size_t bytes, hrx_stream_t stream);

enum class PointerKind { unregistered, device, pinned };
//! What kind of LOOM allocation, if any, contains |pointer|.
PointerKind pointer_kind (const void* pointer);

}  // namespace amrex_loom

#endif  // AMREX_USE_LOOM
#endif  // AMREX_LOOM_H_
