#include <AMReX_GpuSparse.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuError.H>

#include <mutex>
#include <new>
#include <string>
#include <utility>

// Creating a cuSPARSE or rocSPARSE handle allocates device resources and costs
// of the order of a millisecond, far more than a sparse matrix-vector product
// on a coarse multigrid level.  The handles are therefore kept, one per host
// thread, and created on first use.  When a thread exits, its handle becomes
// a spare for the next thread that needs one; amrex::Finalize releases all.
//
// A thread that exits makes no calls into the GPU runtime or library, which
// may have torn down their own per-thread state by then.  Nothing that can
// throw or abort runs while sparse_handles_mutex is held, except the
// registration with amrex::ExecOnFinalize, and a thread that exits while it
// holds the mutex (exit() from an abort or a signal handler) does not take it
// again.  Runtime calls under the lock only return their errors.

namespace amrex::Gpu {

namespace {

#if defined(AMREX_USE_CUDA)
    using SparseHandle = cusparseHandle_t;
    using Stream = cudaStream_t;
#else
    using SparseHandle = rocsparse_handle;
    using Stream = hipStream_t;
#endif

    struct ThreadHandle;

    // The handle of a thread that exited, for the next thread that needs one.
    struct SpareHandle
    {
        SparseHandle handle{};
        Stream stream{};  // the stream it was last bound to
        SpareHandle* next = nullptr;
    };

    // Guards the two lists and finalize_registered, and makes thread exit and
    // amrex::Finalize hand over or release each handle exactly once.
    std::mutex sparse_handles_mutex;
    // the threads that hold a handle
    ThreadHandle* live_head = nullptr;
    // the handles of the threads that exited
    SpareHandle* spare_head = nullptr;
    bool finalize_registered = false;

    // True while this thread holds or takes sparse_handles_mutex.
    thread_local bool this_thread_in_registry = false;

    class RegistryLock
    {
    public:
        RegistryLock () {
            this_thread_in_registry = true;
            sparse_handles_mutex.lock();
        }
        ~RegistryLock () {
            sparse_handles_mutex.unlock();
            this_thread_in_registry = false;
        }
        RegistryLock (RegistryLock const&) = delete;
        RegistryLock (RegistryLock&&) = delete;
        RegistryLock& operator= (RegistryLock const&) = delete;
        RegistryLock& operator= (RegistryLock&&) = delete;
    };

    // Release a handle.  It may still be bound to a stream that no longer
    // exists (e.g., an external stream), so it is bound to the default stream
    // first; the caller has synchronized the device.  Errors are returned.
    bool release (SparseHandle handle) noexcept
    {
#if defined(AMREX_USE_CUDA)
        bool ok = cusparseSetStream(handle, nullptr) == CUSPARSE_STATUS_SUCCESS;
        ok = (cusparseDestroy(handle) == CUSPARSE_STATUS_SUCCESS) && ok;
#else
        bool ok = rocsparse_set_stream(handle, nullptr) == rocsparse_status_success;
        ok = (rocsparse_destroy_handle(handle) == rocsparse_status_success) && ok;
#endif
        return ok;
    }

    struct ThreadHandle
    {
        SparseHandle handle{};
        Stream stream{};  // the stream it was last bound to
        // in the list of live handles if handle is set
        ThreadHandle* prev = nullptr;
        ThreadHandle* next = nullptr;

        ThreadHandle () = default;
        ThreadHandle (ThreadHandle const&) = delete;
        ThreadHandle (ThreadHandle&&) = delete;
        ThreadHandle& operator= (ThreadHandle const&) = delete;
        ThreadHandle& operator= (ThreadHandle&&) = delete;

        // with sparse_handles_mutex held
        void link (SparseHandle h) noexcept
        {
            handle = h;
            prev = nullptr;
            next = live_head;
            if (live_head) { live_head->prev = this; }
            live_head = this;
        }

        // with sparse_handles_mutex held; returns the handle
        SparseHandle unlink () noexcept
        {
            if (prev) { prev->next = next; } else { live_head = next; }
            if (next) { next->prev = prev; }
            prev = next = nullptr;
            return std::exchange(handle, SparseHandle{});
        }

        // The thread exits: its handle becomes a spare, unless amrex::Finalize
        // has released it.  No GPU runtime or library calls here.
        ~ThreadHandle ()
        {
            // exit() while this thread holds the mutex: the process ends
            if (this_thread_in_registry) { return; }
            // allocated outside the lock; without it, the handle is left to
            // the process teardown
            auto* spare = new (std::nothrow) SpareHandle;
            {
                RegistryLock lock;
                if (handle && spare) {
                    spare->stream = stream;
                    spare->handle = unlink();
                    spare->next = spare_head;
                    spare_head = std::exchange(spare, nullptr);
                } else if (handle) {
                    (void) unlink();
                }
            }
            delete spare;  // if not used
        }
    };

    void destroy_sparse_handles ()
    {
        bool ok = true;
        SpareHandle* spares = nullptr;
        {
            // Held while releasing, so that a thread exiting meanwhile waits
            // and then finds its handle released.  Nothing here aborts.
            RegistryLock lock;
            finalize_registered = false;
            spares = std::exchange(spare_head, nullptr);
            if (live_head || spares) {
                // All handles were made on AMReX's device.  Work queued with
                // them may still run, also on streams that AMReX does not
                // know: synchronize the whole device.
                int current = -1;
                int const device = Gpu::Device::deviceId();
#if defined(AMREX_USE_CUDA)
                ok = cudaGetDevice(&current) == cudaSuccess;
                bool const other = ok && current != device;
                if (other) { ok = cudaSetDevice(device) == cudaSuccess; }
                ok = ok && (cudaDeviceSynchronize() == cudaSuccess);
#else
                ok = hipGetDevice(&current) == hipSuccess;
                bool const other = ok && current != device;
                if (other) { ok = hipSetDevice(device) == hipSuccess; }
                ok = ok && (hipDeviceSynchronize() == hipSuccess);
#endif
                // No other thread may make sparse calls during amrex::Finalize,
                // so the owners of live handles do not use them now.
                while (live_head) {
                    ok = release(live_head->unlink()) && ok;
                }
                for (auto* s = spares; s; s = s->next) {
                    ok = release(s->handle) && ok;
                }
#if defined(AMREX_USE_CUDA)
                if (other) { ok = (cudaSetDevice(current) == cudaSuccess) && ok; }
#else
                if (other) { ok = (hipSetDevice(current) == hipSuccess) && ok; }
#endif
            }
        }
        while (spares) {
            delete std::exchange(spares, spares->next);
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ok, "Gpu::destroy_sparse_handles failed");
    }

    SparseHandle get_sparse_handle ()
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::Initialized(),
            "Gpu::getCusparseHandle/getRocsparseHandle: AMReX is not initialized");
        // The handles belong to AMReX's device.
        int device = -1;
#if defined(AMREX_USE_CUDA)
        AMREX_CUDA_SAFE_CALL(cudaGetDevice(&device));
#else
        AMREX_HIP_SAFE_CALL(hipGetDevice(&device));
#endif
        if (device != Gpu::Device::deviceId()) {
            amrex::Abort("Gpu::getCusparseHandle/getRocsparseHandle: this host thread uses"
                         " GPU " + std::to_string(device) + ", but AMReX uses GPU "
                         + std::to_string(Gpu::Device::deviceId()) + ". Set the device"
                         " of the thread to amrex::Gpu::Device::deviceId() first.");
        }

        static thread_local ThreadHandle t;
        Stream const stream = Gpu::gpuStream();
        if (!t.handle) {
            SpareHandle* spare = nullptr;
            {
                RegistryLock lock;
                if (!finalize_registered) {
                    // If this throws, nothing has changed.  If it aborts, the
                    // exit() skips the lock, see ~ThreadHandle.
                    amrex::ExecOnFinalize(destroy_sparse_handles);
                    finalize_registered = true;
                }
                spare = spare_head;
                if (spare) { spare_head = spare->next; }
            }
            SparseHandle handle{};
            if (spare) {
                handle = spare->handle;
                bool const same_stream = spare->stream == stream;
                delete spare;
                // Work that the thread before queued with it on another
                // stream may still be running.
                if (!same_stream) { Gpu::Device::synchronize(); }
            } else {
#if defined(AMREX_USE_CUDA)
                AMREX_CUSPARSE_SAFE_CALL(cusparseCreate(&handle));
#else
                AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
#endif
            }
            RegistryLock lock;
            t.link(handle);
        }
        // The handle is used by this host thread only.  Bind it to the current
        // stream, and reset the pointer mode that AMReX's calls rely on (alpha
        // and beta on the host), in case a caller changed it.
        t.stream = stream;
#if defined(AMREX_USE_CUDA)
        AMREX_CUSPARSE_SAFE_CALL(cusparseSetStream(t.handle, stream));
        AMREX_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(t.handle, CUSPARSE_POINTER_MODE_HOST));
#else
        AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(t.handle, stream));
        AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_pointer_mode(t.handle,
                                                             rocsparse_pointer_mode_host));
#endif
        return t.handle;
    }
}

#if defined(AMREX_USE_CUDA)
cusparseHandle_t getCusparseHandle () { return get_sparse_handle(); }
#else
rocsparse_handle getRocsparseHandle () { return get_sparse_handle(); }
#endif

}

#endif
