#include <AMReX_GpuSparse.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuError.H>

#include <mutex>
#include <string>
#include <thread>
#include <utility>

// Creating a cuSPARSE or rocSPARSE handle allocates device resources and costs
// of the order of a millisecond, far more than a sparse matrix-vector product
// on a coarse multigrid level.  The handles are therefore kept, one per host
// thread: created on first use, released when the thread exits or in
// amrex::Finalize, whichever comes first.
//
// Nothing that can abort or allocate runs while sparse_handles_mutex is held,
// and a thread that exits while it holds the mutex (exit() from a signal
// handler) does not take it again.  The live handles are therefore kept in an
// intrusive list, and runtime calls under the lock only return their errors.

namespace amrex::Gpu {

namespace {

#if defined(AMREX_USE_CUDA)
    using SparseHandle = cusparseHandle_t;
#else
    using SparseHandle = rocsparse_handle;
#endif

    struct ThreadHandle;

    // Guards the list of live handles and finalize_registered, and makes thread
    // exit and amrex::Finalize release each handle exactly once.
    std::mutex sparse_handles_mutex;
    // The threads whose handles were created since amrex::Initialize and are
    // not released yet.
    ThreadHandle* live_head = nullptr;
    bool finalize_registered = false;

    // The thread that loaded AMReX, normally the main thread.  It exits only
    // with the process, also after an abort, and leaves its handle to the
    // process teardown instead of calling into the GPU runtime then.
    std::thread::id const loading_thread = std::this_thread::get_id();

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

    // Run f() with device current, then make the previous device current
    // again.  Errors are returned.
    template <typename F>
    bool on_device (int device, F const& f) noexcept
    {
        int current = -1;
#if defined(AMREX_USE_CUDA)
        bool ok = cudaGetDevice(&current) == cudaSuccess;
        bool const other = ok && current != device;
        if (other) { ok = cudaSetDevice(device) == cudaSuccess; }
        ok = ok && f();
        if (other) { ok = (cudaSetDevice(current) == cudaSuccess) && ok; }
#else
        bool ok = hipGetDevice(&current) == hipSuccess;
        bool const other = ok && current != device;
        if (other) { ok = hipSetDevice(device) == hipSuccess; }
        ok = ok && f();
        if (other) { ok = (hipSetDevice(current) == hipSuccess) && ok; }
#endif
        return ok;
    }

    // Release a handle, with its device current.  It may still be bound to a
    // stream that no longer exists (e.g., an external stream), so it is bound
    // to the default stream first.  Errors are returned.
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
        int device = -1;  // the device it was created on
        // in the list of live handles if handle is set
        ThreadHandle* prev = nullptr;
        ThreadHandle* next = nullptr;

        ThreadHandle () = default;
        ThreadHandle (ThreadHandle const&) = delete;
        ThreadHandle (ThreadHandle&&) = delete;
        ThreadHandle& operator= (ThreadHandle const&) = delete;
        ThreadHandle& operator= (ThreadHandle&&) = delete;

        // with sparse_handles_mutex held
        void link (SparseHandle h, int dev) noexcept
        {
            handle = h;
            device = dev;
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

        // The thread exits: release its handle, unless amrex::Finalize did.
        ~ThreadHandle ()
        {
            // exit() while this thread holds the mutex: the process ends
            if (this_thread_in_registry) { return; }
            SparseHandle h{};
            {
                RegistryLock lock;
                if (handle) { h = unlink(); }
            }
            if (h && std::this_thread::get_id() != loading_thread) {
                (void) on_device(device, [&] () noexcept {
#if defined(AMREX_USE_CUDA)
                    // Work queued with the handle may still run: cuSPARSE
                    // defers the release of a handle's device resources (as
                    // SpMV already relies on), so the device need not wait.
                    return release(h);
#else
                    // rocSPARSE may free the device memory of a handle
                    // right away, so wait for the work queued with it.
                    bool const ok = hipDeviceSynchronize() == hipSuccess;
                    return release(h) && ok;
#endif
                });
            }
        }
    };

    void destroy_sparse_handles ()
    {
        bool ok = true;
        {
            // Held while releasing, so that a thread exiting meanwhile waits
            // and then finds its handle released.  Nothing here aborts.
            RegistryLock lock;
            finalize_registered = false;
            if (live_head) {
                // All handles were made on AMReX's device.  Work queued with
                // them may still run, also on streams that AMReX does not
                // know: synchronize the whole device.
                ok = on_device(Gpu::Device::deviceId(), [&] () noexcept {
#if defined(AMREX_USE_CUDA)
                    bool all = cudaDeviceSynchronize() == cudaSuccess;
#else
                    bool all = hipDeviceSynchronize() == hipSuccess;
#endif
                    while (live_head) {
                        // No other thread may make sparse calls during
                        // amrex::Finalize, so the owners do not use them now.
                        all = release(live_head->unlink()) && all;
                    }
                    return all;
                });
            }
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ok, "Gpu::destroy_sparse_handles failed");
    }

    SparseHandle get_sparse_handle ()
    {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::Initialized(),
            "Gpu::getCusparseHandle/getRocsparseHandle: AMReX is not initialized");
        // The handle belongs to the current device: check that it is AMReX's.
        // This first runtime call also sets up the runtime's own per-thread
        // state before the ThreadHandle below exists.  Thread-local objects
        // are destroyed in the reverse order of their construction, so the
        // ThreadHandle destructor runs while the runtime still works on this
        // thread.
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
        if (!t.handle) {
            SparseHandle handle{};
#if defined(AMREX_USE_CUDA)
            AMREX_CUSPARSE_SAFE_CALL(cusparseCreate(&handle));
#else
            AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
#endif
            bool register_finalize = false;
            {
                RegistryLock lock;
                t.link(handle, device);
                register_finalize = !std::exchange(finalize_registered, true);
            }
            if (register_finalize) {
                try {
                    amrex::ExecOnFinalize(destroy_sparse_handles);
                } catch (...) {
                    // nothing would release the handle: do it now
                    SparseHandle h{};
                    {
                        RegistryLock lock;
                        finalize_registered = false;
                        h = t.unlink();
                    }
                    (void) release(h);
                    throw;
                }
            }
        }
        // The handle is used by this host thread only.  Bind it to the current
        // stream, and reset the pointer mode that AMReX's calls rely on (alpha
        // and beta on the host), in case a caller changed it.
#if defined(AMREX_USE_CUDA)
        AMREX_CUSPARSE_SAFE_CALL(cusparseSetStream(t.handle, Gpu::gpuStream()));
        AMREX_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(t.handle, CUSPARSE_POINTER_MODE_HOST));
#else
        AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(t.handle, Gpu::gpuStream()));
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
