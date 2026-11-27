#include <AMReX_GpuSparse.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuError.H>

#include <mutex>
#include <string>
#include <utility>

// Creating a cuSPARSE or rocSPARSE handle allocates device resources and costs
// of the order of a millisecond, far more than a sparse matrix-vector product
// on a coarse multigrid level.  The handles are therefore kept, one per host
// thread: created on first use, released when the thread exits or in
// amrex::Finalize, whichever comes first.
//
// Nothing that can abort or allocate runs while sparse_handles_mutex is held:
// an abort (also from a failed allocation) can end in exit(), which runs this
// thread's ThreadHandle destructor, and that locks the mutex too.  The live
// handles are therefore kept in an intrusive list, and runtime calls under the
// lock only return their errors.

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

    // Work queued with a handle may still be running, also on streams that
    // AMReX does not know: synchronize the whole device.  Errors are returned.
    bool synchronize_device (int device) noexcept
    {
        // the device the handle was created on, in case the current one differs
        int current = -1;
#if defined(AMREX_USE_CUDA)
        bool ok = cudaGetDevice(&current) == cudaSuccess;
        bool const other = ok && current != device;
        if (other) { ok = cudaSetDevice(device) == cudaSuccess; }
        ok = ok && (cudaDeviceSynchronize() == cudaSuccess);
        if (other) { ok = (cudaSetDevice(current) == cudaSuccess) && ok; }
#else
        bool ok = hipGetDevice(&current) == hipSuccess;
        bool const other = ok && current != device;
        if (other) { ok = hipSetDevice(device) == hipSuccess; }
        ok = ok && (hipDeviceSynchronize() == hipSuccess);
        if (other) { ok = (hipSetDevice(current) == hipSuccess) && ok; }
#endif
        return ok;
    }

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
        int device = -1;
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
            SparseHandle h{};
            {
                std::lock_guard<std::mutex> lock(sparse_handles_mutex);
                if (handle) { h = unlink(); }
            }
            if (h) {
                (void) synchronize_device(device);
                (void) release(h);
            }
        }
    };

    void destroy_sparse_handles ()
    {
        bool ok = true;
        {
            // Held while releasing, so that a thread exiting meanwhile waits
            // and then finds its handle released.  Nothing here aborts.
            std::lock_guard<std::mutex> lock(sparse_handles_mutex);
            finalize_registered = false;
            if (live_head) {
                ok = synchronize_device(Gpu::Device::deviceId());
            }
            while (live_head) {
                // No other thread may make sparse calls during amrex::Finalize,
                // so the owners do not use their handles now.
                ok = release(live_head->unlink()) && ok;
            }
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ok, "Gpu::destroy_sparse_handles failed");
    }

    SparseHandle get_sparse_handle ()
    {
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
                std::lock_guard<std::mutex> lock(sparse_handles_mutex);
                t.link(handle, device);
                register_finalize = !std::exchange(finalize_registered, true);
            }
            if (register_finalize) {
                try {
                    amrex::ExecOnFinalize(destroy_sparse_handles);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(sparse_handles_mutex);
                    finalize_registered = false;
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
