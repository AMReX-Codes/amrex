#include <AMReX_GpuSparse.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuError.H>

#include <mutex>
#include <set>
#include <string>
#include <utility>

// Creating a cuSPARSE or rocSPARSE handle allocates device resources and costs
// of the order of a millisecond, far more than a sparse matrix-vector product
// on a coarse multigrid level.  The handles are therefore kept, one per host
// thread: created on first use, released when the thread exits or in
// amrex::Finalize, whichever comes first.
//
// Nothing that can abort runs while sparse_handles_mutex is held: an abort
// can end in exit(), which runs this thread's ThreadHandle destructor, and
// that locks the mutex too.

namespace amrex::Gpu {

namespace {

#if defined(AMREX_USE_CUDA)
    using SparseHandle = cusparseHandle_t;
#else
    using SparseHandle = rocsparse_handle;
#endif

    struct ThreadHandle;

    // Guards live_handles and finalize_registered, and makes thread exit and
    // amrex::Finalize release each handle exactly once.
    std::mutex sparse_handles_mutex;
    // The handles created since amrex::Initialize that are not released yet.
    std::set<ThreadHandle*> live_handles;
    bool finalize_registered = false;

    // Release a handle.  It may still be bound to a stream that no longer
    // exists (e.g., an external stream), so it is bound to the default stream
    // first; the caller has synchronized the device.  Errors are returned, not
    // aborted on, so that this can run when a thread exits.
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

        ThreadHandle () = default;
        ThreadHandle (ThreadHandle const&) = delete;
        ThreadHandle (ThreadHandle&&) = delete;
        ThreadHandle& operator= (ThreadHandle const&) = delete;
        ThreadHandle& operator= (ThreadHandle&&) = delete;

        // The thread exits: release its handle, unless amrex::Finalize did.
        ~ThreadHandle ()
        {
            SparseHandle h{};
            {
                std::lock_guard<std::mutex> lock(sparse_handles_mutex);
                if (live_handles.erase(this) > 0) {
                    h = std::exchange(handle, SparseHandle{});
                }
            }
            if (h) {
                // work queued with the handle may still be running
#if defined(AMREX_USE_CUDA)
                (void) cudaDeviceSynchronize();
#else
                (void) hipDeviceSynchronize();
#endif
                (void) release(h);
            }
        }
    };

    void destroy_sparse_handles ()
    {
        std::set<ThreadHandle*> handles;
        {
            std::lock_guard<std::mutex> lock(sparse_handles_mutex);
            // Take the handles out first, so that nothing is released twice
            // and amrex::Initialize can start afresh even if a release fails.
            handles = std::exchange(live_handles, {});
            finalize_registered = false;
        }
        if (handles.empty()) { return; }
        // Work queued with a handle may still be running, also on streams that
        // AMReX does not know, so synchronize the whole device.
        Gpu::Device::synchronize();
        bool ok = true;
        for (auto* t : handles) {
            // No other thread may make sparse calls during amrex::Finalize,
            // so the owners do not touch t->handle now.
            ok = release(std::exchange(t->handle, SparseHandle{})) && ok;
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
                t.handle = handle;
                live_handles.insert(&t);
                register_finalize = !std::exchange(finalize_registered, true);
            }
            if (register_finalize) {
                amrex::ExecOnFinalize(destroy_sparse_handles);
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
