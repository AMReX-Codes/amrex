#include <AMReX_GpuSparse.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuError.H>

#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

// Creating a cuSPARSE or rocSPARSE handle allocates device resources and costs
// of the order of a millisecond, far more than a sparse matrix-vector product
// on a coarse multigrid level.  The handles are therefore kept, one per host
// thread and stream of AMReX's pool, which live until amrex::Finalize, plus
// one per host thread for all external streams; they are created on first
// use.  A handle stays on its stream, so the work queued with it stays in
// order; when it is found bound to another stream (a new external stream, or
// a caller bound it), the device is synchronized first.  When a thread exits,
// its handles become spares for the next thread that needs one for the same
// stream; amrex::Finalize releases all.
//
// A thread that exits makes no calls into the GPU runtime or library, which
// may have torn down their own per-thread state by then.  The lists of
// handles hold heap nodes only, never storage of a thread.  Nothing that can
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

    // A handle for one stream of AMReX's pool, or for the external streams
    struct Handle
    {
        SparseHandle handle{};  // none after amrex::Finalize released it
        Stream stream{};  // the pool stream, unless external
        bool external = false;
        // the external stream override it was last used in
        std::uint64_t external_id = 0;
        Handle* prev = nullptr;
        Handle* next = nullptr;

        [[nodiscard]] bool isFor (Stream s, bool ext) const noexcept
        {
            return external == ext && (ext || stream == s);
        }
    };

    // with sparse_handles_mutex held
    struct HandleList
    {
        Handle* head = nullptr;

        void push (Handle* h) noexcept
        {
            h->prev = nullptr;
            h->next = head;
            if (head) { head->prev = h; }
            head = h;
        }

        void remove (Handle* h) noexcept
        {
            if (h->prev) { h->prev->next = h->next; } else { head = h->next; }
            if (h->next) { h->next->prev = h->prev; }
            h->prev = h->next = nullptr;
        }

        Handle* take (Stream stream, bool external) noexcept
        {
            for (Handle* h = head; h; h = h->next) {
                if (h->isFor(stream, external)) {
                    remove(h);
                    return h;
                }
            }
            return nullptr;
        }
    };

    // Guards the lists and finalize_registered, and makes thread exit and
    // amrex::Finalize hand over or release each handle exactly once.
    std::mutex sparse_handles_mutex;
    HandleList live;   // the handles of running threads
    HandleList spare;  // the handles of threads that exited
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

    // The handles of a host thread, one per pool stream it used, and one for
    // the external streams.
    struct ThreadHandles
    {
        std::vector<Handle*> handles;

        ThreadHandles () = default;
        ThreadHandles (ThreadHandles const&) = delete;
        ThreadHandles (ThreadHandles&&) = delete;
        ThreadHandles& operator= (ThreadHandles const&) = delete;
        ThreadHandles& operator= (ThreadHandles&&) = delete;

        // The thread exits: its handles become spares, unless amrex::Finalize
        // has released them.  No GPU runtime or library calls here.
        ~ThreadHandles ()
        {
            // exit() while this thread holds the mutex: the process ends
            if (this_thread_in_registry) { return; }
            RegistryLock lock;
            for (Handle* h : handles) {
                if (h->handle) {
                    live.remove(h);
                    spare.push(h);
                } else {
                    delete h;
                }
            }
        }
    };

    void destroy_sparse_handles ()
    {
        bool ok = true;
        Handle* spares = nullptr;
        {
            // Held while releasing, so that a thread exiting meanwhile waits
            // and then finds its handles released.  Nothing here aborts.
            RegistryLock lock;
            finalize_registered = false;
            spares = std::exchange(spare.head, nullptr);
            if (live.head || spares) {
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
                // No other thread may make sparse calls during amrex::Finalize.
                // Their threads keep the nodes of live handles.  Without the
                // sync, work may still use the handles: leave them.
                bool const synced = ok;
                while (Handle* h = live.head) {
                    live.remove(h);
                    SparseHandle const hh = std::exchange(h->handle, SparseHandle{});
                    if (synced) { ok = release(hh) && ok; }
                }
                for (Handle* h = spares; synced && h; h = h->next) {
                    ok = release(h->handle) && ok;
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
        // a warning, so that the rest of amrex::Finalize still runs
        if (!ok) {
            amrex::Warning("Gpu::destroy_sparse_handles: synchronizing the device or "
                           "releasing the sparse library handles failed");
        }
    }

    SparseHandle get_sparse_handle ()
    {
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

        static thread_local ThreadHandles t;
        Stream const stream = Gpu::gpuStream();
        bool const external = Gpu::Device::usingExternalStream();
        Handle* h = nullptr;
        for (auto it = t.handles.begin(); it != t.handles.end();) {
            Handle* c = *it;
            if (c->isFor(stream, external)) {
                h = c;
            } else if (!c->handle) {
                // released by an earlier amrex::Finalize, in no list
                delete c;
                it = t.handles.erase(it);
                continue;
            }
            ++it;
        }
        if (!h) {
            auto node = std::make_unique<Handle>();
            node->stream = external ? Stream{} : stream;
            node->external = external;
            t.handles.push_back(node.get());
            h = node.release();
        }

        bool created = false;
        if (!h->handle) {
            Handle* reused = nullptr;
            {
                RegistryLock lock;
                if (!finalize_registered) {
                    // If this throws, nothing has changed.  If it aborts, the
                    // exit() skips the lock, see ~ThreadHandles.
                    amrex::ExecOnFinalize(destroy_sparse_handles);
                    finalize_registered = true;
                }
                reused = spare.take(stream, external);
            }
            SparseHandle handle{};
            if (reused) {
                handle = reused->handle;
                delete reused;
            } else {
#if defined(AMREX_USE_CUDA)
                AMREX_CUSPARSE_SAFE_CALL(cusparseCreate(&handle));
#else
                AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
#endif
                created = true;
            }
            RegistryLock lock;
            h->handle = handle;
            live.push(h);
        }

        // The handle may have work queued on another stream: one a caller
        // bound it to (see the header), or the external stream of another
        // override, which is told by its id since a stream address can be
        // reused.  Wait for that work before using the handle here.
        std::uint64_t const external_id = Gpu::Device::externalStreamId();
        if (!created) {
            Stream bound{};
#if defined(AMREX_USE_CUDA)
            AMREX_CUSPARSE_SAFE_CALL(cusparseGetStream(h->handle, &bound));
#else
            AMREX_ROCSPARSE_SAFE_CALL(rocsparse_get_stream(h->handle, &bound));
#endif
            if (bound != stream || h->external_id != external_id) {
                Gpu::Device::synchronize();
            }
        }
        h->external_id = external_id;
        // Bind it to its stream, and reset the pointer mode that AMReX's calls
        // rely on (alpha and beta on the host), in case a caller changed it.
#if defined(AMREX_USE_CUDA)
        AMREX_CUSPARSE_SAFE_CALL(cusparseSetStream(h->handle, stream));
        AMREX_CUSPARSE_SAFE_CALL(cusparseSetPointerMode(h->handle, CUSPARSE_POINTER_MODE_HOST));
#else
        AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(h->handle, stream));
        AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_pointer_mode(h->handle,
                                                             rocsparse_pointer_mode_host));
#endif
        return h->handle;
    }
}

#if defined(AMREX_USE_CUDA)
cusparseHandle_t getCusparseHandle () { return get_sparse_handle(); }
#else
rocsparse_handle getRocsparseHandle () { return get_sparse_handle(); }
#endif

}

#endif
