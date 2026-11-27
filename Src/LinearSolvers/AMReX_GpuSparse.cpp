#include <AMReX_GpuSparse.H>

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)

#include <AMReX.H>
#include <AMReX_GpuDevice.H>
#include <AMReX_GpuError.H>

#include <map>
#include <mutex>
#include <thread>

namespace amrex::Gpu {

namespace {

#if defined(AMREX_USE_CUDA)
    using SparseHandle = cusparseHandle_t;
#else
    using SparseHandle = rocsparse_handle;
#endif

    std::mutex sparse_handles_mutex;
    std::map<std::thread::id, SparseHandle> sparse_handles;
    bool sparse_handles_finalize_registered = false;

    void destroy_sparse_handles ()
    {
        std::lock_guard<std::mutex> lock(sparse_handles_mutex);
        // Work queued with a handle may still be running.
        Gpu::streamSynchronizeAll();
        for (auto const& [tid, handle] : sparse_handles) {
#if defined(AMREX_USE_CUDA)
            AMREX_CUSPARSE_SAFE_CALL(cusparseDestroy(handle));
#else
            AMREX_ROCSPARSE_SAFE_CALL(rocsparse_destroy_handle(handle));
#endif
        }
        sparse_handles.clear();
        // amrex::Initialize may be called again.
        sparse_handles_finalize_registered = false;
    }

    SparseHandle get_sparse_handle ()
    {
        SparseHandle handle;
        {
            std::lock_guard<std::mutex> lock(sparse_handles_mutex);
            auto it = sparse_handles.find(std::this_thread::get_id());
            if (it == sparse_handles.end()) {
#if defined(AMREX_USE_CUDA)
                AMREX_CUSPARSE_SAFE_CALL(cusparseCreate(&handle));
#else
                AMREX_ROCSPARSE_SAFE_CALL(rocsparse_create_handle(&handle));
#endif
                sparse_handles.emplace(std::this_thread::get_id(), handle);
                if (!sparse_handles_finalize_registered) {
                    amrex::ExecOnFinalize(destroy_sparse_handles);
                    sparse_handles_finalize_registered = true;
                }
            } else {
                handle = it->second;
            }
        }
        // The handle is used by this host thread only.
#if defined(AMREX_USE_CUDA)
        AMREX_CUSPARSE_SAFE_CALL(cusparseSetStream(handle, Gpu::gpuStream()));
#else
        AMREX_ROCSPARSE_SAFE_CALL(rocsparse_set_stream(handle, Gpu::gpuStream()));
#endif
        return handle;
    }
}

#if defined(AMREX_USE_CUDA)
cusparseHandle_t getCusparseHandle () { return get_sparse_handle(); }
#else
rocsparse_handle getRocsparseHandle () { return get_sparse_handle(); }
#endif

}

#endif
