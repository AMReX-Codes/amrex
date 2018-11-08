#include <AMReX_CudaAsyncArray.H>
#include <mutex>

#ifdef AMREX_USE_CUDA

namespace {
    std::mutex asyncarray_callback_mutex;
}

extern "C" {
    void CUDART_CB amrex_asyncarray_delete (cudaStream_t stream, cudaError_t error, void* p)
    {
        void** pp = (void**)p;
        void* hp = pp[0];
        void* dp = pp[1];
        std::free(hp);
        std::free(p);
        std::lock_guard<std::mutex> guard(asyncarray_callback_mutex);
        amrex::The_Device_Arena()->free(dp);
    }
}

#endif
