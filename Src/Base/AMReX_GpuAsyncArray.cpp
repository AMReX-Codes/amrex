#include <AMReX_GpuAsyncArray.H>
#include <mutex>

#ifdef AMREX_USE_GPU

extern "C" {
AMREX_HIP_OR_CUDA(
    void  HIPRT_CB amrex_asyncarray_delete ( hipStream_t stream,  hipError_t error, void* p),
    void CUDART_CB amrex_asyncarray_delete (cudaStream_t stream, cudaError_t error, void* p))
    {
        void** pp = (void**)p;
        void* dp = pp[0];
        void* hp = pp[1];
        std::free(hp);
        std::free(p);
        amrex::The_Device_Arena()->free(dp);
    }
}

#endif
