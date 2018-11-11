#include <AMReX_CudaAsyncArray.H>
#include <mutex>

#ifdef AMREX_USE_CUDA

extern "C" {
    void CUDART_CB amrex_asyncarray_delete (cudaStream_t stream, cudaError_t error, void* p)
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
