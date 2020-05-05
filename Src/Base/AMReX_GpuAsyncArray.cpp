#include <AMReX_GpuAsyncArray.H>
#include <mutex>

#ifdef AMREX_USE_GPU

#ifdef __HIP_PLATFORM_HCC__
#define HIPRT_CB 
#endif

#if !defined(AMREX_USE_DPCPP)
extern "C" {
AMREX_HIP_OR_CUDA(
         void HIPRT_CB  amrex_asyncarray_delete ( hipStream_t stream,  hipError_t error, void* p),
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )
         void CUDART_CB amrex_asyncarray_delete (void* p))
#else
         void CUDART_CB amrex_asyncarray_delete (cudaStream_t stream, cudaError_t error, void* p))
#endif
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

#endif
