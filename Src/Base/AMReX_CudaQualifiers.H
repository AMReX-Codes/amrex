#ifndef AMREX_CUDA_QUALIFIERS_H_
#define AMREX_CUDA_QUALIFIERS_H_

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

#define AMREX_CUDA_HOST __host__
#define AMREX_CUDA_DEVICE __device__
#define AMREX_CUDA_GLOBAL __global__
#define AMREX_CUDA_HOST_DEVICE __host__ __device__

#else

#define AMREX_CUDA_HOST
#define AMREX_CUDA_DEVICE
#define AMREX_CUDA_HOST_DEVICE
#define AMREX_CUDA_GLOBAL

#endif

#endif
