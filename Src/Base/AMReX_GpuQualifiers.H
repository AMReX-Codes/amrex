#ifndef AMREX_GPU_QUALIFIERS_H_
#define AMREX_GPU_QUALIFIERS_H_

#if defined(AMREX_USE_CUDA)

#define AMREX_GPU_HOST __host__
#define AMREX_GPU_DEVICE __device__
#define AMREX_GPU_GLOBAL __global__
#define AMREX_GPU_HOST_DEVICE __host__ __device__
#define AMREX_GPU_DEVICE_MANAGED __device__ __managed__

#else

#define AMREX_GPU_HOST
#define AMREX_GPU_DEVICE
#define AMREX_GPU_GLOBAL
#define AMREX_GPU_HOST_DEVICE
#define AMREX_GPU_DEVICE_MANAGED

#endif

#endif
