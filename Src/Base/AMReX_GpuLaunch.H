#ifndef AMREX_GPU_LAUNCH_H_
#define AMREX_GPU_LAUNCH_H_

#include <AMReX_CudaLaunch.H>

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(...) AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(__VA_ARGS__)
#define AMREX_LAUNCH_DEVICE_LAMBDA_BOXIV(...) AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_BOXIV(__VA_ARGS__)


#define AMREX_WRONG_NUM_ARGS(...) static_assert(false,"Wrong number of arguments to macro")
#define AMREX_GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,_9,NAME,...) NAME
#define AMREX_LAUNCH_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)

#endif
