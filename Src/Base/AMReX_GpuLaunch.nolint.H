// Do not include this header anywhere other than AMReX_GpuLaunch.H.
// The purpose of this file is to avoid clang-tidy.

#define AMREX_WRONG_NUM_ARGS(...) static_assert(false,"Wrong number of arguments to macro")
#define AMREX_GET_MACRO(_1,_2,_3,_4,_5,_6,_7,_8,_9,NAME,...) NAME
#define AMREX_LAUNCH_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA(...) AMREX_GET_MACRO(__VA_ARGS__,\
                                                        AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE, \
                                                        AMREX_WRONG_NUM_ARGS, \
                                                        AMREX_WRONG_NUM_ARGS)(__VA_ARGS__)

#if (AMREX_SPACEDIM == 1)
#define AMREX_LAUNCH_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3)      AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE     (a1,a2,a3)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3) AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(a1,a2,a3)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(fl,a1,a2,a3,b1,b2,b3,c1,c2,c3) AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_FLAG(fl,a1,a2,a3)
#elif (AMREX_SPACEDIM == 2)
#define AMREX_LAUNCH_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3)      AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2     (a1,a2,a3,b1,b2,b3)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(a1,a2,a3,b1,b2,b3,c1,c2,c3) AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(a1,a2,a3,b1,b2,b3)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(fl,a1,a2,a3,b1,b2,b3,c1,c2,c3) AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2_FLAG(fl,a1,a2,a3,b1,b2,b3)
#elif (AMREX_SPACEDIM == 3)
#define AMREX_LAUNCH_DEVICE_LAMBDA_DIM(...)      AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_3     (__VA_ARGS__)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM(...) AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(__VA_ARGS__)
#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_DIM_FLAG(...) AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3_FLAG(__VA_ARGS__)
#endif

#define AMREX_FOR_1D(...)      AMREX_GPU_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_FOR_3D(...)      AMREX_GPU_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_FOR_4D(...)      AMREX_GPU_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_PARALLEL_FOR_1D(...)      AMREX_GPU_DEVICE_PARALLEL_FOR_1D(__VA_ARGS__)
#define AMREX_PARALLEL_FOR_3D(...)      AMREX_GPU_DEVICE_PARALLEL_FOR_3D(__VA_ARGS__)
#define AMREX_PARALLEL_FOR_4D(...)      AMREX_GPU_DEVICE_PARALLEL_FOR_4D(__VA_ARGS__)

#define AMREX_HOST_DEVICE_FOR_1D(...) AMREX_GPU_HOST_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_FOR_3D(...) AMREX_GPU_HOST_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_FOR_4D(...) AMREX_GPU_HOST_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_1D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_3D(__VA_ARGS__)
#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_4D(__VA_ARGS__)

#ifdef AMREX_USE_GPU

#ifndef AMREX_USE_SYCL

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept \
            block \
        ); \
    } \
    else { \
        AMREX_PRAGMA_SIMD \
        for (amrex_i_inttype i = 0; i < n; ++i) { \
            block \
        } \
    }}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopConcurrentOnCpu(box, [=] (int i, int j, int k) noexcept \
            block \
        ); \
    }

#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopConcurrentOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
            block \
        ); \
    }

#define AMREX_HOST_DEVICE_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept \
            block \
        ); \
    } \
    else { \
        for (amrex_i_inttype i = 0; i < n; ++i) { \
            block \
        } \
    }}

#define AMREX_HOST_DEVICE_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopOnCpu(box, [=] (int i, int j, int k) noexcept \
            block \
        ); \
    }

#define AMREX_HOST_DEVICE_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::LoopOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
            block \
        ); \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(where_to_run,box,tbox,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(box,tbox,block); \
    } else { \
        auto tbox = box; \
        block; \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_FLAG(where_to_run,bx1,tbx1,block1) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(bx1,tbx1,block1); \
    } else { \
        auto tbx1 = bx1; \
        block1; \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2_FLAG(where_to_run,bx1,tbx1,block1,bx2,tbx2,block2) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(bx1,tbx1,block1,bx2,tbx2,block2); \
    } else { \
        auto tbx1 = bx1; \
        auto tbx2 = bx2; \
        block1; \
        block2; \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3_FLAG(where_to_run,bx1,tbx1,block1,bx2,tbx2,block2,bx3,tbx3,block3) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(bx1,tbx1,block1,bx2,tbx2,block2,bx3,tbx3,block3); \
    } else { \
        auto tbx1 = bx1; \
        auto tbx2 = bx2; \
        auto tbx3 = bx3; \
        block1; \
        block2; \
        block3; \
    }

#else
// xxxxx SYCL todo: host disabled in host device

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#define AMREX_HOST_DEVICE_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(n, [=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }}

#define AMREX_HOST_DEVICE_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#define AMREX_HOST_DEVICE_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        amrex::ParallelFor(box, nc, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept \
            block \
        ); \
    } \
    else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(where_to_run,box,tbox,block) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(box,tbox,block); \
    } else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_FLAG(where_to_run,bx1,tbx1,block1) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(bx1,tbx1,block1); \
    } else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2_FLAG(where_to_run,bx1,tbx1,block1,bx2,tbx2,block2) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(bx1,tbx1,block1,bx2,tbx2,block2); \
    } else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3_FLAG(where_to_run,bx1,tbx1,block1,bx2,tbx2,block2,bx3,tbx3,block3) \
    if ((where_to_run == RunOn::Device) && (Gpu::inLaunchRegion())) \
    { \
        AMREX_LAUNCH_DEVICE_LAMBDA(bx1,tbx1,block1,bx2,tbx2,block2,bx3,tbx3,block3); \
    } else { \
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile"); \
    }

#endif

#else

#define AMREX_HOST_DEVICE_PARALLEL_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    amrex::ignore_unused(where_to_run); \
    AMREX_PRAGMA_SIMD \
    for (amrex_i_inttype i = 0; i < n; ++i) { \
        block \
    }}

#define AMREX_HOST_DEVICE_PARALLEL_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopConcurrentOnCpu(box, [=] (int i, int j, int k) noexcept \
        block \
    );

#define AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopConcurrentOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
        block \
    );

#define AMREX_HOST_DEVICE_FOR_1D_FLAG(where_to_run,n,i,block) \
    {  using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    amrex::ignore_unused(where_to_run); \
    for (amrex_i_inttype i = 0; i < n; ++i) { \
        block \
    }}

#define AMREX_HOST_DEVICE_FOR_3D_FLAG(where_to_run,box,i,j,k,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopOnCpu(box, [=] (int i, int j, int k) noexcept \
        block \
    );

#define AMREX_HOST_DEVICE_FOR_4D_FLAG(where_to_run,box,nc,i,j,k,n,block) \
    amrex::ignore_unused(where_to_run); \
    amrex::LoopOnCpu(box, nc, [=] (int i, int j, int k, int n) noexcept \
        block \
    );

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(where_to_run,box,tbox,block) \
    amrex::ignore_unused(where_to_run); \
    { \
        auto tbox = box; \
        block; \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_FLAG(where_to_run,bx1,tbx1,block1) \
    amrex::ignore_unused(where_to_run); \
    { \
        auto tbx1 = bx1; \
        block1; \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2_FLAG(where_to_run,bx1,tbx1,block1,bx2,tbx2,block2) \
    amrex::ignore_unused(where_to_run); \
    { \
        auto tbx1 = bx1; \
        auto tbx2 = bx2; \
        block1; \
        block2; \
    }

#define AMREX_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3_FLAG(where_to_run,bx1,tbx1,block1,bx2,tbx2,block2,bx3,tbx3,block3) \
    amrex::ignore_unused(where_to_run); \
    { \
        auto tbx1 = bx1; \
        auto tbx2 = bx2; \
        auto tbx3 = bx3; \
        block1; \
        block2; \
        block3; \
    }

#endif
