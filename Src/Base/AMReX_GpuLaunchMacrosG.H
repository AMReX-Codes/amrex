#ifndef AMREX_GPU_LAUNCH_MACROS_G_H_
#define AMREX_GPU_LAUNCH_MACROS_G_H_

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(BX,IV,block) \
    { auto const& amrex_i_bx = BX; \
    if (amrex_i_bx.ok()) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long amrex_i_npts = amrex_i_bx.numPts(); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_npts); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const amrex_i_i : amrex::Gpu::Range(amrex_i_npts)) { \
                const auto IV = amrex_i_bx.atOffset(amrex_i_i); \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        long amrex_i_npts = amrex_i_bx.numPts(); \
        for (auto const amrex_i_i : amrex::Gpu::Range(amrex_i_npts)) { \
            const auto IV = amrex_i_bx.atOffset(amrex_i_i); \
            block \
        } \
    }}}

#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_BOXIV(BX,IV,block) \
    { auto const& amrex_i_bx = BX; \
    if (amrex_i_bx.ok()) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long amrex_i_npts = amrex_i_bx.numPts(); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_npts);   \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const amrex_i_i : amrex::Gpu::Range(amrex_i_npts)) { \
                const auto IV = amrex_i_bx.atOffset(amrex_i_i); \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_BOXIV: cannot call device function from host"); \
    }}}

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    { auto const& amrex_i_tn = TN; \
    if (!amrex::isEmpty(amrex_i_tn)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_tn); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const TI : amrex::Gpu::Range(amrex_i_tn)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI : amrex::Gpu::Range(amrex_i_tn)) { \
            block \
        } \
    }}}

// two fused launches
#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; \
    if (!amrex::isEmpty(amrex_i_tn1) or !amrex::isEmpty(amrex_i_tn2)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        AMREX_LAUNCH_KERNEL(amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
            block2 \
        } \
    }}}

// three fused launches
#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; auto const& amrex_i_tn3 = TN3; \
    if (!amrex::isEmpty(amrex_i_tn1) or !amrex::isEmpty(amrex_i_tn2) or !amrex::isEmpty(amrex_i_tn3)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        const auto amrex_i_ec3 = amrex::Gpu::ExecutionConfig(amrex_i_tn3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        AMREX_LAUNCH_KERNEL(amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
                    break; \
            case 2: for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
            block2 \
        } \
        for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3)) { \
            block3 \
        } \
    }}}

#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    { auto const& amrex_i_tn = TN; \
    if (!amrex::isEmpty(amrex_i_tn)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_tn); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const TI : amrex::Gpu::Range(amrex_i_tn)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE: cannot call device function from host"); \
    }}}

// two fused launches
#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; \
    if (!amrex::isEmpty(amrex_i_tn1) or !amrex::isEmpty(amrex_i_tn2)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        AMREX_LAUNCH_KERNEL(amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }}}

// three fused launches
#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    { auto const& amrex_i_tn1 = TN1; auto const& amrex_i_tn2 = TN2; auto const& amrex_i_tn3 = TN3; \
    if (!amrex::isEmpty(amrex_i_tn1) or !amrex::isEmpty(amrex_i_tn2) or !amrex::isEmpty(amrex_i_tn3)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(amrex_i_tn1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(amrex_i_tn2); \
        const auto amrex_i_ec3 = amrex::Gpu::ExecutionConfig(amrex_i_tn3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        AMREX_LAUNCH_KERNEL(amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
                        block2 \
                    } \
            case 2: for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }}}

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_BOX(bbb,tbb,block) \
    { auto const& amrex_i_bbb = bbb; \
    if (amrex_i_bbb.ok()) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_bbb); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            long amrex_i_numpts = amrex_i_bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / amrex::Gpu::Device::warp_size; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*amrex::Gpu::Device::warp_size; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_GPU_NCELLS_PER_THREAD*amrex::Gpu::Device::warp_size; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += amrex::Gpu::Device::warp_size) \
            { \
                amrex::Box tbb = amrex::Gpu::getThreadBox(amrex_i_bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const amrex::Box& tbb = amrex_i_bbb; \
        block \
    }}}

#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_BOX(bbb,tbb,block) \
    { auto const& amrex_i_bbb = bbb; \
    if (amrex_i_bbb.ok()) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_bbb); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            long amrex_i_numpts = amrex_i_bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / amrex::Gpu::Device::warp_size; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*amrex::Gpu::Device::warp_size; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_GPU_NCELLS_PER_THREAD*amrex::Gpu::Device::warp_size; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += amrex::Gpu::Device::warp_size) \
            { \
                amrex::Box tbb = amrex::Gpu::getThreadBox(amrex_i_bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA: cannot call device function from host"); \
    }}}

#define AMREX_GPU_LAUNCH_HOST_DEVICE_XYZ(bbx,bby,bbz,tbx,tby,tbz,blockx,blocky,blockz) \
    { auto const& amrex_i_bbx = bbx; auto const& amrex_i_bby = bby; auto const& amrex_i_bbz = bbz; \
    if (amrex_i_bbx.ok() or amrex_i_bby.ok() or amrex_i_bbz.ok()) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long max_pts = std::max(amrex_i_bbx.numPts(), std::max(amrex_i_bby.numPts(), amrex_i_bbz.numPts())); \
        const auto amrex_i_ec = Gpu::ExecutionConfig(max_pts); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / amrex::Gpu::Device::warp_size; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*amrex::Gpu::Device::warp_size; \
            const long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_GPU_NCELLS_PER_THREAD*amrex::Gpu::Device::warp_size; \
            long amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += amrex::Gpu::Device::warp_size) \
                { \
                    amrex::Box tbx = amrex::Gpu::getThreadBox(amrex_i_bbx, amrex_i_loc_offset); \
                    if (tbx.ok()) blockx \
                } \
            amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += amrex::Gpu::Device::warp_size) \
                { \
                    amrex::Box tby = amrex::Gpu::getThreadBox(amrex_i_bby, amrex_i_loc_offset); \
                    if (tby.ok()) blocky \
                } \
            amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += amrex::Gpu::Device::warp_size) \
                { \
                    amrex::Box tbz = amrex::Gpu::getThreadBox(amrex_i_bbz, amrex_i_loc_offset); \
                    if (tbz.ok()) blockz \
                } \
         }); \
    } \
    else \
    { \
        const amrex::Box& tbx = amrex_i_bbx; \
        blockx \
        const amrex::Box& tby = amrex_i_bby; \
        blocky \
        const amrex::Box& tbz = amrex_i_bbz; \
        blockz \
    }}}

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_NOBOX(bbb,tbb,block) \
    { auto const& amrex_i_bbb = bbb; \
    if (amrex_i_bbb.ok()) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            const auto amrex_i_lo = amrex_i_bbb.loVect3d(); \
            const auto amrex_i_hi = amrex_i_bbb.hiVect3d(); \
            for (int amrex_i_k = amrex_i_lo[2] + blockIdx.z * blockDim.z + threadIdx.z; amrex_i_k <= amrex_i_hi[2]; amrex_i_k += blockDim.z * gridDim.z) { \
            for (int amrex_i_j = amrex_i_lo[1] + blockIdx.y * blockDim.y + threadIdx.y; amrex_i_j <= amrex_i_hi[1]; amrex_i_j += blockDim.y * gridDim.y) { \
            for (int amrex_i_i = amrex_i_lo[0] + blockIdx.x * blockDim.x + threadIdx.x; amrex_i_i <= amrex_i_hi[0]; amrex_i_i += blockDim.x * gridDim.x) { \
                amrex::Box tamrex_i_bb(IntVect(AMREX_D_DECL(amrex_i_i,amrex_i_j,amrex_i_k)), \
                               IntVect(AMREX_D_DECL(amrex_i_i,amrex_i_j,amrex_i_k)), \
                               amrex_i_bbb.type()); \
                block \
            }}} \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const amrex::Box& tbb = amrex_i_bbb; \
        block \
    }}}

#define AMREX_GPU_LAUNCH_HOST_DEVICE(strategy, ...) \
    { \
      if (amrex::Gpu::inLaunchRegion()) \
      { \
         const auto amrex_i_ec = strategy;                                           \
         AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), __VA_ARGS__); \
         AMREX_GPU_ERROR_CHECK(); \
      } \
      else \
      { \
         amrex::launch_host(__VA_ARGS__); \
      } \
    }

#define AMREX_GPU_LAUNCH_DEVICE(strategy, ...) \
    { \
      if (amrex::Gpu::inLaunchRegion()) \
      { \
         const auto amrex_i_ec = strategy;                                   \
         AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), __VA_ARGS__); \
         AMREX_GPU_ERROR_CHECK(); \
      } \
      else \
      { \
         amrex::Abort("AMREX_GPU_LAUNCH_DEVICE: cannot call device function from host"); \
      } \
    }

// Cannot respect Gpu::inLaunchRegion because function must be __global__.
#ifdef AMREX_USE_HIP

#define AMREX_GPU_LAUNCH_GLOBAL(strategy, function, ...) \
    { \
        const auto amrex_i_ec = strategy; \
        hipLaunchKernelGGL(function, amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), __VA_ARGS__); \
        AMREX_GPU_ERROR_CHECK(); \
    }

#elif defined(AMREX_USE_CUDA)

#define AMREX_GPU_LAUNCH_GLOBAL(strategy, function, ...) \
    { \
        const auto amrex_i_ec = strategy; \
        function<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>(__VA_ARGS__); \
        AMREX_GPU_ERROR_CHECK(); \
    }

#endif
// FOR_1D

#define AMREX_GPU_HOST_DEVICE_FOR_1D(n,i,block) \
    { auto const& amrex_i_n = n; \
      using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if (!amrex::isEmpty(amrex_i_n)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_n); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (amrex_i_inttype i = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     i < amrex_i_n; i += amrex_i_stride) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        AMREX_PRAGMA_SIMD \
        for (amrex_i_inttype i = 0; i < amrex_i_n; ++i) block     \
    }}}

#define AMREX_GPU_DEVICE_FOR_1D(n,i,block) \
    { auto const& amrex_i_n = n; \
      using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if (!amrex::isEmpty(amrex_i_n)) { \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_n); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (amrex_i_inttype i = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     i < amrex_i_n; i += amrex_i_stride) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    }}}

#define AMREX_GPU_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_DEVICE_FOR_1D(__VA_ARGS__)

// FOR_3D

#define AMREX_GPU_HOST_DEVICE_FOR_3D(box,i,j,k,block) \
    { auto const& amrex_i_box = box; \
    if (!amrex::isEmpty(amrex_i_box)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        int amrex_i_ncells = amrex_i_box.numPts();          \
        const auto amrex_i_lo  = amrex::lbound(amrex_i_box); \
        const auto amrex_i_len = amrex::length(amrex_i_box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int amrex_i_icell = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     amrex_i_icell < amrex_i_ncells; amrex_i_icell += amrex_i_stride) { \
                int k =  amrex_i_icell /   (amrex_i_len.x*amrex_i_len.y); \
                int j = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) /   amrex_i_len.x; \
                int i = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) - j*amrex_i_len.x; \
                i += amrex_i_lo.x; \
                j += amrex_i_lo.y; \
                k += amrex_i_lo.z; \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const auto amrex_i_lo = amrex::lbound(amrex_i_box); \
        const auto amrex_i_hi = amrex::ubound(amrex_i_box); \
        for (int k = amrex_i_lo.z; k <= amrex_i_hi.z; ++k) { \
        for (int j = amrex_i_lo.y; j <= amrex_i_hi.y; ++j) { \
        AMREX_PRAGMA_SIMD \
        for (int i = amrex_i_lo.x; i <= amrex_i_hi.x; ++i) { \
            block \
        }}} \
    }}}

#define AMREX_GPU_DEVICE_FOR_3D(box,i,j,k,block) \
    { auto const& amrex_i_box = box; \
    if (!amrex::isEmpty(amrex_i_box)) { \
        int amrex_i_ncells = amrex_i_box.numPts(); \
        const auto amrex_i_lo  = amrex::lbound(amrex_i_box); \
        const auto amrex_i_len = amrex::length(amrex_i_box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int amrex_i_icell = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     amrex_i_icell < amrex_i_ncells; amrex_i_icell += amrex_i_stride) { \
                int k =  amrex_i_icell /   (amrex_i_len.x*amrex_i_len.y); \
                int j = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) /   amrex_i_len.x; \
                int i = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) - j*amrex_i_len.x; \
                i += amrex_i_lo.x; \
                j += amrex_i_lo.y; \
                k += amrex_i_lo.z; \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    }}

#define AMREX_GPU_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_DEVICE_FOR_3D(__VA_ARGS__)

// FOR_4D

#define AMREX_GPU_HOST_DEVICE_FOR_4D(box,ncomp,i,j,k,n,block) \
    { auto const& amrex_i_box = box; auto const& amrex_i_ncomp = ncomp; \
    if (!amrex::isEmpty(amrex_i_box)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        int amrex_i_ncells = amrex_i_box.numPts();          \
        const auto amrex_i_lo  = amrex::lbound(amrex_i_box); \
        const auto amrex_i_len = amrex::length(amrex_i_box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int amrex_i_icell = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     amrex_i_icell < amrex_i_ncells; amrex_i_icell += amrex_i_stride) { \
                int k =  amrex_i_icell /   (amrex_i_len.x*amrex_i_len.y); \
                int j = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) /   amrex_i_len.x; \
                int i = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) - j*amrex_i_len.x; \
                i += amrex_i_lo.x; \
                j += amrex_i_lo.y; \
                k += amrex_i_lo.z; \
                for (int n = 0; n < amrex_i_ncomp; ++n) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const auto amrex_i_lo = amrex::lbound(amrex_i_box); \
        const auto amrex_i_hi = amrex::ubound(amrex_i_box); \
        for (int n = 0; n < amrex_i_ncomp; ++n) { \
        for (int k = amrex_i_lo.z; k <= amrex_i_hi.z; ++k) { \
        for (int j = amrex_i_lo.y; j <= amrex_i_hi.y; ++j) { \
        AMREX_PRAGMA_SIMD \
        for (int i = amrex_i_lo.x; i <= amrex_i_hi.x; ++i) { \
            block \
        }}}} \
    }}}

#define AMREX_GPU_DEVICE_FOR_4D(box,ncomp,i,j,k,n,block) \
    { auto const& amrex_i_box = box; auto const& amrex_i_ncomp = ncomp; \
    if (!amrex::isEmpty(amrex_i_box)) { \
        int amrex_i_ncells = amrex_i_box.numPts(); \
        const auto amrex_i_lo  = amrex::lbound(amrex_i_box); \
        const auto amrex_i_len = amrex::length(amrex_i_box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        AMREX_LAUNCH_KERNEL(amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream(), \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int amrex_i_icell = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     amrex_i_icell < amrex_i_ncells; amrex_i_icell += amrex_i_stride) { \
                int k =  amrex_i_icell /   (amrex_i_len.x*amrex_i_len.y); \
                int j = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) /   amrex_i_len.x; \
                int i = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) - j*amrex_i_len.x; \
                i += amrex_i_lo.x; \
                j += amrex_i_lo.y; \
                k += amrex_i_lo.z; \
                for (int n = 0; n < amrex_i_ncomp; ++n) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    }}

#define AMREX_GPU_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_HOST_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_HOST_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_HOST_DEVICE_FOR_4D(__VA_ARGS__)

#endif
