#ifndef AMREX_GPU_LAUNCH_MACROS_G_H_
#define AMREX_GPU_LAUNCH_MACROS_G_H_

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(BX,IV,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long amrex_i_npts = BX.numPts(); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_npts); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const amrex_i_i : amrex::Gpu::Range(amrex_i_npts)) { \
                const auto IV = BX.atOffset(amrex_i_i); \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        long amrex_i_npts = BX.numPts(); \
        for (auto const amrex_i_i : amrex::Gpu::Range(amrex_i_npts)) { \
            const auto IV = BX.atOffset(amrex_i_i); \
            block \
        } \
    }

#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_BOXIV(BX,IV,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long amrex_i_npts = BX.numPts(); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_npts);   \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const amrex_i_i : amrex::Gpu::Range(amrex_i_npts)) { \
                const auto IV = BX.atOffset(amrex_i_i); \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_BOXIV: cannot call device function from host"); \
    }

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(TN);             \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const TI : amrex::Gpu::Range(TN)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI : amrex::Gpu::Range(TN)) { \
            block \
        } \
    }

// two fused launches
#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(TN2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(TN2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Gpu::Range(TN1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Gpu::Range(TN2)) { \
            block2 \
        } \
    }

// three fused launches
#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(TN2); \
        const auto amrex_i_ec3 = amrex::Gpu::ExecutionConfig(TN3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(TN2)) { \
                        block2 \
                    } \
                    break; \
            case 2: for (auto const TI3 : amrex::Gpu::Range(TN3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Gpu::Range(TN1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Gpu::Range(TN2)) { \
            block2 \
        } \
        for (auto const TI3 : amrex::Gpu::Range(TN3)) { \
            block3 \
        } \
    }

#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        auto amrex_i_ec = amrex::Gpu::ExecutionConfig(TN); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (auto const TI : amrex::Gpu::Range(TN)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE: cannot call device function from host"); \
    }

// two fused launches
#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(TN2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(TN2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }

// three fused launches
#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Gpu::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Gpu::ExecutionConfig(TN2); \
        const auto amrex_i_ec3 = amrex::Gpu::ExecutionConfig(TN3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Gpu::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Gpu::Range(TN2)) { \
                        block2 \
                    } \
            case 2: for (auto const TI3 : amrex::Gpu::Range(TN3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_BOX(bbb,tbb,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(bbb); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / amrex::Gpu::Device::warp_size; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*amrex::Gpu::Device::warp_size; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_GPU_NCELLS_PER_THREAD*amrex::Gpu::Device::warp_size; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += amrex::Gpu::Device::warp_size) \
            { \
                amrex::Box tbb = amrex::Gpu::getThreadBox(bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_ASYNC(bbb,tbb,sync_var,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        auto amrex_i_ec = amrex::Gpu::ExecutionConfig(bbb); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / amrex::Gpu::Device::warp_size; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*amrex::Gpu::Device::warp_size; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_GPU_NCELLS_PER_THREAD*amrex::Gpu::Device::warp_size; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += amrex::Gpu::Device::warp_size) \
            { \
                amrex::Box tbb = amrex::Gpu::getThreadBox(bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        ++sync_var; \
        amrex::Gpu::Device::set_stream_index(sync_var); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const amrex::Box& tbb = bbb; \
        block \
    }  

#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_BOX(bbb,tbb,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(bbb); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / amrex::Gpu::Device::warp_size; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*amrex::Gpu::Device::warp_size; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_GPU_NCELLS_PER_THREAD*amrex::Gpu::Device::warp_size; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += amrex::Gpu::Device::warp_size) \
            { \
                amrex::Box tbb = amrex::Gpu::getThreadBox(bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA: cannot call device function from host"); \
    }

#define AMREX_GPU_LAUNCH_HOST_DEVICE_XYZ(bbx,bby,bbz,tbx,tby,tbz,blockx,blocky,blockz) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long max_pts = std::max(bbx.numPts(), std::max(bby.numPts(), bbz.numPts())); \
        const auto amrex_i_ec = Gpu::ExecutionConfig(max_pts); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / amrex::Gpu::Device::warp_size; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*amrex::Gpu::Device::warp_size; \
            const long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_GPU_NCELLS_PER_THREAD*amrex::Gpu::Device::warp_size; \
            long amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += amrex::Gpu::Device::warp_size) \
                { \
                    amrex::Box tbx = amrex::Gpu::getThreadBox(bbx, amrex_i_loc_offset); \
                    if (tbx.ok()) blockx \
                } \
            amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += amrex::Gpu::Device::warp_size) \
                { \
                    amrex::Box tby = amrex::Gpu::getThreadBox(bby, amrex_i_loc_offset); \
                    if (tby.ok()) blocky \
                } \
            amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_GPU_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += amrex::Gpu::Device::warp_size) \
                { \
                    amrex::Box tbz = amrex::Gpu::getThreadBox(bbz, amrex_i_loc_offset); \
                    if (tbz.ok()) blockz \
                } \
         }); \
    } \
    else \
    { \
        const amrex::Box& tbx = bbx; \
        blockx \
        const amrex::Box& tby = bby; \
        blocky \
        const amrex::Box& tbz = bbz; \
        blockz \
    }

#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_NOBOX(bbb,tbb,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            const auto amrex_i_lo = bbb.loVect3d(); \
            const auto amrex_i_hi = bbb.hiVect3d(); \
            for (int amrex_i_k = amrex_i_lo[2] + blockIdx.z * blockDim.z + threadIdx.z; amrex_i_k <= amrex_i_hi[2]; amrex_i_k += blockDim.z * gridDim.z) { \
            for (int amrex_i_j = amrex_i_lo[1] + blockIdx.y * blockDim.y + threadIdx.y; amrex_i_j <= amrex_i_hi[1]; amrex_i_j += blockDim.y * gridDim.y) { \
            for (int amrex_i_i = amrex_i_lo[0] + blockIdx.x * blockDim.x + threadIdx.x; amrex_i_i <= amrex_i_hi[0]; amrex_i_i += blockDim.x * gridDim.x) { \
                amrex::Box tbb(IntVect(AMREX_D_DECL(amrex_i_i,amrex_i_j,amrex_i_k)), \
                               IntVect(AMREX_D_DECL(amrex_i_i,amrex_i_j,amrex_i_k)), \
                               bbb.type()); \
                block \
            }}} \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_GPU_LAUNCH_HOST_DEVICE(strategy, ...) \
    { \
      if (amrex::Gpu::inLaunchRegion()) \
      { \
         const auto amrex_i_ec = strategy;                                           \
         amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>(__VA_ARGS__); \
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
         amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>(__VA_ARGS__); \
         AMREX_GPU_ERROR_CHECK(); \
      } \
      else \
      { \
         amrex::Abort("AMREX_GPU_LAUNCH_DEVICE: cannot call device function from host"); \
      } \
    }

// Cannot respect Gpu::inLaunchRegion because function must be __global__.
#define AMREX_GPU_LAUNCH_GLOBAL(strategy, function, ...) \
    { \
        const auto amrex_i_ec = strategy;                                             \
        function<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>(__VA_ARGS__); \
        AMREX_GPU_ERROR_CHECK();                                               \
    }

// FOR_1D

#define AMREX_GPU_HOST_DEVICE_FOR_1D(n,i,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(n); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int i = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     i < n; i += amrex_i_stride) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        AMREX_PRAGMA_SIMD \
        for (auto i = decltype(n){0}; i < n; ++i) block     \
    }

#define AMREX_GPU_DEVICE_FOR_1D(n,i,block) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(n); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int i = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     i < n; i += amrex_i_stride) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    }

#define AMREX_GPU_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_DEVICE_FOR_1D(__VA_ARGS__)

// FOR_3D

#define AMREX_GPU_HOST_DEVICE_FOR_3D(box,i,j,k,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        int amrex_i_ncells = box.numPts();          \
        const auto amrex_i_lo  = amrex::lbound(box); \
        const auto amrex_i_len = amrex::length(box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
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
        const auto amrex_i_lo = amrex::lbound(box); \
        const auto amrex_i_hi = amrex::ubound(box); \
        for (int k = amrex_i_lo.z; k <= amrex_i_hi.z; ++k) { \
        for (int j = amrex_i_lo.y; j <= amrex_i_hi.y; ++j) { \
        AMREX_PRAGMA_SIMD \
        for (int i = amrex_i_lo.x; i <= amrex_i_hi.x; ++i) { \
            block \
        }}} \
    }

#define AMREX_GPU_DEVICE_FOR_3D(box,i,j,k,block) \
    { \
        int amrex_i_ncells = box.numPts(); \
        const auto amrex_i_lo  = amrex::lbound(box); \
        const auto amrex_i_len = amrex::length(box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
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
    }

#define AMREX_GPU_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_DEVICE_FOR_3D(__VA_ARGS__)

// FOR_4D

#define AMREX_GPU_HOST_DEVICE_FOR_4D(box,ancomp,i,j,k,n,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        int amrex_i_nc = ancomp; \
        int amrex_i_ncells = box.numPts();          \
        const auto amrex_i_lo  = amrex::lbound(box); \
        const auto amrex_i_len = amrex::length(box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int amrex_i_icell = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     amrex_i_icell < amrex_i_ncells; amrex_i_icell += amrex_i_stride) { \
                int k =  amrex_i_icell /   (amrex_i_len.x*amrex_i_len.y); \
                int j = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) /   amrex_i_len.x; \
                int i = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) - j*amrex_i_len.x; \
                i += amrex_i_lo.x; \
                j += amrex_i_lo.y; \
                k += amrex_i_lo.z; \
                for (int n = 0; n < amrex_i_nc; ++n) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const auto amrex_i_lo = amrex::lbound(box); \
        const auto amrex_i_hi = amrex::ubound(box); \
        for (int n = 0; n < ancomp; ++n) { \
        for (int k = amrex_i_lo.z; k <= amrex_i_hi.z; ++k) { \
        for (int j = amrex_i_lo.y; j <= amrex_i_hi.y; ++j) { \
        AMREX_PRAGMA_SIMD \
        for (int i = amrex_i_lo.x; i <= amrex_i_hi.x; ++i) { \
            block \
        }}}} \
    }

#define AMREX_GPU_DEVICE_FOR_4D(box,ancomp,i,j,k,n,block) \
    { \
        int amrex_i_nc = ancomp; \
        int amrex_i_ncells = box.numPts(); \
        const auto amrex_i_lo  = amrex::lbound(box); \
        const auto amrex_i_len = amrex::length(box); \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_ncells); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Gpu::gpuStream()>>>( \
        [=] AMREX_GPU_DEVICE () noexcept { \
            for (int amrex_i_icell = blockDim.x*blockIdx.x+threadIdx.x, amrex_i_stride = blockDim.x*gridDim.x; \
                     amrex_i_icell < amrex_i_ncells; amrex_i_icell += amrex_i_stride) { \
                int k =  amrex_i_icell /   (amrex_i_len.x*amrex_i_len.y); \
                int j = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) /   amrex_i_len.x; \
                int i = (amrex_i_icell - k*(amrex_i_len.x*amrex_i_len.y)) - j*amrex_i_len.x; \
                i += amrex_i_lo.x; \
                j += amrex_i_lo.y; \
                k += amrex_i_lo.z; \
                for (int n = 0; n < amrex_i_nc; ++n) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    }

#define AMREX_GPU_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_HOST_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_HOST_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_HOST_DEVICE_FOR_4D(__VA_ARGS__)

#endif
