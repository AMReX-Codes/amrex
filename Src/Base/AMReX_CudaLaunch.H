#ifndef AMREX_CUDA_LAUNCH_H_
#define AMREX_CUDA_LAUNCH_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_GpuControl.H>
#include <AMReX_GpuError.H>
#include <AMReX_CudaRange.H>
#include <AMReX_CudaDevice.H>
#include <AMReX_Box.H>
#include <cstddef>

#define AMREX_CUDA_NCELLS_PER_THREAD 2
#define AMREX_CUDA_WARP_SIZE 32
#define AMREX_CUDA_Y_STRIDE 1
#define AMREX_CUDA_Z_STRIDE 1

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)

// ************************************************
// CUDA versions

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(BX,IV,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long amrex_i_npts = BX.numPts(); \
        const auto amrex_i_ec = amrex::Cuda::ExecutionConfig(amrex_i_npts); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            for (auto const amrex_i_i : amrex::Cuda::Range(amrex_i_npts)) { \
                const auto IV = BX.atOffset(amrex_i_i); \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        long amrex_i_npts = BX.numPts(); \
        for (auto const amrex_i_i : amrex::Cuda::Range(amrex_i_npts)) { \
            const auto IV = BX.atOffset(amrex_i_i); \
            block \
        } \
    }

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_BOXIV(BX,IV,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long amrex_i_npts = BX.numPts(); \
        const auto amrex_i_ec = amrex::Cuda::ExecutionConfig(amrex_i_npts);   \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            for (auto const amrex_i_i : amrex::Cuda::Range(amrex_i_npts)) { \
                const auto IV = BX.atOffset(amrex_i_i); \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_BOXIV: cannot call device function from host"); \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Cuda::ExecutionConfig(TN);             \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            for (auto const TI : amrex::Cuda::Range(TN)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI : amrex::Cuda::Range(TN)) { \
            block \
        } \
    }

// two fused launches
#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Cuda::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Cuda::ExecutionConfig(TN2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
            block2 \
        } \
    }

// three fused launches
#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Cuda::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Cuda::ExecutionConfig(TN2); \
        const auto amrex_i_ec3 = amrex::Cuda::ExecutionConfig(TN3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
                        block2 \
                    } \
                    break; \
            case 2: for (auto const TI3 : amrex::Cuda::Range(TN3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
            block2 \
        } \
        for (auto const TI3 : amrex::Cuda::Range(TN3)) { \
            block3 \
        } \
    }

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        auto amrex_i_ec = amrex::Cuda::ExecutionConfig(TN); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            for (auto const TI : amrex::Cuda::Range(TN)) { \
                block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE: cannot call device function from host"); \
    }

// two fused launches
#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Cuda::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Cuda::ExecutionConfig(TN2); \
        dim3 amrex_i_nblocks = amrex::max(amrex_i_ec1.numBlocks.x, \
                                          amrex_i_ec2.numBlocks.x); \
        amrex_i_nblocks.y = 2; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
                        block2 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }

// three fused launches
#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec1 = amrex::Cuda::ExecutionConfig(TN1); \
        const auto amrex_i_ec2 = amrex::Cuda::ExecutionConfig(TN2); \
        const auto amrex_i_ec3 = amrex::Cuda::ExecutionConfig(TN3); \
        dim3 amrex_i_nblocks = amrex::max(amrex::max(amrex_i_ec1.numBlocks.x, \
                                                     amrex_i_ec2.numBlocks.x), \
                                                     amrex_i_ec3.numBlocks.x); \
        amrex_i_nblocks.y = 3; \
        amrex::launch_global<<<amrex_i_nblocks, amrex_i_ec1.numThreads, 0, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            switch (blockIdx.y) { \
            case 0: for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
                        block1 \
                    } \
                    break; \
            case 1: for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
                        block2 \
                    } \
            case 2: for (auto const TI3 : amrex::Cuda::Range(TN3)) { \
                        block3 \
                    } \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOX(bbb,tbb,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Cuda::ExecutionConfig(bbb); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / AMREX_CUDA_WARP_SIZE; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*AMREX_CUDA_WARP_SIZE; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_CUDA_NCELLS_PER_THREAD*AMREX_CUDA_WARP_SIZE; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_CUDA_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += AMREX_CUDA_WARP_SIZE) \
            { \
                amrex::Box tbb = amrex::Cuda::getThreadBox(bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_ASYNC(bbb,tbb,sync_var,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        auto amrex_i_ec = amrex::Cuda::ExecutionConfig(bbb); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / AMREX_CUDA_WARP_SIZE; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*AMREX_CUDA_WARP_SIZE; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_CUDA_NCELLS_PER_THREAD*AMREX_CUDA_WARP_SIZE; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_CUDA_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += AMREX_CUDA_WARP_SIZE) \
            { \
                amrex::Box tbb = amrex::Cuda::getThreadBox(bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        ++sync_var; \
        amrex::Cuda::Device::set_stream_index(sync_var); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        const amrex::Box& tbb = bbb; \
        block \
    }  

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_BOX(bbb,tbb,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Cuda::ExecutionConfig(bbb); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            long amrex_i_numpts = bbb.numPts(); \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / AMREX_CUDA_WARP_SIZE; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*AMREX_CUDA_WARP_SIZE; \
            long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_CUDA_NCELLS_PER_THREAD*AMREX_CUDA_WARP_SIZE; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_CUDA_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_offset += AMREX_CUDA_WARP_SIZE) \
            { \
                amrex::Box tbb = amrex::Cuda::getThreadBox(bbb, amrex_i_offset); \
                if (tbb.ok()) block \
            } \
        }); \
        AMREX_GPU_ERROR_CHECK(); \
    } \
    else { \
        amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE_LAMBDA: cannot call device function from host"); \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_XYZ(bbx,bby,bbz,tbx,tby,tbz,blockx,blocky,blockz) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        long max_pts = std::max(bbx.numPts(), std::max(bby.numPts(), bbz.numPts())); \
        const auto amrex_i_ec = Cuda::ExecutionConfig(max_pts); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
            long amrex_i_tid = blockDim.x*blockIdx.x + threadIdx.x; \
            long amrex_i_wid = amrex_i_tid / AMREX_CUDA_WARP_SIZE; \
            long amrex_i_lid = amrex_i_tid - amrex_i_wid*AMREX_CUDA_WARP_SIZE; \
            const long amrex_i_offset = amrex_i_lid + amrex_i_wid*AMREX_CUDA_NCELLS_PER_THREAD*AMREX_CUDA_WARP_SIZE; \
            long amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_CUDA_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += AMREX_CUDA_WARP_SIZE) \
                { \
                    amrex::Box tbx = amrex::Cuda::getThreadBox(bbx, amrex_i_loc_offset); \
                    if (tbx.ok()) blockx \
                } \
            amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_CUDA_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += AMREX_CUDA_WARP_SIZE) \
                { \
                    amrex::Box tby = amrex::Cuda::getThreadBox(bby, amrex_i_loc_offset); \
                    if (tby.ok()) blocky \
                } \
            amrex_i_loc_offset = amrex_i_offset; \
            for (long amrex_i_i = 0; amrex_i_i < AMREX_CUDA_NCELLS_PER_THREAD; ++amrex_i_i, amrex_i_loc_offset += AMREX_CUDA_WARP_SIZE) \
                { \
                    amrex::Box tbz = amrex::Cuda::getThreadBox(bbz, amrex_i_loc_offset); \
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
 

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_NOBOX(bbb,tbb,block) \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Cuda::ExecutionConfig(); \
        amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>( \
        [=] AMREX_GPU_DEVICE () { \
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

#define AMREX_CUDA_LAUNCH_HOST_DEVICE(strategy, ...) \
    { \
      if (amrex::Gpu::inLaunchRegion()) \
      { \
         const auto amrex_i_ec = strategy;                                           \
         amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>(__VA_ARGS__); \
         AMREX_GPU_ERROR_CHECK(); \
      } \
      else \
      { \
         amrex::launch_host(__VA_ARGS__); \
      } \
    }

#define AMREX_CUDA_LAUNCH_DEVICE(strategy, ...) \
    { \
      if (amrex::Gpu::inLaunchRegion()) \
      { \
         const auto amrex_i_ec = strategy;                                   \
         amrex::launch_global<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>(__VA_ARGS__); \
         AMREX_GPU_ERROR_CHECK(); \
      } \
      else \
      { \
         amrex::Abort("AMREX_CUDA_LAUNCH_DEVICE: cannot call device function from host"); \
      } \
    }

// Cannot respect Gpu::inLaunchRegion because function must be __global__.
#define AMREX_CUDA_LAUNCH_GLOBAL(strategy, function, ...) \
    { \
        const auto amrex_i_ec = strategy;                                             \
        function<<<amrex_i_ec.numBlocks, amrex_i_ec.numThreads, amrex_i_ec.sharedMem, amrex::Cuda::Device::cudaStream()>>>(__VA_ARGS__); \
        AMREX_GPU_ERROR_CHECK();                                               \
    }

#else

// ************************************************
// CPU versions

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(BX,IV,block) \
    { \
        long amrex_i_npts = BX.numPts(); \
        for (auto const amrex_i_i : amrex::Cuda::Range(amrex_i_npts)) { \
            const auto IV = BX.atOffset(amrex_i_i); \
            block \
        } \
    }

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_BOXIV(...) AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOXIV(__VA_ARGS__)

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    { \
        for (auto const TI : amrex::Cuda::Range(TN)) { \
            block \
        } \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(TN1,TI1,block1,TN2,TI2,block2) \
    { \
        for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
            block2 \
        } \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(TN1,TI1,block1,TN2,TI2,block2,TN3,TI3,block3) \
    { \
        for (auto const TI1 : amrex::Cuda::Range(TN1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Cuda::Range(TN2)) { \
            block2 \
        } \
        for (auto const TI3 : amrex::Cuda::Range(TN3)) { \
            block3 \
        } \
    }

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE(...) AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(__VA_ARGS__)
#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_2(...) AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_2(__VA_ARGS__)
#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_RANGE_3(...) AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_RANGE_3(__VA_ARGS__)

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOX(bbb,tbb,block) \
    { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_CUDA_LAUNCH_DEVICE_LAMBDA_BOX(...) AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_BOX(__VA_ARGS__)

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_LAMBDA_ASYNC(bbb,tbb,sync_id,block) \
    { \
        const amrex::Box& tbb = bbb; \
        block \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE_XYZ(bbx,bby,bbz,tbx,tby,tbz,blockx,blocky,blockz) \
    { \
        const amrex::Box& tbx = bbx; \
        blockx \
        const amrex::Box& tby = bby; \
        blocky \
        const amrex::Box& tbz = bbz; \
        blockz \
    }

#define AMREX_CUDA_LAUNCH_HOST_DEVICE(strategy, ...) amrex::launch_host(__VA_ARGS__);
#define AMREX_CUDA_LAUNCH_DEVICE(strategy, ...) amrex::launch_host(__VA_ARGS__);
#define AMREX_CUDA_LAUNCH_GLOBAL(strategy, function, ...) function(__VA_ARGS__);

#endif

namespace amrex {

// ************************************************
//  Variadic lambda function wrappers for C++ CUDA Kernel calls.

    template<class L>
    AMREX_GPU_GLOBAL void launch_global (L f0) { f0(); }

    template<class L, class... Lambdas>
    AMREX_GPU_GLOBAL void launch_global (L f0, Lambdas... fs) { f0(); call_device(fs...); }
    
    template<class L>
    AMREX_GPU_DEVICE void call_device (L f0) { f0(); }
    
    template<class L, class... Lambdas>
    AMREX_GPU_DEVICE void call_device (L f0, Lambdas... fs) { f0(); call_device(fs...); }
    
// CPU variation

    template<class L>
    void launch_host (L f0) { f0(); }
    
    template<class L, class... Lambdas>
    void launch_host (L f0, Lambdas... fs) { f0(); launch_host(fs...); }


    template <class T> class LayoutData;
    class FabArrayBase;

namespace Cuda {

// ************************************************

    struct ComponentBox {
        Box box;
        int ic;
        int nc;
    };

    struct GridSize {
        int numBlocks;
        int numThreads;
        int globalBlockId;
    };

// ************************************************

    void getGridSize (FabArrayBase const& fa, int ngrow, LayoutData<GridSize>& gs, int& ntotblocks);

    AMREX_GPU_HOST_DEVICE
    inline
    Box getThreadBox (const Box& bx, long offset)
    {
#if defined(AMREX_USE_CUDA) && defined(__CUDA_ARCH__)
        const auto len = bx.length3d();
        long k = offset / (len[0]*len[1]);
        long j = (offset - k*(len[0]*len[1])) / len[0];
        long i = (offset - k*(len[0]*len[1])) - j*len[0];
        IntVect iv{AMREX_D_DECL(static_cast<int>(i),
                                static_cast<int>(j),
                                static_cast<int>(k))};
        iv += bx.smallEnd();
        return (bx & Box(iv,iv,bx.type()));
#else
        return bx;
#endif
    }

// ************************************************

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__)
    struct ExecutionConfig {
        ExecutionConfig () {
            Cuda::Device::grid_stride_threads_and_blocks(numBlocks,numThreads);
        }
        ExecutionConfig (const Box& box) {
            // If we change this, we must make sure it doesn't break say getGridSize, which
            // assumes the decomposiiton is 1D.
            Cuda::Device::n_threads_and_blocks( ((box.numPts()+AMREX_CUDA_NCELLS_PER_THREAD-1)/AMREX_CUDA_NCELLS_PER_THREAD), numBlocks, numThreads );
#if 0
            Box b = amrex::surroundingNodes(box);
            b -= box.smallEnd();
            b.coarsen(IntVect(AMREX_D_DECL(1,AMREX_CUDA_Y_STRIDE,AMREX_CUDA_Z_STRIDE)));
            Cuda::Device::c_threads_and_blocks(b.loVect(), b.hiVect(), numBlocks, numThreads);
#endif
        }
        ExecutionConfig (const Box& box, int comps) {
            const Box& b = amrex::surroundingNodes(box);
            Cuda::Device::c_comps_threads_and_blocks(b.loVect(), b.hiVect(), comps, numBlocks, numThreads);
        }
        ExecutionConfig (long N) {
            Cuda::Device::n_threads_and_blocks(N, numBlocks, numThreads);
        }
        ExecutionConfig (dim3 nb, dim3 nt, std::size_t sm=0) : numBlocks(nb), numThreads(nt), sharedMem(sm) {}
        
        dim3 numBlocks;
        dim3 numThreads;
        std::size_t sharedMem = 0;
    };
#endif

}

}

#endif
