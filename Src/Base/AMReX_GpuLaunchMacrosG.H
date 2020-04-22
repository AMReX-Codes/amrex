#ifndef AMREX_GPU_LAUNCH_MACROS_G_H_
#define AMREX_GPU_LAUNCH_MACROS_G_H_

#ifdef AMREX_USE_DPCPP
#define AMREX_GPU_LAUNCH_HOST_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    { auto const& amrex_i_tn = TN; \
    if (!amrex::isEmpty(amrex_i_tn)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        const auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_tn); \
        if (amrex::Gpu::onNullStream()) amrex::Gpu::nonNullStreamSynchronize(); \
        int amrex_i_nthreads_per_block = amrex_i_ec.numThreads.x; \
        int amrex_i_nthreads_total = amrex_i_nthreads_per_block * amrex_i_ec.numBlocks.x; \
        auto& amrex_i_q = amrex::Gpu::Device::streamQueue(); \
        try { \
            amrex_i_q.submit([&] (sycl::handler& amrex_i_h) { \
                amrex_i_h.parallel_for(sycl::nd_range<1>(sycl::range<1>(amrex_i_nthreads_total), \
                                                         sycl::range<1>(amrex_i_nthreads_per_block)), \
                [=] (sycl::nd_item<1> amrex_i_item) \
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size) \
                { \
                    for (auto const TI : amrex::Gpu::Range(amrex_i_tn,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                        block \
                    } \
                }); \
            }); \
        } catch (sycl::exception const& ex) { \
            amrex::Abort(std::string("LAUNCH: ")+ex.what()+"!!!!!"); \
        } \
    } \
    else { \
        for (auto const TI : amrex::Gpu::Range(amrex_i_tn)) { \
            block \
        } \
    }}}
#else
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
#endif

// two fused launches
#ifdef AMREX_USE_DPCPP
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
        if (amrex::Gpu::onNullStream()) amrex::Gpu::nonNullStreamSynchronize(); \
        int amrex_i_nthreads_per_block = amrex_i_ec1.numThreads.x; \
        int amrex_i_nthreads_total = amrex_i_nthreads_per_block * amrex_i_nblocks.x; \
        auto& amrex_i_q = amrex::Gpu::Device::streamQueue(); \
        try { \
            amrex_i_q.submit([&] (sycl::handler& amrex_i_h) { \
                amrex_i_h.parallel_for(sycl::nd_range<2>(sycl::range<2>(amrex_i_nthreads_total,2), \
                                                         sycl::range<2>(amrex_i_nthreads_per_block,1)), \
                [=] (sycl::nd_item<2> amrex_i_item) \
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size) \
                { \
                    switch (amrex_i_item.get_group(1)) { \
                    case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block1 \
                            } \
                    break; \
                    case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block2 \
                            } \
                    } \
                }); \
            }); \
        } catch (sycl::exception const& ex) { \
            amrex::Abort(std::string("LAUNCH: ")+ex.what()+"!!!!!"); \
        } \
    } \
    else { \
        for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1)) { \
            block1 \
        } \
        for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2)) { \
            block2 \
        } \
    }}}
#else
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
#endif

// three fused launches
#ifdef AMREX_USE_DPCPP
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
        if (amrex::Gpu::onNullStream()) amrex::Gpu::nonNullStreamSynchronize(); \
        int amrex_i_nthreads_per_block = amrex_i_ec1.numThreads.x; \
        int amrex_i_nthreads_total = amrex_i_nthreads_per_block * amrex_i_nblocks.x; \
        auto& amrex_i_q = amrex::Gpu::Device::streamQueue(); \
        try { \
            amrex_i_q.submit([&] (sycl::handler& amrex_i_h) { \
                amrex_i_h.parallel_for(sycl::nd_range<2>(sycl::range<2>(amrex_i_nthreads_total,3), \
                                                         sycl::range<2>(amrex_i_nthreads_per_block,1)), \
                [=] (sycl::nd_item<2> amrex_i_item) \
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size) \
                { \
                    switch (amrex_i_item.get_group(1)) { \
                    case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block1 \
                            } \
                            break; \
                    case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block2 \
                            } \
                            break; \
                    case 2: for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block3 \
                            } \
                    } \
                }); \
            }); \
        } catch (sycl::exception const& ex) { \
            amrex::Abort(std::string("LAUNCH: ")+ex.what()+"!!!!!"); \
        } \
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
#else
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
#endif

#ifdef AMREX_USE_DPCPP
#define AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE(TN,TI,block) \
    { auto const& amrex_i_tn = TN; \
    if (!amrex::isEmpty(amrex_i_tn)) { \
    if (amrex::Gpu::inLaunchRegion()) \
    { \
        auto amrex_i_ec = amrex::Gpu::ExecutionConfig(amrex_i_tn); \
        if (amrex::Gpu::onNullStream()) amrex::Gpu::nonNullStreamSynchronize(); \
        int amrex_i_nthreads_per_block = amrex_i_ec.numThreads.x; \
        int amrex_i_nthreads_total = amrex_i_nthreads_per_block * amrex_i_ec.numBlocks.x; \
        auto& amrex_i_q = amrex::Gpu::Device::streamQueue(); \
        try { \
            amrex_i_q.submit([&] (sycl::handler& amrex_i_h) { \
                amrex_i_h.parallel_for(sycl::nd_range<1>(sycl::range<1>(amrex_i_nthreads_total), \
                                                         sycl::range<1>(amrex_i_nthreads_per_block)), \
                [=] (sycl::nd_item<1> amrex_i_item) \
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size) \
                { \
                    for (auto const TI : amrex::Gpu::Range(amrex_i_tn,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                        block \
                    } \
                }); \
            }); \
        } catch (sycl::exception const& ex) { \
            amrex::Abort(std::string("LAUNCH: ")+ex.what()+"!!!!!"); \
        } \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE: cannot call device function from host"); \
    }}}
#else
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
#endif

// two fused launches
#ifdef AMREX_USE_DPCPP
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
        if (amrex::Gpu::onNullStream()) amrex::Gpu::nonNullStreamSynchronize(); \
        int amrex_i_nthreads_per_block = amrex_i_ec1.numThreads.x; \
        int amrex_i_nthreads_total = amrex_i_nthreads_per_block * amrex_i_nblocks.x; \
        auto& amrex_i_q = amrex::Gpu::Device::streamQueue(); \
        try { \
            amrex_i_q.submit([&] (sycl::handler& amrex_i_h) { \
                amrex_i_h.parallel_for(sycl::nd_range<2>(sycl::range<2>(amrex_i_nthreads_total,2), \
                                                         sycl::range<2>(amrex_i_nthreads_per_block,1)), \
                [=] (sycl::nd_item<2> amrex_i_item) \
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size) \
                { \
                    switch (amrex_i_item.get_group(1)) { \
                    case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block1 \
                            } \
                            break; \
                    case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block2 \
                            } \
                    } \
                }); \
            }); \
        } catch (sycl::exception const& ex) { \
            amrex::Abort(std::string("LAUNCH: ")+ex.what()+"!!!!!"); \
        } \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }}}
#else
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
#endif

// three fused launches
#ifdef AMREX_USE_DPCPP
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
        if (amrex::Gpu::onNullStream()) amrex::Gpu::nonNullStreamSynchronize(); \
        int amrex_i_nthreads_per_block = amrex_i_ec1.numThreads.x; \
        int amrex_i_nthreads_total = amrex_i_nthreads_per_block * amrex_i_nblocks.x; \
        auto& amrex_i_q = amrex::Gpu::Device::streamQueue(); \
        try { \
            amrex_i_q.submit([&] (sycl::handler& amrex_i_h) { \
                amrex_i_h.parallel_for(sycl::nd_range<2>(sycl::range<2>(amrex_i_nthreads_total,3), \
                                                         sycl::range<2>(amrex_i_nthreads_per_block,1)), \
                [=] (sycl::nd_item<2> amrex_i_item) \
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size) \
                { \
                    switch (amrex_i_item.get_group(1)) { \
                    case 0: for (auto const TI1 : amrex::Gpu::Range(amrex_i_tn1,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block1 \
                            } \
                            break; \
                    case 1: for (auto const TI2 : amrex::Gpu::Range(amrex_i_tn2,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block2 \
                            } \
                    case 2: for (auto const TI3 : amrex::Gpu::Range(amrex_i_tn3,amrex_i_item.get_global_id(0),amrex_i_item.get_global_range(0))) { \
                                block3 \
                            } \
                    } \
                }); \
            }); \
        } catch (sycl::exception const& ex) { \
            amrex::Abort(std::string("LAUNCH: ")+ex.what()+"!!!!!"); \
        } \
    } \
    else { \
        amrex::Abort("AMREX_GPU_LAUNCH_DEVICE_LAMBDA_RANGE_2: cannot call device function from host"); \
    }}}
#else
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
#endif

// FOR_1D

#define AMREX_GPU_HOST_DEVICE_FOR_1D(n,i,block) \
{ \
    auto const& amrex_i_n = n; \
    using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_n,[=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept block); \
    } else { \
        auto amrex_i_lambda = [=] (amrex_i_inttype i) noexcept block; \
        AMREX_PRAGMA_SIMD \
        for (amrex_i_inttype i = 0; i < amrex_i_n; ++i) amrex_i_lambda(i); \
    } \
}

#define AMREX_GPU_DEVICE_FOR_1D(n,i,block) \
{ \
    using amrex_i_inttype = typename std::remove_const<decltype(n)>::type; \
    amrex::ParallelFor(n,[=] AMREX_GPU_DEVICE (amrex_i_inttype i) noexcept block); \
}

// FOR_3D

#define AMREX_GPU_HOST_DEVICE_FOR_3D(box,i,j,k,block) \
{ \
    auto const& amrex_i_box = box; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept block); \
    } else { \
        amrex::LoopConcurrentOnCpu(amrex_i_box,[=] (int i, int j, int k) noexcept block); \
    } \
}

#define AMREX_GPU_DEVICE_FOR_3D(box,i,j,k,block) \
{ \
    amrex::ParallelFor(box,[=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept block); \
}

// FOR_4D

#define AMREX_GPU_HOST_DEVICE_FOR_4D(box,ncomp,i,j,k,n,block) \
{ \
    auto const& amrex_i_box = box; \
    auto const& amrex_i_ncomp = ncomp; \
    if (amrex::Gpu::inLaunchRegion()) { \
        amrex::ParallelFor(amrex_i_box,amrex_i_ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept block); \
    } else { \
        amrex::LoopConcurrentOnCpu(amrex_i_box,amrex_i_ncomp,[=] (int i, int j, int k, int n) noexcept block); \
    } \
}

#define AMREX_GPU_DEVICE_FOR_4D(box,ncomp,i,j,k,n,block) \
{ \
    amrex::ParallelFor(box,ncomp,[=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept block); \
}

#define AMREX_GPU_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_GPU_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_GPU_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_DEVICE_FOR_4D(__VA_ARGS__)

#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_1D(...) AMREX_GPU_HOST_DEVICE_FOR_1D(__VA_ARGS__)
#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_3D(...) AMREX_GPU_HOST_DEVICE_FOR_3D(__VA_ARGS__)
#define AMREX_GPU_HOST_DEVICE_PARALLEL_FOR_4D(...) AMREX_GPU_HOST_DEVICE_FOR_4D(__VA_ARGS__)

#endif
