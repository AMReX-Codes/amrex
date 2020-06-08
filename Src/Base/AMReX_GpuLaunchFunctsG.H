#ifndef AMREX_GPU_LAUNCH_FUNCTS_G_H_
#define AMREX_GPU_LAUNCH_FUNCTS_G_H_

namespace amrex {

#ifdef AMREX_USE_DPCPP

template <typename L>
void single_task (gpuStream_t stream, L&& f) noexcept
{
    if (Gpu::onNullStream(stream)) Gpu::nonNullStreamSynchronize();
    auto& q = *(stream.queue);
    try {
        q.submit([&] (sycl::handler& h) {
            h.single_task([=] () { f(); });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("single_task: ")+ex.what()+"!!!!!");
    }
}

template<typename L>
void launch (int nblocks, int nthreads_per_block, std::size_t shared_mem_bytes,
             gpuStream_t stream, L&& f) noexcept
{
    if (Gpu::onNullStream(stream)) Gpu::nonNullStreamSynchronize();
    int nthreads_total = nthreads_per_block * nblocks;
    std::size_t shared_mem_numull = (shared_mem_bytes+sizeof(unsigned long long)-1)
        / sizeof(unsigned long long);
    auto& q = *(stream.queue);
    try {
        q.submit([&] (sycl::handler& h) {
            sycl::accessor<unsigned long long, 1, sycl::access::mode::read_write, sycl::access::target::local>
                shared_data(sycl::range<1>(shared_mem_numull), h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                f(Gpu::Handler{item,shared_data.get_pointer()});
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("launch: ")+ex.what()+"!!!!!");
    }
}

template<typename L>
void launch (int nblocks, int nthreads_per_block, gpuStream_t stream, L&& f) noexcept
{
    if (Gpu::onNullStream(stream)) Gpu::nonNullStreamSynchronize();
    int nthreads_total = nthreads_per_block * nblocks;
    auto& q = *(stream.queue);
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                f(item);
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("launch: ")+ex.what()+"!!!!!");
    }
}

template<typename T, typename L>
void launch (T const& n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
    const auto ec = Gpu::ExecutionConfig(n);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (auto const i : Gpu::Range(n,item.get_global_id(0),item.get_global_range(0))) {
                    f(i);
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("launch: ")+ex.what()+"!!!!!");
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void ParallelFor (T n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
    const auto ec = Gpu::ExecutionConfig(n);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (T i = item.get_global_id(0), stride = item.get_global_range(0);
                     i < n; i += stride) {
                    f(i);
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename L>
void ParallelFor (Box const& box, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    f(i,j,k);
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void ParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    for (T n = 0; n < ncomp; ++n) {
                        f(i,j,k,n);
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename L1, typename L2>
void ParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    // xxxxx DPCPP todo: launch separate kernel to reduce kernel size
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
#if 0
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < ncells1) {
                        int k =  icell /   (len1.x*len1.y);
                        int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                        int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                        i += lo1.x;
                        j += lo1.y;
                        k += lo1.z;
                        f1(i,j,k);
                    }
                    if (icell < ncells2) {
                        int k =  icell /   (len2.x*len2.y);
                        int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                        int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                        i += lo2.x;
                        j += lo2.y;
                        k += lo2.z;
                        f2(i,j,k);
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
#endif
}

template <typename L1, typename L2, typename L3>
void ParallelFor (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    // xxxxx DPCPP todo: launch separate kernel to reduce kernel size
    ParallelFor(box1, std::forward<L1>(f1));
    ParallelFor(box2, std::forward<L2>(f2));
    ParallelFor(box3, std::forward<L3>(f3));
#if 0
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2) and amrex::isEmpty(box3)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto lo3  = amrex::lbound(box3);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto len3 = amrex::length(box3);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < ncells1) {
                        int k =  icell /   (len1.x*len1.y);
                        int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                        int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                        i += lo1.x;
                        j += lo1.y;
                        k += lo1.z;
                        f1(i,j,k);
                    }
                    if (icell < ncells2) {
                        int k =  icell /   (len2.x*len2.y);
                        int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                        int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                        i += lo2.x;
                        j += lo2.y;
                        k += lo2.z;
                        f2(i,j,k);
                    }
                    if (icell < ncells3) {
                        int k =  icell /   (len3.x*len3.y);
                        int j = (icell - k*(len3.x*len3.y)) /   len3.x;
                        int i = (icell - k*(len3.x*len3.y)) - j*len3.x;
                        i += lo3.x;
                        j += lo3.y;
                        k += lo3.z;
                        f3(i,j,k);
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
#endif
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    // xxxxx DPCPP todo: launch separate kernel to reduce kernel size
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
#if 0
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < ncells1) {
                        int k =  icell /   (len1.x*len1.y);
                        int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                        int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                        i += lo1.x;
                        j += lo1.y;
                        k += lo1.z;
                        for (T1 n = 0; n < ncomp1; ++n) {
                            f1(i,j,k,n);
                        }
                    }
                    if (icell < ncells2) {
                        int k =  icell /   (len2.x*len2.y);
                        int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                        int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                        i += lo2.x;
                        j += lo2.y;
                        k += lo2.z;
                        for (T2 n = 0; n < ncomp2; ++n) {
                            f2(i,j,k,n);
                        }
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
#endif
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    // xxxxx DPCPP todo: launch separate kernel to reduce kernel size
    ParallelFor(box1, ncomp1, std::forward<L1>(f1));
    ParallelFor(box2, ncomp2, std::forward<L2>(f2));
    ParallelFor(box3, ncomp3, std::forward<L3>(f3));
#if 0
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2) and amrex::isEmpty(box3)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto lo3  = amrex::lbound(box3);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto len3 = amrex::length(box3);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < ncells1) {
                        int k =  icell /   (len1.x*len1.y);
                        int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                        int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                        i += lo1.x;
                        j += lo1.y;
                        k += lo1.z;
                        for (T1 n = 0; n < ncomp1; ++n) {
                            f1(i,j,k,n);
                        }
                    }
                    if (icell < ncells2) {
                        int k =  icell /   (len2.x*len2.y);
                        int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                        int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                        i += lo2.x;
                        j += lo2.y;
                        k += lo2.z;
                        for (T2 n = 0; n < ncomp2; ++n) {
                            f2(i,j,k,n);
                        }
                    }
                    if (icell < ncells3) {
                        int k =  icell /   (len3.x*len3.y);
                        int j = (icell - k*(len3.x*len3.y)) /   len3.x;
                        int i = (icell - k*(len3.x*len3.y)) - j*len3.x;
                        i += lo3.x;
                        j += lo3.y;
                        k += lo3.z;
                        for (T3 n = 0; n < ncomp3; ++n) {
                            f3(i,j,k,n);
                        }
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
#endif
}

template <typename T, typename L1, typename L2>
void FabReduce (Box const& box, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                auto r = init_val;
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    f1(i,j,k,&r);
                }
                f2(r);
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void FabReduce (Box const& box, N ncomp, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                auto r = init_val;
                for (int icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    for (N n = 0; n < ncomp; ++n) {
                        f1(i,j,k,n,&r);
                    }
                }
                f2(r);
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
void VecReduce (N n, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(n)) return;
    auto ec = Gpu::ExecutionConfig(n);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
                sycl::accessor<T, 1, sycl::access::mode::read_write, sycl::access::target::local>
                    shared_data(sycl::range<1>(Gpu::Device::warp_size), h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                auto r = init_val;
                for (N i = item.get_global_id(0), stride = item.get_global_range(0);
                     i < n; i += stride) {
                    f1(i,&r);
                }
                f2(r,Gpu::Handler{item,shared_data.get_pointer()});
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

#else
// CUDA or HIP

template <typename L>
void single_task (gpuStream_t stream, L&& f) noexcept
{
    AMREX_LAUNCH_KERNEL(1, 1, 0, stream,
                        [=] AMREX_GPU_DEVICE () noexcept {f();});
    AMREX_GPU_ERROR_CHECK();
}

template<typename L>
void launch (int nblocks, int nthreads_per_block, std::size_t shared_mem_bytes,
             gpuStream_t stream, L&& f) noexcept
{
    AMREX_LAUNCH_KERNEL(nblocks, nthreads_per_block, shared_mem_bytes,
                        stream, [=] AMREX_GPU_DEVICE () noexcept { f(); });
    AMREX_GPU_ERROR_CHECK();
}

template<typename L>
void launch (int nblocks, int nthreads_per_block, gpuStream_t stream, L&& f) noexcept
{
    launch(nblocks, nthreads_per_block, 0, stream, std::forward<L>(f));
}

template<typename T, typename L>
void launch (T const& n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
    const auto ec = Gpu::ExecutionConfig(n);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (auto const i : Gpu::Range(n)) {
            f(i);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
amrex::EnableIf_t<MaybeDeviceRunnable<L>::value>
ParallelFor (T n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
    const auto ec = Gpu::ExecutionConfig(n);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (T i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             i < n; i += stride) {
            f(i);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L>
amrex::EnableIf_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Box const& box, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            f(i,j,k);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
amrex::EnableIf_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (T n = 0; n < ncomp; ++n) {
                f(i,j,k,n);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L1, typename L2>
amrex::EnableIf_t<MaybeDeviceRunnable<L1>::value and MaybeDeviceRunnable<L2>::value>
ParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                f1(i,j,k);
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                f2(i,j,k);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L1, typename L2, typename L3>
amrex::EnableIf_t<MaybeDeviceRunnable<L1>::value and MaybeDeviceRunnable<L2>::value and MaybeDeviceRunnable<L3>::value>
ParallelFor (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2) and amrex::isEmpty(box3)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto lo3  = amrex::lbound(box3);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto len3 = amrex::length(box3);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                f1(i,j,k);
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                f2(i,j,k);
            }
            if (icell < ncells3) {
                int k =  icell /   (len3.x*len3.y);
                int j = (icell - k*(len3.x*len3.y)) /   len3.x;
                int i = (icell - k*(len3.x*len3.y)) - j*len3.x;
                i += lo3.x;
                j += lo3.y;
                k += lo3.z;
                f3(i,j,k);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
amrex::EnableIf_t<MaybeDeviceRunnable<L1>::value and MaybeDeviceRunnable<L2>::value>
ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
             Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                for (T1 n = 0; n < ncomp1; ++n) {
                    f1(i,j,k,n);
                }
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                for (T2 n = 0; n < ncomp2; ++n) {
                    f2(i,j,k,n);
                }
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
amrex::EnableIf_t<MaybeDeviceRunnable<L1>::value and MaybeDeviceRunnable<L2>::value and MaybeDeviceRunnable<L3>::value>
ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
             Box const& box2, T2 ncomp2, L2&& f2,
             Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) and amrex::isEmpty(box2) and amrex::isEmpty(box3)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
    const auto lo1  = amrex::lbound(box1);
    const auto lo2  = amrex::lbound(box2);
    const auto lo3  = amrex::lbound(box3);
    const auto len1 = amrex::length(box1);
    const auto len2 = amrex::length(box2);
    const auto len3 = amrex::length(box3);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < ncells1) {
                int k =  icell /   (len1.x*len1.y);
                int j = (icell - k*(len1.x*len1.y)) /   len1.x;
                int i = (icell - k*(len1.x*len1.y)) - j*len1.x;
                i += lo1.x;
                j += lo1.y;
                k += lo1.z;
                for (T1 n = 0; n < ncomp1; ++n) {
                    f1(i,j,k,n);
                }
            }
            if (icell < ncells2) {
                int k =  icell /   (len2.x*len2.y);
                int j = (icell - k*(len2.x*len2.y)) /   len2.x;
                int i = (icell - k*(len2.x*len2.y)) - j*len2.x;
                i += lo2.x;
                j += lo2.y;
                k += lo2.z;
                for (T2 n = 0; n < ncomp2; ++n) {
                    f2(i,j,k,n);
                }
            }
            if (icell < ncells3) {
                int k =  icell /   (len3.x*len3.y);
                int j = (icell - k*(len3.x*len3.y)) /   len3.x;
                int i = (icell - k*(len3.x*len3.y)) - j*len3.x;
                i += lo3.x;
                j += lo3.y;
                k += lo3.z;
                for (T3 n = 0; n < ncomp3; ++n) {
                    f3(i,j,k,n);
                }
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L1, typename L2>
amrex::EnableIf_t<MaybeDeviceRunnable<L1>::value and MaybeDeviceRunnable<L2>::value>
FabReduce (Box const& box, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto r = init_val;
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            f1(i,j,k,&r);
        }
        f2(r);
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
amrex::EnableIf_t<MaybeDeviceRunnable<L1>::value and MaybeDeviceRunnable<L2>::value>
FabReduce (Box const& box, N ncomp, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto r = init_val;
        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (N n = 0; n < ncomp; ++n) {
                f1(i,j,k,n,&r);
            }
        }
        f2(r);
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename N, typename T, typename L1, typename L2,
          typename M=amrex::EnableIf_t<std::is_integral<N>::value> >
amrex::EnableIf_t<MaybeDeviceRunnable<L1>::value and MaybeDeviceRunnable<L2>::value>
VecReduce (N n, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(n)) return;
    auto ec = Gpu::ExecutionConfig(n);
    ec.numBlocks.x = std::min(ec.numBlocks.x, static_cast<unsigned int>(Gpu::Device::maxBlocksPerLaunch()));
    AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto r = init_val;
        for (N i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
             i < n; i += stride) {
            f1(i,&r);
        }
        f2(r);
    });
    AMREX_GPU_ERROR_CHECK();
}

#endif

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (T n, L&& f) noexcept
{
    ParallelFor(n,std::forward<L>(f));
}

template <typename L>
void For (Box const& box, L&& f) noexcept
{
    ParallelFor(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void For (Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void For (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void For (Box const& box1, Box const& box2, Box const& box3, L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
amrex::EnableIf_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (T n, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(n,std::forward<L>(f));
    } else {
        AMREX_PRAGMA_SIMD
        for (T i = 0; i < n; ++i) f(i);
    }
}

template <typename L>
amrex::EnableIf_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Box const& box, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(box,std::forward<L>(f));
    } else {
        LoopConcurrentOnCpu(box,std::forward<L>(f));
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
amrex::EnableIf_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(box,ncomp,std::forward<L>(f));
    } else {
        LoopConcurrentOnCpu(box,ncomp,std::forward<L>(f));
    }
}

template <typename L1, typename L2>
amrex::EnableIf_t<MaybeHostDeviceRunnable<L1>::value and MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
    } else {
        LoopConcurrentOnCpu(box1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,std::forward<L2>(f2));
    }
}

template <typename L1, typename L2, typename L3>
amrex::EnableIf_t<MaybeHostDeviceRunnable<L1>::value and MaybeHostDeviceRunnable<L2>::value and MaybeHostDeviceRunnable<L3>::value>
HostDeviceParallelFor (Box const& box1, Box const& box2, Box const& box3,
                       L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
    } else {
        LoopConcurrentOnCpu(box1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,std::forward<L2>(f2));
        LoopConcurrentOnCpu(box3,std::forward<L3>(f3));
    }
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
amrex::EnableIf_t<MaybeHostDeviceRunnable<L1>::value and MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                       Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
    } else {
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
    }
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
amrex::EnableIf_t<MaybeHostDeviceRunnable<L1>::value and MaybeHostDeviceRunnable<L2>::value and MaybeHostDeviceRunnable<L3>::value>
HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                       Box const& box2, T2 ncomp2, L2&& f2,
                       Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(box1,ncomp1,std::forward<L1>(f1),
                    box2,ncomp2,std::forward<L2>(f2),
                    box3,ncomp3,std::forward<L3>(f3));
    } else {
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
        LoopConcurrentOnCpu(box3,ncomp3,std::forward<L3>(f3));
    }
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void HostDeviceFor (T n, L&& f) noexcept
{
    HostDeviceParallelFor(n,std::forward<L>(f));
}

template <typename L>
void HostDeviceFor (Box const& box, L&& f) noexcept
{
    HostDeviceParallelFor(box,std::forward<L>(f));
}

template <typename T, typename L, typename M=amrex::EnableIf_t<std::is_integral<T>::value> >
void HostDeviceFor (Box const& box, T ncomp, L&& f) noexcept
{
    HostDeviceParallelFor(box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    HostDeviceParallelFor(box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceFor (Box const& box1, Box const& box2, Box const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    HostDeviceParallelFor(box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    HostDeviceParallelFor(box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=amrex::EnableIf_t<std::is_integral<T1>::value>,
          typename M2=amrex::EnableIf_t<std::is_integral<T2>::value>,
          typename M3=amrex::EnableIf_t<std::is_integral<T3>::value> >
void HostDeviceFor (Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2,
                    Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    HostDeviceParallelFor(box1,ncomp1,std::forward<L1>(f1),
                          box2,ncomp2,std::forward<L2>(f2),
                          box3,ncomp3,std::forward<L3>(f3));
}

}

#endif
