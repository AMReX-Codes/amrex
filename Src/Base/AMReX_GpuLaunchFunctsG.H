#ifndef AMREX_GPU_LAUNCH_FUNCTS_G_H_
#define AMREX_GPU_LAUNCH_FUNCTS_G_H_
#include <AMReX_Config.H>

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
                f(Gpu::Handler{&item,shared_data.get_pointer()});
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

namespace detail {
    template <typename F, typename N>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, N i, Gpu::Handler const&)
        noexcept -> decltype(f(0))
    {
        f(i);
    }

    template <typename F, typename N>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, N i, Gpu::Handler const& handler)
        noexcept -> decltype(f(0,Gpu::Handler{}))
    {
        f(i,handler);
    }

    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, Gpu::Handler const&)
        noexcept -> decltype(f(0,0,0))
    {
        f(i,j,k);
    }

    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, Gpu::Handler const& handler)
        noexcept -> decltype(f(0,0,0,Gpu::Handler{}))
    {
        f(i,j,k,handler);
    }

    template <typename F, typename T>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, T ncomp, Gpu::Handler const&)
        noexcept -> decltype(f(0,0,0,0))
    {
        for (T n = 0; n < ncomp; ++n) f(i,j,k,n);
    }

    template <typename F, typename T>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, T ncomp, Gpu::Handler const& handler)
        noexcept -> decltype(f(0,0,0,0,Gpu::Handler{}))
    {
        for (T n = 0; n < ncomp; ++n) f(i,j,k,n,handler);
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
    const auto ec = Gpu::ExecutionConfig(n);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        if (info.hasReduction()) {
            q.submit([&] (sycl::handler& h) {
                sycl::accessor<unsigned long long, 1, sycl::access::mode::read_write,
                               sycl::access::target::local>
                    shared_data(sycl::range<1>(Gpu::Device::warp_size), h);
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
                {
                    for (T i = item.get_global_id(0), stride = item.get_global_range(0);
                         i < n; i += stride) {
                        int n_active_threads = amrex::min(n-i+(T)item.get_local_id(0),
                                                          (T)item.get_local_range(0));
                        detail::call_f(f, i, Gpu::Handler{&item, shared_data.get_pointer(),
                                                          n_active_threads});
                    }
                });
            });
        } else {
            q.submit([&] (sycl::handler& h) {
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
                {
                    for (T i = item.get_global_id(0), stride = item.get_global_range(0);
                         i < n; i += stride) {
                        detail::call_f(f, i, Gpu::Handler{&item});
                    }
                });
            });
        }
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename L>
void ParallelFor (Gpu::KernelInfo const& info, Box const& box, L&& f) noexcept
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
        if (info.hasReduction()) {
            q.submit([&] (sycl::handler& h) {
                sycl::accessor<unsigned long long, 1, sycl::access::mode::read_write,
                               sycl::access::target::local>
                    shared_data(sycl::range<1>(Gpu::Device::warp_size), h);
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
                        int n_active_threads = amrex::min(ncells-icell+(int)item.get_local_id(0),
                                                          (int)item.get_local_range(0));
                        detail::call_f(f, i, j, k, Gpu::Handler{&item, shared_data.get_pointer(),
                                                                n_active_threads});
                    }
                });
            });
        } else {
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
                        detail::call_f(f,i,j,k,Gpu::Handler{&item});
                    }
                });
            });
        }
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const& info, Box const& box, T ncomp, L&& f) noexcept
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
        if (info.hasReduction()) {
            q.submit([&] (sycl::handler& h) {
                sycl::accessor<unsigned long long, 1, sycl::access::mode::read_write,
                               sycl::access::target::local>
                    shared_data(sycl::range<1>(Gpu::Device::warp_size), h);
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
                        int n_active_threads = amrex::min(ncells-icell+(int)item.get_local_id(0),
                                                          (int)item.get_local_range(0));
                        detail::call_f(f, i, j, k, ncomp,
                                       Gpu::Handler{&item, shared_data.get_pointer(),
                                                    n_active_threads});
                    }
                });
            });
        } else {
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
                        detail::call_f(f,i,j,k,ncomp,Gpu::Handler{&item});
                    }
                });
            });
        }
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelForRNG (T n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
    const auto ec = Gpu::ExecutionConfig(n);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * amrex::min(ec.numBlocks.x,Gpu::Device::maxBlocksPerLaunch());
    auto& q = Gpu::Device::nullQueue();
    auto& engdescr = *(getRandEngineDescriptor());
    try {
        q.submit([&] (sycl::handler& h) {
            auto engine_acc = engdescr.get_access(h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                int tid = item.get_global_id(0);
                auto engine = engine_acc.load(tid);
                RandomEngine rand_eng{&engine};
                for (T i = tid, stride = item.get_global_range(0); i < n; i += stride) {
                    f(i,rand_eng);
                }
                engine_acc.store(engine, tid);
            });
        });
        q.wait_and_throw(); // because next launch might be on a different queue
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename L>
void ParallelForRNG (Box const& box, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * amrex::min(ec.numBlocks.x,Gpu::Device::maxBlocksPerLaunch());
    auto& q = Gpu::Device::nullQueue();
    auto& engdescr = *(getRandEngineDescriptor());
    try {
        q.submit([&] (sycl::handler& h) {
            auto engine_acc = engdescr.get_access(h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                int tid = item.get_global_id(0);
                auto engine = engine_acc.load(tid);
                RandomEngine rand_eng{&engine};
                for (int icell = tid, stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    f(i,j,k,rand_eng);
                }
                engine_acc.store(engine, tid);
            });
        });
        q.wait_and_throw(); // because next launch might be on a different queue
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelForRNG (Box const& box, T ncomp, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    // If we are on default queue, block all other streams
    if (Gpu::onNullStream()) Gpu::nonNullStreamSynchronize();
    int nthreads_per_block = ec.numThreads.x;
    int nthreads_total = nthreads_per_block * amrex::min(ec.numBlocks.x,Gpu::Device::maxBlocksPerLaunch());
    auto& q = Gpu::Device::streamQueue();
    auto& engdescr = *(getRandEngineDescriptor());
    try {
        q.submit([&] (sycl::handler& h) {
            auto engine_acc = engdescr.get_access(h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
            {
                int tid = item.get_global_id(0);
                auto engine = engine_acc.load(tid);
                RandomEngine rand_eng{&engine};
                for (int icell = tid, stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    for (T n = 0; n < ncomp; ++n) {
                        f(i,j,k,n,rand_eng);
                    }
                }
                engine_acc.store(engine, tid);
            });
        });
        q.wait_and_throw(); // because next launch might be on a different queue
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename L1, typename L2>
void ParallelFor (Gpu::KernelInfo const& /*info*/, Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) return;
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
}

template <typename L1, typename L2, typename L3>
void ParallelFor (Gpu::KernelInfo const& /*info*/,
                  Box const& box1, Box const& box2, Box const& box3,
                  L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) return;
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
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (Gpu::KernelInfo const& /*info*/,
                  Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) return;
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
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (Gpu::KernelInfo const& /*info*/,
                  Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) return;
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
}

template <typename T, typename L1, typename L2>
void FabReduce (Box const& box, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());
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
          typename M=std::enable_if_t<std::is_integral<N>::value> >
void FabReduce (Box const& box, N ncomp, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());
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
          typename M=std::enable_if_t<std::is_integral<N>::value> >
void VecReduce (N n, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(n)) return;
    auto ec = Gpu::ExecutionConfig(n);
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());
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
                f2(r,Gpu::Handler{&item,shared_data.get_pointer()});
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

namespace detail {
    template <typename F, typename N>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, N i, N /*nleft*/)
        noexcept -> decltype(f(0))
    {
        f(i);
    }

    template <typename F, typename N>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, N i, N nleft)
        noexcept -> decltype(f(0,Gpu::Handler{}))
    {
        f(i,Gpu::Handler(amrex::min(nleft,(N)blockDim.x)));
    }

    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, int /*nleft*/)
        noexcept -> decltype(f(0,0,0))
    {
        f(i,j,k);
    }

    template <typename F>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, int nleft)
        noexcept -> decltype(f(0,0,0,Gpu::Handler{}))
    {
        f(i,j,k,Gpu::Handler(amrex::min(nleft,(int)blockDim.x)));
    }

    template <typename F, typename T>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, T ncomp, int /*nleft*/)
        noexcept -> decltype(f(0,0,0,0))
    {
        for (T n = 0; n < ncomp; ++n) f(i,j,k,n);
    }

    template <typename F, typename T>
    AMREX_GPU_DEVICE
    auto call_f (F const& f, int i, int j, int k, T ncomp, int nleft)
        noexcept -> decltype(f(0,0,0,0,Gpu::Handler{}))
    {
        for (T n = 0; n < ncomp; ++n) f(i,j,k,n,Gpu::Handler(amrex::min(nleft,(int)blockDim.x)));
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
#if defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(AMREX_HIP_INDIRECT_FUNCTION))
    if (!info.hasReduction() && Gpu::inFuseRegion() && info.isFusible() && Long(n) <= Gpu::getFuseSizeThreshold()) {
        Gpu::Register(n, f);
    } else
#endif
    {
        amrex::ignore_unused(info);
        const auto ec = Gpu::ExecutionConfig(n);

        AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept {
            for (T i = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 i < n; i += stride) {
                detail::call_f(f, i, (n-i+(T)threadIdx.x));
            }
        });
        AMREX_GPU_ERROR_CHECK();
    }
}

template <typename L>
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const& info, Box const& box, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
#if defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(AMREX_HIP_INDIRECT_FUNCTION))
    if (!info.hasReduction() && Gpu::inFuseRegion() && info.isFusible() && ncells <= Gpu::getFuseSizeThreshold()) {
        Gpu::Register(box, f);
    } else
#endif
    {
        amrex::ignore_unused(info);
        const auto lo  = amrex::lbound(box);
        const auto len = amrex::length(box);
        const auto ec = Gpu::ExecutionConfig(ncells);

        AMREX_LAUNCH_KERNEL(ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept {
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                 icell < ncells; icell += stride)
            {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                detail::call_f(f, i, j, k, (ncells-icell+(int)threadIdx.x));
            }
        });
        AMREX_GPU_ERROR_CHECK();
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const& info, Box const& box, T ncomp, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
#if defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(AMREX_HIP_INDIRECT_FUNCTION))
    if (!info.hasReduction() && Gpu::inFuseRegion() && info.isFusible() && ncells <= Gpu::getFuseSizeThreshold()) {
        Gpu::Register(box, ncomp, f);
    } else
#endif
    {
        amrex::ignore_unused(info);
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
                detail::call_f(f, i, j, k, ncomp, (ncells-icell+(int)threadIdx.x));
            }
        });
        AMREX_GPU_ERROR_CHECK();
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelForRNG (T n, L&& f) noexcept
{
    if (amrex::isEmpty(n)) return;
    randState_t* rand_state = getRandState();
    const auto ec = Gpu::ExecutionConfig(n);
    AMREX_LAUNCH_KERNEL(amrex::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch()),
                        ec.numThreads, 0, Gpu::nullStream(),  // use null stream
    [=] AMREX_GPU_DEVICE () noexcept {
        int tid = blockDim.x*blockIdx.x+threadIdx.x;
        RandomEngine engine{&(rand_state[tid])};
        for (T i = tid, stride = blockDim.x*gridDim.x; i < n; i += stride) {
            f(i,engine);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L>
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelForRNG (Box const& box, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    randState_t* rand_state = getRandState();
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(amrex::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch()),
                        ec.numThreads, 0, Gpu::nullStream(),  // use null stream
    [=] AMREX_GPU_DEVICE () noexcept {
        int tid = blockDim.x*blockIdx.x+threadIdx.x;
        RandomEngine engine{&(rand_state[tid])};
        for (int icell = tid, stride = blockDim.x*gridDim.x; icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            f(i,j,k,engine);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelForRNG (Box const& box, T ncomp, L&& f) noexcept
{
    if (amrex::isEmpty(box)) return;
    randState_t* rand_state = getRandState();
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    const auto ec = Gpu::ExecutionConfig(ncells);
    AMREX_LAUNCH_KERNEL(amrex::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch()),
                        ec.numThreads, 0, Gpu::nullStream(), // use null stream
    [=] AMREX_GPU_DEVICE () noexcept {
        int tid = blockDim.x*blockIdx.x+threadIdx.x;
        RandomEngine engine{&(rand_state[tid])};
        for (int icell = tid, stride = blockDim.x*gridDim.x; icell < ncells; icell += stride) {
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (T n = 0; n < ncomp; ++n) {
                f(i,j,k,n,engine);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename L1, typename L2>
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
ParallelFor (Gpu::KernelInfo const& info,
             Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
#if defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(AMREX_HIP_INDIRECT_FUNCTION))
    if (Gpu::inFuseRegion() && info.isFusible() && ncells <= Gpu::getFuseSizeThreshold()) {
        Gpu::Register(box1, f1);
        Gpu::Register(box2, f2);
    } else
#endif
    {
        amrex::ignore_unused(info);
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
}

template <typename L1, typename L2, typename L3>
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value && MaybeDeviceRunnable<L3>::value>
ParallelFor (Gpu::KernelInfo const& info,
             Box const& box1, Box const& box2, Box const& box3,
             L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
#if defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(AMREX_HIP_INDIRECT_FUNCTION))
    if (Gpu::inFuseRegion() && info.isFusible() && ncells <= Gpu::getFuseSizeThreshold()) {
        Gpu::Register(box1, f1);
        Gpu::Register(box2, f2);
        Gpu::Register(box3, f3);
    } else
#endif
    {
        amrex::ignore_unused(info);
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
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
ParallelFor (Gpu::KernelInfo const& info,
             Box const& box1, T1 ncomp1, L1&& f1,
             Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells = amrex::max(ncells1, ncells2);
#if defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(AMREX_HIP_INDIRECT_FUNCTION))
    if (Gpu::inFuseRegion() && info.isFusible() && ncells <= Gpu::getFuseSizeThreshold()) {
        Gpu::Register(box1, ncomp1, f1);
        Gpu::Register(box2, ncomp2, f2);
    } else
#endif
    {
        amrex::ignore_unused(info);
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
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value && MaybeDeviceRunnable<L3>::value>
ParallelFor (Gpu::KernelInfo const& info,
             Box const& box1, T1 ncomp1, L1&& f1,
             Box const& box2, T2 ncomp2, L2&& f2,
             Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) return;
    int ncells1 = box1.numPts();
    int ncells2 = box2.numPts();
    int ncells3 = box3.numPts();
    int ncells = amrex::max(ncells1, ncells2, ncells3);
#if defined(AMREX_USE_CUDA) || (defined(AMREX_USE_HIP) && defined(AMREX_HIP_INDIRECT_FUNCTION))
    if (Gpu::inFuseRegion() && info.isFusible() && ncells <= Gpu::getFuseSizeThreshold()) {
        Gpu::Register(box1, ncomp1, f1);
        Gpu::Register(box2, ncomp2, f2);
        Gpu::Register(box3, ncomp3, f3);
    } else
#endif
    {
        amrex::ignore_unused(info);
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
}

template <typename T, typename L1, typename L2>
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
FabReduce (Box const& box, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());
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
          typename M=std::enable_if_t<std::is_integral<N>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
FabReduce (Box const& box, N ncomp, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box)) return;
    int ncells = box.numPts();
    const auto lo  = amrex::lbound(box);
    const auto len = amrex::length(box);
    auto ec = Gpu::ExecutionConfig(ncells);
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());
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
          typename M=std::enable_if_t<std::is_integral<N>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
VecReduce (N n, T const& init_val, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(n)) return;
    auto ec = Gpu::ExecutionConfig(n);
    ec.numBlocks.x = std::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch());
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

template <typename L>
void single_task (L&& f) noexcept
{
    single_task(Gpu::gpuStream(), std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    ParallelFor(info, n,std::forward<L>(f));
}

template <typename L>
void For (Gpu::KernelInfo const& info, Box const& box, L&& f) noexcept
{
    ParallelFor(info, box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const& info, Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(info,box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void For (Gpu::KernelInfo const& info,
          Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void For (Gpu::KernelInfo const& info,
          Box const& box1, Box const& box2, Box const& box3,
          L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(info,box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Gpu::KernelInfo const& info,
          Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Gpu::KernelInfo const& info,
          Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(info,
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (T n, L&& f) noexcept
{
    ParallelFor(Gpu::KernelInfo{}, n, std::forward<L>(f));
}

template <typename L>
void ParallelFor (Box const& box, L&& f) noexcept
{
    ParallelFor(Gpu::KernelInfo{}, box, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void ParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void ParallelFor (Box const& box1, Box const& box2, Box const& box3,
                  L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                  Box const& box2, T2 ncomp2, L2&& f2,
                  Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(Gpu::KernelInfo{},
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (T n, L&& f) noexcept
{
    ParallelFor(Gpu::KernelInfo{}, n,std::forward<L>(f));
}

template <typename L>
void For (Box const& box, L&& f) noexcept
{
    ParallelFor(Gpu::KernelInfo{}, box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Box const& box, T ncomp, L&& f) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void For (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void For (Box const& box1, Box const& box2, Box const& box3,
          L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Box const& box1, T1 ncomp1, L1&& f1,
          Box const& box2, T2 ncomp2, L2&& f2,
          Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor(Gpu::KernelInfo{},
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(info,n,std::forward<L>(f));
    } else {
        AMREX_PRAGMA_SIMD
        for (T i = 0; i < n; ++i) f(i);
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (T n, L&& f) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{}, n, std::forward<L>(f));
}

template <typename L>
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, Box const& box, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(info, box,std::forward<L>(f));
    } else {
        LoopConcurrentOnCpu(box,std::forward<L>(f));
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, Box const& box, T ncomp, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(info, box,ncomp,std::forward<L>(f));
    } else {
        LoopConcurrentOnCpu(box,ncomp,std::forward<L>(f));
    }
}

template <typename L1, typename L2>
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
    } else {
        LoopConcurrentOnCpu(box1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,std::forward<L2>(f2));
    }
}

template <typename L1, typename L2, typename L3>
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value && MaybeHostDeviceRunnable<L3>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       Box const& box1, Box const& box2, Box const& box3,
                       L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(info,box1,box2,box3,
                    std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
    } else {
        LoopConcurrentOnCpu(box1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,std::forward<L2>(f2));
        LoopConcurrentOnCpu(box3,std::forward<L3>(f3));
    }
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       Box const& box1, T1 ncomp1, L1&& f1,
                       Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
    } else {
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
    }
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value && MaybeHostDeviceRunnable<L3>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       Box const& box1, T1 ncomp1, L1&& f1,
                       Box const& box2, T2 ncomp2, L2&& f2,
                       Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor(info,
                    box1,ncomp1,std::forward<L1>(f1),
                    box2,ncomp2,std::forward<L2>(f2),
                    box3,ncomp3,std::forward<L3>(f3));
    } else {
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
        LoopConcurrentOnCpu(box3,ncomp3,std::forward<L3>(f3));
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    HostDeviceParallelFor(info,n,std::forward<L>(f));
}

template <typename L>
void HostDeviceFor (Gpu::KernelInfo const& info, Box const& box, L&& f) noexcept
{
    HostDeviceParallelFor(info,box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info, Box const& box, T ncomp, L&& f) noexcept
{
    HostDeviceParallelFor(info,box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceFor (Gpu::KernelInfo const& info,
                    Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    HostDeviceParallelFor(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceFor (Gpu::KernelInfo const& info,
                    Box const& box1, Box const& box2, Box const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    HostDeviceParallelFor(info, box1,box2,box3,
                          std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info,
                    Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    HostDeviceParallelFor(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info,
                    Box const& box1, T1 ncomp1, L1&& f1,
                    Box const& box2, T2 ncomp2, L2&& f2,
                    Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    HostDeviceParallelFor(info,
                          box1,ncomp1,std::forward<L1>(f1),
                          box2,ncomp2,std::forward<L2>(f2),
                          box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{},n,std::forward<L>(f));
}

template <typename L>
void HostDeviceParallelFor (Box const& box, L&& f) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{},box,std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (Box const& box, T ncomp, L&& f) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2>
void HostDeviceParallelFor (Box const& box1, Box const& box2, L1&& f1, L2&& f2) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3>
void HostDeviceParallelFor (Box const& box1, Box const& box2, Box const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{}, box1,box2,box3,
                          std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (Box const& box1, T1 ncomp1, L1&& f1,
                            Box const& box2, T2 ncomp2, L2&& f2,
                            Box const& box3, T3 ncomp3, L3&& f3) noexcept
{
    HostDeviceParallelFor(Gpu::KernelInfo{},
                          box1,ncomp1,std::forward<L1>(f1),
                          box2,ncomp2,std::forward<L2>(f2),
                          box3,ncomp3,std::forward<L3>(f3));
}

}

#endif
