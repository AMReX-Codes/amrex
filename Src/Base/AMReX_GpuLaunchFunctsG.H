#ifndef AMREX_GPU_LAUNCH_FUNCTS_G_H_
#define AMREX_GPU_LAUNCH_FUNCTS_G_H_
#include <AMReX_Config.H>

namespace amrex {

namespace detail {

    // call_f_scalar_handler

    template <typename F, typename N>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_scalar_handler (F const& f, N i, Gpu::Handler const&)
        noexcept -> decltype(f(0))
    {
        f(i);
    }

    template <typename F, typename N>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_scalar_handler (F const& f, N i, Gpu::Handler const& handler)
        noexcept -> decltype(f(0,Gpu::Handler{}))
    {
        f(i, handler);
    }

    // call_f_intvect_inner

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<1> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], 0, 0, args...);
    }

    template <typename F, std::size_t...Ns, class...Args>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<2> iv, Args...args)
        noexcept -> decltype(f(0, 0, 0, args...))
    {
        f(iv[0], iv[1], 0, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv, args...))
    {
        f(iv, args...);
    }

    template <typename F, int dim, std::size_t...Ns, class...Args>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_inner (std::index_sequence<Ns...>, F const& f, IntVectND<dim> iv, Args...args)
        noexcept -> decltype(f(iv[Ns]..., args...))
    {
        f(iv[Ns]..., args...);
    }

    // call_f_intvect

    template <typename F, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect (F const& f, IntVectND<dim> iv)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv);
    }

    // call_f_intvect_engine

    template <typename F, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_engine (F const& f, IntVectND<dim> iv, RandomEngine engine)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, engine))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, engine);
    }

    // call_f_intvect_handler

    template <typename F, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_handler (F const& f, IntVectND<dim> iv, Gpu::Handler const&)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv);
    }

    template <typename F, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_handler (F const& f, IntVectND<dim> iv, Gpu::Handler const& handler)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, Gpu::Handler{}))
    {
        call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, handler);
    }

    // call_f_intvect_ncomp

    template <typename F, typename T, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp (F const& f, IntVectND<dim> iv, T ncomp)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, 0))
    {
        for (T n = 0; n < ncomp; ++n) {
            call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n);
        }
    }

    // call_f_intvect_ncomp_engine

    template <typename F, typename T, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp_engine (F const& f, IntVectND<dim> iv, T ncomp, RandomEngine engine)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, 0, engine))
    {
        for (T n = 0; n < ncomp; ++n) {
            call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n, engine);
        }
    }

    // call_f_intvect_ncomp_handler

    template <typename F, typename T, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp_handler (F const& f, IntVectND<dim> iv, T ncomp, Gpu::Handler const&)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, 0))
    {
        for (T n = 0; n < ncomp; ++n) {
            call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n);
        }
    }

    template <typename F, typename T, int dim>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    auto call_f_intvect_ncomp_handler (F const& f, IntVectND<dim> iv, T ncomp, Gpu::Handler const& handler)
        noexcept -> decltype(call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, 0, Gpu::Handler{}))
    {
        for (T n = 0; n < ncomp; ++n) {
            call_f_intvect_inner(std::make_index_sequence<dim>(), f, iv, n, handler);
        }
    }

}

#ifdef AMREX_USE_SYCL

template <typename L>
void single_task (gpuStream_t stream, L const& f) noexcept
{
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
             gpuStream_t stream, L const& f) noexcept
{
    const auto nthreads_total = std::size_t(nthreads_per_block) * nblocks;
    const std::size_t shared_mem_numull = (shared_mem_bytes+sizeof(unsigned long long)-1)
        / sizeof(unsigned long long);
    auto& q = *(stream.queue);
    try {
        q.submit([&] (sycl::handler& h) {
            sycl::local_accessor<unsigned long long>
                shared_data(sycl::range<1>(shared_mem_numull), h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                f(Gpu::Handler{&item,shared_data.get_multi_ptr<sycl::access::decorated::yes>().get()});
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("launch: ")+ex.what()+"!!!!!");
    }
}

template<typename L>
void launch (int nblocks, int nthreads_per_block, gpuStream_t stream, L const& f) noexcept
{
    const auto nthreads_total = std::size_t(nthreads_per_block) * nblocks;
    auto& q = *(stream.queue);
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                f(item);
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("launch: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename L>
void launch (int nblocks, std::size_t shared_mem_bytes, gpuStream_t stream,
             L const& f) noexcept
{
    const auto nthreads_total = MT * std::size_t(nblocks);
    const std::size_t shared_mem_numull = (shared_mem_bytes+sizeof(unsigned long long)-1)
        / sizeof(unsigned long long);
    auto& q = *(stream.queue);
    try {
        q.submit([&] (sycl::handler& h) {
            sycl::local_accessor<unsigned long long>
                shared_data(sycl::range<1>(shared_mem_numull), h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(MT)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(MT)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                f(Gpu::Handler{&item,shared_data.get_multi_ptr<sycl::access::decorated::yes>().get()});
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("launch: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename L>
void launch (int nblocks, gpuStream_t stream, L const& f) noexcept
{
    const auto nthreads_total = MT * std::size_t(nblocks);
    auto& q = *(stream.queue);
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(MT)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(MT)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                f(item);
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("launch: ")+ex.what()+"!!!!!");
    }
}

template<int MT, typename T, typename L>
void launch (T const& n, L const& f) noexcept
{
    if (amrex::isEmpty(n)) { return; }
    const auto ec = Gpu::makeExecutionConfig<MT>(n);
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(MT)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
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

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const& info, T n, L const& f) noexcept
{
    if (amrex::isEmpty(n)) { return; }
    const auto ec = Gpu::makeExecutionConfig<MT>(n);
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        if (info.hasReduction()) {
            q.submit([&] (sycl::handler& h) {
                sycl::local_accessor<unsigned long long>
                    shared_data(sycl::range<1>(Gpu::Device::warp_size), h);
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                [[sycl::reqd_work_group_size(MT)]]
                [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
                {
                    for (std::size_t i = item.get_global_id(0), stride = item.get_global_range(0);
                         i < std::size_t(n); i += stride) {
                        int n_active_threads = amrex::min(std::size_t(n)-i+item.get_local_id(0),
                                                          item.get_local_range(0));
                        detail::call_f_scalar_handler(f, T(i), Gpu::Handler{&item, shared_data.get_multi_ptr<sycl::access::decorated::yes>().get(),
                                                          n_active_threads});
                    }
                });
            });
        } else {
            q.submit([&] (sycl::handler& h) {
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                [[sycl::reqd_work_group_size(MT)]]
                [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
                {
                    for (std::size_t i = item.get_global_id(0), stride = item.get_global_range(0);
                         i < std::size_t(n); i += stride) {
                        detail::call_f_scalar_handler(f, T(i), Gpu::Handler{&item});
                    }
                });
            });
        }
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename L, int dim>
void ParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::makeExecutionConfig<MT>(box.numPts());
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        if (info.hasReduction()) {
            q.submit([&] (sycl::handler& h) {
                sycl::local_accessor<unsigned long long>
                    shared_data(sycl::range<1>(Gpu::Device::warp_size), h);
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                [[sycl::reqd_work_group_size(MT)]]
                [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
                {
                    for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                         icell < indexer.numPts(); icell += stride) {
                        auto iv = indexer.intVect(icell);
                        int n_active_threads = amrex::min(indexer.numPts()-icell+std::uint64_t(item.get_local_id(0)),
                                                          std::uint64_t(item.get_local_range(0)));
                        detail::call_f_intvect_handler(f, iv, Gpu::Handler{&item, shared_data.get_multi_ptr<sycl::access::decorated::yes>().get(),
                                                                n_active_threads});
                    }
                });
            });
        } else {
            q.submit([&] (sycl::handler& h) {
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                [[sycl::reqd_work_group_size(MT)]]
                [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
                {
                    for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                         icell < indexer.numPts(); icell += stride) {
                        auto iv = indexer.intVect(icell);
                        detail::call_f_intvect_handler(f,iv,Gpu::Handler{&item});
                    }
                });
            });
        }
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::makeExecutionConfig<MT>(box.numPts());
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        if (info.hasReduction()) {
            q.submit([&] (sycl::handler& h) {
                sycl::local_accessor<unsigned long long>
                    shared_data(sycl::range<1>(Gpu::Device::warp_size), h);
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                [[sycl::reqd_work_group_size(MT)]]
                [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
                {
                    for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                         icell < indexer.numPts(); icell += stride) {
                        auto iv = indexer.intVect(icell);
                        int n_active_threads = amrex::min(indexer.numPts()-icell+std::uint64_t(item.get_local_id(0)),
                                                          std::uint64_t(item.get_local_range(0)));
                        detail::call_f_intvect_ncomp_handler(f, iv, ncomp,
                                       Gpu::Handler{&item, shared_data.get_multi_ptr<sycl::access::decorated::yes>().get(),
                                                    n_active_threads});
                    }
                });
            });
        } else {
            q.submit([&] (sycl::handler& h) {
                h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                                 sycl::range<1>(nthreads_per_block)),
                [=] (sycl::nd_item<1> item)
                [[sycl::reqd_work_group_size(MT)]]
                [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
                {
                    for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                         icell < indexer.numPts(); icell += stride) {
                        auto iv = indexer.intVect(icell);
                        detail::call_f_intvect_ncomp_handler(f,iv,ncomp,Gpu::Handler{&item});
                    }
                });
            });
        }
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelForRNG (T n, L const& f) noexcept
{
    if (amrex::isEmpty(n)) { return; }
    const auto ec = Gpu::ExecutionConfig(n);
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * amrex::min(ec.numBlocks.x,Gpu::Device::maxBlocksPerLaunch());
    auto& q = Gpu::Device::streamQueue();
    auto& engdescr = *(getRandEngineDescriptor());
    try {
        q.submit([&] (sycl::handler& h) {
            auto engine_acc = engdescr.get_access(h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(AMREX_GPU_MAX_THREADS)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                auto const tid = item.get_global_id(0);
                auto engine = engine_acc.load(tid);
                RandomEngine rand_eng{&engine};
                for (std::size_t i = tid, stride = item.get_global_range(0); i < std::size_t(n); i += stride) {
                    f(T(i),rand_eng);
                }
                engine_acc.store(engine, tid);
            });
        });
        q.wait_and_throw(); // because next launch might be on a different queue
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename L, int dim>
void ParallelForRNG (BoxND<dim> const& box, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::ExecutionConfig(box.numPts());
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * amrex::min(ec.numBlocks.x,Gpu::Device::maxBlocksPerLaunch());
    auto& q = Gpu::Device::streamQueue();
    auto& engdescr = *(getRandEngineDescriptor());
    try {
        q.submit([&] (sycl::handler& h) {
            auto engine_acc = engdescr.get_access(h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(AMREX_GPU_MAX_THREADS)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                auto const tid = item.get_global_id(0);
                auto engine = engine_acc.load(tid);
                RandomEngine rand_eng{&engine};
                for (std::uint64_t icell = tid, stride = item.get_global_range(0);
                     icell < indexer.numPts(); icell += stride) {
                    auto iv = indexer.intVect(icell);
                    detail::call_f_intvect_engine(f,iv,rand_eng);
                }
                engine_acc.store(engine, tid);
            });
        });
        q.wait_and_throw(); // because next launch might be on a different queue
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelForRNG (BoxND<dim> const& box, T ncomp, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::ExecutionConfig(box.numPts());
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * amrex::min(ec.numBlocks.x,Gpu::Device::maxBlocksPerLaunch());
    auto& q = Gpu::Device::streamQueue();
    auto& engdescr = *(getRandEngineDescriptor());
    try {
        q.submit([&] (sycl::handler& h) {
            auto engine_acc = engdescr.get_access(h);
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(AMREX_GPU_MAX_THREADS)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                auto const tid = item.get_global_id(0);
                auto engine = engine_acc.load(tid);
                RandomEngine rand_eng{&engine};
                for (std::uint64_t icell = tid, stride = item.get_global_range(0);
                     icell < indexer.numPts(); icell += stride) {
                    auto iv = indexer.intVect(icell);
                    detail::call_f_intvect_ncomp_engine(f,iv,ncomp,rand_eng);
                }
                engine_acc.store(engine, tid);
            });
        });
        q.wait_and_throw(); // because next launch might be on a different queue
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename L1, typename L2, int dim>
void ParallelFor (Gpu::KernelInfo const& /*info*/, BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max(box1.numPts(), box2.numPts()));
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(MT)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                auto const ncells = std::max(indexer1.numPts(), indexer2.numPts());
                for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < indexer1.numPts()) {
                        auto iv = indexer1.intVect(icell);
                        detail::call_f_intvect(f1,iv);
                    }
                    if (icell < indexer2.numPts()) {
                        auto iv = indexer2.intVect(icell);
                        detail::call_f_intvect(f2,iv);
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void ParallelFor (Gpu::KernelInfo const& /*info*/,
                  BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                  L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const BoxIndexerND<dim> indexer3(box3);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max({box1.numPts(),box2.numPts(),box3.numPts()}));
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(MT)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                auto const ncells = std::max({indexer1.numPts(), indexer2.numPts(), indexer3.numPts()});
                for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < indexer1.numPts()) {
                        auto iv = indexer1.intVect(icell);
                        detail::call_f_intvect(f1,iv);
                    }
                    if (icell < indexer2.numPts()) {
                        auto iv = indexer2.intVect(icell);
                        detail::call_f_intvect(f2,iv);
                    }
                    if (icell < indexer3.numPts()) {
                        auto iv = indexer3.intVect(icell);
                        detail::call_f_intvect(f3,iv);
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (Gpu::KernelInfo const& /*info*/,
                  BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max(box1.numPts(),box2.numPts()));
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(MT)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                auto const ncells = std::max(indexer1.numPts(), indexer2.numPts());
                for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < indexer1.numPts()) {
                        auto iv = indexer1.intVect(icell);
                        detail::call_f_intvect_ncomp(f1,iv,ncomp1);
                    }
                    if (icell < indexer2.numPts()) {
                        auto iv = indexer2.intVect(icell);
                        detail::call_f_intvect_ncomp(f2,iv,ncomp2);
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (Gpu::KernelInfo const& /*info*/,
                  BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                  BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const BoxIndexerND<dim> indexer3(box3);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max({box1.numPts(),box2.numPts(),box3.numPts()}));
    const auto nthreads_per_block = ec.numThreads.x;
    const auto nthreads_total = std::size_t(nthreads_per_block) * ec.numBlocks.x;
    auto& q = Gpu::Device::streamQueue();
    try {
        q.submit([&] (sycl::handler& h) {
            h.parallel_for(sycl::nd_range<1>(sycl::range<1>(nthreads_total),
                                             sycl::range<1>(nthreads_per_block)),
            [=] (sycl::nd_item<1> item)
            [[sycl::reqd_work_group_size(MT)]]
            [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
            {
                auto const ncells = std::max({indexer1.numPts(), indexer2.numPts(), indexer3.numPts()});
                for (std::uint64_t icell = item.get_global_id(0), stride = item.get_global_range(0);
                     icell < ncells; icell += stride) {
                    if (icell < indexer1.numPts()) {
                        auto iv = indexer1.intVect(icell);
                        detail::call_f_intvect_ncomp(f1,iv,ncomp1);
                    }
                    if (icell < indexer2.numPts()) {
                        auto iv = indexer2.intVect(icell);
                        detail::call_f_intvect_ncomp(f2,iv,ncomp2);
                    }
                    if (icell < indexer3.numPts()) {
                        auto iv = indexer3.intVect(icell);
                        detail::call_f_intvect_ncomp(f3,iv,ncomp3);
                    }
                }
            });
        });
    } catch (sycl::exception const& ex) {
        amrex::Abort(std::string("ParallelFor: ")+ex.what()+"!!!!!");
    }
}

#else
// CUDA or HIP

template <typename L>
void single_task (gpuStream_t stream, L const& f) noexcept
{
    AMREX_LAUNCH_KERNEL(Gpu::Device::warp_size, 1, 1, 0, stream,
                        [=] AMREX_GPU_DEVICE () noexcept {f();});
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename L>
void launch (int nblocks, std::size_t shared_mem_bytes, gpuStream_t stream,
             L const& f) noexcept
{
    AMREX_LAUNCH_KERNEL(MT, nblocks, MT, shared_mem_bytes, stream,
                        [=] AMREX_GPU_DEVICE () noexcept { f(); });
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename L>
void launch (int nblocks, gpuStream_t stream, L const& f) noexcept
{
    AMREX_LAUNCH_KERNEL(MT, nblocks, MT, 0, stream,
                        [=] AMREX_GPU_DEVICE () noexcept { f(); });
    AMREX_GPU_ERROR_CHECK();
}

template<typename L>
void launch (int nblocks, int nthreads_per_block, std::size_t shared_mem_bytes,
             gpuStream_t stream, L const& f) noexcept
{
    AMREX_ASSERT(nthreads_per_block <= AMREX_GPU_MAX_THREADS);
    AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS, nblocks, nthreads_per_block, shared_mem_bytes,
                        stream, [=] AMREX_GPU_DEVICE () noexcept { f(); });
    AMREX_GPU_ERROR_CHECK();
}

template<typename L>
void launch (int nblocks, int nthreads_per_block, gpuStream_t stream, L&& f) noexcept
{
    launch(nblocks, nthreads_per_block, 0, stream, std::forward<L>(f));
}

template<int MT, typename T, typename L>
void launch (T const& n, L const& f) noexcept
{
    if (amrex::isEmpty(n)) { return; }
    const auto ec = Gpu::makeExecutionConfig<MT>(n);
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (auto const i : Gpu::Range(n)) {
            f(i);
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const&, T n, L const& f) noexcept
{
    if (amrex::isEmpty(n)) { return; }
    const auto ec = Gpu::makeExecutionConfig<MT>(n);
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (Long i = Long(blockDim.x)*blockIdx.x+threadIdx.x, stride = Long(blockDim.x)*gridDim.x;
             i < Long(n); i += stride) {
            detail::call_f_scalar_handler(f, T(i),
                Gpu::Handler(amrex::min((std::uint64_t(n)-i+(std::uint64_t)threadIdx.x),
                (std::uint64_t)blockDim.x)));
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename L, int dim>
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::makeExecutionConfig<MT>(box.numPts());
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (std::uint64_t icell = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x, stride = std::uint64_t(blockDim.x)*gridDim.x;
             icell < indexer.numPts(); icell += stride)
        {
            auto iv = indexer.intVect(icell);
            detail::call_f_intvect_handler(f, iv,
                Gpu::Handler(amrex::min((indexer.numPts()-icell+(std::uint64_t)threadIdx.x),
                (std::uint64_t)blockDim.x)));
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const&, BoxND<dim> const& box, T ncomp, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::makeExecutionConfig<MT>(box.numPts());
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        for (std::uint64_t icell = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x, stride = std::uint64_t(blockDim.x)*gridDim.x;
             icell < indexer.numPts(); icell += stride) {
            auto iv = indexer.intVect(icell);
            detail::call_f_intvect_ncomp_handler(f, iv, ncomp,
                Gpu::Handler(amrex::min((indexer.numPts()-icell+(std::uint64_t)threadIdx.x),
                (std::uint64_t)blockDim.x)));
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelForRNG (T n, L const& f) noexcept
{
    if (amrex::isEmpty(n)) { return; }
    randState_t* rand_state = getRandState();
    const auto ec = Gpu::ExecutionConfig(n);
    AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS,
                        amrex::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch()),
                        ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        Long tid = Long(blockDim.x)*blockIdx.x+threadIdx.x;
        RandomEngine engine{&(rand_state[tid])};
        for (Long i = tid, stride = Long(blockDim.x)*gridDim.x; i < Long(n); i += stride) {
            f(T(i),engine);
        }
    });
    Gpu::streamSynchronize(); // To avoid multiple streams using RNG
    AMREX_GPU_ERROR_CHECK();
}

template <typename L, int dim>
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelForRNG (BoxND<dim> const& box, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    randState_t* rand_state = getRandState();
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::ExecutionConfig(box.numPts());
    AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS,
                        amrex::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch()),
                        ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto const tid = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x;
        RandomEngine engine{&(rand_state[tid])};
        for (std::uint64_t icell = tid, stride = std::uint64_t(blockDim.x)*gridDim.x; icell < indexer.numPts(); icell += stride) {
            auto iv = indexer.intVect(icell);
            detail::call_f_intvect_engine(f, iv, engine);
        }
    });
    Gpu::streamSynchronize(); // To avoid multiple streams using RNG
    AMREX_GPU_ERROR_CHECK();
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelForRNG (BoxND<dim> const& box, T ncomp, L const& f) noexcept
{
    if (amrex::isEmpty(box)) { return; }
    randState_t* rand_state = getRandState();
    const BoxIndexerND<dim> indexer(box);
    const auto ec = Gpu::ExecutionConfig(box.numPts());
    AMREX_LAUNCH_KERNEL(AMREX_GPU_MAX_THREADS,
                        amrex::min(ec.numBlocks.x, Gpu::Device::maxBlocksPerLaunch()),
                        ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto const tid = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x;
        RandomEngine engine{&(rand_state[tid])};
        for (std::uint64_t icell = tid, stride = std::uint64_t(blockDim.x)*gridDim.x; icell < indexer.numPts(); icell += stride) {
            auto iv = indexer.intVect(icell);
            detail::call_f_intvect_ncomp_engine(f, iv, ncomp, engine);
        }
    });
    Gpu::streamSynchronize(); // To avoid multiple streams using RNG
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename L1, typename L2, int dim>
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
ParallelFor (Gpu::KernelInfo const&,
             BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max(box1.numPts(),box2.numPts()));
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto const ncells = std::max(indexer1.numPts(), indexer2.numPts());
        for (std::uint64_t icell = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x, stride = std::uint64_t(blockDim.x)*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < indexer1.numPts()) {
                auto iv = indexer1.intVect(icell);
                detail::call_f_intvect(f1, iv);
            }
            if (icell < indexer2.numPts()) {
                auto iv = indexer2.intVect(icell);
                detail::call_f_intvect(f2, iv);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename L1, typename L2, typename L3, int dim>
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value && MaybeDeviceRunnable<L3>::value>
ParallelFor (Gpu::KernelInfo const&,
             BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
             L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const BoxIndexerND<dim> indexer3(box3);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max({box1.numPts(),box2.numPts(),box3.numPts()}));
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto const ncells = std::max({indexer1.numPts(), indexer2.numPts(), indexer3.numPts()});
        for (std::uint64_t icell = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x, stride = std::uint64_t(blockDim.x)*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < indexer1.numPts()) {
                auto iv = indexer1.intVect(icell);
                detail::call_f_intvect(f1, iv);
            }
            if (icell < indexer2.numPts()) {
                auto iv = indexer2.intVect(icell);
                detail::call_f_intvect(f2, iv);
            }
            if (icell < indexer3.numPts()) {
                auto iv = indexer3.intVect(icell);
                detail::call_f_intvect(f3, iv);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
ParallelFor (Gpu::KernelInfo const&,
             BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
             BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max(box1.numPts(),box2.numPts()));
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto const ncells = std::max(indexer1.numPts(), indexer2.numPts());
        for (std::uint64_t icell = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x, stride = std::uint64_t(blockDim.x)*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < indexer1.numPts()) {
                auto iv = indexer1.intVect(icell);
                detail::call_f_intvect_ncomp(f1, iv, ncomp1);
            }
            if (icell < indexer2.numPts()) {
                auto iv = indexer2.intVect(icell);
                detail::call_f_intvect_ncomp(f2, iv, ncomp2);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value && MaybeDeviceRunnable<L3>::value>
ParallelFor (Gpu::KernelInfo const&,
             BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
             BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
             BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (amrex::isEmpty(box1) && amrex::isEmpty(box2) && amrex::isEmpty(box3)) { return; }
    const BoxIndexerND<dim> indexer1(box1);
    const BoxIndexerND<dim> indexer2(box2);
    const BoxIndexerND<dim> indexer3(box3);
    const auto ec = Gpu::makeExecutionConfig<MT>(std::max({box1.numPts(),box2.numPts(),box3.numPts()}));
    AMREX_LAUNCH_KERNEL(MT, ec.numBlocks, ec.numThreads, 0, Gpu::gpuStream(),
    [=] AMREX_GPU_DEVICE () noexcept {
        auto const ncells = std::max({indexer1.numPts(), indexer2.numPts(), indexer3.numPts()});
        for (std::uint64_t icell = std::uint64_t(blockDim.x)*blockIdx.x+threadIdx.x, stride = std::uint64_t(blockDim.x)*gridDim.x;
             icell < ncells; icell += stride) {
            if (icell < indexer1.numPts()) {
                auto iv = indexer1.intVect(icell);
                detail::call_f_intvect_ncomp(f1, iv, ncomp1);
            }
            if (icell < indexer2.numPts()) {
                auto iv = indexer2.intVect(icell);
                detail::call_f_intvect_ncomp(f2, iv, ncomp2);
            }
            if (icell < indexer3.numPts()) {
                auto iv = indexer3.intVect(icell);
                detail::call_f_intvect_ncomp(f3, iv, ncomp3);
            }
        }
    });
    AMREX_GPU_ERROR_CHECK();
}

#endif

template <typename L>
void single_task (L&& f) noexcept
{
    single_task(Gpu::gpuStream(), std::forward<L>(f));
}

template<typename T, typename L>
void launch (T const& n, L&& f) noexcept
{
    launch<AMREX_GPU_MAX_THREADS>(n, std::forward<L>(f));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, n, std::forward<L>(f));
}

template <typename L, int dim>
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, box, std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeDeviceRunnable<L>::value>
ParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, box, ncomp, std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
ParallelFor (Gpu::KernelInfo const& info,
             BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, box1, box2, std::forward<L1>(f1),
                                       std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value && MaybeDeviceRunnable<L3>::value>
ParallelFor (Gpu::KernelInfo const& info,
             BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
             L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, box1, box2, box3, std::forward<L1>(f1),
                                       std::forward<L2>(f2), std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value>
ParallelFor (Gpu::KernelInfo const& info,
             BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
             BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, box1, ncomp1, std::forward<L1>(f1),
                                             box2, ncomp2, std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
std::enable_if_t<MaybeDeviceRunnable<L1>::value && MaybeDeviceRunnable<L2>::value && MaybeDeviceRunnable<L3>::value>
ParallelFor (Gpu::KernelInfo const& info,
             BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
             BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
             BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, box1, ncomp1, std::forward<L1>(f1),
                                             box2, ncomp2, std::forward<L2>(f2),
                                             box3, ncomp3, std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    ParallelFor<MT>(info, n,std::forward<L>(f));
}

template <typename L, int dim>
void For (Gpu::KernelInfo const& info, BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info, box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void For (Gpu::KernelInfo const& info, BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor<MT>(info, box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info,box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor<MT>(info,box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor<MT>(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
          L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info,box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
          L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor<MT>(info,box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor<MT>(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(info,
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (Gpu::KernelInfo const& info,
          BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor<MT>(info,
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (T n, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{}, n, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (T n, L&& f) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{}, n, std::forward<L>(f));
}

template <typename L, int dim>
void ParallelFor (BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{}, box, std::forward<L>(f));
}

template <int MT, typename L, int dim>
void ParallelFor (BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{}, box, std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void ParallelFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                  L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void ParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                  L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                  BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void ParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                  BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                  BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (T n, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{}, n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (T n, L&& f) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{}, n,std::forward<L>(f));
}

template <typename L, int dim>
void For (BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{}, box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void For (BoxND<dim> const& box, L&& f) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{}, box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void For (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
          L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void For (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
          L1&& f1, L2&& f2, L3&& f3) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box1,box2,box3,std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void For (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
          BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
          BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    ParallelFor<MT>(Gpu::KernelInfo{},
                box1,ncomp1,std::forward<L1>(f1),
                box2,ncomp2,std::forward<L2>(f2),
                box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<AMREX_GPU_MAX_THREADS>(info,n,std::forward<L>(f));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        AMREX_PRAGMA_SIMD
        for (T i = 0; i < n; ++i) { f(i); }
#endif
    }
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<MT>(info,n,std::forward<L>(f));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        AMREX_PRAGMA_SIMD
        for (T i = 0; i < n; ++i) { f(i); }
#endif
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (T n, L&& f) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{}, n, std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (T n, L&& f) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{}, n, std::forward<L>(f));
}

template <typename L, int dim>
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<AMREX_GPU_MAX_THREADS>(info, box,std::forward<L>(f));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box,std::forward<L>(f));
#endif
    }
}

template <int MT, typename L, int dim>
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<MT>(info, box,std::forward<L>(f));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box,std::forward<L>(f));
#endif
    }
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<AMREX_GPU_MAX_THREADS>(info, box,ncomp,std::forward<L>(f));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box,ncomp,std::forward<L>(f));
#endif
    }
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<MT>(info, box,ncomp,std::forward<L>(f));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box,ncomp,std::forward<L>(f));
#endif
    }
}

template <typename L1, typename L2, int dim>
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<AMREX_GPU_MAX_THREADS>(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,std::forward<L2>(f2));
#endif
    }
}

template <int MT, typename L1, typename L2, int dim>
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<MT>(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,std::forward<L2>(f2));
#endif
    }
}

template <int MT, typename L1, typename L2, typename L3, int dim>
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value && MaybeHostDeviceRunnable<L3>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                       L1&& f1, L2&& f2, L3&& f3) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<MT>(info,box1,box2,box3,
                    std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,std::forward<L2>(f2));
        LoopConcurrentOnCpu(box3,std::forward<L3>(f3));
#endif
    }
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                       BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<AMREX_GPU_MAX_THREADS>(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
#endif
    }
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                       BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<MT>(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
#endif
    }
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value && MaybeHostDeviceRunnable<L3>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                       BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                       BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<AMREX_GPU_MAX_THREADS>(info,
                    box1,ncomp1,std::forward<L1>(f1),
                    box2,ncomp2,std::forward<L2>(f2),
                    box3,ncomp3,std::forward<L3>(f3));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
        LoopConcurrentOnCpu(box3,ncomp3,std::forward<L3>(f3));
#endif
    }
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
std::enable_if_t<MaybeHostDeviceRunnable<L1>::value && MaybeHostDeviceRunnable<L2>::value && MaybeHostDeviceRunnable<L3>::value>
HostDeviceParallelFor (Gpu::KernelInfo const& info,
                       BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                       BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                       BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    if (Gpu::inLaunchRegion()) {
        ParallelFor<MT>(info,
                    box1,ncomp1,std::forward<L1>(f1),
                    box2,ncomp2,std::forward<L2>(f2),
                    box3,ncomp3,std::forward<L3>(f3));
    } else {
#ifdef AMREX_USE_SYCL
        amrex::Abort("amrex:: HOST_DEVICE disabled for Intel.  It takes too long to compile");
#else
        LoopConcurrentOnCpu(box1,ncomp1,std::forward<L1>(f1));
        LoopConcurrentOnCpu(box2,ncomp2,std::forward<L2>(f2));
        LoopConcurrentOnCpu(box3,ncomp3,std::forward<L3>(f3));
#endif
    }
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(info,n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info, T n, L&& f) noexcept
{
    HostDeviceParallelFor<MT>(info,n,std::forward<L>(f));
}

template <typename L, int dim>
void HostDeviceFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, L&& f) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(info,box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void HostDeviceFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, L&& f) noexcept
{
    HostDeviceParallelFor<MT>(info,box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(info,box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info, BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    HostDeviceParallelFor<MT>(info,box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    HostDeviceParallelFor<MT>(info,box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(info, box1,box2,box3,
                          std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                    L1&& f1, L2&& f2, L3&& f3) noexcept
{
    HostDeviceParallelFor<MT>(info, box1,box2,box3,
                          std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    HostDeviceParallelFor<MT>(info,box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                    BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(info,
                          box1,ncomp1,std::forward<L1>(f1),
                          box2,ncomp2,std::forward<L2>(f2),
                          box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceFor (Gpu::KernelInfo const& info,
                    BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                    BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                    BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    HostDeviceParallelFor<MT>(info,
                          box1,ncomp1,std::forward<L1>(f1),
                          box2,ncomp2,std::forward<L2>(f2),
                          box3,ncomp3,std::forward<L3>(f3));
}

template <typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},n,std::forward<L>(f));
}

template <int MT, typename T, typename L, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (T n, L&& f) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{},n,std::forward<L>(f));
}

template <typename L, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box, L&& f) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box,std::forward<L>(f));
}

template <int MT, typename L, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box, L&& f) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{},box,std::forward<L>(f));
}

template <typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <int MT, typename T, typename L, int dim, typename M=std::enable_if_t<std::is_integral<T>::value> >
void HostDeviceParallelFor (BoxND<dim> const& box, T ncomp, L&& f) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{},box,ncomp,std::forward<L>(f));
}

template <typename L1, typename L2, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <int MT, typename L1, typename L2, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, L1&& f1, L2&& f2) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{},box1,box2,std::forward<L1>(f1),std::forward<L2>(f2));
}

template <typename L1, typename L2, typename L3, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{}, box1,box2,box3,
                          std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <int MT, typename L1, typename L2, typename L3, int dim>
void HostDeviceParallelFor (BoxND<dim> const& box1, BoxND<dim> const& box2, BoxND<dim> const& box3,
                            L1&& f1, L2&& f2, L3&& f3) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{}, box1,box2,box3,
                          std::forward<L1>(f1),std::forward<L2>(f2),std::forward<L3>(f3));
}

template <typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <int MT, typename T1, typename T2, typename L1, typename L2, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{},box1,ncomp1,std::forward<L1>(f1),box2,ncomp2,std::forward<L2>(f2));
}

template <typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                            BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    HostDeviceParallelFor<AMREX_GPU_MAX_THREADS>(Gpu::KernelInfo{},
                          box1,ncomp1,std::forward<L1>(f1),
                          box2,ncomp2,std::forward<L2>(f2),
                          box3,ncomp3,std::forward<L3>(f3));
}

template <int MT, typename T1, typename T2, typename T3, typename L1, typename L2, typename L3, int dim,
          typename M1=std::enable_if_t<std::is_integral<T1>::value>,
          typename M2=std::enable_if_t<std::is_integral<T2>::value>,
          typename M3=std::enable_if_t<std::is_integral<T3>::value> >
void HostDeviceParallelFor (BoxND<dim> const& box1, T1 ncomp1, L1&& f1,
                            BoxND<dim> const& box2, T2 ncomp2, L2&& f2,
                            BoxND<dim> const& box3, T3 ncomp3, L3&& f3) noexcept
{
    HostDeviceParallelFor<MT>(Gpu::KernelInfo{},
                          box1,ncomp1,std::forward<L1>(f1),
                          box2,ncomp2,std::forward<L2>(f2),
                          box3,ncomp3,std::forward<L3>(f3));
}

}

#endif
