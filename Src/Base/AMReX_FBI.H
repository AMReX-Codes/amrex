#ifndef AMREX_FBI_H_
#define AMREX_FBI_H_

template <class FAB>
struct FabCopyTag {
    FAB const* sfab;
    Box dbox;
    IntVect offset; // sbox.smallEnd() - dbox.smallEnd()
};

struct VoidCopyTag {
    char const* p;
    Box dbox;
};

namespace detail {

#ifdef AMREX_USE_GPU

template <class T0, class T1>
struct CellStore
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator() (T0* d, T1 s) const noexcept
    {
      *d = static_cast<T0>(s);
    }
};

template <class T0, class T1>
struct CellAdd
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator() (T0* d, T1 s) const noexcept
    {
        *d += static_cast<T0>(s);
    }
};

template <class T0, class T1>
struct CellAtomicAdd
{
    template<class U0=T0, std::enable_if_t<amrex::HasAtomicAdd<U0>::value,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator() (U0* d, T1 s) const noexcept
    {
        Gpu::Atomic::AddNoRet(d, static_cast<U0>(s));
    }
};

template <class T0, class T1, class F>
void
fab_to_fab (Vector<Array4CopyTag<T0, T1> > const& copy_tags, int scomp, int dcomp, int ncomp,
            F && f)
{
    TagVector<Array4CopyTag<T0, T1>> tv{copy_tags};

    detail::ParallelFor_doit(tv,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, int i, int j, int k, Array4CopyTag<T0, T1> const& tag) noexcept
    {
        if (icell < ncells) {
            for (int n = 0; n < ncomp; ++n) {
                f(&(tag.dfab(i,j,k,n+dcomp)),
                  tag.sfab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z,n+scomp));
            }
        }
    });
}

template <class T0, class T1, class F>
void
fab_to_fab (Vector<Array4CopyTag<T0, T1> > const& copy_tags, int scomp, int dcomp,
            int ncomp, F && f, Vector<Array4Tag<int> > const& masks)
{
    using TagType = Array4MaskCopyTag<T0, T1>;
    Vector<TagType> tags;
    const int N = copy_tags.size();
    tags.reserve(N);
    for (int i = 0; i < N; ++i) {
        tags.push_back(TagType{copy_tags[i].dfab, copy_tags[i].sfab, masks[i].dfab,
                               copy_tags[i].dbox, copy_tags[i].offset});
    }

    amrex::Abort("xxxxx TODO This function still has a bug.  Even if we fix the bug, it should still be avoided because it is slow due to the lack of atomic operations for this type.");

    TagVector<TagType> tv{tags};

    detail::ParallelFor_doit(tv,
    [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
        sycl::nd_item<1> const& item,
#endif
        int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
    {
#ifdef AMREX_USE_SYCL
        int g_tid = item.get_global_id(0);
        int g_wid = g_tid / Gpu::Device::warp_size;

        int* m = (icell < ncells) ? tag.mask.ptr(i,j,k) : nullptr;
        int mypriority = g_wid+1;
        int to_try  = 1;
        while (true) {
            int msk = (m && to_try) ? Gpu::Atomic::CAS(m, 0, mypriority) : 0;
            if (sycl::all_of_group(item.get_sub_group(), msk == 0)) {  // 0 means lock acquired
                break; // all threads have acquired.
            } else {
                if (sycl::any_of_group(item.get_sub_group(), msk > mypriority)) {
                    if (m) { *m = 0; } // yield
                    sycl::atomic_fence(sycl::memory_order::acq_rel, sycl::memory_scope::device);
                    to_try = 1;
                } else {
                    to_try = (msk > 0); // hold on to my lock
                }
            }
        };

        if (icell < ncells) {
            for (int n = 0; n < ncomp; ++n) {
                f(&(tag.dfab(i,j,k,n+dcomp)),
                  tag.sfab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z,n+scomp));
            }
        }

        if (m) *m = 0;

#else

        int g_tid = blockDim.x*blockIdx.x + threadIdx.x;
        int g_wid = g_tid / Gpu::Device::warp_size;

        int* m = (icell < ncells) ? tag.mask.ptr(i,j,k) : nullptr;
        int mypriority = g_wid+1;
        int to_try  = 1;
        while (true) {
            int msk = (m && to_try) ? atomicCAS(m, 0, mypriority) : 0;
#ifdef AMREX_USE_CUDA
            if (__all_sync(0xffffffff, msk == 0)) {  // 0 means lock acquired
#elif defined(AMREX_USE_HIP)
            if (__all(msk == 0)) {
#endif
                break; // all threads have acquired.
            } else {
#ifdef AMREX_USE_CUDA
                if (__any_sync(0xffffffff, msk > mypriority)) {
#elif defined(AMREX_USE_HIP)
                if (__any(msk > mypriority)) {
#endif
                    if (m) *m = 0; // yield
                    __threadfence();
                    to_try = 1;
                } else {
                    to_try = (msk > 0); // hold on to my lock
                }
            }
        };

        if (icell < ncells) {
            for (int n = 0; n < ncomp; ++n) {
                f(&(tag.dfab(i,j,k,n+dcomp)),
                  tag.sfab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z,n+scomp));
            }
        }

        if (m) *m = 0;
#endif
    });
}

template <typename T0, typename T1,
          std::enable_if_t<amrex::IsStoreAtomic<T0>::value,int> = 0>
void
fab_to_fab_atomic_cpy (Vector<Array4CopyTag<T0, T1> > const& copy_tags, int scomp,
                       int dcomp, int ncomp, Vector<Array4Tag<int> > const&)
{
    fab_to_fab<T0, T1>(copy_tags, scomp, dcomp, ncomp, CellStore<T0, T1>());
}

template <typename T0, typename T1,
          std::enable_if_t<!amrex::IsStoreAtomic<T0>::value,int> = 0>
void
fab_to_fab_atomic_cpy (Vector<Array4CopyTag<T0, T1> > const& copy_tags, int scomp,
                       int dcomp, int ncomp, Vector<Array4Tag<int> > const& masks)
{
    fab_to_fab(copy_tags, scomp, dcomp, ncomp, CellStore<T0, T1>(), masks);
}

template <typename T0, typename T1,
          std::enable_if_t<amrex::HasAtomicAdd<T0>::value,int> = 0>
void
fab_to_fab_atomic_add (Vector<Array4CopyTag<T0, T1> > const& copy_tags, int scomp,
                       int dcomp, int ncomp, Vector<Array4Tag<int> > const&)
{
    fab_to_fab(copy_tags, scomp, dcomp, ncomp, CellAtomicAdd<T0, T1>());
}

template <typename T0, typename T1,
          std::enable_if_t<!amrex::HasAtomicAdd<T0>::value,int> = 0>
void
fab_to_fab_atomic_add (Vector<Array4CopyTag<T0, T1> > const& copy_tags, int scomp,
                       int dcomp, int ncomp, Vector<Array4Tag<int> > const& masks)
{
    fab_to_fab(copy_tags, scomp, dcomp, ncomp, CellAdd<T0, T1>(), masks);
}

#endif /* AMREX_USE_GPU */

}

template <class FAB>
void
FabArray<FAB>::FB_local_copy_cpu (const FB& TheFB, int scomp, int ncomp)
{
    auto const& LocTags = *(TheFB.m_LocTags);
    auto N_locs = static_cast<int>(LocTags.size());
    if (N_locs == 0) { return; }
    bool is_thread_safe = TheFB.m_threadsafe_loc;
    if (is_thread_safe)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = LocTags[i];

            BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            const FAB* sfab = &(get(tag.srcIndex));
                  FAB* dfab = &(get(tag.dstIndex));
            dfab->template copy<RunOn::Host>(*sfab, tag.sbox, scomp, tag.dbox, scomp, ncomp);
        }
    }
    else
    {
        LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = LocTags[i];

            BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            loc_copy_tags[tag.dstIndex].push_back
                ({this->fabPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
        }
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            const auto& tags = loc_copy_tags[mfi];
            auto dfab = this->array(mfi);
            for (auto const & tag : tags)
            {
                auto const sfab = tag.sfab->array();
                const auto offset = tag.offset.dim3();
                amrex::LoopConcurrentOnCpu(tag.dbox, ncomp,
                [=] (int i, int j, int k, int n) noexcept
                {
                    dfab(i,j,k,n+scomp) = sfab(i+offset.x,j+offset.y,k+offset.z,n+scomp);
                });
            }
        }
    }
}

#ifdef AMREX_USE_GPU

template <class FAB>
auto FabArray<FAB>::FB_get_local_copy_tag_vector (const FB& TheFB)
    -> TagVector<Array4CopyTag<value_type>> const*
{
    auto const& LocTags = *(TheFB.m_LocTags);
    int N_locs = LocTags.size();

    using TagType = Array4CopyTag<value_type>;

    TagVector<TagType>* tv;
    if (auto it = m_local_copy_handler.find(TheFB.m_LocTags.get());
        it != m_local_copy_handler.end())
    {
        tv = it->second;
    } else {
        Vector<TagType> loc_copy_tags;
        loc_copy_tags.reserve(N_locs);

        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = LocTags[i];

            BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            int li = this->localindex(tag.dstIndex);
            loc_copy_tags.push_back
                ({this->atLocalIdx(li).array(),
                  this->fabPtr(tag.srcIndex)->const_array(),
                  tag.dbox,
                  (tag.sbox.smallEnd()-tag.dbox.smallEnd()).dim3()});
        }

        tv = new TagVector<TagType>(loc_copy_tags);
        // Note that during the lifetime of this FabArray LocTags's pointer
        // is stable.
        m_local_copy_handler[TheFB.m_LocTags.get()] = tv;
    }

    return tv;
}

template <class FAB>
void
FabArray<FAB>::FB_local_copy_gpu (const FB& TheFB, int scomp, int ncomp)
{
    auto const& LocTags = *(TheFB.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) { return; }
    bool is_thread_safe = TheFB.m_threadsafe_loc;

    if (!amrex::IsStoreAtomic<value_type>::value && !is_thread_safe) {
        amrex::Abort("xxxxx TODO This function still has a bug for non-atomic types. Even if we fix the bug, it should still be avoided because it is slow due to the lack of atomic operations for this type.");
    }

    using TagType = Array4CopyTag<value_type>;

    auto* tv = FB_get_local_copy_tag_vector(TheFB);

    detail::ParallelFor_doit(*tv,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
    {
        if (icell < ncells) {
            for (int n = 0; n < ncomp; ++n) {
                tag.dfab(i,j,k,n+scomp) = tag.sfab(i+tag.offset.x,
                                                   j+tag.offset.y,
                                                   k+tag.offset.z,n+scomp);
            }
        }
    });
}

template <class FAB>
void
FabArray<FAB>::CMD_local_setVal_gpu (typename FabArray<FAB>::value_type x,
                                    const CommMetaData& thecmd, int scomp, int ncomp)
{
    auto const& LocTags = *(thecmd.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) { return; }
    bool is_thread_safe = thecmd.m_threadsafe_loc;

    using TagType = Array4BoxTag<value_type>;
    Vector<TagType> loc_setval_tags;
    loc_setval_tags.reserve(N_locs);

    AMREX_ALWAYS_ASSERT(amrex::IsStoreAtomic<value_type>::value || is_thread_safe);

    for (int i = 0; i < N_locs; ++i)
    {
        const CopyComTag& tag = LocTags[i];
        BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
        loc_setval_tags.push_back({this->array(tag.dstIndex), tag.dbox});
    }

    amrex::ParallelFor(loc_setval_tags, ncomp,
    [x,scomp] AMREX_GPU_DEVICE (int i, int j, int k, int n, TagType const& tag) noexcept
    {
        tag.dfab(i,j,k,n+scomp) = x;
    });
}

template <class FAB>
void
FabArray<FAB>::CMD_remote_setVal_gpu (typename FabArray<FAB>::value_type x,
                                    const CommMetaData& thecmd, int scomp, int ncomp)
{
    auto const& RcvTags = *(thecmd.m_RcvTags);
    bool is_thread_safe = thecmd.m_threadsafe_rcv;

    using TagType = Array4BoxTag<value_type>;
    Vector<TagType> rcv_setval_tags;

    for (auto it = RcvTags.begin(); it != RcvTags.end(); ++it) {
        for (auto const& tag: it->second) {
            rcv_setval_tags.push_back({this->array(tag.dstIndex), tag.dbox});
        }
    }

    if (rcv_setval_tags.empty()) { return; }

    AMREX_ALWAYS_ASSERT(amrex::IsStoreAtomic<value_type>::value || is_thread_safe);

    amrex::ParallelFor(rcv_setval_tags, ncomp,
    [x,scomp] AMREX_GPU_DEVICE (int i, int j, int k, int n, TagType const& tag) noexcept
    {
        tag.dfab(i,j,k,n+scomp) = x;
    });
}

#if defined(__CUDACC__) && defined (AMREX_USE_CUDA)
template <class FAB>
void
FabArray<FAB>::FB_local_copy_cuda_graph_1 (const FB& TheFB, int scomp, int ncomp)
{
    const int N_locs = (*TheFB.m_LocTags).size();
    LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
    for (int i = 0; i < N_locs; ++i)
    {
        const CopyComTag& tag = (*TheFB.m_LocTags)[i];

        BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
        BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

        loc_copy_tags[tag.dstIndex].push_back
            ({this->fabPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
    }

    // Create Graph if one is needed.
    if ( !(TheFB.m_localCopy.ready()) )
    {
        const_cast<FB&>(TheFB).m_localCopy.resize(N_locs);

        int idx = 0;
        // Record the graph.
        for (MFIter mfi(*this, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi)
        {
            amrex::Gpu::Device::startGraphRecording( (mfi.LocalIndex() == 0),
                                                     const_cast<FB&>(TheFB).m_localCopy.getHostPtr(0),
                                                     (TheFB).m_localCopy.getDevicePtr(0),
                                                     std::size_t(sizeof(CopyMemory)*N_locs) );

            const auto& tags = loc_copy_tags[mfi];
            for (auto const & tag : tags)
            {
                const auto offset = tag.offset.dim3();
                CopyMemory* cmem = TheFB.m_localCopy.getDevicePtr(idx++);
                AMREX_HOST_DEVICE_FOR_3D (tag.dbox, i, j, k,
                {
                    // Build the Array4's.
                    auto const dst = cmem->getDst<value_type>();
                    auto const src = cmem->getSrc<value_type>();
                    for (int n = 0; n < cmem->ncomp; ++n) {
                        dst(i,j,k,(cmem->scomp)+n) = src(i+offset.x,j+offset.y,k+offset.z,(cmem->scomp)+n);
                    }
                });
            }

            bool last_iter = mfi.LocalIndex() == (this->local_size()-1);
            cudaGraphExec_t graphExec = amrex::Gpu::Device::stopGraphRecording(last_iter);
            if (last_iter) { const_cast<FB&>(TheFB).m_localCopy.setGraph( graphExec ); }
        }
    }

    // Setup Launch Parameters
    // This is perfectly threadable, right?
    // Additional optimization -> Check to see whether values need to be reset?
    // Can then remove this setup and memcpy from CudaGraph::executeGraph.
    int idx = 0;
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        auto const dst_array = this->array(mfi);
        const auto& tags = loc_copy_tags[mfi];
        for (auto const & tag : tags)
        {
            const_cast<FB&>(TheFB).m_localCopy.setParams(idx++, makeCopyMemory(tag.sfab->array(),
                                                                               dst_array,
                                                                               scomp, ncomp));
        }
    }

    // Launch Graph
    TheFB.m_localCopy.executeGraph();
}

#ifdef AMREX_USE_MPI
template <class FAB>
void
FabArray<FAB>::FB_local_copy_cuda_graph_n (const FB& TheFB, int scomp, int ncomp)
{
    const int N_locs = TheFB.m_LocTags->size();

    int launches = 0; // Used for graphs only.
    LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
    for (int i = 0; i < N_locs; ++i)
    {
        const CopyComTag& tag = (*TheFB.m_LocTags)[i];

        BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.dstIndex]));
        BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.srcIndex]));

        if (distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc())
        {
            loc_copy_tags[tag.dstIndex].push_back
                ({this->fabPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
            launches++;
        }
    }

    FillBoundary_test();

    if ( !(TheFB.m_localCopy.ready()) )
    {
        const_cast<FB&>(TheFB).m_localCopy.resize(launches);

        int idx = 0;
        int cuda_stream = 0;
        for (MFIter mfi(*this, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi)
        {
            const auto& tags = loc_copy_tags[mfi];
            for (int t = 0; t<tags.size(); ++t)
            {
                Gpu::Device::setStreamIndex(cuda_stream++);
                amrex::Gpu::Device::startGraphRecording( (idx == 0),
                                                         const_cast<FB&>(TheFB).m_localCopy.getHostPtr(0),
                                                         (TheFB).m_localCopy.getDevicePtr(0),
                                                         std::size_t(sizeof(CopyMemory)*launches) );

                const auto& tag = tags[t];
                const Dim3 offset = tag.offset.dim3();

                CopyMemory* cmem = TheFB.m_localCopy.getDevicePtr(idx++);
                AMREX_HOST_DEVICE_FOR_3D(tag.dbox, i, j, k,
                {
                    auto const dst = cmem->getDst<value_type>();
                    auto const src = cmem->getSrc<value_type>();
                    for (int n = 0; n < cmem->ncomp; ++n) {
                        dst(i,j,k,(cmem->scomp)+n) = src(i+offset.x,j+offset.y,k+offset.z,(cmem->scomp)+n);
                    }
                });

                bool last_iter = idx == launches;
                cudaGraphExec_t graphExec = Gpu::Device::stopGraphRecording(last_iter);
                if (last_iter) { const_cast<FB&>(TheFB).m_localCopy.setGraph( graphExec ); }
            }
        }
    }

    // Setup Launch Parameters
    // This is perfectly threadable, right?
    int idx = 0;
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        const auto& dst_array = this->array(mfi);
        const auto& tags = loc_copy_tags[mfi];
        for (auto const & tag : tags)
        {
            const_cast<FB&>(TheFB).m_localCopy.setParams(idx++, makeCopyMemory(tag.sfab->array(),
                                                                               dst_array,
                                                                               scomp, ncomp));
        }
    }

    // Launch Graph without synch. Local work is entirely independent.
    TheFB.m_localCopy.executeGraph(false);
}
#endif /* AMREX_USE_MPI */

#endif /* __CUDACC__ */

#endif /* AMREX_USE_GPU */

#ifdef AMREX_USE_MPI

#ifdef AMREX_USE_GPU

#if defined(__CUDACC__) && defined(AMREX_USE_CUDA)

template <class FAB>
void
FabArray<FAB>::FB_pack_send_buffer_cuda_graph (const FB& TheFB, int scomp, int ncomp,
                                               Vector<char*>& send_data,
                                               Vector<std::size_t> const& send_size,
                                               Vector<typename FabArray<FAB>::CopyComTagsContainer const*> const& send_cctc)
{
    const int N_snds = send_data.size();
    if (N_snds == 0) { return; }

    if ( !(TheFB.m_copyToBuffer.ready()) )
    {
        // Set size of CudaGraph buffer.
        // Is the conditional ever expected false?
        int launches = 0;
        for (int send = 0; send < N_snds; ++send) {
            if (send_size[send] > 0) {
                launches += send_cctc[send]->size();
            }
        }
        const_cast<FB&>(TheFB).m_copyToBuffer.resize(launches);

        // Record the graph.
        int idx = 0;
        for (Gpu::StreamIter sit(N_snds,Gpu::StreamItInfo().DisableDeviceSync());
             sit.isValid(); ++sit)
        {
            amrex::Gpu::Device::startGraphRecording( (sit() == 0),
                                                     const_cast<FB&>(TheFB).m_copyToBuffer.getHostPtr(0),
                                                     (TheFB).m_copyToBuffer.getDevicePtr(0),
                                                     std::size_t(sizeof(CopyMemory)*launches) );

            const int j = sit();
            if (send_size[j] > 0)
            {
                auto const& cctc = *send_cctc[j];
                for (auto const& tag : cctc)
                {
                    const Box& bx = tag.sbox;
                    CopyMemory* cmem = TheFB.m_copyToBuffer.getDevicePtr(idx++);
                    AMREX_HOST_DEVICE_FOR_3D (bx, ii, jj, kk,
                    {
                        auto const pfab = cmem->getDst<value_type>();
                        auto const sfab = cmem->getSrc<value_type>();
                        for (int n = 0; n < cmem->ncomp; ++n)
                        {
                            pfab(ii,jj,kk,n) = sfab(ii,jj,kk,n+(cmem->scomp));
                        }
                    });
                }
            }

            bool last_iter = sit() == (N_snds-1);
            cudaGraphExec_t graphExec = amrex::Gpu::Device::stopGraphRecording(last_iter);
            if (last_iter) { const_cast<FB&>(TheFB).m_copyToBuffer.setGraph( graphExec ); }
        }
    }

    // Setup Launch Parameters
    int idx = 0;
    for (int send = 0; send < N_snds; ++send)
    {
        const int j = send;
        if (send_size[j] > 0)
        {
            char* dptr = send_data[j];
            auto const& cctc = *send_cctc[j];
            for (auto const& tag : cctc)
            {
                const_cast<FB&>(TheFB).m_copyToBuffer.setParams(idx++, makeCopyMemory(this->array(tag.srcIndex),
                                                                                       amrex::makeArray4((value_type*)(dptr),
                                                                                                         tag.sbox,
                                                                                                         ncomp),
                                                                                       scomp, ncomp));

                dptr += (tag.sbox.numPts() * ncomp * sizeof(value_type));
            }
            amrex::ignore_unused(send_size);
            BL_ASSERT(dptr <= send_data[j] + send_size[j]);
        }
    }

    // Launch Graph synched, so copyToBuffer is complete prior to posting sends.
    TheFB.m_copyToBuffer.executeGraph();
}

template <class FAB>
void
FabArray<FAB>::FB_unpack_recv_buffer_cuda_graph (const FB& TheFB, int dcomp, int ncomp,
                                                 Vector<char*> const& recv_data,
                                                 Vector<std::size_t> const& recv_size,
                                                 Vector<CopyComTagsContainer const*> const& recv_cctc,
                                                 bool /*is_thread_safe*/)
{
    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) { return; }

    int launches = 0;
    LayoutData<Vector<VoidCopyTag> > recv_copy_tags(boxArray(),DistributionMap());
    for (int k = 0; k < N_rcvs; ++k)
    {
        if (recv_size[k] > 0)
        {
            const char* dptr = recv_data[k];
            auto const& cctc = *recv_cctc[k];
            for (auto const& tag : cctc)
            {
                recv_copy_tags[tag.dstIndex].push_back({dptr,tag.dbox});
                dptr += tag.dbox.numPts() * ncomp * sizeof(value_type);
                launches++;
            }
            amrex::ignore_unused(recv_size);
            BL_ASSERT(dptr <= recv_data[k] + recv_size[k]);
        }
    }

    if ( !(TheFB.m_copyFromBuffer.ready()) )
    {
        const_cast<FB&>(TheFB).m_copyFromBuffer.resize(launches);

        int idx = 0;
        for (MFIter mfi(*this, MFItInfo().DisableDeviceSync()); mfi.isValid(); ++mfi)
        {
            amrex::Gpu::Device::startGraphRecording( (mfi.LocalIndex() == 0),
                                                     const_cast<FB&>(TheFB).m_copyFromBuffer.getHostPtr(0),
                                                     (TheFB).m_copyFromBuffer.getDevicePtr(0),
                                                     std::size_t(sizeof(CopyMemory)*launches) );

            const auto& tags = recv_copy_tags[mfi];
            for (auto const & tag : tags)
            {
                CopyMemory* cmem = TheFB.m_copyFromBuffer.getDevicePtr(idx++);
                AMREX_HOST_DEVICE_FOR_3D (tag.dbox, i, j, k,
                {
                    auto const pfab = cmem->getSrc<value_type>();
                    auto const dfab = cmem->getDst<value_type>();
                    for (int n = 0; n < cmem->ncomp; ++n)
                    {
                        dfab(i,j,k,n+(cmem->scomp)) = pfab(i,j,k,n);
                    }
                });
            }

            bool last_iter = mfi.LocalIndex() == (this->local_size()-1);
            cudaGraphExec_t graphExec = amrex::Gpu::Device::stopGraphRecording(last_iter);
            if (last_iter) { const_cast<FB&>(TheFB).m_copyFromBuffer.setGraph( graphExec ); }
        }
    }

    // Setup graph.
    int idx = 0;
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        auto dst_array = this->array(mfi);
        const auto & tags = recv_copy_tags[mfi];
        for (auto const & tag : tags)
        {
            const_cast<FB&>(TheFB).m_copyFromBuffer.setParams(idx++, makeCopyMemory(amrex::makeArray4((value_type*)(tag.p),
                                                                                                      tag.dbox,
                                                                                                      ncomp),
                                                                                    dst_array,
                                                                                    dcomp, ncomp));
        }
    }

    // Launch Graph - synced because next action is freeing recv buffer.
    TheFB.m_copyFromBuffer.executeGraph();
}

#endif /* __CUDACC__ */

template <class FAB>
template <typename BUF>
auto
FabArray<FAB>::get_send_copy_tag_vector (Vector<char*> const& send_data,
                                         Vector<std::size_t> const& send_size,
                                         Vector<CopyComTagsContainer const*> const& send_cctc,
                                         int ncomp) const
    -> TagVector<Array4CopyTag<char,value_type>> const*
{
    using TagType = Array4CopyTag<char,value_type>;

    auto kit = std::find_if(send_cctc.begin(), send_cctc.end(),
                            [] (CopyComTagsContainer const* p) { return p != nullptr; });
    if (kit == send_cctc.end()) {
        return nullptr;
    }

    auto key = std::make_pair(*kit, std::make_pair(sizeof(BUF),ncomp));

    TagVector<TagType>* tv;
    if (auto it = m_send_copy_handler.find(key); it != m_send_copy_handler.end()) {
        tv = it->second;
    } else {
        const int N_snds = send_data.size();
        Vector<TagType> snd_copy_tags;
        snd_copy_tags.reserve(N_snds);

        for (int j = 0; j < N_snds; ++j)
        {
            if (send_size[j] > 0)
            {
                std::intptr_t offset = send_data[j]-send_data[0];
                auto* dptr = reinterpret_cast<char*>(offset);
                auto const& cctc = *send_cctc[j];
                for (auto const& tag : cctc)
                {
                    snd_copy_tags.emplace_back(TagType{
                            amrex::makeArray4(dptr, tag.sbox, ncomp),
                            this->array(tag.srcIndex),
                            tag.sbox,
                            Dim3{0,0,0}
                        });
                    dptr += (tag.sbox.numPts() * ncomp * sizeof(BUF));
                }
            }
        }

        if (m_send_copy_handler.size() > 64) {
            // Just in case. If this is used in ParallelCopy, it's possible
            // that the size of this map may increase indefinitely.
            Gpu::streamSynchronize();
            for (auto& [k,v] : m_send_copy_handler) {
                delete v;
            }
            m_send_copy_handler.clear();
        }

        tv = new TagVector<TagType>(snd_copy_tags);
        m_send_copy_handler[key] = tv;
    }

    return tv;
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::pack_send_buffer_gpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                     Vector<char*> const& send_data,
                                     Vector<std::size_t> const& send_size,
                                     Vector<CopyComTagsContainer const*> const& send_cctc)
{
    const int N_snds = send_data.size();
    if (N_snds == 0) { return; }

    char* pbuffer = send_data[0];
    using TagType = Array4CopyTag<char, value_type>;

    auto* tv = src.template get_send_copy_tag_vector<BUF>(send_data, send_size,
                                                          send_cctc, ncomp);
    if (tv == nullptr) { return; }

    detail::ParallelFor_doit(*tv,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
    {
        if (icell < ncells) {
            Array4<BUF> dfab{reinterpret_cast<BUF*>
                (pbuffer + reinterpret_cast<std::intptr_t>(tag.dfab.p)),
                tag.dfab.begin, tag.dfab.end, tag.dfab.ncomp};
            for (int n = 0; n < dfab.ncomp; ++n) {
                dfab(i,j,k,n) = (BUF)tag.sfab(i+tag.offset.x,
                                              j+tag.offset.y,
                                              k+tag.offset.z,n+scomp);
            }
        }
    });

    Gpu::streamSynchronize();
}

template <class FAB>
template <typename BUF>
auto
FabArray<FAB>::get_recv_copy_tag_vector (Vector<char*> const& recv_data,
                                         Vector<std::size_t> const& recv_size,
                                         Vector<CopyComTagsContainer const*> const& recv_cctc,
                                         int ncomp)
    -> TagVector<Array4CopyTag<value_type,char>> const*
{
    using TagType = Array4CopyTag<value_type,char>;

    auto kit = std::find_if(recv_cctc.begin(), recv_cctc.end(),
                            [] (CopyComTagsContainer const* p) { return p != nullptr; });
    if (kit == recv_cctc.end()) {
        return nullptr;
    }

    auto key = std::make_pair(*kit, std::make_pair(sizeof(BUF),ncomp));

    TagVector<TagType>* tv;
    if (auto it = m_recv_copy_handler.find(key); it != m_recv_copy_handler.end()) {
        tv = it->second;
    } else {
        const int N_rcvs = recv_cctc.size();
        Vector<TagType> recv_copy_tags;
        recv_copy_tags.reserve(N_rcvs);

        for (int k = 0; k < N_rcvs; ++k)
        {
            if (recv_size[k] > 0)
            {
                std::intptr_t offset = recv_data[k]-recv_data[0];
                auto* dptr = reinterpret_cast<char const*>(offset);
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    const int li = this->localindex(tag.dstIndex);
                    recv_copy_tags.emplace_back(TagType{
                            this->atLocalIdx(li).array(),
                            amrex::makeArray4(dptr, tag.dbox, ncomp),
                            tag.dbox,
                            Dim3{0,0,0}
                        });
                    dptr += tag.dbox.numPts() * ncomp * sizeof(BUF);
                }
            }
        }

        if (m_recv_copy_handler.size() > 64) {
            // Just in case. If this is used in ParallelCopy, it's possible
            // that the size of this map may increase indefinitely.
            Gpu::streamSynchronize();
            for (auto& [k,v] : m_recv_copy_handler) {
                delete v;
            }
            m_recv_copy_handler.clear();
        }

        tv = new TagVector<TagType>(recv_copy_tags);
        m_recv_copy_handler[key] = tv;
    }

    return tv;
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::unpack_recv_buffer_gpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                       Vector<char*> const& recv_data,
                                       Vector<std::size_t> const& recv_size,
                                       Vector<CopyComTagsContainer const*> const& recv_cctc,
                                       CpOp op, bool is_thread_safe)
{
    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) { return; }

    if (!is_thread_safe)
    {
        if ((op == FabArrayBase::COPY && !amrex::IsStoreAtomic<value_type>::value) ||
            (op == FabArrayBase::ADD  && !amrex::HasAtomicAdd <value_type>::value))
        {
            amrex::Abort("xxxxx TODO This function still has a bug for non-atomic types. Even if we fix the bug, it should still be avoided because it is slow due to the lack of atomic operations for this type.");
        }
    }

    char* pbuffer = recv_data[0];
    using TagType = Array4CopyTag<value_type, char>;

    auto* tv = dst.template get_recv_copy_tag_vector<BUF>(recv_data, recv_size,
                                                          recv_cctc, ncomp);
    if (tv == nullptr) { return; }

    if (op == FabArrayBase::COPY)
    {
        detail::ParallelFor_doit(*tv,
            [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
                sycl::nd_item<1> const& /*item*/,
#endif
                int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
        {
            if (icell < ncells) {
                Array4<BUF const> sfab{reinterpret_cast<BUF const*>
                    (pbuffer + reinterpret_cast<std::intptr_t>(tag.sfab.p)),
                    tag.sfab.begin, tag.sfab.end, tag.sfab.ncomp};
                for (int n = 0; n < sfab.ncomp; ++n) {
                    tag.dfab(i,j,k,n+dcomp) = (value_type)sfab(i,j,k,n);
                }
            }
        });
    }
    else
    {
        if (is_thread_safe) {
            detail::ParallelFor_doit(*tv,
                [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
                    sycl::nd_item<1> const& /*item*/,
#endif
                    int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
            {
                if (icell < ncells) {
                    Array4<BUF const> sfab{reinterpret_cast<BUF const*>
                        (pbuffer + reinterpret_cast<std::intptr_t>(tag.sfab.p)),
                        tag.sfab.begin, tag.sfab.end, tag.sfab.ncomp};
                    for (int n = 0; n < sfab.ncomp; ++n) {
                        tag.dfab(i,j,k,n+dcomp) += (value_type)sfab(i,j,k,n);
                    }
                }
            });
        } else {
            if constexpr (amrex::HasAtomicAdd<value_type>::value) {
                detail::ParallelFor_doit(*tv,
                    [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_SYCL
                        sycl::nd_item<1> const& /*item*/,
#endif
                        int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
                {
                    if (icell < ncells) {
                        Array4<BUF const> sfab{reinterpret_cast<BUF const*>
                            (pbuffer + reinterpret_cast<std::intptr_t>(tag.sfab.p)),
                            tag.sfab.begin, tag.sfab.end, tag.sfab.ncomp};
                        for (int n = 0; n < sfab.ncomp; ++n) {
                            Gpu::Atomic::AddNoRet(tag.dfab.ptr(i,j,k,n+dcomp),
                                                  (value_type)sfab(i,j,k,n));
                        }
                    }
                });
            } else {
                amrex::Abort("unpack_recv_buffer_gpu: should NOT get here");
            }
        }
    }

    Gpu::streamSynchronize();
}

#endif /* AMREX_USE_GPU */

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::pack_send_buffer_cpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                     Vector<char*> const& send_data,
                                     Vector<std::size_t> const& send_size,
                                     Vector<CopyComTagsContainer const*> const& send_cctc)
{
    amrex::ignore_unused(send_size);

    auto const N_snds = static_cast<int>(send_data.size());
    if (N_snds == 0) { return; }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int j = 0; j < N_snds; ++j)
    {
        if (send_size[j] > 0)
        {
            char* dptr = send_data[j];
            auto const& cctc = *send_cctc[j];
            for (auto const& tag : cctc)
            {
                const Box& bx = tag.sbox;
                auto const sfab = src.array(tag.srcIndex);
                auto pfab = amrex::makeArray4((BUF*)(dptr),bx,ncomp);
                amrex::LoopConcurrentOnCpu( bx, ncomp,
                [=] (int ii, int jj, int kk, int n) noexcept
                {
                    pfab(ii,jj,kk,n) = static_cast<BUF>(sfab(ii,jj,kk,n+scomp));
                });
                dptr += (bx.numPts() * ncomp * sizeof(BUF));
            }
            BL_ASSERT(dptr <= send_data[j] + send_size[j]);
        }
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::unpack_recv_buffer_cpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                       Vector<char*> const& recv_data,
                                       Vector<std::size_t> const& recv_size,
                                       Vector<CopyComTagsContainer const*> const& recv_cctc,
                                       CpOp op, bool is_thread_safe)
{
    amrex::ignore_unused(recv_size);

    auto const N_rcvs = static_cast<int>(recv_cctc.size());
    if (N_rcvs == 0) { return; }

    if (is_thread_safe)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
        for (int k = 0; k < N_rcvs; ++k)
        {
            if (recv_size[k] > 0)
            {
                const char* dptr = recv_data[k];
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    const Box& bx  = tag.dbox;
                    FAB& dfab = dst[tag.dstIndex];
                    if (op == FabArrayBase::COPY)
                    {
                        dfab.template copyFromMem<RunOn::Host, BUF>(bx, dcomp, ncomp, dptr);
                    }
                    else
                    {
                        dfab.template addFromMem<RunOn::Host, BUF>(tag.dbox, dcomp, ncomp, dptr);
                    }
                    dptr += bx.numPts() * ncomp * sizeof(BUF);
                }
                BL_ASSERT(dptr <= recv_data[k] + recv_size[k]);
            }
        }
    }
    else
    {
        LayoutData<Vector<VoidCopyTag> > recv_copy_tags;
        recv_copy_tags.define(dst.boxArray(),dst.DistributionMap());
        for (int k = 0; k < N_rcvs; ++k)
        {
            if (recv_size[k] > 0)
            {
                const char* dptr = recv_data[k];
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    recv_copy_tags[tag.dstIndex].push_back({dptr,tag.dbox});
                    dptr += tag.dbox.numPts() * ncomp * sizeof(BUF);
                }
                BL_ASSERT(dptr <= recv_data[k] + recv_size[k]);
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi(dst); mfi.isValid(); ++mfi)
        {
            const auto& tags = recv_copy_tags[mfi];
            auto dfab = dst.array(mfi);
            for (auto const & tag : tags)
            {
                auto pfab = amrex::makeArray4((BUF*)(tag.p), tag.dbox, ncomp);
                if (op == FabArrayBase::COPY)
                {
                    amrex::LoopConcurrentOnCpu(tag.dbox, ncomp,
                    [=] (int i, int j, int k, int n) noexcept
                    {
                        dfab(i,j,k,n+dcomp) = pfab(i,j,k,n);
                    });
                }
                else
                {
                    amrex::LoopConcurrentOnCpu(tag.dbox, ncomp,
                    [=] (int i, int j, int k, int n) noexcept
                    {
                        dfab(i,j,k,n+dcomp) += pfab(i,j,k,n);
                    });
                }
            }
        }
    }
}

#endif /* AMREX_USE_MPI */

#endif
