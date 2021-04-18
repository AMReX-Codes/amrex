#ifndef AMREX_FBI_H_
#define AMREX_FBI_H_

template <class FAB>
struct FabCopyTag {
    FAB const* sfab;
    Box dbox;
    IntVect offset; // sbox.smallEnd() - dbox.smallEnd()
};

template <class T>
struct Array4CopyTag {
    Array4<T      > dfab;
    Array4<T const> sfab;
    Box dbox;
    Dim3 offset; // sbox.smallEnd() - dbox.smallEnd()

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class T>
struct Array4MaskCopyTag {
    Array4<T      > dfab;
    Array4<T const> sfab;
    Array4<int    > mask;
    Box dbox;
    Dim3 offset; // sbox.smallEnd() - dbox.smallEnd()

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

struct VoidCopyTag {
    char const* p;
    Box dbox;
};

template <class T>
struct Array4BoxTag {
    Array4<T> dfab;
    Box       dbox;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

#ifdef AMREX_USE_GPU
namespace detail {
template <class TagType, class F>
void
ParallelFor_doit (Vector<TagType> const& tags, F && f)
{
    const int ntags = tags.size();
    if (ntags == 0) return;

    int ntotwarps = 0;
    Vector<int> nwarps;
    nwarps.reserve(ntags+1);
    for (int i = 0; i < ntags; ++i)
    {
        auto& tag = tags[i];
        nwarps.push_back(ntotwarps);
        ntotwarps += static_cast<int>(tag.box().numPts()+Gpu::Device::warp_size-1)/Gpu::Device::warp_size;
    }
    nwarps.push_back(ntotwarps);

    std::size_t sizeof_tags = ntags*sizeof(TagType);
    std::size_t offset_nwarps = Arena::align(sizeof_tags);
    std::size_t sizeof_nwarps = (ntags+1)*sizeof(int);
    std::size_t total_buf_size = offset_nwarps + sizeof_nwarps;

    char* h_buffer = (char*)The_Pinned_Arena()->alloc(total_buf_size);
    char* d_buffer = (char*)The_Arena()->alloc(total_buf_size);

    std::memcpy(h_buffer, tags.data(), sizeof_tags);
    std::memcpy(h_buffer+offset_nwarps, nwarps.data(), sizeof_nwarps);
    Gpu::htod_memcpy_async(d_buffer, h_buffer, total_buf_size);

    auto d_tags = reinterpret_cast<TagType*>(d_buffer);
    auto d_nwarps = reinterpret_cast<int*>(d_buffer+offset_nwarps);

    constexpr int nthreads = 256;
    constexpr int nwarps_per_block = nthreads/Gpu::Device::warp_size;
    int nblocks = (ntotwarps + nwarps_per_block-1) / nwarps_per_block;

    amrex::launch(nblocks, nthreads, Gpu::gpuStream(),
#ifdef AMREX_USE_DPCPP
    [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
    AMREX_REQUIRE_SUBGROUP_SIZE(Gpu::Device::warp_size)
#else
    [=] AMREX_GPU_DEVICE () noexcept
#endif
    {
#ifdef AMREX_USE_DPCPP
        int g_tid = item.get_global_id(0);
#else
        int g_tid = blockDim.x*blockIdx.x + threadIdx.x;
#endif
        int g_wid = g_tid / Gpu::Device::warp_size;
        if (g_wid >= ntotwarps) return;

        int tag_id = -10000;
        {
            int lo = 0;
            int hi = ntags;
            while (lo <= hi)
            {
                int mid = (lo+hi)/2;
                if (g_wid >= d_nwarps[mid] && g_wid < d_nwarps[mid+1]) {
                    tag_id = mid;
                    break;
                } else if (g_wid < d_nwarps[mid]) {
                    hi = mid-1;
                } else {
                    lo = mid+1;
                }
            };
        }

        auto tag = d_tags[tag_id];
        int ncells = tag.box().numPts();
        int b_wid = g_wid - d_nwarps[tag_id]; // b_wid'th warp on this box
#ifdef AMREX_USE_DPCPP
        int lane = item.get_local_id(0) % Gpu::Device::warp_size;
#else
        int lane = threadIdx.x % Gpu::Device::warp_size;
#endif
        int icell = b_wid*Gpu::Device::warp_size + lane;

        const auto len = amrex::length(tag.box());
        const auto lo  = amrex::lbound(tag.box());
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;

#ifdef AMREX_USE_DPCPP
        f(item, icell, ncells, i, j, k, tag);
#else
        f(icell, ncells, i, j, k, tag);
#endif
    });

    Gpu::synchronize();
    The_Pinned_Arena()->free(h_buffer);
    The_Arena()->free(d_buffer);
}
}

template <class TagType, class F>
void
ParallelFor (Vector<TagType> const& tags, int ncomp, F && f)
{
    detail::ParallelFor_doit(tags,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_DPCPP
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
        {
            if (icell < ncells) {
                for (int n = 0; n < ncomp; ++n) {
                    f(i,j,k,n,tag);
                }
            }
        });
}

#endif

namespace detail {

#ifdef AMREX_USE_GPU

template <class T>
struct CellStore
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator() (T* d, T s) const noexcept
    {
        *d = s;
    }
};

template <class T>
struct CellAdd
{
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator() (T* d, T s) const noexcept
    {
        *d += s;
    }
};

template <class T>
struct CellAtomicAdd
{
    template<class U=T,
             std::enable_if_t<amrex::HasAtomicAdd<U>::value,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator() (U* d, U s) const noexcept
    {
        Gpu::Atomic::AddNoRet(d,s);
    }
};

template <class T, class F>
void
fab_to_fab (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
            F && f)
{
    detail::ParallelFor_doit(copy_tags,
        [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_DPCPP
            sycl::nd_item<1> const& /*item*/,
#endif
            int icell, int ncells, int i, int j, int k, Array4CopyTag<T> const tag) noexcept
        {
            if (icell < ncells) {
                for (int n = 0; n < ncomp; ++n) {
                    f(&(tag.dfab(i,j,k,n+dcomp)),
                      tag.sfab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z,n+scomp));
                }
            }
        });
}

template <class T, class F>
void
fab_to_fab (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
            F && f, Vector<Array4<int> > const& masks)
{
    typedef Array4MaskCopyTag<T> TagType;
    Vector<TagType> tags;
    const int N = copy_tags.size();
    tags.reserve(N);
    for (int i = 0; i < N; ++i) {
        tags.push_back(TagType{copy_tags[i].dfab, copy_tags[i].sfab, masks[i],
                               copy_tags[i].dbox, copy_tags[i].offset});
    }

    amrex::Abort("xxxxx TODO This function still has a bug.  Even if we fix the bug, it should still be avoided because it is slow due to the lack of atomic operations for this type.");

    detail::ParallelFor_doit(tags,
    [=] AMREX_GPU_DEVICE (
#ifdef AMREX_USE_DPCPP
        sycl::nd_item<1> const& item,
#endif
        int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
    {
#ifdef AMREX_USE_DPCPP
        int g_tid = item.get_global_id(0);
        int g_wid = g_tid / Gpu::Device::warp_size;

        int* m = (icell < ncells) ? tag.mask.ptr(i,j,k) : nullptr;
        int mypriority = g_wid+1;
        int to_try  = 1;
        while (true) {
            int msk = (m && to_try) ? Gpu::Atomic::CAS(m, 0, mypriority) : 0;
            if (sycl::ONEAPI::all_of(item.get_sub_group(), msk == 0)) {  // 0 means lock acquired
                break; // all threads have acquired.
            } else {
                if (sycl::ONEAPI::any_of(item.get_sub_group(), msk > mypriority)) {
                    if (m) *m = 0; // yield
                    item.mem_fence(); // xxxxx DPCPP todo: This is block level, but needs to be device level fence, which is currently a PR in intel/llvm
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

template <typename T, std::enable_if_t<amrex::IsStoreAtomic<T>::value,int> = 0>
void
fab_to_fab_atomic_cpy (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
                       Vector<Array4<int> > const&)
{
    fab_to_fab<T>(copy_tags, scomp, dcomp, ncomp, CellStore<T>());
}

template <typename T, std::enable_if_t<!amrex::IsStoreAtomic<T>::value,int> = 0>
void
fab_to_fab_atomic_cpy (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
                       Vector<Array4<int> > const& masks)
{
    fab_to_fab<T>(copy_tags, scomp, dcomp, ncomp, CellStore<T>(), masks);
}

template <typename T, std::enable_if_t<amrex::HasAtomicAdd<T>::value,int> = 0>
void
fab_to_fab_atomic_add (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
                       Vector<Array4<int> > const&)
{
    fab_to_fab<T>(copy_tags, scomp, dcomp, ncomp, CellAtomicAdd<T>());
}

template <typename T, std::enable_if_t<!amrex::HasAtomicAdd<T>::value,int> = 0>
void
fab_to_fab_atomic_add (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
                       Vector<Array4<int> > const& masks)
{
    fab_to_fab<T>(copy_tags, scomp, dcomp, ncomp, CellAdd<T>(), masks);
}

#endif /* AMREX_USE_GPU */

}

template <class FAB>
void
FabArray<FAB>::FB_local_copy_cpu (const FB& TheFB, int scomp, int ncomp)
{
    auto const& LocTags = *(TheFB.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) return;
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
void
FabArray<FAB>::FB_local_copy_gpu (const FB& TheFB, int scomp, int ncomp)
{
    auto const& LocTags = *(TheFB.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) return;
    bool is_thread_safe = TheFB.m_threadsafe_loc;

    typedef Array4CopyTag<value_type> TagType;
    Vector<TagType> loc_copy_tags;
    loc_copy_tags.reserve(N_locs);

    Vector<BaseFab<int> > maskfabs;
    Vector<Array4<int> > masks;
    if (!amrex::IsStoreAtomic<value_type>::value && !is_thread_safe)
    {
        maskfabs.resize(this->local_size());
        masks.reserve(N_locs);
    }

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

        if (maskfabs.size() > 0) {
            if (!maskfabs[li].isAllocated()) {
                maskfabs[li].resize(this->atLocalIdx(li).box());
            }
            masks.push_back(maskfabs[li].array());
        }
    }

    if (maskfabs.size() > 0) {
        Gpu::FuseSafeGuard fsg(maskfabs.size() >= Gpu::getFuseNumKernelsThreshold());
        for (Gpu::StreamIter sit(maskfabs.size()); sit.isValid(); ++sit) {
            BaseFab<int>& mskfab = maskfabs[sit()];
            const Array4<int>& msk = mskfab.array();
            const Box& bx = mskfab.box();
            amrex::ParallelFor(Gpu::KernelInfo{}.setFusible(true), bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                msk(i,j,k) = 0;
            });
        }
        Gpu::LaunchFusedKernels();
    }

    if (is_thread_safe) {
        detail::fab_to_fab<value_type>(loc_copy_tags, scomp, scomp, ncomp,
                                       detail::CellStore<value_type>());
    } else {
        detail::fab_to_fab_atomic_cpy<value_type>(loc_copy_tags, scomp, scomp, ncomp, masks);
    }
}

template <class FAB>
void
FabArray<FAB>::CMD_local_setVal_gpu (typename FabArray<FAB>::value_type x,
                                    const CommMetaData& thecmd, int scomp, int ncomp)
{
    auto const& LocTags = *(thecmd.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) return;
    bool is_thread_safe = thecmd.m_threadsafe_loc;

    typedef Array4BoxTag<value_type> TagType;
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

    typedef Array4BoxTag<value_type> TagType;
    Vector<TagType> rcv_setval_tags;

    for (auto it = RcvTags.begin(); it != RcvTags.end(); ++it) {
        for (auto const& tag: it->second) {
            rcv_setval_tags.push_back({this->array(tag.dstIndex), tag.dbox});
        }
    }

    if (rcv_setval_tags.empty()) return;

    AMREX_ALWAYS_ASSERT(amrex::IsStoreAtomic<value_type>::value || is_thread_safe);

    amrex::ParallelFor(rcv_setval_tags, ncomp,
    [x,scomp] AMREX_GPU_DEVICE (int i, int j, int k, int n, TagType const& tag) noexcept
    {
        tag.dfab(i,j,k,n+scomp) = x;
    });
}

#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10))
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

#endif /* CUDA >= 10 */

#endif /* AMREX_USE_GPU */

#ifdef AMREX_USE_MPI

#ifdef AMREX_USE_GPU

#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )

template <class FAB>
void
FabArray<FAB>::FB_pack_send_buffer_cuda_graph (const FB& TheFB, int scomp, int ncomp,
                                               Vector<char*>& send_data,
                                               Vector<std::size_t> const& send_size,
                                               Vector<typename FabArray<FAB>::CopyComTagsContainer const*> const& send_cctc)
{
    const int N_snds = send_data.size();
    if (N_snds == 0) return;

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
    if (N_rcvs == 0) return;

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

#endif /* CUDA >= 10 */

template <class FAB>
void
FabArray<FAB>::pack_send_buffer_gpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                     Vector<char*> const& send_data,
                                     Vector<std::size_t> const& send_size,
                                     Vector<CopyComTagsContainer const*> const& send_cctc)
{
    amrex::ignore_unused(send_size);

    const int N_snds = send_data.size();
    if (N_snds == 0) return;

    char* pbuffer = send_data[0];
    std::size_t szbuffer = 0;
#if 0
    // For linear solver test on summit, this is slower than writing to
    // pinned memory directly on device.
    if (! ParallelDescriptor::UseGpuAwareMpi()) {
        // Memory in send_data is pinned.
        szbuffer = (send_data[N_snds-1]-send_data[0]) + send_size[N_snds-1];
        pbuffer = (char*)The_Arena()->alloc(szbuffer);
    }
#endif

    typedef Array4CopyTag<value_type> TagType;
    Vector<TagType> snd_copy_tags;
    for (int j = 0; j < N_snds; ++j)
    {
        if (send_size[j] > 0)
        {
            std::size_t offset = send_data[j]-send_data[0];
            char* dptr = pbuffer + offset;
            auto const& cctc = *send_cctc[j];
            for (auto const& tag : cctc)
            {
                snd_copy_tags.emplace_back(TagType{
                    amrex::makeArray4((value_type*)(dptr), tag.sbox, ncomp),
                    src.array(tag.srcIndex),
                    tag.sbox,
                    Dim3{0,0,0}
                });
                dptr += (tag.sbox.numPts() * ncomp * sizeof(value_type));
            }
            BL_ASSERT(dptr <= pbuffer + offset + send_size[j]);
        }
    }

    detail::fab_to_fab<value_type>(snd_copy_tags, scomp, 0, ncomp,
                                   detail::CellStore<value_type>());

    // There is Gpu::synchronize in fab_to_fab.

    if (pbuffer != send_data[0]) {
        Gpu::copyAsync(Gpu::deviceToHost,pbuffer,pbuffer+szbuffer,send_data[0]);
        Gpu::synchronize();
        The_Arena()->free(pbuffer);
    }
}

template <class FAB>
void
FabArray<FAB>::unpack_recv_buffer_gpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                       Vector<char*> const& recv_data,
                                       Vector<std::size_t> const& recv_size,
                                       Vector<CopyComTagsContainer const*> const& recv_cctc,
                                       CpOp op, bool is_thread_safe)
{
    amrex::ignore_unused(recv_size);

    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    char* pbuffer = recv_data[0];
#if 0
    std::size_t szbuffer = 0;
    // For linear solver test on summit, this is slower than writing to
    // pinned memory directly on device.
    if (! ParallelDescriptor::UseGpuAwareMpi()) {
        // Memory in recv_data is pinned.
        szbuffer = (recv_data[N_rcvs-1]-recv_data[0]) + recv_size[N_rcvs-1];
        pbuffer = (char*)The_Arena()->alloc(szbuffer);
        Gpu::copyAsync(Gpu::hostToDevice,recv_data[0],recv_data[0]+szbuffer,pbuffer);
        Gpu::synchronize();
    }
#endif

    typedef Array4CopyTag<value_type> TagType;
    Vector<TagType> recv_copy_tags;
    recv_copy_tags.reserve(N_rcvs);

    Vector<BaseFab<int> > maskfabs;
    Vector<Array4<int> > masks;
    if (!is_thread_safe)
    {
        if ((op == FabArrayBase::COPY && !amrex::IsStoreAtomic<value_type>::value) ||
            (op == FabArrayBase::ADD  && !amrex::HasAtomicAdd <value_type>::value))
        {
            maskfabs.resize(dst.local_size());
        }
    }

    for (int k = 0; k < N_rcvs; ++k)
    {
        if (recv_size[k] > 0)
        {
            std::size_t offset = recv_data[k]-recv_data[0];
            const char* dptr = pbuffer + offset;
            auto const& cctc = *recv_cctc[k];
            for (auto const& tag : cctc)
            {
                const int li = dst.localindex(tag.dstIndex);
                recv_copy_tags.emplace_back(TagType{
                    dst.atLocalIdx(li).array(),
                    amrex::makeArray4((value_type const*)(dptr), tag.dbox, ncomp),
                    tag.dbox,
                    Dim3{0,0,0}
                });
                dptr += tag.dbox.numPts() * ncomp * sizeof(value_type);

                if (maskfabs.size() > 0) {
                    if (!maskfabs[li].isAllocated()) {
                        maskfabs[li].resize(dst.atLocalIdx(li).box());
                    }
                    masks.push_back(maskfabs[li].array());
                }
            }
            BL_ASSERT(dptr <= pbuffer + offset + recv_size[k]);
        }
    }

    if (maskfabs.size() > 0) {
        Gpu::FuseSafeGuard fsg(maskfabs.size() >= Gpu::getFuseNumKernelsThreshold());
        for (Gpu::StreamIter sit(maskfabs.size()); sit.isValid(); ++sit) {
            BaseFab<int>& mskfab = maskfabs[sit()];
            const Array4<int>& msk = mskfab.array();
            const Box& bx = mskfab.box();
            amrex::ParallelFor(Gpu::KernelInfo().setFusible(true), bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                msk(i,j,k) = 0;
            });
        }
        Gpu::LaunchFusedKernels();
    }

    if (op == FabArrayBase::COPY)
    {
        if (is_thread_safe) {
            detail::fab_to_fab<value_type>(recv_copy_tags, 0, dcomp, ncomp,
                                           detail::CellStore<value_type>());
        } else {
            detail::fab_to_fab_atomic_cpy<value_type>(recv_copy_tags, 0, dcomp, ncomp, masks);
        }
    }
    else
    {
        if (is_thread_safe) {
            detail::fab_to_fab<value_type>(recv_copy_tags, 0, dcomp, ncomp,
                                           detail::CellAdd<value_type>());
        } else {
            detail::fab_to_fab_atomic_add<value_type>(recv_copy_tags, 0, dcomp, ncomp, masks);
        }
    }

    // There is Gpu::synchronize in fab_to_fab.

    if (pbuffer != recv_data[0]) {
        The_Arena()->free(pbuffer);
    }
}

#endif /* AMREX_USE_GPU */

template <class FAB>
void
FabArray<FAB>::pack_send_buffer_cpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                     Vector<char*> const& send_data,
                                     Vector<std::size_t> const& send_size,
                                     Vector<CopyComTagsContainer const*> const& send_cctc)
{
    amrex::ignore_unused(send_size);

    const int N_snds = send_data.size();
    if (N_snds == 0) return;

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
                auto pfab = amrex::makeArray4((value_type*)(dptr),bx,ncomp);
                amrex::LoopConcurrentOnCpu( bx, ncomp,
                [=] (int ii, int jj, int kk, int n) noexcept
                {
                    pfab(ii,jj,kk,n) = sfab(ii,jj,kk,n+scomp);
                });
                dptr += (bx.numPts() * ncomp * sizeof(value_type));
            }
            BL_ASSERT(dptr <= send_data[j] + send_size[j]);
        }
    }
}

template <class FAB>
void
FabArray<FAB>::unpack_recv_buffer_cpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                       Vector<char*> const& recv_data,
                                       Vector<std::size_t> const& recv_size,
                                       Vector<CopyComTagsContainer const*> const& recv_cctc,
                                       CpOp op, bool is_thread_safe)
{
    amrex::ignore_unused(recv_size);

    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

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
                        dfab.template copyFromMem<RunOn::Host>(bx, dcomp, ncomp, dptr);
                    }
                    else
                    {
                        dfab.template addFromMem<RunOn::Host>(tag.dbox, dcomp, ncomp, dptr);
                    }
                    dptr += bx.numPts() * ncomp * sizeof(value_type);
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
                    dptr += tag.dbox.numPts() * ncomp * sizeof(value_type);
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
                auto pfab = amrex::makeArray4((value_type*)(tag.p), tag.dbox, ncomp);
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
