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
};

struct VoidCopyTag {
    char const* p;
    Box dbox;
};

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
             amrex::EnableIf_t<amrex::HasAtomicAdd<U>::value,int> = 0>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
    operator() (U* d, U s) const noexcept
    {
        atomicAdd(d,s);
    }
};

template <class T, class F>
void
fab_to_fab (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
            F && f)
{
    typedef Array4CopyTag<T> TagType;

    const int N_locs = copy_tags.size();
    if (N_locs == 0) return;

    int ntotwarps = 0;
    Vector<int> nwarps;
    nwarps.reserve(N_locs+1);
    for (int i = 0; i < N_locs; ++i)
    {
        auto& tag = copy_tags[i];
        nwarps.push_back(ntotwarps);
        ntotwarps += static_cast<int>(tag.dbox.numPts()+Gpu::Device::warp_size-1)/Gpu::Device::warp_size;
    }
    nwarps.push_back(ntotwarps);

    std::size_t nbytes = N_locs*sizeof(TagType);
    auto d_tags = static_cast<TagType*>(The_Device_Arena()->alloc(nbytes));
    Gpu::htod_memcpy(d_tags, copy_tags.data(), nbytes);

    nbytes = (N_locs+1)*sizeof(int);
    auto d_nwarps = static_cast<int*>(The_Device_Arena()->alloc(nbytes));
    Gpu::htod_memcpy(d_nwarps, nwarps.data(), nbytes);

    constexpr int nthreads = 128;
    constexpr int nwarps_per_block = nthreads/Gpu::Device::warp_size;
    int nblocks = (ntotwarps + nwarps_per_block-1) / nwarps_per_block;
    amrex::launch_global<<<nblocks,nthreads>>>(
    [=] AMREX_GPU_DEVICE () noexcept
    {
        int g_tid = blockDim.x*blockIdx.x + threadIdx.x;
        int g_wid = g_tid / Gpu::Device::warp_size;
        if (g_wid >= ntotwarps) return;

        int tag_id;
        {
            int lo = 0;
            int hi = N_locs;
            while (lo <= hi)
            {
                int mid = (lo+hi)/2;
                if (g_wid >= d_nwarps[mid] and g_wid < d_nwarps[mid+1]) {
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
        int ncells = tag.dbox.numPts();
        int b_wid = g_wid - d_nwarps[tag_id]; // b_wid'th warp on this box
        int lane = threadIdx.x % Gpu::Device::warp_size;
        int icell = b_wid*Gpu::Device::warp_size + lane;
        if (icell < ncells) {
            const auto len = amrex::length(tag.dbox);
            const auto lo  = amrex::lbound(tag.dbox);
            int k =  icell /   (len.x*len.y);
            int j = (icell - k*(len.x*len.y)) /   len.x;
            int i = (icell - k*(len.x*len.y)) - j*len.x;
            i += lo.x;
            j += lo.y;
            k += lo.z;
            for (int n = 0; n < ncomp; ++n) {
                f(&(tag.dfab(i,j,k,n+dcomp)),
                  tag.sfab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z,n+scomp));
            }
        }
    });

    Gpu::synchronize();
    The_Device_Arena()->free(d_nwarps);
    The_Device_Arena()->free(d_tags);
}

template <class T, class F>
void
fab_to_fab (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
            F && f, Vector<Array4<int> > const& masks)
{
    typedef Array4CopyTag<T> TagType;

    const int N_locs = copy_tags.size();
    if (N_locs == 0) return;

    int ntotwarps = 0;
    Vector<int> nwarps;
    nwarps.reserve(N_locs+1);
    for (int i = 0; i < N_locs; ++i)
    {
        auto& tag = copy_tags[i];
        nwarps.push_back(ntotwarps);
        ntotwarps += static_cast<int>(tag.dbox.numPts()+Gpu::Device::warp_size-1)/Gpu::Device::warp_size;
    }
    nwarps.push_back(ntotwarps);

    std::size_t nbytes = N_locs*sizeof(TagType);
    auto d_tags = static_cast<TagType*>(The_Device_Arena()->alloc(nbytes));
    Gpu::htod_memcpy(d_tags, copy_tags.data(), nbytes);

    nbytes = (N_locs+1)*sizeof(int);
    auto d_nwarps = static_cast<int*>(The_Device_Arena()->alloc(nbytes));
    Gpu::htod_memcpy(d_nwarps, nwarps.data(), nbytes);

    nbytes = masks.size()*sizeof(Array4<int>);
    auto d_masks = static_cast<Array4<int>*>(The_Device_Arena()->alloc(nbytes));
    Gpu::htod_memcpy(d_masks, masks.data(), nbytes);

    constexpr int nthreads = 128;
    constexpr int nwarps_per_block = nthreads/Gpu::Device::warp_size;
    int nblocks = (ntotwarps + nwarps_per_block-1) / nwarps_per_block;
    amrex::launch_global<<<nblocks,nthreads>>>(
    [=] AMREX_GPU_DEVICE () noexcept
    {
        int g_tid = blockDim.x*blockIdx.x + threadIdx.x;
        int g_wid = g_tid / Gpu::Device::warp_size;
        if (g_wid >= ntotwarps) return;

        int tag_id;
        {
            int lo = 0;
            int hi = N_locs;
            while (lo <= hi)
            {
                int mid = (lo+hi)/2;
                if (g_wid >= d_nwarps[mid] and g_wid < d_nwarps[mid+1]) {
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
        int ncells = tag.dbox.numPts();
        int b_wid = g_wid - d_nwarps[tag_id]; // b_wid'th warp on this box
        int lane = threadIdx.x % Gpu::Device::warp_size;
        int icell = b_wid*Gpu::Device::warp_size + lane;

        const auto len = amrex::length(tag.dbox);
        const auto lo  = amrex::lbound(tag.dbox);
        int k =  icell /   (len.x*len.y);
        int j = (icell - k*(len.x*len.y)) /   len.x;
        int i = (icell - k*(len.x*len.y)) - j*len.x;
        i += lo.x;
        j += lo.y;
        k += lo.z;

        int* m = (icell < ncells) ? d_masks[tag_id].ptr(i,j,k) : nullptr;
        int mypriority = g_wid+1;
        int to_try  = 1;
        while (true) {
            int msk = (m && to_try) ? atomicCAS(m, 0, mypriority) : 0;
            if (__all_sync(0xffffffff, msk == 0)) {  // 0 means lock acquired
                break; // all threads have acquired.
            } else {
                if (__any_sync(0xffffffff, msk > mypriority)) {
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
    });

    Gpu::synchronize();
    The_Device_Arena()->free(d_masks);
    The_Device_Arena()->free(d_nwarps);
    The_Device_Arena()->free(d_tags);
}

template <typename T, amrex::EnableIf_t<amrex::IsStoreAtomic<T>::value,int> = 0>
void
fab_to_fab_atomic_cpy (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
                       Vector<Array4<int> > const&)
{
    fab_to_fab<T>(copy_tags, scomp, dcomp, ncomp, CellStore<T>());
}

template <typename T, amrex::EnableIf_t<!amrex::IsStoreAtomic<T>::value,int> = 0>
void
fab_to_fab_atomic_cpy (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
                       Vector<Array4<int> > const& masks)
{
    fab_to_fab<T>(copy_tags, scomp, dcomp, ncomp, CellStore<T>(), masks);
}

template <typename T, amrex::EnableIf_t<amrex::HasAtomicAdd<T>::value,int> = 0>
void
fab_to_fab_atomic_add (Vector<Array4CopyTag<T> > const& copy_tags, int scomp, int dcomp, int ncomp,
                       Vector<Array4<int> > const&)
{
    fab_to_fab<T>(copy_tags, scomp, dcomp, ncomp, CellAtomicAdd<T>());
}

template <typename T, amrex::EnableIf_t<!amrex::HasAtomicAdd<T>::value,int> = 0>
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
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = LocTags[i];

            BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            const FAB* sfab = &(get(tag.srcIndex));
                  FAB* dfab = &(get(tag.dstIndex));
            dfab->copy(*sfab, tag.sbox, scomp, tag.dbox, scomp, ncomp);
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
#ifdef _OPENMP
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
                const Box sbox = tag.dbox + tag.offset;
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
    if (!amrex::IsStoreAtomic<value_type>::value and !is_thread_safe)
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
        for (Gpu::StreamIter sit(maskfabs.size()); sit.isValid(); ++sit) {
            BaseFab<int>& mskfab = maskfabs[sit()];
            const Array4<int>& msk = mskfab.array();
            const Box& bx = mskfab.box();
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                msk(i,j,k) = 0;
            });
        }
    }

    if (is_thread_safe) {
        detail::fab_to_fab<value_type>(loc_copy_tags, scomp, scomp, ncomp,
                                       detail::CellStore<value_type>());
    } else {
        detail::fab_to_fab_atomic_cpy<value_type>(loc_copy_tags, scomp, scomp, ncomp, masks);
    }
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
                                               Vector<int> const& send_size,
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
            if (send_data[send] != nullptr) {
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
            char* dptr = send_data[j];
            if (dptr != nullptr)
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
        char* dptr = send_data[j];
        if (dptr != nullptr)
        {
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
            BL_ASSERT(dptr == send_data[j] + send_size[j]); 
        }
    }
    
    // Launch Graph synched, so copyToBuffer is complete prior to posting sends.
    TheFB.m_copyToBuffer.executeGraph();
}

template <class FAB>
void
FabArray<FAB>::FB_unpack_recv_buffer_cuda_graph (const FB& TheFB, int dcomp, int ncomp,
                                                 Vector<char*> const& recv_data,
                                                 Vector<int> const& recv_size,
                                                 Vector<CopyComTagsContainer const*> const& recv_cctc,
                                                 bool is_thread_safe)
{
    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    int launches = 0;
    LayoutData<Vector<VoidCopyTag> > recv_copy_tags(boxArray(),DistributionMap());
    for (int k = 0; k < N_rcvs; ++k)
    {
        const char* dptr = recv_data[k];
        if (dptr != nullptr)
        {
            auto const& cctc = *recv_cctc[k];
            for (auto const& tag : cctc)
            {
                recv_copy_tags[tag.dstIndex].push_back({dptr,tag.dbox});
                dptr += tag.dbox.numPts() * ncomp * sizeof(value_type);
                launches++;
            }
            BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
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
                                     Vector<char*>& send_data,
                                     Vector<int> const& send_size,
                                     Vector<CopyComTagsContainer const*> const& send_cctc)
{
    const int N_snds = send_data.size();
    if (N_snds == 0) return;

    typedef Array4CopyTag<value_type> TagType;
    Vector<TagType> snd_copy_tags;
    for (int j = 0; j < N_snds; ++j)
    {
        char* dptr = send_data[j];
        if (dptr != nullptr)
        {
            auto const& cctc = *send_cctc[j];
            for (auto const& tag : cctc)
            {
                snd_copy_tags.push_back
                    ({amrex::makeArray4((value_type*)(dptr),tag.sbox,ncomp),
                      src.array(tag.srcIndex), tag.sbox, Dim3{0,0,0}});
                dptr += (tag.sbox.numPts() * ncomp * sizeof(value_type));
            }
            BL_ASSERT(dptr == send_data[j] + send_size[j]); 
        }
    }

    detail::fab_to_fab<value_type>(snd_copy_tags, scomp, 0, ncomp,
                                   detail::CellStore<value_type>());
}

template <class FAB>
void
FabArray<FAB>::unpack_recv_buffer_gpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                       Vector<char*> const& recv_data,
                                       Vector<int> const& recv_size,
                                       Vector<CopyComTagsContainer const*> const& recv_cctc,
                                       CpOp op, bool is_thread_safe)
{
    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    typedef Array4CopyTag<value_type> TagType;
    Vector<TagType> recv_copy_tags;

    Vector<BaseFab<int> > maskfabs;
    Vector<Array4<int> > masks;
    if (!is_thread_safe)
    {
        if ((op == FabArrayBase::COPY and !amrex::IsStoreAtomic<value_type>::value) or
            (op == FabArrayBase::ADD  and !amrex::HasAtomicAdd <value_type>::value))
        {
            maskfabs.resize(dst.local_size());
        }
    }

    for (int k = 0; k < N_rcvs; ++k)
    {
        const char* dptr = recv_data[k];
        if (dptr != nullptr)
        {
            auto const& cctc = *recv_cctc[k];
            for (auto const& tag : cctc)
            {
                const int li = dst.localindex(tag.dstIndex);
                recv_copy_tags.push_back
                    ({dst.atLocalIdx(li).array(),
                      amrex::makeArray4((value_type*)(dptr),tag.dbox,ncomp),
                      tag.dbox, Dim3{0,0,0}});
                dptr += tag.dbox.numPts() * ncomp * sizeof(value_type);

                if (maskfabs.size() > 0) {
                    if (!maskfabs[li].isAllocated()) {
                        maskfabs[li].resize(dst.atLocalIdx(li).box());
                    }
                    masks.push_back(maskfabs[li].array());
                }
            }
            BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
        }
    }

    if (maskfabs.size() > 0) {
        for (Gpu::StreamIter sit(maskfabs.size()); sit.isValid(); ++sit) {
            BaseFab<int>& mskfab = maskfabs[sit()];
            const Array4<int>& msk = mskfab.array();
            const Box& bx = mskfab.box();
            amrex::ParallelFor(bx,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                msk(i,j,k) = 0;
            });
        }
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
}

#endif /* AMREX_USE_GPU */

template <class FAB>
void
FabArray<FAB>::pack_send_buffer_cpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                     Vector<char*>& send_data,
                                     Vector<int> const& send_size,
                                     Vector<CopyComTagsContainer const*> const& send_cctc)
{
    const int N_snds = send_data.size();
    if (N_snds == 0) return;

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int j = 0; j < N_snds; ++j)
    {
        char* dptr = send_data[j];
        if (dptr != nullptr)
        {
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
            BL_ASSERT(dptr == send_data[j] + send_size[j]); 
        }
    }
}

template <class FAB>
void
FabArray<FAB>::unpack_recv_buffer_cpu (FabArray<FAB>& dst, int dcomp, int ncomp,
                                       Vector<char*> const& recv_data,
                                       Vector<int> const& recv_size,
                                       Vector<CopyComTagsContainer const*> const& recv_cctc,
                                       CpOp op, bool is_thread_safe)
{
    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    if (is_thread_safe)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int k = 0; k < N_rcvs; ++k)
        {
            const char* dptr = recv_data[k];
            if (dptr != nullptr)
            {
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    const Box& bx  = tag.dbox;
                    FAB& dfab = dst[tag.dstIndex];
                    if (op == FabArrayBase::COPY)
                    {
                        dfab.copyFromMem(bx, dcomp, ncomp, dptr);
                    }
                    else
                    {
                        dfab.addFromMem(tag.dbox, dcomp, ncomp, dptr);
                    }
                    dptr += bx.numPts() * ncomp * sizeof(value_type);
                }
                BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
            }
        }
    }
    else
    {
        LayoutData<Vector<VoidCopyTag> > recv_copy_tags;
        recv_copy_tags.define(dst.boxArray(),dst.DistributionMap());
        for (int k = 0; k < N_rcvs; ++k)
        {
            const char* dptr = recv_data[k];
            if (dptr != nullptr)
            {
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    recv_copy_tags[tag.dstIndex].push_back({dptr,tag.dbox});
                    dptr += tag.dbox.numPts() * ncomp * sizeof(value_type);
                }
                BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
            }
        }

#ifdef _OPENMP
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
