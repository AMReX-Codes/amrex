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

template <class TagType, class F>
void
fab_to_fab_store (Vector<TagType> const& tags, int scomp, int dcomp, int ncomp, F&&f)
{
    amrex::ParallelFor(tags,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, TagType const& tag) noexcept
    {
        int m = Gpu::Atomic::Add(&(tag.mask(i,j,k)), 1);
        if (m == 0) {
            for (int n = 0; n < ncomp; ++n) {
                f(&(tag.dfab(i,j,k,n+dcomp)),
                  tag.sfab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z,n+scomp));
            }
        }
    });
}

template <class TagType, class F>
void
fab_to_fab_other (Vector<TagType> const& tags, int scomp, int dcomp, int ncomp, F&&f)
{
    amrex::ParallelFor(tags,
        [=] AMREX_GPU_DEVICE (int i, int j, int k, TagType const& tag) noexcept
    {
        int* m = &(tag.mask(i,j,k));
        bool my_turn = false;
        do {
#if defined(AMREX_USE_SYCL)
            my_turn = (Gpu::Atomic::Exch(m, 1) == 0);
#else
            my_turn = (Gpu::Atomic::CAS(m, 0, 1) == 0);
#endif
            if (my_turn) {
#if defined(AMREX_USE_SYCL)
                sycl::atomic_fence(sycl::memory_order::acq_rel, sycl::memory_scope::device);
#else
                __threadfence();
#endif
                for (int n = 0; n < ncomp; ++n) {
                    f(&(tag.dfab(i,j,k,n+dcomp)),
                      tag.sfab(i+tag.offset.x,j+tag.offset.y,k+tag.offset.z,n+scomp));
                }
#if defined(AMREX_USE_SYCL)
                sycl::atomic_fence(sycl::memory_order::acq_rel, sycl::memory_scope::device);
#else
                __threadfence(); // It appears that we need this fence too if a GPU is shared.
#endif
                Gpu::Atomic::Exch(m, 0);
            }
            else {
#if defined(AMREX_USE_CUDA)

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 700)
#if defined(_WIN32)
                volatile int tmp; // prevent optimization
                for (int c = 0; c < 2; ++c) {
                    ++tmp;
                }
#else
                for (int c = 0; c < 2; ++c) {
                    __asm__ volatile(""); // prevent optimization
                }
#endif
#else
                __nanosleep(1);
#endif

#elif defined(AMREX_USE_HIP)

                __builtin_amdgcn_s_sleep(1);

#elif defined(AMREX_USE_SYCL)

                for (int c = 0; c < 2; ++c) {
                    __asm__ volatile(""); // prevent optimization
                }

#endif
            }
        } while (!my_turn);
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

    if constexpr (std::is_same_v<F, CellStore<T0,T1>>)
    {
        fab_to_fab_store(tags, scomp, dcomp, ncomp, std::forward<F>(f));
    }
    else
    {
        fab_to_fab_other(tags, scomp, dcomp, ncomp, std::forward<F>(f));
    }
    // Note that Tag ParalleFor has a stream sync call in the end.
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

template <typename T0, typename T1, class F>
void deterministic_fab_to_fab (Vector<Array4CopyTag<T0,T1>> const& a_tags, int scomp,
                               int dcomp, int ncomp, F const& f)
{
    if (a_tags.empty()) { return; }

    using TagType = Array4CopyTag<T0,T1>;

    struct TiledTag {
        int tag_index;
        std::pair<int,Box> dindex_tilebox;
        bool operator< (TiledTag const& rhs) const noexcept {
            return this->dindex_tilebox < rhs.dindex_tilebox;
        }
        bool operator!= (TiledTag const& rhs) const noexcept {
            return this->dindex_tilebox != rhs.dindex_tilebox;
        }
    };
    Vector<TiledTag> tiled_tags;

    auto const ixtype = a_tags[0].dbox.ixType();

    constexpr int tile_size = 64;
    for (int itag = 0; itag < a_tags.size(); ++itag) {
        auto const& tag = a_tags[itag];
        auto const& dlo = tag.dbox.smallEnd();
        auto const& dhi = tag.dbox.bigEnd();
        IntVect tlo(AMREX_D_DECL(amrex::coarsen<tile_size>(dlo[0]),
                                 amrex::coarsen<tile_size>(dlo[1]),
                                 amrex::coarsen<tile_size>(dlo[2])));
        IntVect thi(AMREX_D_DECL(amrex::coarsen<tile_size>(dhi[0]),
                                 amrex::coarsen<tile_size>(dhi[1]),
                                 amrex::coarsen<tile_size>(dhi[2])));
#if (AMREX_SPACEDIM == 3)
        for (int kt = tlo[2]; kt <= thi[2]; ++kt)
#endif
        {
#if (AMREX_SPACEDIM >= 2)
            for (int jt = tlo[1]; jt <= thi[1]; ++jt)
#endif
            {
                for (int it = tlo[0]; it <= thi[0]; ++it)
                {
                    IntVect lo(AMREX_D_DECL(it*tile_size,
                                            jt*tile_size,
                                            kt*tile_size));
                    tiled_tags.push_back(TiledTag{itag, std::make_pair
                            (tag.dindex, Box(lo, lo+(tile_size-1), ixtype))});
                }
            }
        }
    }

    std::sort(tiled_tags.begin(), tiled_tags.end());

    Gpu::HostVector<unsigned int> h_ntags;
    Gpu::HostVector<TagType> h_tags;
    h_tags.reserve(tiled_tags.size());

    for (unsigned int itag = 0; itag < tiled_tags.size(); ++itag) {
        if (itag == 0) {
            h_ntags.push_back(0);
        } else if (tiled_tags[itag-1] != tiled_tags[itag]) {
            h_ntags.push_back(itag);
        }
        auto const& ttag = tiled_tags[itag];
        auto const& btag = a_tags[ttag.tag_index];
        h_tags.push_back({btag.dfab, btag.dindex, btag.sfab,
                          btag.dbox & ttag.dindex_tilebox.second, btag.offset});
    }
    h_ntags.push_back((unsigned int)tiled_tags.size());

    Gpu::DeviceVector<TagType> d_tags(h_tags.size());
    Gpu::DeviceVector<unsigned int> d_ntags(h_ntags.size());
    Gpu::copyAsync(Gpu::hostToDevice,h_tags.begin(),h_tags.end(),d_tags.begin());
    Gpu::copyAsync(Gpu::hostToDevice,h_ntags.begin(),h_ntags.end(),d_ntags.begin());
    auto const* ptag = d_tags.data();
    auto const* pntags = d_ntags.data();
    auto const nblocks = int(h_ntags.size()-1);
    constexpr auto nthreads = 256;
    amrex::launch<nthreads>(nblocks, Gpu::gpuStream(),
#ifdef AMREX_USE_SYCL
        [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
        [[sycl::reqd_work_group_size(nthreads)]]
#else
        [=] AMREX_GPU_DEVICE () noexcept
#endif
    {
#ifdef AMREX_USE_SYCL
        Dim1 blockIdx{item.get_group_linear_id()};
        Dim1 threadIdx{item.get_local_linear_id()};
#endif

        for (unsigned int itag = pntags[blockIdx.x]; itag < pntags[blockIdx.x+1]; ++itag) {
            auto const tag = ptag[itag];
            auto ncells = int(tag.dbox.numPts());
            const auto len = amrex::length(tag.dbox);
            const auto lo  = amrex::lbound(tag.dbox);
            for (int icell = int(threadIdx.x); icell < ncells; icell += nthreads) {
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                for (int n = 0; n < ncomp; ++n) {
                    f(tag.dfab.ptr(i,j,k,n+dcomp),
                      tag.sfab(i + tag.offset.x,
                               j + tag.offset.y,
                               k + tag.offset.z, n+scomp));
                }
            }

            if (itag+1 < pntags[blockIdx.x+1]) {
#ifdef AMREX_USE_SYCL
                sycl::group_barrier(item.get_group());
#else
                __syncthreads();
#endif
            }
        }
    });
    Gpu::streamSynchronize();
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

template <class FAB>
void
FabArray<FAB>::FB_local_add_cpu (const FB& TheFB, int scomp, int ncomp)
{
    auto const& LocTags = *(TheFB.m_LocTags);
    auto N_locs = static_cast<int>(LocTags.size());
    if (N_locs == 0) { return; }

    LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
    // We must make a temporary copy of the data first to avoid race condition.
    std::vector<FAB> src_fabs(N_locs);
    for (int itag = 0; itag < N_locs; ++itag) {
        const CopyComTag& tag = LocTags[itag];
        src_fabs[itag].resize(tag.sbox,ncomp);
        loc_copy_tags[tag.dstIndex].push_back
            (FabCopyTag<FAB>{&(src_fabs[itag]),
                 tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int itag = 0; itag < N_locs; ++itag) {
        const CopyComTag& tag = LocTags[itag];
        src_fabs[itag].template copy<RunOn::Host>(this->operator[](tag.srcIndex), scomp, 0, ncomp);
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(*this); mfi.isValid(); ++mfi)
    {
        const auto& tags = loc_copy_tags[mfi];
        const auto& dfab = this->array(mfi);
        for (auto const & tag : tags)
        {
            auto const sfab = tag.sfab->array();
            const auto offset = tag.offset.dim3();
            amrex::LoopConcurrentOnCpu(tag.dbox, ncomp,
            [&] (int i, int j, int k, int n) noexcept
            {
                dfab(i,j,k,n+scomp) += sfab(i+offset.x,j+offset.y,k+offset.z,n);
            });
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
    if (auto it = m_fb_local_copy_handler.find(TheFB.m_id);
        it != m_fb_local_copy_handler.end())
    {
        tv = it->second.get();
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
                ({this->atLocalIdx(li).array(), tag.dstIndex,
                  this->fabPtr(tag.srcIndex)->const_array(),
                  tag.dbox,
                  (tag.sbox.smallEnd()-tag.dbox.smallEnd()).dim3()});
        }

        auto utv = std::make_unique<TagVector<TagType>>(loc_copy_tags);
        tv = utv.get();
        m_fb_local_copy_handler[TheFB.m_id] = std::move(utv);
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

    using TagType = Array4CopyTag<value_type>;

    if (is_thread_safe || amrex::IsStoreAtomic<value_type>::value)
    {
        auto* tv = FB_get_local_copy_tag_vector(TheFB);

        detail::ParallelFor_doit(*tv,
            [=] AMREX_GPU_DEVICE (
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
    else
    {
        Vector<TagType> loc_copy_tags;
        loc_copy_tags.reserve(N_locs);

        Vector<BaseFab<int>> maskfabs(this->local_size());
        Vector<Array4Tag<int> > masks_unique;
        masks_unique.reserve(this->local_size());
        Vector<Array4Tag<int> > masks;
        masks.reserve(N_locs);

        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = LocTags[i];

            BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            int li = this->localindex(tag.dstIndex);
            loc_copy_tags.push_back
                ({this->atLocalIdx(li).array(), tag.dstIndex,
                  this->fabPtr(tag.srcIndex)->const_array(),
                  tag.dbox,
                  (tag.sbox.smallEnd()-tag.dbox.smallEnd()).dim3()});

            if (!maskfabs[li].isAllocated()) {
                maskfabs[li].resize(this->atLocalIdx(li).box());
                masks_unique.emplace_back(Array4Tag<int>{maskfabs[li].array()});
            }
            masks.emplace_back(Array4Tag<int>{maskfabs[li].array()});
        }

        amrex::ParallelFor(masks_unique,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Array4Tag<int> const& msk) noexcept
            {
                msk.dfab(i,j,k) = 0;
            });

        detail::fab_to_fab_atomic_cpy<value_type, value_type>(
            loc_copy_tags, scomp, scomp, ncomp, masks);
    }
}

template <class FAB>
void
FabArray<FAB>::FB_local_add_gpu (const FB& TheFB, int scomp, int ncomp, bool deterministic)
{
    auto const& LocTags = *(TheFB.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) { return; }

    using TagType = Array4CopyTag<value_type>;

    Vector<TagType> loc_copy_tags_1, loc_copy_tags_2;
    loc_copy_tags_1.reserve(N_locs);
    loc_copy_tags_2.reserve(N_locs);

    Vector<FAB> src_fabs(N_locs);
    for (int itag = 0; itag < N_locs; ++itag) {
        const CopyComTag& tag = LocTags[itag];
        src_fabs[itag].resize(tag.sbox,ncomp);
        loc_copy_tags_1.push_back(
            TagType{src_fabs[itag].array(), -1,
                    this->const_array(tag.srcIndex,scomp), tag.sbox,
                    Dim3{0,0,0}});
        loc_copy_tags_2.push_back(
            TagType{this->array(tag.dstIndex,scomp), tag.dstIndex,
                    src_fabs[itag].const_array(), tag.dbox,
                    (tag.sbox.smallEnd()-tag.dbox.smallEnd()).dim3()});
    }

    // Note that we have shifted the starting component to zero in the code above.

    // TODO: We could try to cache the tags like in FillBoundary.

    detail::fab_to_fab(loc_copy_tags_1, 0, 0, ncomp,
                       detail::CellStore<value_type, value_type>{});
    if (deterministic || ! amrex::HasAtomicAdd<value_type>::value) {
        detail::deterministic_fab_to_fab(loc_copy_tags_2, 0, 0, ncomp,
                                         detail::CellAdd<value_type,value_type>{});
    } else {
        if constexpr(amrex::HasAtomicAdd<value_type>::value) {
            detail::fab_to_fab(loc_copy_tags_2, 0, 0, ncomp,
                               detail::CellAtomicAdd<value_type, value_type>{});
        }
        ((void)0);
    }

    // Note that fab_to_fab is synchronous.
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
                                         int ncomp, std::uint64_t id) const
    -> TagVector<CommSendBufTag<value_type>> const*
{
    using TagType = CommSendBufTag<value_type>;

    auto kit = std::find_if(send_cctc.begin(), send_cctc.end(),
                            [] (CopyComTagsContainer const* p) { return p != nullptr; });
    if (kit == send_cctc.end()) {
        return nullptr;
    }

    auto get_tags = [&] () -> Vector<TagType>
    {
        Vector<TagType> snd_copy_tags;
        char* pbuf = send_data[0];
        const int N_snds = send_data.size();
        for (int j = 0; j < N_snds; ++j)
        {
            if (send_size[j] > 0)
            {
                char* dptr = send_data[j];
                auto const& cctc = *send_cctc[j];
                for (auto const& tag : cctc)
                {
                    snd_copy_tags.emplace_back
                        (TagType{this->const_array(tag.srcIndex), dptr-pbuf, tag.sbox});
                    dptr += (tag.sbox.numPts() * ncomp * sizeof(BUF));
                }
            }
        }
        return snd_copy_tags;
    };

    TagVector<TagType>* tv;
    std::tuple<std::uint64_t,std::size_t,int> key{id, sizeof(BUF), ncomp};

    if (auto it = m_send_copy_handler.find(key); it != m_send_copy_handler.end()) {
        tv = it->second.get();
    } else {
        if (m_send_copy_handler.size() > 32) {
            // Just in case. If this is used in ParallelCopy, it's possible
            // that the sending FabArray is the same, but the receiving
            // FabArray is different every time. Then the size of this map
            // could increase indefinitely.
            m_send_copy_handler.clear();
        }
        auto snd_copy_tags = get_tags();
        auto utv = std::make_unique<TagVector<TagType>>(snd_copy_tags);
        tv = utv.get();
        m_send_copy_handler[key] = std::move(utv);
    }

    return tv;
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::pack_send_buffer_gpu (FabArray<FAB> const& src, int scomp, int ncomp,
                                     Vector<char*> const& send_data,
                                     Vector<std::size_t> const& send_size,
                                     Vector<CopyComTagsContainer const*> const& send_cctc,
                                     std::uint64_t id)
{
    const int N_snds = send_data.size();
    if (N_snds == 0) { return; }

    using TagType = CommSendBufTag<value_type>;

    auto* tv = src.template get_send_copy_tag_vector<BUF>
        (send_data, send_size, send_cctc, ncomp, id);
    if (tv == nullptr) { return; }

    char* pbuffer = send_data[0];

    detail::ParallelFor_doit(*tv,
        [=] AMREX_GPU_DEVICE (
            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
    {
        if (icell < ncells) {
            Array4<BUF> dfab{(BUF*)(pbuffer+tag.poff),
                amrex::begin(tag.bx), amrex::end(tag.bx), ncomp};
            for (int n = 0; n < ncomp; ++n) {
                dfab(i,j,k,n) = (BUF)tag.sfab(i,j,k,n+scomp);
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
                                         int ncomp, std::uint64_t id)
    -> TagVector<CommRecvBufTag<value_type>> const*
{
    using TagType = CommRecvBufTag<value_type>;

    auto kit = std::find_if(recv_cctc.begin(), recv_cctc.end(),
                            [] (CopyComTagsContainer const* p) { return p != nullptr; });
    if (kit == recv_cctc.end()) {
        return nullptr;
    }

    auto get_tags = [&] () -> Vector<TagType>
    {
        Vector<TagType> recv_copy_tags;
        char* pbuf = recv_data[0];
        const int N_rcvs = recv_cctc.size();
        for (int k = 0; k < N_rcvs; ++k)
        {
            if (recv_size[k] > 0)
            {
                char* dptr = recv_data[k];
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    const int li = this->localindex(tag.dstIndex);
                    recv_copy_tags.emplace_back
                        (TagType{this->atLocalIdx(li).array(), dptr-pbuf, tag.dbox});
                    dptr += tag.dbox.numPts() * ncomp * sizeof(BUF);
                }
            }
        }
        return recv_copy_tags;
    };

    TagVector<TagType>* tv;
    std::tuple<std::uint64_t,std::size_t,int> key{id, sizeof(BUF), ncomp};

    if (auto it = m_recv_copy_handler.find(key); it != m_recv_copy_handler.end()) {
        tv = it->second.get();
    } else {
        if (m_recv_copy_handler.size() > 32) {
            // Just in case. If this is used in ParallelCopy, it's possible
            // that the receiving FabArray is the same, but the sending
            // FabArray is different every time. Then the size of this map
            // could increase indefinitely.
            m_recv_copy_handler.clear();
        }
        auto recv_copy_tags = get_tags();
        auto utv = std::make_unique<TagVector<TagType>>(recv_copy_tags);
        tv = utv.get();
        m_recv_copy_handler[key] = std::move(utv);
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
                                       CpOp op, bool is_thread_safe, std::uint64_t id,
                                       bool deterministic)
{
    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) { return; }

    bool use_mask = false;
    if (!is_thread_safe)
    {
        if ((op == FabArrayBase::COPY && !amrex::IsStoreAtomic<value_type>::value) ||
            (op == FabArrayBase::ADD  && !amrex::HasAtomicAdd <value_type>::value))
        {
            use_mask = true;
        }
    }

    if (deterministic)
    {
        AMREX_ALWAYS_ASSERT(op == FabArrayBase::ADD); // Only ADD for now
        using TagType = Array4CopyTag<value_type,BUF>;
        Vector<TagType> tags;
        tags.reserve(N_rcvs);
        for (int k = 0; k < N_rcvs; ++k) {
            if (recv_size[k] > 0) {
                char const* dptr = recv_data[k];
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc) {
                    tags.emplace_back(
                        TagType{dst.array(tag.dstIndex), tag.dstIndex,
                                Array4<BUF const>((BUF const*)dptr,
                                                  amrex::begin(tag.dbox),
                                                  amrex::end(tag.dbox), ncomp),
                                tag.dbox, Dim3{0,0,0}});
                    dptr += tag.dbox.numPts() * ncomp * sizeof(BUF);
                }
            }
        }
        if constexpr (amrex::IsAddAssignable<value_type>::value) {
            detail::deterministic_fab_to_fab<value_type,BUF>
                (tags, 0, dcomp, ncomp, detail::CellAdd<value_type,BUF>{});
        } else {
            amrex::Abort("SumBoundary requires operator+=");
        }
    }
    else if (!use_mask)
    {
        using TagType = CommRecvBufTag<value_type>;
        auto* tv = dst.template get_recv_copy_tag_vector<BUF>
            (recv_data, recv_size, recv_cctc, ncomp, id);
        if (tv == nullptr) { return; }

        char* pbuffer = recv_data[0];

        if (op == FabArrayBase::COPY)
        {
            detail::ParallelFor_doit(*tv,
                [=] AMREX_GPU_DEVICE (
                    int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
            {
                if (icell < ncells) {
                    Array4<BUF const> sfab{(BUF const*)(pbuffer+tag.poff),
                        amrex::begin(tag.bx), amrex::end(tag.bx), ncomp};
                    for (int n = 0; n < ncomp; ++n) {
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
                        int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
                {
                    if (icell < ncells) {
                        Array4<BUF const> sfab{(BUF const*)(pbuffer+tag.poff),
                            amrex::begin(tag.bx), amrex::end(tag.bx), ncomp};
                        for (int n = 0; n < ncomp; ++n) {
                            tag.dfab(i,j,k,n+dcomp) += (value_type)sfab(i,j,k,n);
                        }
                    }
                });
            } else {
                if constexpr (amrex::HasAtomicAdd<value_type>::value) {
                    detail::ParallelFor_doit(*tv,
                        [=] AMREX_GPU_DEVICE (
                            int icell, int ncells, int i, int j, int k, TagType const& tag) noexcept
                    {
                        if (icell < ncells) {
                            Array4<BUF const> sfab{(BUF const*)(pbuffer+tag.poff),
                                amrex::begin(tag.bx), amrex::end(tag.bx), ncomp};
                            for (int n = 0; n < ncomp; ++n) {
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
    else
    {
        char* pbuffer = recv_data[0];

        using TagType = Array4CopyTag<value_type, BUF>;
        Vector<TagType> recv_copy_tags;
        recv_copy_tags.reserve(N_rcvs);

        Vector<BaseFab<int> > maskfabs(dst.local_size());
        Vector<Array4Tag<int> > masks_unique;
        masks_unique.reserve(dst.local_size());
        Vector<Array4Tag<int> > masks;

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
                            dst.atLocalIdx(li).array(), tag.dstIndex,
                            amrex::makeArray4((BUF const*)(dptr), tag.dbox, ncomp),
                            tag.dbox,
                            Dim3{0,0,0}
                        });
                    dptr += tag.dbox.numPts() * ncomp * sizeof(BUF);

                    if (!maskfabs[li].isAllocated()) {
                        maskfabs[li].resize(dst.atLocalIdx(li).box());
                        masks_unique.emplace_back(Array4Tag<int>{maskfabs[li].array()});
                    }
                    masks.emplace_back(Array4Tag<int>{maskfabs[li].array()});
                }
                BL_ASSERT(dptr <= pbuffer + offset + recv_size[k]);
            }
        }

        amrex::ParallelFor(masks_unique,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, Array4Tag<int> const& msk) noexcept
            {
                msk.dfab(i,j,k) = 0;
            });

        if (op == FabArrayBase::COPY)
        {
            detail::fab_to_fab_atomic_cpy<value_type, BUF>(
                recv_copy_tags, 0, dcomp, ncomp, masks);
        }
        else
        {
            detail::fab_to_fab_atomic_add<value_type, BUF>(
                recv_copy_tags, 0, dcomp, ncomp, masks);
        }

        // There is Gpu::streamSynchronize in fab_to_fab.
    }
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
