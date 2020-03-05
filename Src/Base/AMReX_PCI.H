#ifndef AMREX_PCI_H_
#define AMREX_PCI_H_

template <class FAB>
void
FabArray<FAB>::PC_local_cpu (const CPC& thecpc, FabArray<FAB> const& src,
                             int scomp, int dcomp, int ncomp, CpOp op)
{
    int N_locs = thecpc.m_LocTags->size();
    if (N_locs == 0) return;
    bool is_thread_safe = thecpc.m_threadsafe_loc;

    if (is_thread_safe)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];
            if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                // avoid self copy or plus
                const FAB* sfab = &(src[tag.srcIndex]);
                      FAB* dfab = &(get(tag.dstIndex));
                if (op == FabArrayBase::COPY)
                {
                    dfab->template copy<RunOn::Host>(*sfab, tag.sbox, scomp, tag.dbox, dcomp, ncomp);
                }
                else
                {
                    dfab->template plus<RunOn::Host>(*sfab, tag.sbox, tag.dbox, scomp, dcomp, ncomp);
                }
            }
        }
    }
    else
    {
        LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
        for (int i = 0; i < N_locs; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];
            if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                loc_copy_tags[tag.dstIndex].push_back
                    ({src.fabPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
            }
        }

#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            const auto& tags = loc_copy_tags[mfi];
            auto dfab = this->array(mfi);
            if (op == FabArrayBase::COPY)
            {
                for (auto const & tag : tags)
                {
                    auto const sfab = tag.sfab->array();
                    Dim3 offset = tag.offset.dim3();
                    amrex::LoopConcurrentOnCpu (tag.dbox, ncomp,
                    [=] (int i, int j, int k, int n) noexcept
                    {
                        dfab(i,j,k,dcomp+n) = sfab(i+offset.x,j+offset.y,k+offset.z,scomp+n);
                    });
                }
            }
            else
            {
                for (auto const & tag : tags)
                {
                    auto const sfab = tag.sfab->array();
                    Dim3 offset = tag.offset.dim3();
                    amrex::LoopConcurrentOnCpu (tag.dbox, ncomp,
                    [=] (int i, int j, int k, int n) noexcept
                    {
                        dfab(i,j,k,dcomp+n) += sfab(i+offset.x,j+offset.y,k+offset.z,scomp+n);
                    });
                }
            }
        }
    }
}

#ifdef AMREX_USE_GPU
template <class FAB>
void
FabArray<FAB>::PC_local_gpu (const CPC& thecpc, FabArray<FAB> const& src,
                             int scomp, int dcomp, int ncomp, CpOp op)
{
    int N_locs = thecpc.m_LocTags->size();
    if (N_locs == 0) return;
    bool is_thread_safe = thecpc.m_threadsafe_loc;

    typedef Array4CopyTag<value_type> TagType;
    Vector<TagType> loc_copy_tags;
    loc_copy_tags.reserve(N_locs);

    Vector<BaseFab<int> > maskfabs;
    Vector<Array4<int> > masks;
    if (!is_thread_safe)
    {
        if ((op == FabArrayBase::COPY and !amrex::IsStoreAtomic<value_type>::value) or
            (op == FabArrayBase::ADD  and !amrex::HasAtomicAdd <value_type>::value))
        {
            maskfabs.resize(this->local_size());
            masks.reserve(N_locs);
        }
    }

    for (int i = 0; i < N_locs; ++i)
    {
        const CopyComTag& tag = (*thecpc.m_LocTags)[i];
        if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
            int li = this->localindex(tag.dstIndex);
            loc_copy_tags.push_back
                ({this->atLocalIdx(li).array(),
                  src.fabPtr(tag.srcIndex)->const_array(),
                  tag.dbox,
                  (tag.sbox.smallEnd()-tag.dbox.smallEnd()).dim3()});

            if (maskfabs.size() > 0) {
                if (!maskfabs[li].isAllocated()) {
                    maskfabs[li].resize(this->atLocalIdx(li).box());
                }
                masks.push_back(maskfabs[li].array());
            }
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
            detail::fab_to_fab<value_type>(loc_copy_tags, scomp, dcomp, ncomp,
                                           detail::CellStore<value_type>());
        } else {
            detail::fab_to_fab_atomic_cpy<value_type>(loc_copy_tags, scomp, dcomp, ncomp, masks);
        }
    }
    else
    {
        if (is_thread_safe) {
            detail::fab_to_fab<value_type>(loc_copy_tags, scomp, dcomp, ncomp,
                                           detail::CellAdd<value_type>());
        } else {
            detail::fab_to_fab_atomic_add<value_type>(loc_copy_tags, scomp, dcomp, ncomp, masks);
        }
    }
}
#endif

#endif
