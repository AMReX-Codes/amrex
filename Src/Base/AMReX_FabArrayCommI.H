
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

template <class FAB>
template <class FOO, class BAR>  // FOO fools nvcc
void
FabArray<FAB>::FBEP_nowait (int scomp, int ncomp, const IntVect& nghost,
                            const Periodicity& period, bool cross,
			    bool enforce_periodicity_only)
{
    fb_cross = cross;
    fb_epo   = enforce_periodicity_only;
    fb_scomp = scomp;
    fb_ncomp = ncomp;
    fb_nghost = nghost;
    fb_period = period;

    fb_recv_reqs.clear();

    bool work_to_do;
    if (enforce_periodicity_only) {
	work_to_do = period.isAnyPeriodic();
    } else {
	work_to_do = nghost.max() > 0;
    }
    if (!work_to_do) return;

    const FB& TheFB = getFB(nghost, period, cross, enforce_periodicity_only);

    if (ParallelContext::NProcsSub() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = (*TheFB.m_LocTags).size();
        bool is_thread_safe = FAB::isCopyOMPSafe() && TheFB.m_threadsafe_loc;
        if (Gpu::inLaunchRegion() || !is_thread_safe)
        {
            LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
            for (int i = 0; i < N_loc; ++i)
            {
                const CopyComTag& tag = (*TheFB.m_LocTags)[i];

                BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
                BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

                loc_copy_tags[tag.dstIndex].push_back
                    ({this->fabHostPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
            }
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*this); mfi.isValid(); ++mfi)
            {
                const auto& tags = loc_copy_tags[mfi];
                auto dfab = this->array(mfi);
                for (auto const & tag : tags)
                {
                    auto const sfab = tag.sfab->array();
                    const auto offset = tag.offset.dim3();
                    AMREX_HOST_DEVICE_FOR_4D ( tag.dbox, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,n+scomp) = sfab(i+offset.x,j+offset.y,k+offset.z,n+scomp);
                    });
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < N_loc; ++i)
            {
                const CopyComTag& tag = (*TheFB.m_LocTags)[i];

                BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
                BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());
                
                const FAB* sfab = &(get(tag.srcIndex));
                      FAB* dfab = &(get(tag.dstIndex));
                dfab->copy(*sfab, tag.sbox, scomp, tag.dbox, scomp, ncomp);
            }
        }

        return;
    }

#ifdef BL_USE_MPI

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int preSeqNum;
    if (!FAB::preAllocatable()) {
        preSeqNum = ParallelDescriptor::SeqNum();
    }
    int SeqNum = ParallelDescriptor::SeqNum();

    const int N_locs = TheFB.m_LocTags->size();
    const int N_rcvs = TheFB.m_RcvTags->size();
    const int N_snds = TheFB.m_SndTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0)
        // No work to do.
        return;

    //
    // Before we post recv, let's preprocess sends in case FAB is not preAllocatable
    //
    char*&                          the_send_data = fb_the_send_data;
    Vector<char*> &                     send_data = fb_send_data;
    Vector<int>                         send_size;
    Vector<int>                         send_rank;
    Vector<MPI_Request>&                send_reqs = fb_send_reqs;
    Vector<const CopyComTagsContainer*> send_cctc;
    Vector<Vector<int> >                indv_send_size;
    Vector<MPI_Request>                 pre_reqs;

    fb_tag = SeqNum;

    if (N_snds > 0)
    {
        fb_send_data.clear();
        fb_send_reqs.clear();

	send_data.reserve(N_snds);
	send_size.reserve(N_snds);
	send_rank.reserve(N_snds);
        send_reqs.reserve(N_snds);
	send_cctc.reserve(N_snds);
        indv_send_size.reserve(N_snds);

        std::size_t total_volume = 0;
        for (auto const& kv : *TheFB.m_SndTags)
        {
            Vector<int> iss;                
            auto const& cctc = kv.second;

            std::size_t nbytes = 0;
            if (FAB::preAllocatable())
            {
                for (auto const& cct : kv.second)
                {
                    nbytes += (*this)[cct.srcIndex].nBytes(cct.sbox,scomp,ncomp);
                }
            }
            else
            {
                for (auto const& tag : cctc)
                {
                    std::size_t b = (*this)[tag.srcIndex].nBytes(tag.sbox,scomp,ncomp);
                    nbytes += b;
                    iss.push_back(static_cast<int>(b));
                }
            }
            
            BL_ASSERT(nbytes < std::numeric_limits<int>::max());
            
            total_volume += nbytes;

            send_data.push_back(nullptr);
            send_size.push_back(static_cast<int>(nbytes));
            send_rank.push_back(kv.first);
            send_reqs.push_back(MPI_REQUEST_NULL);
            send_cctc.push_back(&cctc);
            indv_send_size.push_back(std::move(iss));
        }

        if (total_volume > 0)
        {
            the_send_data = static_cast<char*>(amrex::The_FA_Arena()->alloc(total_volume));
            char* p = the_send_data;
            for (int i = 0, N = send_size.size(); i < N; ++i) {
                if (send_size[i] > 0) {
                    send_data[i] = p;
                    p += send_size[i];
                }
            }
        }
    }

    if (!FAB::preAllocatable())
    {
        pre_reqs.resize(N_snds,MPI_REQUEST_NULL);
        for (int j=0; j<N_snds; ++j)
        {
            pre_reqs[j] = ParallelDescriptor::Asend
                (indv_send_size[j].data(), indv_send_size[j].size(),
                 ParallelContext::global_to_local_rank(send_rank[j]),
                 preSeqNum,
                 ParallelContext::CommunicatorSub()).req();
        }
    }

    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    fb_the_recv_data = nullptr;

    if (N_rcvs > 0) {
        PostRcvs(*TheFB.m_RcvTags, fb_the_recv_data,
                 fb_recv_data, fb_recv_size, fb_recv_from, fb_recv_reqs,
                 scomp, ncomp, SeqNum, preSeqNum);
        fb_recv_stat.resize(N_rcvs);
    }

    //
    // Post send's
    //
    if (N_snds > 0)
    {
        bool is_thread_safe = FAB::isCopyOMPSafe();
#ifdef _OPENMP
#pragma omp parallel if (is_thread_safe && Gpu::notInLaunchRegion())
#endif
        for (Gpu::StreamIter sit(N_snds,is_thread_safe); sit.isValid(); ++sit)
        {
            const int j = sit();
            char* dptr = send_data[j];
            if (dptr != nullptr)
            {
                auto const& cctc = *send_cctc[j];
                for (auto const& tag : cctc)
                {
                    const Box& bx = tag.sbox;
                    auto const sfab = this->array(tag.srcIndex);
                    auto pfab = amrex::makeArray4((value_type*)(dptr),bx);
                    AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, ii, jj, kk, n,
                    {
                        pfab(ii,jj,kk,n) = sfab(ii,jj,kk,n+scomp);
                    });
                    dptr += (bx.numPts() * ncomp * sizeof(value_type));
                }
                BL_ASSERT(dptr == send_data[j] + send_size[j]); 
            }
        }

        int send_counter = 0;
        while (send_counter < N_snds)
        {
            int j;
            if (FAB::preAllocatable())
            {
                j = send_counter;
            }
            else
            {
                MPI_Status status;
                ParallelDescriptor::Waitany(pre_reqs, j, status);
            }

            BL_ASSERT(send_size[j] > 0);

            send_reqs[j] = ParallelDescriptor::Asend
                (send_data[j], send_size[j],
                 ParallelContext::global_to_local_rank(send_rank[j]),
                 SeqNum,
                 ParallelContext::CommunicatorSub()).req();

            ++send_counter;
	}
    }

    FillBoundary_test();

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    {
        bool is_thread_safe = FAB::isCopyOMPSafe() && TheFB.m_threadsafe_loc;
        if (Gpu::inLaunchRegion() || !is_thread_safe)
        {
            LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
            for (int i = 0; i < N_locs; ++i)
            {
                const CopyComTag& tag = (*TheFB.m_LocTags)[i];

                BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.dstIndex]));
                BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.srcIndex]));
	    
                if (distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc())
                {
                    loc_copy_tags[tag.dstIndex].push_back
                        ({this->fabHostPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
                }
            }

            FillBoundary_test();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            {
                int cuda_stream = 0;
                for (MFIter mfi(*this); mfi.isValid(); ++mfi)
                {
                    const auto& tags = loc_copy_tags[mfi];
                    const auto& dfab = this->array(mfi);
                    for (const auto& tag : tags ) {
                        const Dim3 offset = tag.offset.dim3();
                        const auto sfab = tag.sfab->array();
#ifdef AMREX_USE_GPU
                        Gpu::Device::setStreamIndex(cuda_stream++);
#endif
                        AMREX_HOST_DEVICE_FOR_4D ( tag.dbox, ncomp, i, j, k, n,
                        {
                            dfab(i,j,k,scomp+n) = sfab(i+offset.x,j+offset.y,k+offset.z,scomp+n);
                        });
                    }
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < N_locs; ++i)
            {
                const CopyComTag& tag = (*TheFB.m_LocTags)[i];

                BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.dstIndex]));
                BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.srcIndex]));
	    
                if (distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc())
                {
                    const FAB* sfab = &(get(tag.srcIndex));
                          FAB* dfab = &(get(tag.dstIndex));
                    dfab->copy(*sfab, tag.sbox, scomp, tag.dbox, scomp, ncomp);
                }
	    }
	}
    }

    FillBoundary_test();
#endif /*BL_USE_MPI*/
}

template <class FAB>
template <class FOO, class BAR>  // FOO fools nvcc
void
FabArray<FAB>::FillBoundary_finish ()
{
    BL_PROFILE("FillBoundary_finish()");

    if ( n_grow.allLE(IntVect::TheZeroVector()) && !fb_epo ) return; // For epo (Enforce Periodicity Only), there may be no ghost cells.

    if (ParallelContext::NProcsSub() == 1) return;

#ifdef BL_USE_MPI

    const FB& TheFB = getFB(fb_nghost,fb_period,fb_cross,fb_epo);

    const int N_rcvs = TheFB.m_RcvTags->size();
    const int N_snds = TheFB.m_SndTags->size();

    Vector<const CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);
    LayoutData<Vector<VoidCopyTag> > recv_copy_tags;
    if (N_rcvs > 0)
    {
	for (int k = 0; k < N_rcvs; k++) 
	{
            if (fb_recv_size[k] > 0)
            {
                auto const& cctc = TheFB.m_RcvTags->at(fb_recv_from[k]);
                recv_cctc[k] = &cctc;
            }
	}
        bool is_thread_safe = FAB::isCopyOMPSafe() && TheFB.m_threadsafe_rcv;
        if (Gpu::inLaunchRegion() || !is_thread_safe)
        {
            recv_copy_tags.define(boxArray(),DistributionMap());
            for (int k = 0; k < N_rcvs; ++k)
            {
                const char* dptr = fb_recv_data[k];
                if (dptr != nullptr)
                {
                    auto const& cctc = *recv_cctc[k];
                    for (auto const& tag : cctc)
                    {
                        recv_copy_tags[tag.dstIndex].push_back({dptr,tag.dbox});
                        dptr += tag.dbox.numPts() * fb_ncomp * sizeof(value_type);
                    }
                    BL_ASSERT(dptr == fb_recv_data[k] + fb_recv_size[k]);
                }
            }
        }
    }

    int actual_n_rcvs = N_rcvs - std::count(fb_recv_data.begin(), fb_recv_data.end(), nullptr);

    if (actual_n_rcvs > 0) {
        ParallelDescriptor::Waitall(fb_recv_reqs, fb_recv_stat);
#ifdef AMREX_DEBUG
        if (!CheckRcvStats(fb_recv_stat, fb_recv_size, MPI_CHAR, fb_tag))
        {
            amrex::Abort("FillBoundary_finish failed with wrong message size");
        }
#endif
    }

    if (N_rcvs > 0)
    {
        if (!recv_copy_tags.empty())
        {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(*this); mfi.isValid(); ++mfi)
            {
                const auto& tags = recv_copy_tags[mfi];
                auto dfab = this->array(mfi);
                const int scomp = fb_scomp;
                const int ncomp = fb_ncomp;
                for (auto const & tag : tags)
                {
                    auto pfab = amrex::makeArray4((value_type*)(tag.p), tag.dbox);
                    AMREX_HOST_DEVICE_FOR_4D ( tag.dbox, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,n+scomp) = pfab(i,j,k,n);
                    });
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int k = 0; k < N_rcvs; ++k)
            {
                const char* dptr = fb_recv_data[k];
                if (dptr != nullptr)
                {
                    auto const& cctc = *recv_cctc[k];
                    for (auto const& tag : cctc)
                    {
                        const Box& bx  = tag.dbox;
                        FAB* dfab = &(get(tag.dstIndex));
                        const int scomp = fb_scomp;
                        const int ncomp = fb_ncomp;
                        dfab->copyFromMem(bx, scomp, ncomp, dptr);
                        dptr += bx.numPts() * ncomp * sizeof(value_type);
                    }
                    BL_ASSERT(dptr == fb_recv_data[k] + fb_recv_size[k]);
                }
            }
        }

        if (fb_the_recv_data)
        {
	    amrex::The_FA_Arena()->free(fb_the_recv_data);
            fb_the_recv_data = nullptr;
	}
    }

    if (N_snds > 0) {
        Vector<MPI_Status> stats;
        FabArrayBase::WaitForAsyncSends(N_snds,fb_send_reqs,fb_send_data,stats);
        amrex::The_FA_Arena()->free(fb_the_send_data);
        fb_the_send_data = nullptr;
    }

#endif // MPI
}

template <class FAB>
void
FabArray<FAB>::ParallelCopy (const FabArray<FAB>& src,
                             int                  scomp,
                             int                  dcomp,
                             int                  ncomp,
                             const IntVect&       snghost,
                             const IntVect&       dnghost,
                             const Periodicity&   period,
                             CpOp                 op,
                             const FabArrayBase::CPC * a_cpc)
{
    BL_PROFILE("FabArray::ParallelCopy()");

    if (size() == 0 || src.size() == 0) return;

    BL_ASSERT(op == FabArrayBase::COPY || op == FabArrayBase::ADD);
    BL_ASSERT(boxArray().ixType() == src.boxArray().ixType());

    BL_ASSERT(src.nGrowVect().allGE(snghost));
    BL_ASSERT(    nGrowVect().allGE(dnghost));

    if ((src.boxArray().ixType().cellCentered() || op == FabArrayBase::COPY) &&
        (boxarray == src.boxarray && distributionMap == src.distributionMap)
	&& snghost == IntVect::TheZeroVector() && dnghost == IntVect::TheZeroVector()
        && !period.isAnyPeriodic())
    {
        //
        // Short-circuit full intersection code if we're doing copy()s or if
        // we're doing plus()s on cell-centered data.  Don't do plus()s on
        // non-cell-centered data this simplistic way.
        //
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe() && Gpu::notInLaunchRegion())
#endif
        for (MFIter fai(*this,TilingIfNotGPU()); fai.isValid(); ++fai)
        {
            const Box& bx = fai.tilebox();

            // avoid self copy or plus
	    if (this != &src) {
                auto const sfab = src.array(fai);
                auto       dfab = this->array(fai);
		if (op == FabArrayBase::COPY) {
                    AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,dcomp+n) = sfab(i,j,k,scomp+n);
                    });
		} else {
                    AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,dcomp+n) += sfab(i,j,k,scomp+n);
                    });
                }
	    }
        }

        return;
    }

    const CPC& thecpc = (a_cpc) ? *a_cpc : getCPC(dnghost, src, snghost, period);

    if (ParallelContext::NProcsSub() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = (*thecpc.m_LocTags).size();
        bool is_thread_safe = FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc;
        if (Gpu::inLaunchRegion() || !is_thread_safe)
        {
            LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(boxArray(),DistributionMap());
            for (int i = 0; i < N_loc; ++i)
            {
                const CopyComTag& tag = (*thecpc.m_LocTags)[i];
                if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                    loc_copy_tags[tag.dstIndex].push_back
                        ({src.fabHostPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
                }
            }
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
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
                        AMREX_HOST_DEVICE_FOR_4D ( tag.dbox, ncomp, i, j, k, n,
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
                        AMREX_HOST_DEVICE_FOR_4D ( tag.dbox, ncomp, i, j, k, n,
                        {
                            dfab(i,j,k,dcomp+n) += sfab(i+offset.x,j+offset.y,k+offset.z,scomp+n);
                        });
                    }
                }
            }
        }
        else
        {
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int i = 0; i < N_loc; ++i)
            {
                const CopyComTag& tag = (*thecpc.m_LocTags)[i];
                if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                    // avoid self copy or plus
                    const FAB* sfab = &(src[tag.srcIndex]);
                          FAB* dfab = &(get(tag.dstIndex));
                    if (op == FabArrayBase::COPY)
                    {
                        dfab->copy(*sfab, tag.sbox, scomp, tag.dbox, dcomp, ncomp);
                    }
                    else
                    {
                        dfab->plus(*sfab, tag.sbox, tag.dbox, scomp, dcomp, ncomp);
                    }
		}
	    }
        }

        return;
    }

#ifdef BL_USE_MPI

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int preSeqNum;
    if (!FAB::preAllocatable()) {
        preSeqNum = ParallelDescriptor::SeqNum();
    }
    int SeqNum  = ParallelDescriptor::SeqNum();

    const int N_snds = thecpc.m_SndTags->size();
    const int N_rcvs = thecpc.m_RcvTags->size();
    const int N_locs = thecpc.m_LocTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0) {
        //
        // No work to do.
        //
        return;
    }

    //
    // Send/Recv at most MaxComp components at a time to cut down memory usage.
    //
    int NCompLeft = ncomp;

    for (int ipass = 0, SC = scomp, DC = dcomp; ipass < ncomp; )
    {
        const int NC = std::min(NCompLeft,FabArrayBase::MaxComp);

        //
        // Before we post recv, let's preprocess sends in case FAB is not preAllocatable
        //
        char*                               the_send_data;
	Vector<char*>                       send_data;
	Vector<int>                         send_size;
	Vector<int>                         send_rank;
	Vector<MPI_Request>                 send_reqs;
	Vector<const CopyComTagsContainer*> send_cctc;
        Vector<Vector<int> >                indv_send_size;
        Vector<MPI_Request>                 pre_reqs;

	if (N_snds > 0)
	{
	    send_data.reserve(N_snds);
	    send_size.reserve(N_snds);
	    send_rank.reserve(N_snds);
            send_reqs.reserve(N_snds);
	    send_cctc.reserve(N_snds);
            indv_send_size.reserve(N_snds);

            std::size_t total_volume = 0;
            for (auto const& kv : *thecpc.m_SndTags)
	    {
                Vector<int> iss;                
                auto const& cctc = kv.second;

                std::size_t nbytes = 0;
                if (FAB::preAllocatable())
                {
                    for (auto const& cct : kv.second)
                    {
                        nbytes += src[cct.srcIndex].nBytes(cct.sbox,SC,NC);
                    }
                }
                else
                {
                    for (auto const& tag : cctc)
                    {
                        std::size_t b = src[tag.srcIndex].nBytes(tag.sbox,SC,NC);
                        nbytes += b;
                        iss.push_back(static_cast<int>(b));
                    }
                }
		
		BL_ASSERT(nbytes < std::numeric_limits<int>::max());

                total_volume += nbytes;

		send_data.push_back(nullptr);
                send_size.push_back(static_cast<int>(nbytes));
                send_rank.push_back(kv.first);
                send_reqs.push_back(MPI_REQUEST_NULL);
                send_cctc.push_back(&cctc);
                indv_send_size.push_back(std::move(iss));
	    }

            if (total_volume > 0)
            {
                the_send_data = static_cast<char*>(amrex::The_FA_Arena()->alloc(total_volume));
                char* p = the_send_data;
                for (int i = 0, N = send_size.size(); i < N; ++i) {
                    if (send_size[i] > 0) {
                        send_data[i] = p;
                        p += send_size[i];
                    }
                }
            }
        }

        if (!FAB::preAllocatable())
        {
            pre_reqs.resize(N_snds,MPI_REQUEST_NULL);
            for (int j=0; j<N_snds; ++j)
            {
                pre_reqs[j] = ParallelDescriptor::Asend
                    (indv_send_size[j].data(), indv_send_size[j].size(),
                     ParallelContext::global_to_local_rank(send_rank[j]),
                     preSeqNum,
                     ParallelContext::CommunicatorSub()).req();
            }
        }

        Vector<int>         recv_from;
        Vector<char*>       recv_data;
        Vector<int>         recv_size;
        Vector<MPI_Request> recv_reqs;
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        char* the_recv_data = nullptr;

        int actual_n_rcvs = 0;
	if (N_rcvs > 0) {
            PostRcvs(*thecpc.m_RcvTags, the_recv_data,
                     recv_data, recv_size, recv_from, recv_reqs, SC, NC, SeqNum, preSeqNum);
            actual_n_rcvs = N_rcvs - std::count(recv_size.begin(), recv_size.end(), 0);
	}

	//
	// Post send's
	// 
	if (N_snds > 0)
	{
            bool is_thread_safe = FAB::isCopyOMPSafe();
#ifdef _OPENMP
#pragma omp parallel if (is_thread_safe && Gpu::notInLaunchRegion())
#endif
            for (Gpu::StreamIter sit(N_snds,is_thread_safe); sit.isValid(); ++sit)
	    {
                const int j = sit();
		char* dptr = send_data[j];
                if (dptr != nullptr)
                {
                    auto const& cctc = *send_cctc[j];
                    for (auto const& tag : cctc)
                    {
                        const Box& bx = tag.sbox;
                        auto const sfab = src.array(tag.srcIndex);
                        auto pfab = amrex::makeArray4((value_type*)(dptr),bx);
                        AMREX_HOST_DEVICE_FOR_4D ( bx, NC, ii, jj, kk, n,
                        {
                            pfab(ii,jj,kk,n) = sfab(ii,jj,kk,SC+n);
                        });

                        dptr += (bx.numPts() * NC * sizeof(value_type));
                    }
                    BL_ASSERT(dptr == send_data[j] + send_size[j]);
		}
	    }

            int send_counter = 0;
            while (send_counter < N_snds)
            {
                int j;

                if (FAB::preAllocatable())
                {
                    j = send_counter;
                }
                else
                {
                    MPI_Status status;
                    ParallelDescriptor::Waitany(pre_reqs, j, status);
                }

                if (send_size[j] > 0)
                {
                    send_reqs[j] = ParallelDescriptor::Asend
                        (send_data[j], send_size[j],
                         ParallelContext::global_to_local_rank(send_rank[j]),
                         SeqNum,
                         ParallelContext::CommunicatorSub()).req();
                }

                ++send_counter;
	    }
	}

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
	{
            bool is_thread_safe = FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc;
            if (Gpu::inLaunchRegion() || !is_thread_safe)
            {
                // gpu version or omp over dest fabs

                LayoutData<Vector<FabCopyTag<FAB> > > copy_tags(boxArray(),DistributionMap());
                for (int j = 0; j < N_locs; ++j) {
                    const CopyComTag& tag = (*thecpc.m_LocTags)[j];
                    if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                        copy_tags[tag.dstIndex].push_back
                            ({src.fabHostPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
                    }
                }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
                    auto dfab = this->array(mfi);
                    const auto& fab_copy_tags = copy_tags[mfi];
                    if (op == FabArrayBase::COPY) {
                        for (auto const& fab_tag : fab_copy_tags) {
                            Dim3 offset = fab_tag.offset.dim3();
                            auto const sfab = fab_tag.sfab->array();
                            AMREX_HOST_DEVICE_FOR_4D ( fab_tag.dbox, NC, i, j, k, n,
                            {
                                dfab(i,j,k,DC+n) = sfab(i+offset.x,j+offset.y,k+offset.z,SC+n);
                            });
                        }
                    } else {
                        for (auto const& fab_tag : fab_copy_tags) {
                            Dim3 offset = fab_tag.offset.dim3();
                            auto const sfab = fab_tag.sfab->array();
                            AMREX_HOST_DEVICE_FOR_4D ( fab_tag.dbox, NC, i, j, k, n,
                            {
                                dfab(i,j,k,DC+n) += sfab(i+offset.x,j+offset.y,k+offset.z,SC+n);
                            });
                        }
                    }
                }
            }
            else
            {
                // cpu version, omp over dest tile boxes
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int j = 0; j < N_locs; ++j)
                {
                    const CopyComTag& tag = (*thecpc.m_LocTags)[j];
                    // avoid self copy or plus
                    if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                        const FAB* sfab = &(src[tag.srcIndex]);
                              FAB* dfab = &(get(tag.dstIndex));
                        if (op == FabArrayBase::COPY) {
                            dfab->copy(*sfab, tag.sbox, SC, tag.dbox, DC, NC);
                        } else {
                            dfab->plus(*sfab, tag.sbox, tag.dbox, SC, DC, NC);
                        }
                    }
                }
            }
        }

        // prepare for receive
        Vector<const CopyComTagsContainer*> recv_cctc;
        LayoutData<Vector<VoidCopyTag> > recv_copy_tags;
        if (N_rcvs > 0)
        {
            recv_cctc.resize(N_rcvs,nullptr);
	    for (int k = 0; k < N_rcvs; ++k)
	    {
                if (recv_size[k] > 0)
                {
                    auto const& cctc = thecpc.m_RcvTags->at(recv_from[k]);
                    recv_cctc[k] = &cctc;
                }
	    }

            bool is_thread_safe = FAB::isCopyOMPSafe() && thecpc.m_threadsafe_rcv;
            if (Gpu::inLaunchRegion() || !is_thread_safe)
            {
                recv_copy_tags.define(boxArray(),DistributionMap());
                for (int k = 0; k < N_rcvs; ++k) {
                    const char* dptr = recv_data[k];
                    if (dptr != nullptr)
                    {
                        auto const& cctc = *recv_cctc[k];
                        for (auto const& tag : cctc)
                        {
                            recv_copy_tags[tag.dstIndex].push_back({dptr,tag.dbox});
                            dptr += tag.dbox.numPts() * NC * sizeof(value_type);
                        }
                        BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
                    }
                }
            }
        }

	//
	//  wait and unpack
	//
        if (actual_n_rcvs > 0) {
            Vector<MPI_Status> stats(N_rcvs);
            ParallelDescriptor::Waitall(recv_reqs, stats);
            if (!CheckRcvStats(stats, recv_size, MPI_CHAR, SeqNum))
            {
                amrex::Abort("ParallelCopy failed with wrong message size");
            }
	}

	if (N_rcvs > 0)
	{
            if (!recv_copy_tags.empty())
            {
                // gpu version or omp over dest fabs
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
                    auto dfab = this->array(mfi);
                    const auto& void_copy_tags = recv_copy_tags[mfi];
                    if (op == FabArrayBase::COPY)
                    {
                        for (auto const& void_tag : void_copy_tags) {
                            auto pfab = amrex::makeArray4((value_type*)(void_tag.p), void_tag.dbox);
                            AMREX_HOST_DEVICE_FOR_4D ( void_tag.dbox, NC, i, j, k, n,
                            {
                                dfab(i,j,k,DC+n) = pfab(i,j,k,n);
                            });
                        }
                    } else {
                        for (auto const& void_tag : void_copy_tags) {
                            auto pfab = amrex::makeArray4((value_type*)(void_tag.p), void_tag.dbox);
                            AMREX_HOST_DEVICE_FOR_4D ( void_tag.dbox, NC, i, j, k, n,
                            {
                                dfab(i,j,k,DC+n) += pfab(i,j,k,n);
                            });
                        }
                    }
                }
            }
            else
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
                            FAB* dfab = &(get(tag.dstIndex));
                            if (op == FabArrayBase::COPY)
                            {
                                dfab->copyFromMem(tag.dbox,DC,NC,dptr);
                            }
                            else
                            {
                                dfab->addFromMem(tag.dbox, DC, NC, dptr);
                            }
                        
                            dptr += tag.dbox.numPts() * NC * sizeof(value_type);
                        }
                        BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
                    }
                }
            }

            if (the_recv_data)
            {
                amrex::The_FA_Arena()->free(the_recv_data);
                the_recv_data = nullptr;
            }
	}
	
        if (N_snds > 0) {
            if (! thecpc.m_SndTags->empty()) {
                Vector<MPI_Status> stats;
                FabArrayBase::WaitForAsyncSends(N_snds,send_reqs,send_data,stats);
	    }
            amrex::The_FA_Arena()->free(the_send_data);
            the_send_data = nullptr;
        }

        ipass     += NC;
        SC        += NC;
        DC        += NC;
        NCompLeft -= NC;
    }

    return;

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::copyTo (FAB&       dest,
		       const Box& subbox,
		       int        scomp,
		       int        dcomp,
		       int        ncomp,
		       int        nghost) const
{
    BL_PROFILE("FabArray::copy(fab)");

    BL_ASSERT(dcomp + ncomp <= dest.nComp());
    BL_ASSERT(nghost <= nGrow());

    if (ParallelContext::NProcsSub() == 1)
    {
        for (int j = 0, N = size(); j < N; ++j)
        {
	    const Box& bx = amrex::grow(boxarray[j],nghost);
	    const Box& destbox = bx & subbox;
	    if (destbox.ok())
            {
                dest.copy(get(j),destbox,scomp,destbox,dcomp,ncomp);
            }
        }

        return;
    }

    //
    //  Note that subbox must be identical on each process!!
    //
#ifdef AMREX_DEBUG
    {
	BoxCommHelper bch(subbox);	
	ParallelDescriptor::Bcast(bch.data(), bch.size(), 0, ParallelContext::CommunicatorSub());
	const Box& bx0 = bch.make_box();
	BL_ASSERT(subbox == bx0);
    }
#endif

    FAB ovlp;

    std::vector< std::pair<int,Box> > isects;
    boxarray.intersections(subbox, isects, false, nghost);

    for (int j = 0, M = isects.size(); j < M; ++j)
    {
	const int  k  = isects[j].first;
	const Box& bx = isects[j].second;

	ovlp.resize(bx,ncomp);

	if (ParallelDescriptor::MyProc() == distributionMap[k])
	{
	    ovlp.copy(get(k),bx,scomp,bx,0,ncomp);
	}

	const int N = bx.numPts()*ncomp;

	ParallelDescriptor::Bcast(ovlp.dataPtr(),N,
                                  ParallelContext::global_to_local_rank(distributionMap[k]),
                                  ParallelContext::CommunicatorSub());

	dest.copy(ovlp,bx,0,bx,dcomp,ncomp);
    }
}


#ifdef BL_USE_MPI
template <class FAB>
void
FabArray<FAB>::PostRcvs (const MapOfCopyComTagContainers&  m_RcvTags,
                         char*&                            the_recv_data,
                         Vector<char*>&                    recv_data,
                         Vector<int>&                      recv_size,
                         Vector<int>&                      recv_from,
                         Vector<MPI_Request>&              recv_reqs,
                         int                               icomp,
                         int                               ncomp,
                         int                               SeqNum,
                         int                               preSeqNum)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();
    recv_reqs.clear();

    std::size_t TotalRcvsVolume = 0;
    for (const auto& kv : m_RcvTags) // loop over senders
    {
        std::size_t nbytes = 0;
        if (FAB::preAllocatable())
        {
            for (auto const& cct : kv.second)
            {
                nbytes += (*this)[cct.dstIndex].nBytes(cct.dbox,icomp,ncomp);
            }
        }

        BL_ASSERT(nbytes < std::numeric_limits<int>::max());

        TotalRcvsVolume += nbytes;

        recv_data.push_back(nullptr);
        recv_size.push_back(static_cast<int>(nbytes));
        recv_from.push_back(kv.first);
        recv_reqs.push_back(MPI_REQUEST_NULL);
    }

    const int nrecv = recv_from.size();
    
    Vector<Vector<int> > indv_recv_size;
    Vector<MPI_Request> pre_reqs(nrecv,MPI_REQUEST_NULL);

    MPI_Comm comm = ParallelContext::CommunicatorSub();

    if (!FAB::preAllocatable())
    {
        BL_ASSERT(preSeqNum >= 0);
        indv_recv_size.resize(nrecv);
        for (int k = 0; k < nrecv; ++k)
        {
            auto n = m_RcvTags.at(recv_from[k]).size();
            indv_recv_size[k].resize(n);
            pre_reqs[k] = ParallelDescriptor::Arecv(indv_recv_size[k].data(), n,
                                                    ParallelContext::global_to_local_rank(recv_from[k]),
                                                    preSeqNum, comm).req();
        }
    }

    if (TotalRcvsVolume == 0)
    {
        the_recv_data = nullptr;
    }
    else
    {
        the_recv_data = static_cast<char*>(amrex::The_FA_Arena()->alloc(TotalRcvsVolume));
    }

    int recv_counter = 0;
    char* p = the_recv_data;
    while (recv_counter < nrecv)
    {
        int i;

        if (FAB::preAllocatable())
        {
            i = recv_counter;
        }
        else
        {
            MPI_Status status;
            ParallelDescriptor::Waitany(pre_reqs, i, status);

            for (auto x : indv_recv_size[i])
            {
                recv_size[i] += x;
            }
        }

        if (recv_size[i] > 0)
        {
            recv_data[i] = p;
            p += recv_size[i];
            recv_reqs[i] = ParallelDescriptor::Arecv(recv_data[i], recv_size[i],
                                                     ParallelContext::global_to_local_rank(recv_from[i]),
                                                     SeqNum, comm).req();
        }
        
        ++recv_counter;
    }
}
#endif

template <class FAB>
void
FabArray<FAB>::Redistribute (const FabArray<FAB>& src,
                             int                  scomp,
                             int                  dcomp,
                             int                  ncomp,
                             const IntVect&       nghost)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(boxArray() == src.boxArray(),
                                     "FabArray::Redistribute: must have the same BoxArray");

    if (ParallelContext::NProcsSub() == 1)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter fai(*this,true); fai.isValid(); ++fai)
        {
            const Box& bx = fai.growntilebox(nghost);
            auto const sfab = src.array(fai);
            auto       dfab = this->array(fai);
            AMREX_HOST_DEVICE_FOR_4D ( bx, ncomp, i, j, k, n,
            {
                dfab(i,j,k,n+dcomp) = sfab(i,j,k,n+scomp);
            });
        }

        return;
    }

#ifdef BL_USE_MPI

    FabArrayBase::CPC cpc(boxArray(), nghost, DistributionMap(), src.DistributionMap());

    ParallelCopy(src, scomp, dcomp, ncomp, nghost, nghost, Periodicity::NonPeriodic(),
                 FabArrayBase::COPY, &cpc);

#endif
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_test ()
{
#ifdef BL_USE_MPI
#ifndef AMREX_DEBUG
    if (!fb_recv_reqs.empty()) {
        int flag;
        MPI_Testall(fb_recv_reqs.size(), fb_recv_reqs.data(), &flag,
                    fb_recv_stat.data());
    }
#endif
#endif
}

#ifdef BL_USE_MPI
struct CommInfo {
    CommInfo () : req(MPI_REQUEST_NULL), data(nullptr), size(0) {}
    Vector<FabArrayBase::CopyComTagsContainer const*> cctcs;
    Vector<int> imfs;
    MPI_Request req;
    char*       data;
    std::size_t size;
};
#endif

template <class FAB>
void
FillBoundary (Vector<FabArray<FAB>*> const& mf, const Periodicity& period)
{
    BL_PROFILE("FillBoundary(Vector)");

    const int nummfs = mf.size();
    if (ParallelContext::NProcsSub() == 1 || !FAB::preAllocatable())
    {
        for (int imf = 0; imf < nummfs; ++imf) {
            mf[imf]->FillBoundary(period);
        }
    }
    else
    {
#ifdef BL_USE_MPI
        int SeqNum = ParallelDescriptor::SeqNum();
        MPI_Comm comm = ParallelContext::CommunicatorSub();
        int myproc = ParallelDescriptor::MyProc();
        using value_type = typename FAB::value_type;

        Vector<int> scomp(nummfs,0);
        Vector<int> ncomp;
        Vector<IntVect> nghost;
        for (auto pmf : mf) {
            ncomp.push_back(pmf->nComp());
            nghost.push_back(pmf->nGrowVect());
        }

        Vector<FabArrayBase::FB const*> TheFB;
        int N_locs_tot = 0, N_rcvs_tot = 0, N_snds_tot = 0;
        for (int imf = 0; imf < nummfs; ++imf) {
            TheFB.push_back(&(mf[imf]->getFB(nghost[imf], period, false, false)));
            N_locs_tot += TheFB[imf]->m_LocTags->size();
            N_rcvs_tot += TheFB[imf]->m_RcvTags->size();
            N_snds_tot += TheFB[imf]->m_SndTags->size();
        }

        if (N_locs_tot == 0 && N_rcvs_tot == 0 && N_snds_tot == 0) {
            return;
        }

        char* the_send_data = nullptr;
        char* the_recv_data = nullptr;

        std::map<int,CommInfo> send_info;
        std::map<int,CommInfo> recv_info;

        if (N_rcvs_tot > 0)
        {
            std::size_t total_volume = 0;
            for (int imf = 0; imf < nummfs; ++imf)
            {
                for (auto const& kv : *(TheFB[imf]->m_RcvTags))
                {
                    const int recv_from = kv.first;
                    auto const& cctc = kv.second;

                    std::size_t nbytes = 0;
                    for (auto const& cct : cctc) {
                        nbytes += (*mf[imf])[cct.dstIndex].nBytes(cct.dbox,scomp[imf],ncomp[imf]);
                    }
                    total_volume += nbytes;

                    CommInfo& comm_info = recv_info[recv_from];
                    comm_info.cctcs.push_back(&cctc);
                    comm_info.imfs.push_back(imf);
                    comm_info.size += nbytes;
                }
            }

            the_recv_data = static_cast<char*>(amrex::The_FA_Arena()->alloc(total_volume));
            char* p = the_recv_data;
            for (auto it = recv_info.begin(); it != recv_info.end(); ++it)
            {
                const int recv_from = it->first;
                CommInfo& comm_info = it->second;
                comm_info.data = p;
                AMREX_ASSERT(comm_info.size < std::numeric_limits<int>::max());
                comm_info.req = ParallelDescriptor::Arecv(p, comm_info.size,
                                                          ParallelContext::global_to_local_rank(recv_from),
                                                          SeqNum, comm).req();
                p += comm_info.size;
            }
        }

        if (N_snds_tot > 0)
        {
            std::size_t total_volume = 0;
            for (int imf = 0; imf < nummfs; ++imf)
            {
                for (auto const& kv : *(TheFB[imf]->m_SndTags))
                {
                    const int send_to = kv.first;
                    auto const& cctc = kv.second;

                    std::size_t nbytes = 0;
                    for (auto const& cct : cctc) {
                        nbytes += (*mf[imf])[cct.srcIndex].nBytes(cct.sbox,scomp[imf],ncomp[imf]);
                    }
                    total_volume += nbytes;

                    CommInfo& comm_info = send_info[send_to];
                    comm_info.cctcs.push_back(&cctc);
                    comm_info.imfs.push_back(imf);
                    comm_info.size += nbytes;
                }
            }

            const int N_snds = send_info.size();

            Vector<int> send_to_ranks;
            send_to_ranks.reserve(send_info.size());

            the_send_data = static_cast<char*>(amrex::The_FA_Arena()->alloc(total_volume));
            char* ptmp = the_send_data;
            for (auto it = send_info.begin(); it != send_info.end(); ++it) {
                send_to_ranks.push_back(it->first);
                it->second.data = ptmp;
                ptmp += it->second.size;
            }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (Gpu::StreamIter sit(N_snds); sit.isValid(); ++sit)
            {
                const int send_to = send_to_ranks[sit()];
                CommInfo& comm_info = send_info[send_to];
                char* dptr = comm_info.data;
                auto const& cctcs = comm_info.cctcs;
                auto const& imfs = comm_info.imfs;
                for (int ifa = 0, nfa = imfs.size(); ifa < nfa; ++ifa)
                {
                    auto const& cctc = *cctcs[ifa];
                    const int imf = imfs[ifa];
                    auto& fabarray = *mf[imf];
                    const int sc = scomp[imf];
                    const int nc = ncomp[imf];
                    for (auto const& tag : cctc)
                    {
                        const Box& bx = tag.sbox;
                        auto const sfab = fabarray.array(tag.srcIndex);
                        auto pfab = amrex::makeArray4((value_type*)(dptr),bx);
                        AMREX_HOST_DEVICE_FOR_4D ( bx, nc, i, j, k, n,
                        {
                            pfab(i,j,k,n) = sfab(i,j,k,n+sc);
                        });
                        dptr += bx.numPts()*nc*sizeof(value_type);
                    }
                }
                AMREX_ASSERT(dptr == comm_info.data + comm_info.size);
            }

            for (auto it = send_info.begin(); it != send_info.end(); ++it)
            {
                const int send_to = it->first;
                CommInfo& comm_info = it->second;
                comm_info.req = ParallelDescriptor::Asend(comm_info.data, comm_info.size,
                                                          ParallelContext::global_to_local_rank(send_to),
                                                          SeqNum, comm).req();
            }
        }

        if (N_locs_tot > 0)
        {
            for (int imf = 0; imf < nummfs; ++imf)
            {
                auto& fa = *mf[imf];
                const int sc = scomp[imf];
                const int nc = ncomp[imf];
                const auto& ng = nghost[imf];
                const DistributionMapping& dm = fa.DistributionMap();
                const auto& tags = *(TheFB[imf]->m_LocTags);
                const int N_locs = tags.size();
                bool is_thread_safe = FAB::isCopyOMPSafe() && TheFB[imf]->m_threadsafe_loc;
                if (Gpu::inLaunchRegion() || !is_thread_safe)
                {
                    LayoutData<Vector<FabCopyTag<FAB> > > loc_copy_tags(fa.boxArray(),dm);
                    for (int i = 0; i < N_locs; ++i)
                    {
                        const auto& tag = tags[i];
                        if (dm[tag.dstIndex] == myproc)
                        {
                            loc_copy_tags[tag.dstIndex].push_back
                                ({fa.fabHostPtr(tag.srcIndex), tag.dbox, tag.sbox.smallEnd()-tag.dbox.smallEnd()});
                        }
                    }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                    {
                        int cuda_stream = 0;
                        for (MFIter mfi(fa); mfi.isValid(); ++mfi)
                        {
                            const auto& lc_tags = loc_copy_tags[mfi];
                            const auto& dfab = fa.array(mfi);
                            for (const auto& tag : lc_tags) {
                                const Dim3 offset = tag.offset.dim3();
                                const auto sfab = tag.sfab->array();
#ifdef AMREX_USE_GPU
                                Gpu::Device::setStreamIndex(cuda_stream++);
#endif
                                AMREX_HOST_DEVICE_FOR_4D ( tag.dbox, nc, i, j, k, n,
                                {
                                    dfab(i,j,k,sc+n) = sfab(i+offset.x,j+offset.y,k+offset.z,sc+n);
                                });
                            }
                        }
                    }
                }
                else
                {
#ifdef _OPENMP
#pragma omp parallel for
#endif
                    for (int i = 0; i < N_locs; ++i)
                    {
                        const auto& tag = tags[i];
                        if (myproc == dm[tag.dstIndex])
                        {
                            const FAB* sfab = &(fa[tag.srcIndex]);
                                  FAB* dfab = &(fa[tag.dstIndex]);
                            dfab->copy(*sfab, tag.sbox, sc, tag.dbox, sc, nc);
                        }
                    }
                }
            }
        }

        if (!recv_info.empty())
        {
            const int N_rcvs = recv_info.size();
            Vector<MPI_Request> recv_reqs;
            Vector<int> recv_size;
            Vector<CommInfo*> recv_info_v;
            for (auto it = recv_info.begin(); it != recv_info.end(); ++it)
            {
                CommInfo& comm_info = it->second;
                recv_reqs.push_back(comm_info.req);
                recv_size.push_back(static_cast<int>(comm_info.size));
                recv_info_v.push_back(&comm_info);
            }
            Vector<MPI_Status> recv_stat(N_rcvs);

            ParallelDescriptor::Waitall(recv_reqs, recv_stat);
#ifdef AMREX_DEBUG
            if (!FabArrayBase::CheckRcvStats(recv_stat, recv_size, MPI_CHAR, SeqNum))
            {
                amrex::Abort("amrex::FillBoundary failed with wrong message size");
            }
#endif

            bool is_thread_safe = FAB::isCopyOMPSafe();
            for (int imf = 0; imf < nummfs; ++imf) {
                is_thread_safe = is_thread_safe && TheFB[imf]->m_threadsafe_rcv;
            }
            if (Gpu::inLaunchRegion() || !is_thread_safe)
            {
                Vector<LayoutData<Vector<VoidCopyTag> > > recv_copy_tags_all(nummfs);
                for (int imf = 0; imf < nummfs; ++imf) {
                    recv_copy_tags_all[imf].define(mf[imf]->boxArray(),mf[imf]->DistributionMap());
                }
                for (int k = 0; k < N_rcvs; ++k)
                {
                    CommInfo& comm_info = *recv_info_v[k];
                    const char* dptr = comm_info.data;
                    auto const& cctcs = comm_info.cctcs;
                    auto const& imfs = comm_info.imfs;
                    for (int ifa = 0, nfa = imfs.size(); ifa < nfa; ++ifa)
                    {
                        auto const& cctc = *cctcs[ifa];
                        const int imf = imfs[ifa];
                        auto & recv_copy_tags = recv_copy_tags_all[imf];
                        const int sc = scomp[imf];
                        const int nc = ncomp[imf];
                        for (auto const& tag : cctc)
                        {
                            recv_copy_tags[tag.dstIndex].push_back({dptr,tag.dbox});
                            dptr += tag.dbox.numPts()*nc*sizeof(value_type);
                        }
                    }
                    AMREX_ASSERT(dptr == comm_info.data + comm_info.size);
                }
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (int imf = 0; imf < nummfs; ++imf)
                {
                    auto& fabarray = *mf[imf];
                    auto const& recv_copy_tags = recv_copy_tags_all[imf];
                    const int sc = scomp[imf];
                    const int nc = ncomp[imf];
                    for (MFIter mfi(fabarray); mfi.isValid(); ++mfi)
                    {
                        const auto& tags = recv_copy_tags[mfi];
                        auto dfab = fabarray.array(mfi);
                        for (auto const& tag : tags)
                        {
                            auto pfab = amrex::makeArray4((value_type*)(tag.p), tag.dbox);
                            AMREX_HOST_DEVICE_FOR_4D ( tag.dbox, nc, i, j, k, n,
                            {
                                dfab(i,j,k,n+sc) = pfab(i,j,k,n);
                            });
                        }
                    }
                }
            }
            else
            {
#ifdef _OPENMP
#pragma omp parallel for
#endif
                for (int k = 0; k < N_rcvs; ++k)
                {
                    CommInfo& comm_info = *recv_info_v[k];
                    const char* dptr = comm_info.data;
                    auto const& cctcs = comm_info.cctcs;
                    auto const& imfs = comm_info.imfs;
                    for (int ifa = 0, nfa = imfs.size(); ifa < nfa; ++ifa)
                    {
                        auto const& cctc = *cctcs[ifa];
                        const int imf = imfs[ifa];
                        auto& fabarray = *mf[imf];
                        const int sc = scomp[imf];
                        const int nc = ncomp[imf];
                        for (auto const& tag : cctc)
                        {
                            const Box& bx = tag.dbox;
                            FAB* dfab = &(fabarray[tag.dstIndex]);
                            dfab->copyFromMem(bx, sc, nc, dptr);
                            dptr += bx.numPts()*nc*sizeof(value_type);
                        }
                    }
                    AMREX_ASSERT(dptr == comm_info.data + comm_info.size);
                }
            }
        }

        if (!send_info.empty())
        {
            const int N_snds = send_info.size();
            Vector<MPI_Request> send_reqs;
            for (auto it = send_info.begin(); it != send_info.end(); ++it)
            {
                send_reqs.push_back(it->second.req);
            }
            Vector<MPI_Status> send_stat(N_snds);
            ParallelDescriptor::Waitall(send_reqs, send_stat);
        }

        amrex::The_FA_Arena()->free(the_send_data);
        amrex::The_FA_Arena()->free(the_recv_data);
#endif
    }
}
