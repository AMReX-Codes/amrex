
#include <AMReX_FBI.H>
#include <AMReX_PCI.H>

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
	int N_locs = (*TheFB.m_LocTags).size();
        if (N_locs == 0) return;
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10))
            if (Gpu::inGraphRegion())
            {
                FB_local_copy_cuda_graph_1(TheFB, scomp, ncomp);
            }
            else
#endif
            {
                FB_local_copy_gpu(TheFB, scomp, ncomp);
            }
        }
        else
#endif
        {
            FB_local_copy_cpu(TheFB, scomp, ncomp);
        }

        return;
    }

#ifdef BL_USE_MPI

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum = ParallelDescriptor::SeqNum();
    fb_tag = SeqNum;

    const int N_locs = TheFB.m_LocTags->size();
    const int N_rcvs = TheFB.m_RcvTags->size();
    const int N_snds = TheFB.m_SndTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0)
        // No work to do.
        return;

    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    fb_the_recv_data = nullptr;

    if (N_rcvs > 0) {
        PostRcvs(*TheFB.m_RcvTags, fb_the_recv_data,
                 fb_recv_data, fb_recv_size, fb_recv_from, fb_recv_reqs,
                 scomp, ncomp, SeqNum);
        fb_recv_stat.resize(N_rcvs);
    }

    //
    // Post send's
    //
    char*&                          the_send_data = fb_the_send_data;
    Vector<char*> &                     send_data = fb_send_data;
    Vector<int>                         send_size;
    Vector<int>                         send_rank;
    Vector<MPI_Request>&                send_reqs = fb_send_reqs;
    Vector<const CopyComTagsContainer*> send_cctc;

    if (N_snds > 0)
    {
        fb_send_data.clear();
        fb_send_reqs.clear();

	send_data.reserve(N_snds);
	send_size.reserve(N_snds);
	send_rank.reserve(N_snds);
        send_reqs.reserve(N_snds);
	send_cctc.reserve(N_snds);

        std::size_t total_volume = 0;
        for (auto const& kv : *TheFB.m_SndTags)
        {
            Vector<int> iss;                
            auto const& cctc = kv.second;

            std::size_t nbytes = 0;
            for (auto const& cct : kv.second)
            {
                nbytes += (*this)[cct.srcIndex].nBytes(cct.sbox,scomp,ncomp);
            }
            
            BL_ASSERT(nbytes < std::numeric_limits<int>::max());
            
            total_volume += nbytes;

            send_data.push_back(nullptr);
            send_size.push_back(static_cast<int>(nbytes));
            send_rank.push_back(kv.first);
            send_reqs.push_back(MPI_REQUEST_NULL);
            send_cctc.push_back(&cctc);
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

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10))
            if (Gpu::inGraphRegion()) {
                FB_pack_send_buffer_cuda_graph(TheFB, scomp, ncomp, send_data, send_size, send_cctc);
            }
            else
#endif
            {
                pack_send_buffer_gpu(*this, scomp, ncomp, send_data, send_size, send_cctc);
            }
        }
        else
#endif
        {
            pack_send_buffer_cpu(*this, scomp, ncomp, send_data, send_size, send_cctc);
        }

        for (int j = 0; j < N_snds; ++j)
        {
            if (send_size[j] > 0) {
                send_reqs[j] = ParallelDescriptor::Asend
                    (send_data[j], send_size[j],
                     ParallelContext::global_to_local_rank(send_rank[j]),
                     SeqNum,
                     ParallelContext::CommunicatorSub()).req();
            }
	}
    }

    FillBoundary_test();

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    if (N_locs > 0)
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )
            if (Gpu::inGraphRegion()) {
                FB_local_copy_cuda_graph_n(TheFB, scomp, ncomp);
            }
            else
#endif
            {
                FB_local_copy_gpu(TheFB, scomp, ncomp);
            }
        }
        else
#endif
        {
            FB_local_copy_cpu(TheFB, scomp, ncomp);
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

#ifdef AMREX_USE_MPI

    const FB& TheFB = getFB(fb_nghost,fb_period,fb_cross,fb_epo);
    const int N_rcvs = TheFB.m_RcvTags->size();
    if (N_rcvs > 0)
    {
        Vector<const CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);
        for (int k = 0; k < N_rcvs; k++) 
        {
            if (fb_recv_size[k] > 0)
            {
                auto const& cctc = TheFB.m_RcvTags->at(fb_recv_from[k]);
                recv_cctc[k] = &cctc;
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

        bool is_thread_safe = TheFB.m_threadsafe_rcv;

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
#if ( defined(__CUDACC__) && (__CUDACC_VER_MAJOR__ >= 10) )
            if (Gpu::inGraphRegion())
            {
                FB_unpack_recv_buffer_cuda_graph(TheFB, fb_scomp, fb_ncomp,
                                                 fb_recv_data, fb_recv_size,
                                                 recv_cctc, is_thread_safe);
            }
            else
#endif
            {
                unpack_recv_buffer_gpu(*this, fb_scomp, fb_ncomp, fb_recv_data, fb_recv_size,
                                       recv_cctc, FabArrayBase::COPY, is_thread_safe);
            }
        }
        else
#endif
        {
            unpack_recv_buffer_cpu(*this, fb_scomp, fb_ncomp, fb_recv_data, fb_recv_size,
                                   recv_cctc, FabArrayBase::COPY, is_thread_safe);
        }

        if (fb_the_recv_data)
        {
            amrex::The_FA_Arena()->free(fb_the_recv_data);
            fb_the_recv_data = nullptr;
        }
    }

    const int N_snds = TheFB.m_SndTags->size();
    if (N_snds > 0) {
        Vector<MPI_Status> stats;
        FabArrayBase::WaitForAsyncSends(N_snds,fb_send_reqs,fb_send_data,stats);
        amrex::The_FA_Arena()->free(fb_the_send_data);
        fb_the_send_data = nullptr;
    }
#endif
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter fai(*this,TilingIfNotGPU()); fai.isValid(); ++fai)
        {
            const Box& bx = fai.tilebox();

            // avoid self copy or plus
	    if (this != &src) {
                auto const sfab = src.array(fai);
                auto       dfab = this->array(fai);
		if (op == FabArrayBase::COPY) {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
                    {
                        dfab(i,j,k,dcomp+n) = sfab(i,j,k,scomp+n);
                    });
		} else {
                    AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
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
	int N_locs = (*thecpc.m_LocTags).size();
        if (N_locs == 0) return;
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            PC_local_gpu(thecpc, src, scomp, dcomp, ncomp, op);
        }
        else
#endif
        {
            PC_local_cpu(thecpc, src, scomp, dcomp, ncomp, op);
        }

        return;
    }

#ifdef BL_USE_MPI

    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
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
                     recv_data, recv_size, recv_from, recv_reqs, SC, NC, SeqNum);
            actual_n_rcvs = N_rcvs - std::count(recv_size.begin(), recv_size.end(), 0);
	}

	//
	// Post send's
	//
        char*                               the_send_data;
	Vector<char*>                       send_data;
	Vector<int>                         send_size;
	Vector<int>                         send_rank;
	Vector<MPI_Request>                 send_reqs;
	Vector<const CopyComTagsContainer*> send_cctc;

	if (N_snds > 0)
	{
	    send_data.reserve(N_snds);
	    send_size.reserve(N_snds);
	    send_rank.reserve(N_snds);
            send_reqs.reserve(N_snds);
	    send_cctc.reserve(N_snds);

            std::size_t total_volume = 0;
            for (auto const& kv : *thecpc.m_SndTags)
	    {
                Vector<int> iss;                
                auto const& cctc = kv.second;

                std::size_t nbytes = 0;
                for (auto const& cct : kv.second)
                {
                    nbytes += src[cct.srcIndex].nBytes(cct.sbox,SC,NC);
                }
		
		BL_ASSERT(nbytes < std::numeric_limits<int>::max());

                total_volume += nbytes;

		send_data.push_back(nullptr);
                send_size.push_back(static_cast<int>(nbytes));
                send_rank.push_back(kv.first);
                send_reqs.push_back(MPI_REQUEST_NULL);
                send_cctc.push_back(&cctc);
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

#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion())
            {
                pack_send_buffer_gpu(src, SC, NC, send_data, send_size, send_cctc);
            }
            else
#endif
            {
                pack_send_buffer_cpu(src, SC, NC, send_data, send_size, send_cctc);
            }

            for (int j = 0; j < N_snds; ++j)
            {
                if (send_size[j] > 0) {
                    send_reqs[j] = ParallelDescriptor::Asend
                        (send_data[j], send_size[j],
                         ParallelContext::global_to_local_rank(send_rank[j]),
                         SeqNum,
                         ParallelContext::CommunicatorSub()).req();
                }
	    }
	}

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
        if (N_locs > 0)
	{
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion())
            {
                PC_local_gpu(thecpc, src, SC, DC, NC, op);
            }
            else
#endif
            {
                PC_local_cpu(thecpc, src, SC, DC, NC, op);
            }
        }

        if (N_rcvs > 0)
        {
            Vector<const CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);
	    for (int k = 0; k < N_rcvs; ++k)
	    {
                if (recv_size[k] > 0)
                {
                    auto const& cctc = thecpc.m_RcvTags->at(recv_from[k]);
                    recv_cctc[k] = &cctc;
                }
	    }

            if (actual_n_rcvs > 0) {
                Vector<MPI_Status> stats(N_rcvs);
                ParallelDescriptor::Waitall(recv_reqs, stats);
#ifdef AMREX_DEBUG
                if (!CheckRcvStats(stats, recv_size, MPI_CHAR, SeqNum))
                {
                    amrex::Abort("ParallelCopy failed with wrong message size");
                }
#endif
            }

            bool is_thread_safe = thecpc.m_threadsafe_rcv;

#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion())
            {
                unpack_recv_buffer_gpu(*this, DC, NC, recv_data, recv_size, recv_cctc,
                                       op, is_thread_safe);
            }
            else
#endif
            {
                unpack_recv_buffer_cpu(*this, DC, NC, recv_data, recv_size, recv_cctc,
                                       op, is_thread_safe);
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
                         int                               SeqNum)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();
    recv_reqs.clear();

    std::size_t TotalRcvsVolume = 0;
    for (const auto& kv : m_RcvTags) // loop over senders
    {
        std::size_t nbytes = 0;
        for (auto const& cct : kv.second)
        {
            nbytes += (*this)[cct.dstIndex].nBytes(cct.dbox,icomp,ncomp);
        }

        BL_ASSERT(nbytes < std::numeric_limits<int>::max());

        TotalRcvsVolume += nbytes;

        recv_data.push_back(nullptr);
        recv_size.push_back(static_cast<int>(nbytes));
        recv_from.push_back(kv.first);
        recv_reqs.push_back(MPI_REQUEST_NULL);
    }

    const int nrecv = recv_from.size();
    
    MPI_Comm comm = ParallelContext::CommunicatorSub();

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
        int i = recv_counter;

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
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D ( bx, ncomp, i, j, k, n,
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

template <class FAB>
void
FillBoundary (Vector<FabArray<FAB>*> const& mf, const Periodicity& period)
{
    BL_PROFILE("FillBoundary(Vector)");
    const int nummfs = mf.size();
#if 1
    for (int imf = 0; imf < nummfs; ++imf) {
        mf[imf]->FillBoundary(period);
    }
#else
    for (int imf = 0; imf < nummfs; ++imf) {
        mf[imf]->FillBoundary_nowait(period);
    }
    for (int imf = 0; imf < nummfs; ++imf) {
        mf[imf]->FillBoundary_finish();
    }
#endif
}
