

template <class FAB>
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
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && TheFB.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*TheFB.m_LocTags)[i];

            BL_ASSERT(distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(distributionMap[tag.srcIndex] == ParallelDescriptor::MyProc());

            get(tag.dstIndex).copy(get(tag.srcIndex),tag.sbox,scomp,tag.dbox,scomp,ncomp);
        }

        return;
    }
   
#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
    ParallelDescriptor::Mode.incr_upcxx();
#endif

#if defined(BL_USE_MPI3)
    BL_ASSERT(FAB::preAllocatable());
#else
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

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
    Vector<char*> &                     send_data = fb_send_data;
    Vector<int>                         send_size;
    Vector<int>                         send_rank;
    Vector<MPI_Request>&                send_reqs = fb_send_reqs;
    Vector<const CopyComTagsContainer*> send_cctc;
    Vector<Vector<int> >                indv_send_size;
    Vector<MPI_Request>                 pre_reqs;

    fb_tag = SeqNum;

#if defined BL_USE_UPCXX || defined BL_USE_MPI3 
    int actual_n_snds = 0;
#endif
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

        for (auto const& kv : *TheFB.m_SndVols)
        {
            Vector<int> iss;                
            auto const& cctc = TheFB.m_SndTags->at(kv.first);

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
            
            char* data = nullptr;
            if (nbytes > 0)
            {
                data = static_cast<char*>
#ifdef BL_USE_UPCXX
                    (BLPgas::alloc(nbytes));
#else
                    (amrex::The_Arena()->alloc(nbytes));
#endif
            }
                    
            send_data.push_back(data);
            send_size.push_back(static_cast<int>(nbytes));
            send_rank.push_back(kv.first);
            send_reqs.push_back(MPI_REQUEST_NULL);
            send_cctc.push_back(&cctc);
            indv_send_size.push_back(std::move(iss));
        }

#if defined BL_USE_UPCXX || defined BL_USE_MPI3         
        actual_n_snds = N_snds - std::count(send_size.begin(), send_size.end(), 0);
#endif
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
#ifdef BL_USE_MPI3
    MPI_Group tgroup, rgroup, sgroup;
    if (ParallelDescriptor::MPIOneSided()) {
	MPI_Comm_group(ParallelDescriptor::Communicator(), &tgroup);
    }
#endif

    fb_the_recv_data = nullptr;

    if (N_rcvs > 0) {
#ifdef BL_USE_UPCXX
	PostRcvs_PGAS(*TheFB.m_RcvVols, fb_the_recv_data, fb_recv_data,
                      fb_recv_size, fb_recv_from,
                      scomp, ncomp, SeqNum, &BLPgas::fb_recv_event);
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    PostRcvs_MPI_Onesided(*TheFB.m_RcvVols, fb_the_recv_data, fb_recv_data,
                                  fb_recv_size, fb_recv_from, fb_recv_reqs, fb_recv_disp,
                                  scomp, ncomp, SeqNum, ParallelDescriptor::fb_win);
	    MPI_Group_incl(tgroup, fb_recv_from.size(), fb_recv_from.dataPtr(), &rgroup);
	    MPI_Win_post(rgroup, 0, ParallelDescriptor::fb_win);
#endif
	} else {
	    PostRcvs(*TheFB.m_RcvVols, *TheFB.m_RcvTags,
                     fb_recv_data, fb_recv_size, fb_recv_from, fb_recv_reqs,
                     scomp, ncomp, SeqNum, preSeqNum);
	}
#endif
    }

    //
    // Post send's
    //
    if (N_snds > 0)
    {
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe())
#endif
	for (int j=0; j<N_snds; ++j)
	{
            char* dptr = send_data[j];
            if (dptr != nullptr)
            {
                auto const& cctc = *send_cctc[j];
                for (auto const& tag : cctc)
                {
                    const Box& bx = tag.sbox;
                    auto n = (*this)[tag.srcIndex].copyToMem(bx,scomp,ncomp,dptr);
                    dptr += n;
                }
                BL_ASSERT(dptr == send_data[j] + send_size[j]);
            }
	}

#ifdef BL_USE_UPCXX

	BLPgas::fb_send_counter = 0;

	for (int i=0; i<N_snds; ++i) {
            if (send_size[i] > 0) {
                BLPgas::Send(upcxx::global_ptr<void>((void *)send_data[i], upcxx::myrank()),
                             send_rank[i], send_size[i], SeqNum,
                             &BLPgas::fb_send_event, &BLPgas::fb_send_counter);
            }
        }

	// Need to make sure at least half of the sends have been started
	while (BLPgas::fb_send_counter < actual_n_snds) {
	    upcxx::advance();
        }

#else  // MPI

	if (ParallelDescriptor::MPIOneSided())
	{
#if defined(BL_USE_MPI3)
	    Vector<MPI_Aint> send_disp(N_snds,0);
            send_reqs.clear();
	    
	    for (int i=0; i<N_snds; ++i) {
                if (send_size[i] > 0) {
                    send_reqs.push_back(ParallelDescriptor::Arecv
                                        (&send_disp[i],1,send_rank[i],SeqNum).req());
                }
	    }
		
	    MPI_Group_incl(tgroup, N_snds, send_rank.dataPtr(), &sgroup);
	    MPI_Win_start(sgroup,0,ParallelDescriptor::fb_win);
	    
	    int send_counter = 0;	
	    while (send_counter < actual_n_snds) {
		MPI_Status status;
		int index;
	
                ParallelDescriptor::Waitany(send_reqs, index, status);	
		
		BL_ASSERT(status.MPI_TAG == SeqNum);
		BL_ASSERT(status.MPI_SOURCE == send_rank[index]);

                MPI_Put(send_data[index], send_size[index], MPI_CHAR, send_rank[index],
                        send_disp[index], send_size[index], MPI_CHAR, ParallelDescriptor::fb_win);
		
		++send_counter;
	    }
#endif // BL_USE_MPI3
	} 
	else 
	{
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
#endif // #ifdef BL_USE_UPCXX #else 
    }

    if (ParallelDescriptor::MPIOneSided()) {
#ifdef BL_USE_MPI3
	MPI_Group_free(&tgroup);
	if (N_rcvs > 0) MPI_Group_free(&rgroup);
	if (N_snds > 0) MPI_Group_free(&sgroup);
#endif
    }

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    if (ParallelDescriptor::TeamSize() > 1 && TheFB.m_threadsafe_loc)
    {
#ifdef BL_USE_TEAM
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
	ParallelDescriptor::team_for(0, N_locs, [&] (int i) 
        {
	    const auto& tag = (*TheFB.m_LocTags)[i];

	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.dstIndex]));
	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.srcIndex]));

	    get(tag.dstIndex).copy(get(tag.srcIndex),tag.sbox,scomp,tag.dbox,scomp,ncomp);
	});
#endif
    }
    else
    {
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && TheFB.m_threadsafe_loc)
#endif
	for (int i=0; i<N_locs; ++i)
	{
	    const CopyComTag& tag = (*TheFB.m_LocTags)[i];

	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.dstIndex]));
	    BL_ASSERT(ParallelDescriptor::sameTeam(distributionMap[tag.srcIndex]));
	    
	    if (distributionMap[tag.dstIndex] == ParallelDescriptor::MyProc()) {
		get(tag.dstIndex).copy(get(tag.srcIndex),tag.sbox,scomp,tag.dbox,scomp,ncomp);
	    }
	}
    }
#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::FillBoundary_finish ()
{
    if ( n_grow.allLE(IntVect::TheZeroVector()) && !fb_epo ) return; // For epo (Enforce Periodicity Only), there may be no ghost cells.

    if (ParallelContext::NProcsSub() == 1) return;

#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
    ParallelDescriptor::Mode.decr_upcxx();
#endif

#if !defined(BL_USE_MPI3)
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

    const FB& TheFB = getFB(fb_nghost,fb_period,fb_cross,fb_epo);

    const int N_rcvs = TheFB.m_RcvTags->size();
    const int N_snds = TheFB.m_SndTags->size();

    int actual_n_rcvs = N_rcvs - std::count(fb_recv_data.begin(), fb_recv_data.end(), nullptr);

#ifdef BL_USE_UPCXX
    if (actual_n_rcvs > 0) BLPgas::fb_recv_event.wait();
#else 
    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	if (N_snds > 0) MPI_Win_complete(ParallelDescriptor::fb_win);
	if (N_rcvs > 0) MPI_Win_wait    (ParallelDescriptor::fb_win);
#endif
    } else {
	if (actual_n_rcvs > 0) {
	    Vector<MPI_Status> stats(N_rcvs);
            ParallelDescriptor::Waitall(fb_recv_reqs, stats);
	    if (!CheckRcvStats(stats, fb_recv_size, MPI_CHAR, fb_tag))
            {
                amrex::Abort("FillBoundary_finish failed with wrong message size");
            }
	}
    }
#endif

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

#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && TheFB.m_threadsafe_rcv)
#endif

	for (int k = 0; k < N_rcvs; k++) 
	{
	    const char* dptr = fb_recv_data[k];
            if (dptr != nullptr)
            {
                auto const& cctc = *recv_cctc[k];
                for (auto const& tag : cctc)
                {
                    const Box& bx  = tag.dbox;
                    std::size_t n = (*this)[tag.dstIndex].copyFromMem(bx,fb_scomp,fb_ncomp,dptr);
                    dptr += n;
                }

                BL_ASSERT(dptr == fb_recv_data[k] + fb_recv_size[k]);
            }
	}

        if (fb_the_recv_data)
        {
#ifdef BL_USE_UPCXX
            BLPgas::free(fb_the_recv_data);
#else
            if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
                MPI_Win_detach(ParallelDescriptor::fb_win, fb_the_recv_data);
#endif
            }
	    amrex::The_Arena()->free(fb_the_recv_data);
#endif
	}
        else
        {
            for (auto p : fb_recv_data) {
                amrex::The_Arena()->free(p);
            }
        }        
    }

    if (N_snds > 0) {
#ifdef BL_USE_UPCXX
        FabArrayBase::WaitForAsyncSends_PGAS(N_snds,fb_send_data,
					     &BLPgas::fb_send_event,
					     &BLPgas::fb_send_counter);
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    for (int i = 0; i < N_snds; ++i) {
		if (fb_send_data[i]) amrex::The_Arena()->free(fb_send_data[i]);
            }
#endif
        } else {
            Vector<MPI_Status> stats;
            FabArrayBase::WaitForAsyncSends(N_snds,fb_send_reqs,fb_send_data,stats);
        }
#endif
    }

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

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
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
        for (MFIter fai(*this,true); fai.isValid(); ++fai)
        {
            const Box& bx = fai.tilebox();

	    if (this != &src) {
		// avoid self copy or plus
		if (op == FabArrayBase::COPY) {
		    get(fai).copy(src[fai],bx,scomp,bx,dcomp,ncomp);
		} else {
		    get(fai).plus(src[fai],bx,bx,scomp,dcomp,ncomp);
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
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];

	    if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		// avoid self copy or plus
		if (op == FabArrayBase::COPY) {
		    get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,scomp,tag.dbox,dcomp,ncomp);
		} else {
		    get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,scomp,dcomp,ncomp);
		}
	    }
        }

        return;
    }

#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
#endif

#if defined(BL_USE_MPI3)
    BL_ASSERT(FAB::preAllocatable());
#else
    BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

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

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0)
        //
        // No work to do.
        //
        return;

#ifdef BL_USE_MPI3
    MPI_Group tgroup, rgroup, sgroup;
    if (ParallelDescriptor::MPIOneSided()) {
	MPI_Comm_group(ParallelDescriptor::Communicator(), &tgroup);
    }
#endif

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

	Vector<char*>                       send_data;
	Vector<int>                         send_size;
	Vector<int>                         send_rank;
	Vector<MPI_Request>                 send_reqs;
	Vector<const CopyComTagsContainer*> send_cctc;
        Vector<Vector<int> >                 indv_send_size;
        Vector<MPI_Request>                 pre_reqs;

#if defined BL_USE_UPCXX || defined BL_USE_MPI3 
        int actual_n_snds = 0;
#endif
	if (N_snds > 0)
	{
	    send_data.reserve(N_snds);
	    send_size.reserve(N_snds);
	    send_rank.reserve(N_snds);
            send_reqs.reserve(N_snds);
	    send_cctc.reserve(N_snds);
            indv_send_size.reserve(N_snds);

            for (auto const& kv : *thecpc.m_SndVols)
	    {
                Vector<int> iss;                
                auto const& cctc = thecpc.m_SndTags->at(kv.first);

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

                char* data = nullptr;
                if (nbytes > 0)
                {
                    data = static_cast<char*>
#ifdef BL_USE_UPCXX
                        (BLPgas::alloc(nbytes));
#else
                        (amrex::The_Arena()->alloc(nbytes));
#endif
                }
                    
		send_data.push_back(data);
                send_size.push_back(static_cast<int>(nbytes));
                send_rank.push_back(kv.first);
                send_reqs.push_back(MPI_REQUEST_NULL);
                send_cctc.push_back(&cctc);
                indv_send_size.push_back(std::move(iss));
	    }

#if defined BL_USE_UPCXX || defined BL_USE_MPI3
            actual_n_snds = N_snds - std::count(send_size.begin(), send_size.end(), 0);
#endif
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
#ifdef BL_USE_MPI3
	Vector<MPI_Aint>    recv_disp;
#endif
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        char* the_recv_data = nullptr;

        int actual_n_rcvs = 0;
	if (N_rcvs > 0) {
#ifdef BL_USE_UPCXX
	    PostRcvs_PGAS(*thecpc.m_RcvVols, the_recv_data, recv_data,
                          recv_size, recv_from, 
                          SC, NC, SeqNum, &BLPgas::cp_recv_event);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
                PostRcvs_MPI_Onesided(*thecpc.m_RcvVols, the_recv_data, recv_data, 
                                      recv_size, recv_from, recv_reqs, recv_disp,
                                      SC, NC, SeqNum, ParallelDescriptor::cp_win);
		MPI_Group_incl(tgroup, recv_from.size(), recv_from.dataPtr(), &rgroup);
		MPI_Win_post(rgroup, 0, ParallelDescriptor::cp_win);
#endif
	    } else {
                PostRcvs(*thecpc.m_RcvVols, *thecpc.m_RcvTags,
                         recv_data, recv_size, recv_from, recv_reqs, SC, NC, SeqNum, preSeqNum);
	    }
#endif
            actual_n_rcvs = N_rcvs - std::count(recv_size.begin(), recv_size.end(), 0);
	}

	//
	// Post send's
	// 
	if (N_snds > 0)
	{
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe())
#endif
	    for (int j=0; j<N_snds; ++j)
	    {
		char* dptr = send_data[j];
                if (dptr != nullptr)
                {
                    auto const& cctc = *send_cctc[j];
                    for (auto const& tag : cctc)
                    {
                        const Box& bx = tag.sbox;
                        auto n = src[tag.srcIndex].copyToMem(bx,SC,NC,dptr);
                        dptr += n;
                    }
                    BL_ASSERT(dptr == send_data[j] + send_size[j]);
		}
	    }

#ifdef BL_USE_UPCXX
	    
	    BLPgas::cp_send_counter = 0;

	    for (int i=0; i<N_snds; ++i) {
                if (send_size[i] > 0) {
                    BLPgas::Send(upcxx::global_ptr<void>((void *)send_data[i], upcxx::myrank()),
                                 send_rank[i], send_size[i], SeqNum,
                                 &BLPgas::cp_send_event, &BLPgas::cp_send_counter);
                }
	    }
	    
	    // Need to make sure at least half of the sends have been started
	    while (BLPgas::cp_send_counter < actual_n_snds) {
		upcxx::advance();
            }

#else
	    
	    if (ParallelDescriptor::MPIOneSided())
	    {
#if defined(BL_USE_MPI3)
		Vector<MPI_Aint> send_disp(N_snds,0);
                send_reqs.clear();

		for (int i=0; i<N_snds; ++i) {
                    if (send_size[i] > 0) {
                        send_reqs[i] = ParallelDescriptor::Arecv
                            (&send_disp[i],1,send_rank[i],SeqNum).req();
                    }
		}
		
		MPI_Group_incl(tgroup, N_snds, send_rank.data(), &sgroup);
		MPI_Win_start(sgroup,0,ParallelDescriptor::cp_win);
		
		int send_counter = 0;	
		while (send_counter < actual_n_snds) {
		    MPI_Status status;
		    int index;
		   
                    ParallelDescriptor::Waitany(send_reqs, index, status); 
		    
		    BL_ASSERT(status.MPI_TAG == SeqNum);
		    BL_ASSERT(status.MPI_SOURCE == send_rank[index]);
		    
		    MPI_Put(send_data[index], send_size[index], MPI_CHAR, send_rank[index],
			    send_disp[index], send_size[index], MPI_CHAR, ParallelDescriptor::cp_win);
		    
		    ++send_counter;
		}
#endif
	    }
            else
            {
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
                        if (FabArrayBase::do_async_sends)
                        {
                            send_reqs[j] = ParallelDescriptor::Asend
                                (send_data[j], send_size[j],
                                 ParallelContext::global_to_local_rank(send_rank[j]),
                                 SeqNum,
                                 ParallelContext::CommunicatorSub()).req();
                        }
                        else
                        {
                            ParallelDescriptor::Send
                                (send_data[j], send_size[j],
                                 ParallelContext::global_to_local_rank(send_rank[j]),
                                 SeqNum,
                                 ParallelContext::CommunicatorSub());
                            amrex::The_Arena()->free(send_data[j]);
                        }
                    }

                    ++send_counter;
		}
	    }
#endif
	}

#ifdef BL_USE_MPI3
	if (ParallelDescriptor::MPIOneSided()) {
	    if (N_rcvs > 0) MPI_Group_free(&rgroup);
	    if (N_snds > 0) MPI_Group_free(&sgroup);
	}
#endif

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
	if (ParallelDescriptor::TeamSize() > 1 && thecpc.m_threadsafe_loc)
	{
#ifdef BL_USE_TEAM
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
	    ParallelDescriptor::team_for(0, N_locs, [&] (int j) 
            {
		const CopyComTag& tag = (*thecpc.m_LocTags)[j];
		
		if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		    // avoid self copy or plus
		    if (op == FabArrayBase::COPY) {
			get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,SC,tag.dbox,DC,NC);
		    } else {
			get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,SC,DC,NC);
		    }
		}
	    });
#endif	    
	}
	else 
	{
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc)
#endif
	    for (int j=0; j<N_locs; ++j)
	    {
		const CopyComTag& tag = (*thecpc.m_LocTags)[j];

		if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		    // avoid self copy or plus
		    if (op == FabArrayBase::COPY) {
			get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,SC,tag.dbox,DC,NC);
		    } else {
			get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,SC,DC,NC);
		    }
		}
	    }
	}

	//
	//  wait and unpack
	//
#ifdef BL_USE_UPCXX
        if (actual_n_rcvs > 0) BLPgas::cp_recv_event.wait();
#else
	if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
	    if (N_snds > 0) MPI_Win_complete(ParallelDescriptor::cp_win);
	    if (N_rcvs > 0) MPI_Win_wait    (ParallelDescriptor::cp_win);
#endif
	} else {
	    if (actual_n_rcvs > 0) {
		Vector<MPI_Status> stats(N_rcvs);
                ParallelDescriptor::Waitall(recv_reqs, stats);
		if (!CheckRcvStats(stats, recv_size, MPI_CHAR, SeqNum))
                {
                    amrex::Abort("ParallelCopy failed with wrong message size");
                }
	    }
	}
#endif	

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
	    
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_rcv)
#endif
	    {
		FAB fab;

#ifdef _OPENMP
#pragma omp for
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
                            std::size_t n;
                            if (op == FabArrayBase::COPY)
                            {
                                n = get(tag.dstIndex).copyFromMem(bx,DC,NC,dptr);
                            }
                            else
                            {
                                fab.resize(bx,NC);
                                n = fab.copyFromMem(bx,0,NC,dptr);
                                get(tag.dstIndex).plus(fab,bx,bx,0,DC,NC);
                            }
                            dptr += n;
                        }
                        BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
		    }
		}
	    }

            if (the_recv_data)
            {
#ifdef BL_USE_UPCXX
                BLPgas::free(the_recv_data);
#else
                if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
                    MPI_Win_detach(ParallelDescriptor::cp_win, the_recv_data);
#endif
                }
                amrex::The_Arena()->free(the_recv_data);
#endif
            }
            else
            {
                for (auto p : recv_data) {
                    amrex::The_Arena()->free(p);
                }
            }
	}
	
        if (N_snds > 0) {
#ifdef  BL_USE_UPCXX
	    FabArrayBase::WaitForAsyncSends_PGAS(N_snds,send_data,
					         &BLPgas::cp_send_event,
					         &BLPgas::cp_send_counter);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		for (int i = 0; i < N_snds; ++i) {
		    if (send_data[i]) amrex::The_Arena()->free(send_data[i]);
                }
#endif
	    } else {
		if (FabArrayBase::do_async_sends && ! thecpc.m_SndTags->empty()) {
		    Vector<MPI_Status> stats;
		    FabArrayBase::WaitForAsyncSends(N_snds,send_reqs,send_data,stats);
		}
	    }
#endif
        }

        ipass     += NC;
        SC        += NC;
        DC        += NC;
        NCompLeft -= NC;
    }

#ifdef BL_USE_MPI3
    if (ParallelDescriptor::MPIOneSided()) {
	MPI_Group_free(&tgroup);
    }
#endif

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

    return;

#endif /*BL_USE_MPI*/
}

template <class FAB>
typename FabArray<FAB>::CopierHandle
FabArray<FAB>::ParallelCopy_nowait (const FabArray<FAB>& src,
                                    int                  scomp,
                                    int                  dcomp,
                                    int                  ncomp,
                                    const IntVect&       snghost,
                                    const IntVect&       dnghost,
                                    const Periodicity&   period,
                                    CpOp                 op,
                                    const FabArrayBase::CPC * a_cpc)
{
    BL_PROFILE("FabArray::ParallelCopy_nowait()");

    if (size() == 0 || src.size() == 0) return CopierHandle();

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
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
        for (MFIter fai(*this,true); fai.isValid(); ++fai)
        {
            const Box& bx = fai.tilebox();

	    if (this != &src) {
		// avoid self copy or plus
		if (op == FabArrayBase::COPY) {
		    get(fai).copy(src[fai],bx,scomp,bx,dcomp,ncomp);
		} else {
		    get(fai).plus(src[fai],bx,bx,scomp,dcomp,ncomp);
		}
	    }
        }

        return CopierHandle();
    }

    const CPC& thecpc = (a_cpc) ? *a_cpc : getCPC(dnghost, src, snghost, period);

    if (ParallelContext::NProcsSub() == 1)
    {
        //
        // There can only be local work to do.
        //
	int N_loc = (*thecpc.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[i];

	    if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		// avoid self copy or plus
		if (op == FabArrayBase::COPY) {
		    get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,scomp,tag.dbox,dcomp,ncomp);
		} else {
		    get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,scomp,dcomp,ncomp);
		}
	    }
        }

        return CopierHandle();
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
    const int SeqNum  = ParallelDescriptor::SeqNum();

    const int N_snds = thecpc.m_SndTags->size();
    const int N_rcvs = thecpc.m_RcvTags->size();
    const int N_locs = thecpc.m_LocTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0) {
        // No work to do.
        return CopierHandle();
    }

    CopierHandle return_handle(*this, thecpc);
    CopierHandleImpl& handle = *(return_handle.m_impl);

    handle.op = op;
    handle.dcomp = dcomp;
    handle.ncomp = ncomp;
    handle.N_snds = N_snds;
    handle.N_rcvs = N_rcvs;
    handle.SeqNum = SeqNum;

    Vector<int>                         send_size;
    Vector<int>                         send_rank;
    Vector<const CopyComTagsContainer*> send_cctc;
    Vector<Vector<int> >                indv_send_size;
    Vector<MPI_Request>                 pre_reqs;

    Vector<MPI_Request>& send_reqs = handle.send_reqs;
    Vector<char*>      & send_data = handle.send_data;
    Vector<int>        & recv_from = handle.recv_from;
    Vector<char*>      & recv_data = handle.recv_data;
    Vector<int>        & recv_size = handle.recv_size;
    Vector<MPI_Request>& recv_reqs = handle.recv_reqs;

    //
    // Before we post recv, let's preprocess sends in case FAB is not preAllocatable
    //
    if (N_snds > 0)
    {
        send_data.reserve(N_snds);
        send_size.reserve(N_snds);
        send_rank.reserve(N_snds);
        send_reqs.reserve(N_snds);
        send_cctc.reserve(N_snds);
        indv_send_size.reserve(N_snds);
        
        for (auto const& kv : *thecpc.m_SndVols)
        {
            Vector<int> iss;                
            auto const& cctc = thecpc.m_SndTags->at(kv.first);
            
            std::size_t nbytes = 0;
            if (FAB::preAllocatable())
            {
                for (auto const& cct : kv.second)
                {
                    nbytes += src[cct.srcIndex].nBytes(cct.sbox,scomp,ncomp);
                }
            }
            else
            {
                for (auto const& tag : cctc)
                {
                    std::size_t b = src[tag.srcIndex].nBytes(tag.sbox,scomp,ncomp);
                    nbytes += b;
                    iss.push_back(static_cast<int>(b));
                }
            }
		
            BL_ASSERT(nbytes < std::numeric_limits<int>::max());

            char* data = nullptr;
            if (nbytes > 0)
            {
                data = static_cast<char*>(amrex::The_Arena()->alloc(nbytes));
            }
                    
            send_data.push_back(data);
            send_size.push_back(static_cast<int>(nbytes));
            send_rank.push_back(kv.first);
            send_reqs.push_back(MPI_REQUEST_NULL);
            send_cctc.push_back(&cctc);
            indv_send_size.push_back(std::move(iss));
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
    // Post rcvs
    //
    handle.actual_n_rcvs = 0;
    if (N_rcvs > 0) {
        PostRcvs(*thecpc.m_RcvVols, *thecpc.m_RcvTags,
                 recv_data, recv_size, recv_from, recv_reqs, scomp, ncomp, SeqNum, preSeqNum);
        handle.actual_n_rcvs = N_rcvs - std::count(recv_size.begin(), recv_size.end(), 0);
    }

    //
    // Post send's
    // 
    if (N_snds > 0)
    {
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe())
#endif
        for (int j=0; j<N_snds; ++j)
        {
            char* dptr = send_data[j];
            if (dptr != nullptr)
            {
                auto const& cctc = *send_cctc[j];
                for (auto const& tag : cctc)
                {
                    const Box& bx = tag.sbox;
                    auto n = src[tag.srcIndex].copyToMem(bx,scomp,ncomp,dptr);
                    dptr += n;
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
                if (FabArrayBase::do_async_sends)
                {
                    send_reqs[j] = ParallelDescriptor::Asend
                        (send_data[j], send_size[j],
                         ParallelContext::global_to_local_rank(send_rank[j]),
                         SeqNum,
                         ParallelContext::CommunicatorSub()).req();
                }
                else
                {
                    ParallelDescriptor::Send
                        (send_data[j], send_size[j],
                         ParallelContext::global_to_local_rank(send_rank[j]),
                         SeqNum,
                         ParallelContext::CommunicatorSub());
                    amrex::The_Arena()->free(send_data[j]);
                }
            }
            
            ++send_counter;
        }
    }

    //
    // Do the local work.  Hope for a bit of communication/computation overlap.
    //
    if (ParallelDescriptor::TeamSize() > 1 && thecpc.m_threadsafe_loc)
    {
#ifdef BL_USE_TEAM
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe())
#endif
        ParallelDescriptor::team_for(0, N_locs, [&] (int j) 
            {
		const CopyComTag& tag = (*thecpc.m_LocTags)[j];
		
		if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
		    // avoid self copy or plus
		    if (op == FabArrayBase::COPY) {
			get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,scomp,tag.dbox,dcomp,ncomp);
		    } else {
			get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,scomp,dcomp,ncomp);
		    }
		}
	    });
#endif	    
    }
    else 
    {
#ifdef _OPENMP
#pragma omp parallel for if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_loc)
#endif
        for (int j=0; j<N_locs; ++j)
        {
            const CopyComTag& tag = (*thecpc.m_LocTags)[j];
            
            if (this != &src || tag.dstIndex != tag.srcIndex || tag.sbox != tag.dbox) {
                // avoid self copy or plus
                if (op == FabArrayBase::COPY) {
                    get(tag.dstIndex).copy(src[tag.srcIndex],tag.sbox,scomp,tag.dbox,dcomp,ncomp);
                } else {
                    get(tag.dstIndex).plus(src[tag.srcIndex],tag.sbox,tag.dbox,scomp,dcomp,ncomp);
                }
            }
        }
    }

    return return_handle;

#endif /*BL_USE_MPI*/
}

template <class FAB>
void
FabArray<FAB>::CopierHandleImpl::finish ()
{
#ifdef BL_USE_MPI
    BL_PROFILE("FabArray::ParallelCopy_finish()");

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
        Vector<const CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);
        
        for (int k = 0; k < N_rcvs; ++k)
        {
            if (recv_size[k] > 0)
            {
                auto const& cctc = thecpc.m_RcvTags->at(recv_from[k]);
                recv_cctc[k] = &cctc;
            }
        }
	
#ifdef _OPENMP
#pragma omp parallel if (FAB::isCopyOMPSafe() && thecpc.m_threadsafe_rcv)
#endif
        {
            FAB fab;
            
#ifdef _OPENMP
#pragma omp for
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
                        std::size_t n;
                        if (op == FabArrayBase::COPY)
                        {
                            n = dstfa[tag.dstIndex].copyFromMem(bx,dcomp,ncomp,dptr);
                        }
                        else
                        {
                            fab.resize(bx,ncomp);
                            n = fab.copyFromMem(bx,0,ncomp,dptr);
                            dstfa[tag.dstIndex].plus(fab,bx,bx,0,dcomp,ncomp);
                        }
                        dptr += n;
                    }
                    BL_ASSERT(dptr == recv_data[k] + recv_size[k]);
                }
            }
        }
        
        for (auto p : recv_data) {
            amrex::The_Arena()->free(p);
        }
    }
	
    if (N_snds > 0) {
        if (FabArrayBase::do_async_sends && ! thecpc.m_SndTags->empty()) {
            Vector<MPI_Status> stats;
            FabArrayBase::WaitForAsyncSends(N_snds,send_reqs,send_data,stats);
        }
    }

#ifdef BL_USE_TEAM
    ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#endif
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
FabArray<FAB>::PostRcvs (const MapOfCopyComTagContainers&  m_RcvVols,
                         const MapOfCopyComTagContainers&  m_RcvTags,
                         Vector<char*>&                     recv_data,
                         Vector<int>&                       recv_size,
                         Vector<int>&                       recv_from,
                         Vector<MPI_Request>&               recv_reqs,
                         int                               icomp,
                         int                               ncomp,
                         int                               SeqNum,
                         int                               preSeqNum)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();
    recv_reqs.clear();

    for (const auto& kv : m_RcvVols) // loop over senders
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

    int recv_counter = 0;
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
            recv_data[i] = static_cast<char*>(amrex::The_Arena()->alloc(recv_size[i]));
            recv_reqs[i] = ParallelDescriptor::Arecv(recv_data[i], recv_size[i],
                                                     ParallelContext::global_to_local_rank(recv_from[i]),
                                                     SeqNum, comm).req();
        }
        
        ++recv_counter;
    }
}
#endif


#ifdef BL_USE_MPI3
template <class FAB>
void
FabArray<FAB>::PostRcvs_MPI_Onesided (const MapOfCopyComTagContainers&  m_RcvVols,
                                      char*&                            the_recv_data,
                                      Vector<char*>&                     recv_data,
                                      Vector<int>&                       recv_size,
                                      Vector<int>&                       recv_from,
                                      Vector<MPI_Request>&               recv_reqs,
                                      Vector<MPI_Aint> &                 recv_disp,
                                      int                               icomp,
                                      int                               ncomp,
                                      int                               SeqNum,
                                      MPI_Win &                         win)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();
    recv_reqs.clear();
    recv_disp.clear();

    std::size_t TotalRcvsVolume = 0;
    for (const auto& kv : m_RcvVols) // loop over senders
    {
        std::size_t nbytes = 0;
        for (const auto& cct : kv.second)
        {
            nbytes += (*this)[cct.dstIndex].nBytes(cct.dbox,icomp,ncomp);
        }

        BL_ASSERT(nbytes < std::numeric_limits<int>::max());

        TotalRcvsVolume += nbytes;

        recv_data.push_back(nullptr);
        recv_size.push_back(static_cast<int>(nbytes));
        recv_from.push_back(kv.first);
        recv_reqs.push_back(MPI_REQUEST_NULL);
        recv_disp.push_back(0);
    }

    if (TotalRcvsVolume == 0)
    {
        the_recv_data = nullptr;
    }
    else
    {
        the_recv_data = static_cast<char*>(amrex::The_Arena()->alloc(TotalRcvsVolume));

        MPI_Win_attach(win, the_recv_data, TotalRcvsVolume);

        int nrecv = recv_from.size();
        
        std::size_t offset = 0;
        for (int i = 0; i < nrecv; ++i)
        {
            if (recv_size[i] > 0)
            {
                recv_data[i] = the_recv_data + offset;
                MPI_Get_address(recv_data[i], &recv_disp[i]);
                recv_reqs[i] = ParallelDescriptor::Asend(&recv_disp[i],1,recv_from[i],SeqNum).req();
                offset += recv_size[i];
            }
        }
    }
}
#endif   // BL_USE_MPI3


#ifdef BL_USE_UPCXX
template <class FAB>
void
FabArray<FAB>::PostRcvs_PGAS (const MapOfCopyComTagContainers&  m_RcvVols,
                              char*&                            the_recv_data,
                              Vector<char*>&                     recv_data,
                              Vector<int>&                       recv_size,
                              Vector<int>&                       recv_from,
                              int                               icomp,
                              int                               ncomp,
                              int                               SeqNum,
                              upcxx::event*                     recv_event)
{
    recv_data.clear();
    recv_size.clear();
    recv_from.clear();

    std::size_t TotalRcvsVolume = 0;
    for (const auto& kv : m_RcvVols) // loop over senders
    {
        std::size_t nbytes = 0;
        for (const auto& cct : kv.second)
        {
            nbytes += (*this)[cct.dstIndex].nBytes(cct.dbox,icomp,ncomp);
        }

        BL_ASSERT(nbytes < std::numeric_limits<int>::max());

        TotalRcvsVolume += nbytes;

        recv_data.push_back(nullptr);
        recv_size.push_back(static_cast<int>(nbytes));
        recv_from.push_back(kv.first);
    }

    if (TotalRcvsVolume == 0)
    {
        the_recv_data = nullptr;
    }
    else
    {
        the_recv_data = static_cast<char*>(BLPgas::alloc(TotalRcvsVolume));

        int nrecv = recv_from.size();

        std::size_t offset = 0;
        for (int i = 0; i < nrecv; ++i)
        {
            if (recv_size[i] > 0)
            {
                recv_data[i] = the_recv_data + offset;

                /// Send an AM to the sender with the recv pointer
                upcxx::global_ptr<void> dst_ptr =
                    upcxx::global_ptr<void>(recv_data[i], upcxx::myrank());
                // Increment the event reference before launching the remote task
                recv_event->incref();
                // Launch a remote task on the sender to send me data
                BLPgas::Request(recv_from[i], dst_ptr, recv_size[i], SeqNum, recv_event);
                upcxx::advance(); // poll the UPC++ progress engine and the network
                offset += recv_size[i];
            }
        }
    }
}
#endif /* BL_USE_UPCXX */

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
#pragma omp parallel
#endif
        for (MFIter fai(*this,true); fai.isValid(); ++fai)
        {
            const Box& bx = fai.growntilebox(nghost);
            get(fai).copy(src[fai],bx,scomp,bx,dcomp,ncomp);
        }

        return;
    }

#ifdef BL_USE_MPI

    FabArrayBase::CPC cpc(boxArray(), nghost, DistributionMap(), src.DistributionMap());

    ParallelCopy(src, scomp, dcomp, ncomp, nghost, nghost, Periodicity::NonPeriodic(),
                 FabArrayBase::COPY, &cpc);

#endif
}
