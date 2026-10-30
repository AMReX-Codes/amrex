#ifndef AMREX_FB_UNPACK_ON_ARRIVAL_H_
#define AMREX_FB_UNPACK_ON_ARRIVAL_H_
#include <AMReX_Config.H>

#include <algorithm>

/**
 * \file AMReX_FBUnpackOnArrival.H
 * \brief GPU FillBoundary and SumBoundary with GPU-aware MPI. Messages are
 *        unpacked as they arrive when at least `amrex.fb.on_arrival_min_bytes`
 *        are received per call, otherwise together after all have arrived.
 *
 * Everything here runs on the current GPU stream. The pack kernel relies on
 * that to see what earlier kernels wrote into the FabArray. On return the
 * stream is synchronized unless inside a Gpu::NoSyncRegion, where only work
 * on the current stream may still be running. Code that packs or unpacks on
 * another stream must first wait on the current stream and must not leave
 * work on the other stream behind.
 */

#ifdef AMREX_FB_UNPACK_ON_ARRIVAL

namespace amrex {

/// \cond DOXYGEN_IGNORE
namespace detail {

// Unpack functor: ADD with plain += when the tags are thread safe, atomics
// otherwise; COPY stores (benign races store the same value).
template <typename BUF, typename T, bool ADD, bool ATOMIC>
struct FBOnArrivalUnpack {
    char const* dbase;
    int dcomp;
    int ncomp;

    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void operator() (int icell, int ncells, int i, int j, int k,
                     CommRecvBufTag<T> const& tag) const noexcept
    {
        if (icell < ncells) {
            Array4<BUF const> sfab{(BUF const*)(dbase + tag.poff),
                                   amrex::begin(tag.bx), amrex::end(tag.bx), ncomp};
            for (int n = 0; n < ncomp; ++n) {
                if constexpr (ADD && ATOMIC) {
                    Gpu::Atomic::AddNoRet(&tag.dfab(i,j,k,n+dcomp), (T)sfab(i,j,k,n));
                } else if constexpr (ADD) {
                    tag.dfab(i,j,k,n+dcomp) += (T)sfab(i,j,k,n);
                } else {
                    tag.dfab(i,j,k,n+dcomp) = (T)sfab(i,j,k,n);
                }
            }
        }
    }
};

template <typename BUF, typename T>
void fb_oa_pack (TagVector<CommSendBufTag<T>> const& tv, char* dbase, int scomp, int ncomp)
{
    detail::ParallelFor_doit(tv,
        [=] AMREX_GPU_DEVICE (
            int icell, int ncells, int i, int j, int k, CommSendBufTag<T> const& tag) noexcept
    {
        if (icell < ncells) {
            Array4<BUF> dfab{(BUF*)(dbase + tag.poff), amrex::begin(tag.bx),
                             amrex::end(tag.bx), ncomp};
            for (int n = 0; n < ncomp; ++n) {
                dfab(i,j,k,n) = (BUF)tag.sfab(i,j,k,n+scomp);
            }
        }
    });
}

// Ranges of warps of a TagVector, passed to the kernel by value. The
// messages completed by one Waitsome are unpacked by one launch over them.
struct FBOnArrivalWarpRanges {
    static constexpr int max_ranges = 16;
    int begin[max_ranges] = {0};    // first warp of each range in the TagVector
    int prefix[max_ranges+1] = {0}; // prefix sum of the range sizes
    int n = 0;

    // Returns false if the table is full and the range cannot be merged.
    bool add (int b, int nw) noexcept {
        if (n > 0 && begin[n-1] + (prefix[n]-prefix[n-1]) == b) {
            prefix[n] += nw;
        } else if (n < max_ranges) {
            begin[n] = b;
            prefix[n+1] = prefix[n] + nw;
            ++n;
        } else {
            return false;
        }
        return true;
    }
};

template <typename T, typename F>
void fb_oa_unpack_ranges (TagVector<CommRecvBufTag<T>> const& tv,
                          FBOnArrivalWarpRanges const& r, F const& f)
{
    const int nwarps = r.prefix[r.n];
    if (nwarps == 0) { return; }

    constexpr int nthreads = TagVector<CommRecvBufTag<T>>::nthreads;
    constexpr int nwarps_per_block = nthreads / Gpu::Device::warp_size;
    const int nblocks = (nwarps + nwarps_per_block-1) / nwarps_per_block;
    const auto d_tags = tv.d_tags;
    const auto d_nwarps = tv.d_nwarps;
    const auto ntags = tv.ntags;

    amrex::launch<nthreads>(nblocks, Gpu::gpuStream(),
#ifdef AMREX_USE_SYCL
    [=] AMREX_GPU_DEVICE (sycl::nd_item<1> const& item) noexcept
    [[sycl::reqd_work_group_size(nthreads)]]
    [[sycl::reqd_sub_group_size(Gpu::Device::warp_size)]]
#else
    [=] AMREX_GPU_DEVICE () noexcept
#endif
    {
#ifdef AMREX_USE_SYCL
        std::size_t g_tid = item.get_global_id(0);
#else
        auto g_tid = std::size_t(blockDim.x)*blockIdx.x + threadIdx.x;
#endif
        auto l_wid = int(g_tid / Gpu::Device::warp_size);
        if (l_wid >= nwarps) { return; }
        int ir = amrex::bisect(r.prefix, 0, r.n, l_wid);
        int g_wid = r.begin[ir] + (l_wid - r.prefix[ir]);
        int tag_id = amrex::bisect(d_nwarps, 0, ntags, g_wid);
        int b_wid = g_wid - d_nwarps[tag_id];
#ifdef AMREX_USE_SYCL
        int lane = item.get_local_id(0) % Gpu::Device::warp_size;
#else
        int lane = threadIdx.x % Gpu::Device::warp_size;
#endif
        int icell = b_wid*Gpu::Device::warp_size + lane;
        tagparfor_call_f(icell, d_tags[tag_id], f);
    });
}

}
/// \endcond

template <class FAB>
template <typename BUF>
FBOnArrivalMeta<FAB>*
FabArray<FAB>::FB_on_arrival_meta (const FB& TheFB, int ncomp)
{
    RecvSendCopyHandlerKey key{TheFB.m_id, sizeof(BUF), ncomp};
    if (auto it = m_fb_on_arrival_meta.find(key); it != m_fb_on_arrival_meta.end()) {
        return it->second.get();
    }

    if (m_fb_on_arrival_meta.size() >= 32) { m_fb_on_arrival_meta.clear(); }

    auto e = std::make_unique<FBOnArrivalMeta<FAB>>();

    auto make_layout = [&] (MapOfCopyComTagContainers const& tags, bool is_send,
                            typename FBOnArrivalMeta<FAB>::Layout& lo)
    {
        for (auto const& [rank, cctc] : tags) {
            std::size_t nbytes = 0;
            for (auto const& cct : cctc) {
                nbytes += (is_send ? cct.sbox : cct.dbox).numPts() * ncomp * sizeof(BUF);
            }
            std::size_t acd = ParallelDescriptor::sizeof_selected_comm_data_type(nbytes);
            nbytes = amrex::aligned_size(acd, nbytes);
            lo.nbytes = amrex::aligned_size(std::max(alignof(BUF), acd), lo.nbytes);
            lo.rank.push_back(rank);
            lo.size.push_back(nbytes);
            lo.offset.push_back(lo.nbytes);
            lo.nbytes += nbytes;
        }
    };

    make_layout(*TheFB.m_SndTags, true, e->snd);
    make_layout(*TheFB.m_RcvTags, false, e->rcv);

    {
        Vector<CommSendBufTag<value_type>> tags;
        int m = 0;
        for (auto const& kv : *TheFB.m_SndTags) {
            auto off = static_cast<std::ptrdiff_t>(e->snd.offset[m]);
            for (auto const& cct : kv.second) {
                tags.push_back(CommSendBufTag<value_type>{
                    .sfab = this->const_array(cct.srcIndex),
                    .poff = off, .bx = cct.sbox});
                off += cct.sbox.numPts() * ncomp * sizeof(BUF);
            }
            ++m;
        }
        e->send_tv.define(tags);
    }

    {
        // Tags are ordered by message, so each message owns a contiguous
        // range of warps in recv_tv (same warp count as TagVector::define).
        Vector<CommRecvBufTag<value_type>> tags;
        e->rcv_warp_begin.push_back(0);
        int m = 0;
        for (auto const& kv : *TheFB.m_RcvTags) {
            auto off = static_cast<std::ptrdiff_t>(e->rcv.offset[m]);
            int nwarps = e->rcv_warp_begin.back();
            for (auto const& cct : kv.second) {
                const int li = this->localindex(cct.dstIndex);
                tags.push_back(CommRecvBufTag<value_type>{
                    .dfab = this->atLocalIdx(li).array(),
                    .poff = off, .bx = cct.dbox});
                off += cct.dbox.numPts() * ncomp * sizeof(BUF);
                nwarps += static_cast<int>((detail::get_tag_size(tags.back())
                                            + Gpu::Device::warp_size-1)
                                           / Gpu::Device::warp_size);
            }
            e->rcv_warp_begin.push_back(nwarps);
            ++m;
        }
        e->recv_tv.define(tags);
        AMREX_ALWAYS_ASSERT(e->rcv_warp_begin.back() == e->recv_tv.ntotwarps);
    }

    auto* r = e.get();
    m_fb_on_arrival_meta[key] = std::move(e);
    return r;
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FB_on_arrival_nowait (const FB& TheFB, int scomp, int ncomp, int SeqNum,
                                     bool sumboundary)
{
    auto* e = FB_on_arrival_meta<BUF>(TheFB, ncomp);

    const int N_locs = TheFB.m_LocTags->size();
    const auto N_rcvs = static_cast<int>(e->rcv.rank.size());
    const auto N_snds = static_cast<int>(e->snd.rank.size());

    fbd = std::make_unique<FBData<FAB>>();
    fbd->fb    = &TheFB;
    fbd->scomp = scomp;
    fbd->ncomp = ncomp;
    fbd->tag   = SeqNum;
    fbd->oa_meta = e;

    auto alloc = [] (std::size_t n) {
        return n > 0 ? static_cast<char*>(The_Comms_Arena()->alloc(n)) : nullptr;
    };
    fbd->the_recv_data = alloc(e->rcv.nbytes);
    fbd->the_send_data = alloc(e->snd.nbytes);

    MPI_Comm comm = ParallelContext::CommunicatorSub();

    if (N_rcvs > 0) {
        BL_PROFILE("FB::post_recvs");
        fbd->recv_reqs.assign(N_rcvs, MPI_REQUEST_NULL);
        fbd->recv_stat.resize(N_rcvs);
        for (int i = 0; i < N_rcvs; ++i) {
            if (e->rcv.size[i] > 0) {
                const int rank = ParallelContext::global_to_local_rank(e->rcv.rank[i]);
                fbd->recv_reqs[i] = ParallelDescriptor::Arecv
                    (fbd->the_recv_data + e->rcv.offset[i], e->rcv.size[i],
                     rank, SeqNum, comm).req();
            }
        }
    }

    if (N_snds > 0) {
        {
            BL_PROFILE("FB::pack");
            detail::fb_oa_pack<BUF>(e->send_tv, fbd->the_send_data, scomp, ncomp);
        }
        {
            BL_PROFILE("FB::pack_sync");
            Gpu::streamSynchronize();
        }
        BL_PROFILE("FB::post_sends");
        fbd->send_reqs.assign(N_snds, MPI_REQUEST_NULL);
        for (int j = 0; j < N_snds; ++j) {
            if (e->snd.size[j] > 0) {
                const int rank = ParallelContext::global_to_local_rank(e->snd.rank[j]);
                fbd->send_reqs[j] = ParallelDescriptor::Asend
                    (fbd->the_send_data + e->snd.offset[j], e->snd.size[j],
                     rank, SeqNum, comm).req();
            }
        }
    }

    if (N_rcvs > 0) { FillBoundary_test(); }

    if (N_locs > 0) {
        if (sumboundary) {
            if constexpr (amrex::IsAddAssignable<value_type>::value) {
                BL_PROFILE("FB::local_add");
                FB_local_add_gpu(TheFB, scomp, ncomp, false);
            }
        } else {
            BL_PROFILE("FB::local_copy");
            FB_local_copy_gpu(TheFB, scomp, ncomp);
        }
        if (N_rcvs > 0) { FillBoundary_test(); }
    }
}

template <class FAB>
template <typename BUF>
void
FabArray<FAB>::FB_on_arrival_finish ()
{
    auto* e = fbd->oa_meta;
    const FB* TheFB = fbd->fb;
    const auto N_rcvs = static_cast<int>(e->rcv.rank.size());
    const auto N_snds = static_cast<int>(e->snd.rank.size());

    if (N_rcvs > 0) {
        const bool sumboundary = TheFB->m_sb_snghost.allGE(0);
        // Below the threshold, wait for all messages and unpack them in one
        // launch (the Waitsome loop costs more than it hides for tiny messages).
        const bool on_arrival = e->rcv.nbytes >=
            static_cast<std::size_t>(std::max(Long(0), FabArrayBase::m_fb_on_arrival_min_bytes));
        auto unpack = [&] (auto const& f) {
            if (!on_arrival) {
                {
                    BL_PROFILE("FB::wait_recvs");
                    ParallelDescriptor::Waitall(fbd->recv_reqs, fbd->recv_stat);
                }
                BL_PROFILE("FB::unpack");
                detail::ParallelFor_doit(e->recv_tv, f);
                return;
            }
            // The messages completed by each Waitsome are unpacked by one
            // launch. Messages already completed by FillBoundary_test have
            // null requests, which Waitsome never returns, so they go first.
            detail::FBOnArrivalWarpRanges ranges;
            auto flush = [&] () {
                if (ranges.n > 0) {
                    BL_PROFILE("FB::unpack");
                    detail::fb_oa_unpack_ranges(e->recv_tv, ranges, f);
                    ranges.n = 0;
                }
            };
            auto add = [&] (int m) {
                const int b = e->rcv_warp_begin[m];
                const int nw = e->rcv_warp_begin[m+1] - b;
                if (nw > 0 && !ranges.add(b, nw)) {
                    flush();
                    ranges.add(b, nw);
                }
            };
            for (int m = 0; m < N_rcvs; ++m) {
                if (fbd->recv_reqs[m] == MPI_REQUEST_NULL) { add(m); }
            }
            flush();
            Vector<int> indx(N_rcvs);
            for (;;) {
                int completed = 0;
                {
                    BL_PROFILE("FB::wait_recvs");
                    ParallelDescriptor::Waitsome(fbd->recv_reqs, completed, indx, fbd->recv_stat);
                }
                if (completed == MPI_UNDEFINED || completed == 0) { break; }
                std::sort(indx.begin(), indx.begin()+completed);
                for (int c = 0; c < completed; ++c) { add(indx[c]); }
                flush();
            }
        };

        using T = value_type;
        if (sumboundary) {
            if constexpr (amrex::IsAddAssignable<T>::value) {
                if (TheFB->m_threadsafe_rcv) {
                    unpack(detail::FBOnArrivalUnpack<BUF,T,true,false>
                           {fbd->the_recv_data, fbd->scomp, fbd->ncomp});
                } else if constexpr (amrex::HasAtomicAdd<T>::value) {
                    unpack(detail::FBOnArrivalUnpack<BUF,T,true,true>
                           {fbd->the_recv_data, fbd->scomp, fbd->ncomp});
                } else {
                    amrex::Abort("FB_on_arrival_finish: SumBoundary needs thread-safe tags or atomic add");
                }
            }
        } else {
            unpack(detail::FBOnArrivalUnpack<BUF,T,false,false>
                   {fbd->the_recv_data, fbd->scomp, fbd->ncomp});
        }
    }

    // In a no-sync region the recv buffer is released in stream order.
    const bool sync = !Gpu::inNoSyncRegion();
    if (sync) {
        BL_PROFILE("FB::unpack_sync");
        Gpu::streamSynchronize();
    }
    if (fbd->the_recv_data) {
        if (sync) {
            The_Comms_Arena()->free(fbd->the_recv_data);
        } else {
            Gpu::freeAsync(The_Comms_Arena(), fbd->the_recv_data);
        }
    }

    if (N_snds > 0) {
        BL_PROFILE("FB::wait_sends");
        Vector<MPI_Status> stats(fbd->send_reqs.size());
        ParallelDescriptor::Waitall(fbd->send_reqs, stats);
    }
    // The pack kernel finished before the sends were posted.
    if (fbd->the_send_data) { The_Comms_Arena()->free(fbd->the_send_data); }

    fbd.reset();
}

}

#endif
#endif
