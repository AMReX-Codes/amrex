#ifndef AMREX_NONLOCAL_BC_IMPL_H_
#define AMREX_NONLOCAL_BC_IMPL_H_
#include <AMReX_Config.H>

namespace amrex { namespace NonLocalBC {

struct Rotate90ClockWise {
    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(iv[1], -1-iv[0], iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{a.y, -1-a.x, a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.bigEnd  (2))}));
    }
};

struct Rotate90CounterClockWise {
    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(-1-iv[1], iv[0], iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{-1-a.y, a.x, a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.bigEnd  (2))}));
    }
};

struct Rotate90DstToSrc
{
    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        if (a.x < 0) {
            return Rotate90ClockWise()(a);
        } else {
            return Rotate90CounterClockWise()(a);
        }
    }
};

struct Rotate180Fn {
    int Ly;

    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(-1-iv[0], Ly-1-iv[1], iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{-1-a.x, Ly-1-a.y, a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.bigEnd  (1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.smallEnd(1),
                                                   box.bigEnd  (2))}));
    }
};

struct PolarFn {
    int Lx, Ly;

    AMREX_GPU_HOST_DEVICE
    int i_index (int i) const noexcept {
        return (i < Lx/2) ? -1-i : 2*Lx-1-i;
    }

    AMREX_GPU_HOST_DEVICE
    int j_index (int j) const noexcept {
        return (j < Ly/2) ? j+Ly/2 : j-Ly/2;
    }

    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(i_index(iv[0]), j_index(iv[1]), iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{i_index(a.x), j_index(a.y), a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.bigEnd  (2))}));
    }
};

struct PolarFn2 { // for the x-y corners
    int Lx, Ly;

    AMREX_GPU_HOST_DEVICE
    int i_index (int i) const noexcept {
        return (i < Lx/2) ? -1-i : 2*Lx-1-i;
    }

    AMREX_GPU_HOST_DEVICE
    int j_index (int j) const noexcept {
        if (j < 0) {
            return j+Ly/2;
        } else if (j >= Ly) {
            return j-Ly/2;
        } else if (j < Ly/2) {
            return j-Ly/2;
        } else {
            return j+Ly/2;
        }
    }

    AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(i_index(iv[0]), j_index(iv[1]), iv[2])};
    }

    AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{i_index(a.x), j_index(a.y), a.z};
    }

    Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.bigEnd  (2))}));
    }
};

template <class FAB, class DTOS>
void
local_copy_cpu (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
                DTOS dtos)
{
    auto const& LocTags = *(cmd.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) return;
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int itag = 0; itag < N_locs; ++itag) {
        const auto& tag = LocTags[itag];
        auto const& sfab = mf.const_array(tag.srcIndex);
        auto const& dfab = mf.array      (tag.dstIndex);
        amrex::LoopConcurrentOnCpu(tag.dbox, ncomp, [=] (int i, int j, int k, int n) noexcept
        {
            auto const si = dtos(Dim3{i,j,k});
            dfab(i,j,k,scomp+n) = sfab(si.x,si.y,si.z,scomp+n);
        });
    }
}

template <class FAB, class DTOS>
void
unpack_recv_buffer_cpu (FabArray<FAB>& mf, int scomp, int ncomp,
                        Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos)
{
    amrex::ignore_unused(recv_size);

    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    using T = typename FAB::value_type;
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int ircv = 0; ircv < N_rcvs; ++ircv) {
        const char* dptr = recv_data[ircv];
        auto const& cctc = *recv_cctc[ircv];
        for (auto const& tag : cctc) {
            auto const& dfab = mf.array(tag.dstIndex);
            auto const& sfab = amrex::makeArray4((T const*)(dptr), tag.sbox, ncomp);
            amrex::LoopConcurrentOnCpu(tag.dbox, ncomp, [=] (int i, int j, int k, int n) noexcept
            {
                auto const si = dtos(Dim3{i,j,k});
                dfab(i,j,k,scomp+n) = sfab(si.x,si.y,si.z,n);
            });
            dptr += tag.sbox.numPts() * ncomp * sizeof(T);
            AMREX_ASSERT(dptr <= recv_data[ircv] + recv_size[ircv]);
        }
    }
}

#ifdef AMREX_USE_GPU
template <class T>
struct Array4Array4Box {
    Array4<T      > dfab;
    Array4<T const> sfab;
    Box dbox;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Box const& box () const noexcept { return dbox; }
};

template <class FAB, class DTOS>
void
local_copy_gpu (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
                DTOS dtos)
{
    auto const& LocTags = *(cmd.m_LocTags);
    int N_locs = LocTags.size();
    if (N_locs == 0) return;

    using T = typename FAB::value_type;
    Vector<Array4Array4Box<T> > loc_copy_tags;
    loc_copy_tags.reserve(N_locs);
    for (int itag = 0; itag < N_locs; ++itag) {
        const auto& tag = LocTags[itag];
        loc_copy_tags.push_back({mf.array(tag.dstIndex), mf.const_array(tag.srcIndex), tag.dbox});
    }

    ParallelFor(loc_copy_tags, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                                            Array4Array4Box<T> const& tag) noexcept
    {
        auto const si = dtos(Dim3{i,j,k});
        tag.dfab(i,j,k,scomp+n) = tag.sfab(si.x,si.y,si.z,scomp+n);
    });
}

template <class FAB, class DTOS>
void
unpack_recv_buffer_gpu (FabArray<FAB>& mf, int scomp, int ncomp,
                        Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos)
{
    amrex::ignore_unused(recv_size);

    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) return;

    char* pbuffer = recv_data[0];
#if 0
    std::size_t szbuffer = 0;
    // For linear solver test on summit, this is slower than writing to
    // pinned memory directly on device.
    if (not ParallelDescriptor::UseGpuAwareMpi()) {
        // Memory in recv_data is pinned.
        szbuffer = (recv_data[N_rcvs-1]-recv_data[0]) + recv_size[N_rcvs-1];
        pbuffer = (char*)The_Arena()->alloc(szbuffer);
        Gpu::copyAsync(Gpu::hostToDevice,recv_data[0],recv_data[0]+szbuffer,pbuffer);
        Gpu::synchronize();
    }
#endif

    using T = typename FAB::value_type;
    using TagType = Array4Array4Box<T>;
    Vector<TagType> tags;
    tags.reserve(N_rcvs);

    for (int k = 0; k < N_rcvs; ++k)
    {
        std::size_t offset = recv_data[k]-recv_data[0];
        const char* dptr = pbuffer + offset;
        auto const& cctc = *recv_cctc[k];
        for (auto const& tag : cctc)
        {
            tags.emplace_back(TagType{mf.array(tag.dstIndex),
                                      amrex::makeArray4((T const*)dptr, tag.sbox, ncomp),
                                      tag.dbox});
            dptr += tag.dbox.numPts() * ncomp * sizeof(T);
            BL_ASSERT(dptr <= pbuffer + offset + recv_size[k]);
        }
    }

    ParallelFor(tags, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                                   Array4Array4Box<T> const& tag) noexcept
    {
        auto const si = dtos(Dim3{i,j,k});
        tag.dfab(i,j,k,scomp+n) = tag.sfab(si.x,si.y,si.z,n);
    });

    // There is Gpu::synchronize in ParalleFor above

    if (pbuffer != recv_data[0]) {
        The_Arena()->free(pbuffer);
    }
}
#endif

struct CommHandler
{
#ifdef AMREX_USE_MPI
    int mpi_tag;
    char* the_recv_data = nullptr;
    char* the_send_data = nullptr;
    //
    Vector<int>         recv_from;
    Vector<char*>       recv_data;
    Vector<std::size_t> recv_size;
    Vector<MPI_Request> recv_reqs;
    //
    Vector<MPI_Request> send_reqs;
#endif
};

template <class FAB, class DTOS>
CommHandler
Comm_nowait (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             DTOS dtos)
{
#ifdef AMREX_USE_MPI
    if (ParallelContext::NProcsSub() == 1)
#endif
    {
        int N_locs = (*cmd.m_LocTags).size();
        if (N_locs == 0) return CommHandler{};
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            local_copy_gpu(mf, scomp, ncomp, cmd, dtos);
        } else
#endif
        {
            local_copy_cpu(mf, scomp, ncomp, cmd, dtos);
        }
        return CommHandler{};
    }

#ifdef AMREX_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum = ParallelDescriptor::SeqNum();

    const int N_locs = cmd.m_LocTags->size();
    const int N_rcvs = cmd.m_RcvTags->size();
    const int N_snds = cmd.m_SndTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0) {
        // No work to do.
        return CommHandler{};
    }

    CommHandler handler{};
    handler.mpi_tag = SeqNum;

    if (N_rcvs > 0) {
        mf.PostRcvs(*cmd.m_RcvTags, handler.the_recv_data, handler.recv_data, handler.recv_size,
                    handler.recv_from, handler.recv_reqs, ncomp, SeqNum);
    }

    if (N_snds > 0)
    {
        Vector<std::size_t> send_size;
        Vector<int> send_rank;
        Vector<char*> send_data;
        Vector<const FabArrayBase::CopyComTagsContainer*> send_cctc;
        mf.PrepareSendBuffers(*cmd.m_SndTags, handler.the_send_data, send_data,
                              send_size, send_rank, handler.send_reqs, send_cctc, ncomp);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            FabArray<FAB>::pack_send_buffer_gpu(mf, scomp, ncomp, send_data,
                                                send_size, send_cctc);
        } else
#endif
        {
            FabArray<FAB>::pack_send_buffer_cpu(mf, scomp, ncomp, send_data,
                                                send_size, send_cctc);
        }

        FabArray<FAB>::PostSnds(send_data, send_size, send_rank, handler.send_reqs, SeqNum);
    }

    if (N_locs > 0)
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            local_copy_gpu(mf, scomp, ncomp, cmd, dtos);
        } else
#endif
        {
            local_copy_cpu(mf, scomp, ncomp, cmd, dtos);
        }
    }

    return handler;
#endif
}

#ifdef AMREX_USE_MPI
template <class FAB, class DTOS>
void
Comm_finish (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             CommHandler& handler, DTOS dtos)
{
    if (ParallelContext::NProcsSub() == 1) return;

    const int N_rcvs = cmd.m_RcvTags->size();
    if (N_rcvs > 0)
    {
        Vector<const FabArrayBase::CopyComTagsContainer*> recv_cctc(N_rcvs,nullptr);
        for (int k = 0; k < N_rcvs; ++k) {
            auto const& cctc = cmd.m_RcvTags->at(handler.recv_from[k]);
            recv_cctc[k] = &cctc;
        }

        Vector<MPI_Status> stats(N_rcvs);
        ParallelDescriptor::Waitall(handler.recv_reqs, stats);
#ifdef AMREX_D_DEBUG
        if (!CheckRcvStats(stats, handler.recv_size, handler.mpi_tag)) {
            amrex::Abort("NonLocalBC::Comm_finish failed with wrong message size");
        }
#endif

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            unpack_recv_buffer_gpu(mf, scomp, ncomp, handler.recv_data,
                                   handler.recv_size, recv_cctc, dtos);
        } else
#endif
        {
            unpack_recv_buffer_cpu(mf, scomp, ncomp, handler.recv_data,
                                   handler.recv_size, recv_cctc, dtos);
        }

        amrex::The_FA_Arena()->free(handler.the_recv_data);
    }

    const int N_snds = cmd.m_SndTags->size();
    if (N_snds > 0) {
        Vector<MPI_Status> stats(handler.send_reqs.size());
        ParallelDescriptor::Waitall(handler.send_reqs, stats);
        amrex::The_FA_Arena()->free(handler.the_send_data);
    }
}
#endif

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain)
{
    BL_PROFILE("Rotate90");

    AMREX_ASSERT(domain.cellCentered());
    AMREX_ASSERT(domain.smallEnd() == 0);
    AMREX_ASSERT(domain.length(0) == domain.length(1));
    AMREX_ASSERT(mf.is_cell_centered());
    AMREX_ASSERT(scomp < mf.nComp() && scomp+ncomp <= mf.nComp());
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()) && nghost[0] == nghost[1]);

    if (nghost[0] <= 0) return;

    const FabArrayBase::RB90& TheRB90 = mf.getRB90(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, TheRB90,Rotate90DstToSrc{});

    Box corner(-nghost, IntVect{AMREX_D_DECL(-1,-1,domain.bigEnd(2)+nghost[2])});
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
        Box const& bx = corner & mfi.fabbox();
        if (bx.ok()) {
            auto const& fab = mf.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx,ncomp,i,j,k,n,
            {
                fab(i,j,k,n) = fab(-i-1,-j-1,k,n);
            });
        }
    }

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, TheRB90, handler, Rotate90DstToSrc{});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, Box const& domain)
{
    Rotate90(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain)
{
    BL_PROFILE("Rotate180");

    AMREX_ASSERT(domain.cellCentered());
    AMREX_ASSERT(domain.smallEnd() == 0);
    AMREX_ASSERT(domain.length(1) % 2 == 0);
    AMREX_ASSERT(mf.is_cell_centered());
    AMREX_ASSERT(scomp < mf.nComp() && scomp+ncomp <= mf.nComp());
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()));

    if (nghost[0] <= 0) return;

    const FabArrayBase::RB180& TheRB180 = mf.getRB180(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, TheRB180, Rotate180Fn{domain.length(1)});

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, TheRB180, handler, Rotate180Fn{domain.length(1)});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, Box const& domain)
{
    Rotate180(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, int scomp, int ncomp, IntVect const& nghost, Box const& domain)
{
    BL_PROFILE("FillPolar");

    AMREX_ASSERT(domain.cellCentered());
    AMREX_ASSERT(domain.smallEnd() == 0);
    AMREX_ASSERT(domain.length(1) % 2 == 0);
    AMREX_ASSERT(mf.is_cell_centered());
    AMREX_ASSERT(scomp < mf.nComp() && scomp+ncomp <= mf.nComp());
    AMREX_ASSERT(nghost.allLE(mf.nGrowVect()));

    if (nghost[0] <= 0) return;

    const FabArrayBase::PolarB& ThePolarB = mf.getPolarB(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, ThePolarB,
                               PolarFn{domain.length(0), domain.length(1)});

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, ThePolarB, handler,
                PolarFn{domain.length(0), domain.length(1)});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
amrex::EnableIf_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, Box const& domain)
{
    FillPolar(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

}}

#endif
