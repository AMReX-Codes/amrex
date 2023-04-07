#ifndef AMREX_NONLOCAL_BC_H_
#include "AMReX_NonLocalBC.H"
#endif

#ifndef AMREX_NONLOCAL_BC_IMPL_H_
#define AMREX_NONLOCAL_BC_IMPL_H_
#include <AMReX_Config.H>
#include <AMReX_TypeTraits.H>

namespace amrex::NonLocalBC {
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

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int i_index (int i) const noexcept {
        return (i < Lx/2) ? -1-i : 2*Lx-1-i;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int j_index (int j) const noexcept {
        return (j < Ly/2) ? j+Ly/2 : j-Ly/2;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(i_index(iv[0]), j_index(iv[1]), iv[2])};
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{i_index(a.x), j_index(a.y), a.z};
    }

    [[nodiscard]] Box operator() (Box const& box) const noexcept {
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

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int i_index (int i) const noexcept {
        return (i < Lx/2) ? -1-i : 2*Lx-1-i;
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    int j_index (int j) const noexcept {
        if (j < 0) { // NOLINT
            return j+Ly/2;
        } else if (j >= Ly) { // NOLINT
            return j-Ly/2;
        } else if (j < Ly/2) {
            return j-Ly/2;
        } else {
            return j+Ly/2;
        }
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    IntVect operator() (IntVect const& iv) const noexcept {
        return IntVect{AMREX_D_DECL(i_index(iv[0]), j_index(iv[1]), iv[2])};
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& a) const noexcept {
        return Dim3{i_index(a.x), j_index(a.y), a.z};
    }

    [[nodiscard]] Box operator() (Box const& box) const noexcept {
        return Box(operator()(IntVect{AMREX_D_DECL(box.bigEnd  (0),
                                                   box.smallEnd(1),
                                                   box.smallEnd(2))}),
                   operator()(IntVect{AMREX_D_DECL(box.smallEnd(0),
                                                   box.bigEnd  (1),
                                                   box.bigEnd  (2))}));
    }
};

//! \brief Perform the local copies from src to dest without doing any MPI communication.
//!
//! This function assumes that all destination and source boxes stored in the local copy comm tags
//! are related by the DTOS. If this is not the case the behaviour will at best be caught as an
//! assertion.
//!
//! \param[out] dest  The destination FabArray that will be filled with data.
//!
//! \param[in]  src   The source FabArray where data will be taken from.
//!
//! \param[in]  dcomp The first component at the destination.
//!
//! \param[in]  scomp The first component at the source.
//!
//! \param[in]  ncomp The number of components that will be copied.
//!
//! \param[in]  local_tags The vector of copy com tags that describes each local copy
//!                        transaction.
//!
//! \param[in]  dtos  The dest to source index mapping that will be used in the copy function.
//!
//! \param[in]  proj  A FAB projection that might transform the data.
//!
//! \return Nothing.
template <class FAB, class DTOS, class Proj>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
local_copy_cpu (FabArray<FAB>& dest, const FabArray<FAB>& src, int dcomp, int scomp, int ncomp,
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS dtos, Proj proj) noexcept {
    const auto N_locs = static_cast<int>(local_tags.size());
    if (N_locs == 0) return;
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
    for (int itag = 0; itag < N_locs; ++itag) {
        const auto& tag = local_tags[itag];
        auto const& sfab =  src.const_array(tag.srcIndex);
        auto const& dfab = dest.array      (tag.dstIndex);
        amrex::LoopConcurrentOnCpu(tag.dbox, ncomp, [=] (int i, int j, int k, int n) noexcept
        {
            auto const si = dtos(Dim3{i,j,k});
            dfab(i,j,k,dcomp+n) = proj(sfab,si,scomp+n);
        });
    }
}

//! \brief Unpack the received data into the local FABs.
//!
//! This function assumes that all destination and source boxes stored in the local copy comm tags
//! are related by the DTOS. If this is not the case the behaviour will at best be caught as an
//! assertion.
template <class FAB, class DTOS, class Proj>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
unpack_recv_buffer_cpu (FabArray<FAB>& mf, int dcomp, int ncomp, Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos, Proj proj) noexcept {
    amrex::ignore_unused(recv_size);

    const auto N_rcvs = static_cast<int>(recv_cctc.size());
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
            amrex::LoopConcurrentOnCpu(tag.dbox, ncomp, [=](int i, int j, int k, int n) noexcept {
                auto const si = dtos(Dim3{i, j, k});
                dfab(i, j, k, dcomp + n) = proj(sfab, si, n);
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

template <class FAB, class DTOS, class Proj>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
local_copy_gpu (FabArray<FAB>& dest, const FabArray<FAB>& src, int dcomp, int scomp, int ncomp,
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS dtos, Proj proj) noexcept {
    int N_locs = local_tags.size();
    if (N_locs == 0) return;

    using T = typename FAB::value_type;
    Vector<Array4Array4Box<T> > loc_copy_tags;
    loc_copy_tags.reserve(N_locs);
    for (auto const& tag : local_tags) {
        loc_copy_tags.push_back({dest.array(tag.dstIndex), src.const_array(tag.srcIndex), tag.dbox});
    }

    ParallelFor(loc_copy_tags, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                                            Array4Array4Box<T> const& tag) noexcept
    {
        auto const si = dtos(Dim3{i,j,k});
        tag.dfab(i,j,k,dcomp+n) = proj(tag.sfab, si, scomp+n);
    });
}

template <class FAB, class DTOS, class Proj>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3, DTOS, Dim3>() && IsFabProjection<Proj, FAB>()>
unpack_recv_buffer_gpu (FabArray<FAB>& mf, int scomp, int ncomp,
                        Vector<char*> const& recv_data,
                        Vector<std::size_t> const& recv_size,
                        Vector<FabArrayBase::CopyComTagsContainer const*> const& recv_cctc,
                        DTOS dtos, Proj proj)
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
        Gpu::streamSynchronize();
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
        tag.dfab(i,j,k,scomp+n) = proj(tag.sfab, si ,n);
    });

    // There is Gpu::streamSynchronize in ParallelFor above

    if (pbuffer != recv_data[0]) {
        The_Arena()->free(pbuffer);
    }
}
#endif

template <typename DTOS, typename>
MultiBlockCommMetaData::MultiBlockCommMetaData (const FabArrayBase& dst, const Box& dstbox, const FabArrayBase& src,
                                                const IntVect& ngrow, DTOS dtos)
    : MultiBlockCommMetaData(dst.boxArray(), dst.DistributionMap(), dstbox, src.boxArray(),
                                src.DistributionMap(), ngrow, dtos) {}

template <typename DTOS, typename>
MultiBlockCommMetaData::MultiBlockCommMetaData (const BoxArray& dstba, const DistributionMapping& dstdm,
                                                const Box& dstbox, const BoxArray& srcba,
                                                const DistributionMapping& srcdm, const IntVect& ngrow, DTOS dtos) {
    define(dstba, dstdm, dstbox, srcba, srcdm, ngrow, dtos);
}

template <typename DTOS>
std::enable_if_t<IsIndexMapping<DTOS>::value>
MultiBlockCommMetaData::define (const BoxArray& dstba, const DistributionMapping& dstdm, const Box& dstbox,
                                const BoxArray& srcba, const DistributionMapping& srcdm, const IntVect& ngrow,
                                DTOS dtos) {
    m_LocTags = std::make_unique<FabArrayBase::CopyComTagsContainer>();
    m_SndTags = std::make_unique<FabArrayBase::MapOfCopyComTagContainers>();
    m_RcvTags = std::make_unique<FabArrayBase::MapOfCopyComTagContainers>();
    const int myproc = ParallelDescriptor::MyProc();
    for (int i = 0, N = static_cast<int>(dstba.size()); i < N; ++i) {
        const int dest_owner = dstdm[i];
        const Box partial_dstbox = grow(dstba[i], ngrow) & dstbox;
        if (partial_dstbox.isEmpty()) {
            continue;
        }
        const Box partial_dstbox_mapped_in_src = Image(dtos, partial_dstbox).setType(srcba.ixType());
        enum { not_first_only = 0, first_only = 1 };
        std::vector<std::pair<int, Box>> boxes_from_src =
            srcba.intersections(partial_dstbox_mapped_in_src, not_first_only, ngrow);
        for (std::pair<int, Box> counted_box : boxes_from_src) {
            const int k = counted_box.first;
            const Box src_box = counted_box.second;
            AMREX_ASSERT(k < srcba.size());
            const int src_owner = srcdm[k];
            if (dest_owner == myproc || src_owner == myproc) {
                if (src_owner == dest_owner) {
                    const BoxList tilelist(src_box, FabArrayBase::comm_tile_size);
                    for (const Box& tilebox : tilelist) {
                        const Box inverse_image = InverseImage(dtos, tilebox).setType(dstba.ixType());
                        if ((inverse_image & partial_dstbox).ok()) {
                            m_LocTags->emplace_back(inverse_image, tilebox, i, k);
                        }
                    }
                } else {
                    const Box inverse_image = InverseImage(dtos, src_box).setType(dstba.ixType());
                    if ((inverse_image & partial_dstbox).ok()) {
                        FabArrayBase::CopyComTagsContainer& copy_tags =
                            (src_owner == myproc) ? (*m_SndTags)[dest_owner]
                                                    : (*m_RcvTags)[src_owner];
                        copy_tags.emplace_back(inverse_image, src_box, i, k);
                    }
                }
            }
        }
    }
}

template <class FAB, class DTOS>
#ifdef AMREX_USE_MPI
AMREX_NODISCARD
#endif
CommHandler
Comm_nowait (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             DTOS dtos)
{
#ifdef AMREX_USE_MPI
    if (ParallelContext::NProcsSub() == 1)
#endif
    {
        if (cmd.m_LocTags->empty()) { return CommHandler{}; }
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            local_copy_gpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        } else
#endif
        {
            local_copy_cpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        }
        return CommHandler{};
    }

#ifdef AMREX_USE_MPI
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    int SeqNum = ParallelDescriptor::SeqNum();

    const auto N_locs = cmd.m_LocTags->size();
    const auto N_rcvs = cmd.m_RcvTags->size();
    const auto N_snds = cmd.m_SndTags->size();

    if (N_locs == 0 && N_rcvs == 0 && N_snds == 0) {
        // No work to do.
        return CommHandler{};
    }

    CommHandler handler{};
    handler.mpi_tag = SeqNum;

    if (N_rcvs > 0) {
        handler.recv.the_data = mf.PostRcvs(*cmd.m_RcvTags, handler.recv.data, handler.recv.size,
                                            handler.recv.rank, handler.recv.request, ncomp, SeqNum);
    }

    if (N_snds > 0) {
        handler.send.the_data =
            mf.PrepareSendBuffers(*cmd.m_SndTags, handler.send.data, handler.send.size,
                                  handler.send.rank, handler.send.request, handler.send.cctc, ncomp);

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            FabArray<FAB>::pack_send_buffer_gpu(mf, scomp, ncomp, handler.send.data,
                                                handler.send.size, handler.send.cctc);
        } else
#endif
        {
            FabArray<FAB>::pack_send_buffer_cpu(mf, scomp, ncomp, handler.send.data,
                                                handler.send.size, handler.send.cctc);
        }

        FabArray<FAB>::PostSnds(handler.send.data, handler.send.size, handler.send.rank, handler.send.request, SeqNum);
    }

    if (N_locs > 0)
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            local_copy_gpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        } else
#endif
        {
            local_copy_cpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos);
        }
    }

    return handler;
#endif
}

#ifdef AMREX_USE_MPI
template <class FAB, class DTOS>
void
Comm_finish (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             CommHandler handler, DTOS dtos)
{
    if (ParallelContext::NProcsSub() == 1) return;

    const auto N_rcvs = static_cast<int>(cmd.m_RcvTags->size());
    if (N_rcvs > 0)
    {
        handler.recv.cctc.resize(N_rcvs, nullptr);
        for (int k = 0; k < N_rcvs; ++k) {
            auto const& cctc = cmd.m_RcvTags->at(handler.recv.rank[k]);
            handler.recv.cctc[k] = &cctc;
        }
        handler.recv.stats.resize(handler.recv.request.size());
        ParallelDescriptor::Waitall(handler.recv.request, handler.recv.stats);
#ifdef AMREX_DEBUG
        if (!CheckRcvStats(handler.recv.stats, handler.recv.size, handler.mpi_tag)) {
            amrex::Abort("NonLocalBC::Comm_finish failed with wrong message size");
        }
#endif

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion())
        {
            unpack_recv_buffer_gpu(mf, scomp, ncomp, handler.recv.data,
                                   handler.recv.size, handler.recv.cctc, dtos);
        } else
#endif
        {
            unpack_recv_buffer_cpu(mf, scomp, ncomp, handler.recv.data,
                                   handler.recv.size, handler.recv.cctc, dtos);
        }
    }

    if ( ! cmd.m_SndTags->empty() ) {
        handler.send.stats.resize(handler.send.request.size());
        ParallelDescriptor::Waitall(handler.send.request, handler.send.stats);
    }
}
#endif

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
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
    Comm_finish(mf, scomp, ncomp, TheRB90, std::move(handler), Rotate90DstToSrc{});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
Rotate90 (FabArray<FAB>& mf, Box const& domain)
{
    Rotate90(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
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
    Comm_finish(mf, scomp, ncomp, TheRB180, std::move(handler), Rotate180Fn{domain.length(1)});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
Rotate180 (FabArray<FAB>& mf, Box const& domain)
{
    Rotate180(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
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
    Comm_finish(mf, scomp, ncomp, ThePolarB, std::move(handler),
                PolarFn{domain.length(0), domain.length(1)});
#else
    amrex::ignore_unused(handler);
#endif
}

template <class FAB>
std::enable_if_t<IsBaseFab<FAB>::value>
FillPolar (FabArray<FAB>& mf, Box const& domain)
{
    FillPolar(mf, 0, mf.nComp(), mf.nGrowVect(), domain);
}

extern template MultiBlockCommMetaData ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox,
                                                    const FabArray<FArrayBox>& src, int destcomp,
                                                    int srccomp, int numcomp, const IntVect& ngrow,
                                                    MultiBlockIndexMapping, Identity);
}

#endif
