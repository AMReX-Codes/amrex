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
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS const& dtos, Proj const& proj) noexcept {
    const auto N_locs = static_cast<int>(local_tags.size());
    if (N_locs == 0) { return; }
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
                        DTOS const& dtos, Proj const& proj) noexcept {
    amrex::ignore_unused(recv_size);

    const auto N_rcvs = static_cast<int>(recv_cctc.size());
    if (N_rcvs == 0) { return; }

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
                FabArrayBase::CopyComTagsContainer const& local_tags, DTOS const& dtos, Proj const& proj) noexcept {
    int N_locs = local_tags.size();
    if (N_locs == 0) { return; }

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
                        DTOS const& dtos, Proj const& proj)
{
    amrex::ignore_unused(recv_size);

    const int N_rcvs = recv_cctc.size();
    if (N_rcvs == 0) { return; }

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
                                                const IntVect& ngrow, DTOS const& dtos)
    : MultiBlockCommMetaData(dst.boxArray(), dst.DistributionMap(), dstbox, src.boxArray(),
                                src.DistributionMap(), ngrow, dtos) {}

template <typename DTOS, typename>
MultiBlockCommMetaData::MultiBlockCommMetaData (const BoxArray& dstba, const DistributionMapping& dstdm,
                                                const Box& dstbox, const BoxArray& srcba,
                                                const DistributionMapping& srcdm, const IntVect& ngrow, DTOS const& dtos) {
    define(dstba, dstdm, dstbox, srcba, srcdm, ngrow, dtos);
}

template <typename DTOS>
std::enable_if_t<IsIndexMapping<DTOS>::value>
MultiBlockCommMetaData::define (const BoxArray& dstba, const DistributionMapping& dstdm, const Box& dstbox,
                                const BoxArray& srcba, const DistributionMapping& srcdm, const IntVect& ngrow,
                                DTOS const& dtos) {
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

template <class FAB, class DTOS, class Proj>
#ifdef AMREX_USE_MPI
AMREX_NODISCARD
#endif
CommHandler
Comm_nowait (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             DTOS const& dtos, Proj const& proj)
{
#ifdef AMREX_USE_MPI
    if (ParallelContext::NProcsSub() == 1)
#endif
    {
        if (cmd.m_LocTags->empty()) { return CommHandler{}; }
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            local_copy_gpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos, proj);
        } else
#endif
        {
            local_copy_cpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos, proj);
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
        handler.recv.the_data = FabArray<FAB>::PostRcvs(*cmd.m_RcvTags, handler.recv.data, handler.recv.size,
                                            handler.recv.rank, handler.recv.request, ncomp, SeqNum);
    }

    if (N_snds > 0) {
        handler.send.the_data =
            FabArray<FAB>::PrepareSendBuffers(*cmd.m_SndTags, handler.send.data, handler.send.size,
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
            local_copy_gpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos, proj);
        } else
#endif
        {
            local_copy_cpu(mf, mf, scomp, scomp, ncomp, *cmd.m_LocTags, dtos, proj);
        }
    }

    return handler;
#endif
}

#ifdef AMREX_USE_MPI
template <class FAB, class DTOS, class Proj>
void
Comm_finish (FabArray<FAB>& mf, int scomp, int ncomp, FabArrayBase::CommMetaData const& cmd,
             CommHandler handler, DTOS const& dtos, Proj const& proj)
{
    if (ParallelContext::NProcsSub() == 1) { return; }

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
                                   handler.recv.size, handler.recv.cctc, dtos, proj);
        } else
#endif
        {
            unpack_recv_buffer_cpu(mf, scomp, ncomp, handler.recv.data,
                                   handler.recv.size, handler.recv.cctc, dtos, proj);
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

    if (nghost[0] <= 0) { return; }

    const FabArrayBase::RB90& TheRB90 = mf.getRB90(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, TheRB90,Rotate90DstToSrc{},
                               Identity{});

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
    Comm_finish(mf, scomp, ncomp, TheRB90, std::move(handler), Rotate90DstToSrc{},
                Identity{});
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

    if (nghost[0] <= 0) { return; }

    const FabArrayBase::RB180& TheRB180 = mf.getRB180(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, TheRB180,
                               Rotate180Fn{domain.length(1)}, Identity{});

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, TheRB180, std::move(handler),
                Rotate180Fn{domain.length(1)}, Identity{});
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

    if (nghost[0] <= 0) { return; }

    const FabArrayBase::PolarB& ThePolarB = mf.getPolarB(nghost, domain);

    auto handler = Comm_nowait(mf, scomp, ncomp, ThePolarB,
                               PolarFn{domain.length(0), domain.length(1)},
                               Identity{});

#ifdef AMREX_USE_MPI
    Comm_finish(mf, scomp, ncomp, ThePolarB, std::move(handler),
                PolarFn{domain.length(0), domain.length(1)}, Identity{});
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

template <typename FAB, typename DTOS, typename Proj>
std::enable_if_t<IsBaseFab<FAB>() &&
                 IsCallableR<Dim3,DTOS,Dim3>() &&
                 IsFabProjection<Proj,FAB>(),
                 CommHandler>
FillBoundary_nowait (FabArray<FAB>& mf, const FabArrayBase::CommMetaData& cmd,
                     int scomp, int ncomp, DTOS const& dtos, Proj const& proj)
{
    BL_PROFILE("FillBoundary_nowait(cmd)");
    AMREX_ASSERT(scomp < mf.nComp() && scomp+ncomp <= mf.nComp());
    return Comm_nowait(mf, scomp, ncomp, cmd, dtos, proj);
}

template <typename FAB, typename DTOS, typename Proj>
std::enable_if_t<IsBaseFab<FAB>() &&
                 IsCallableR<Dim3,DTOS,Dim3>() &&
                 IsFabProjection<Proj,FAB>()>
FillBoundary_finish (CommHandler handler,
                     FabArray<FAB>& mf, const FabArrayBase::CommMetaData& cmd,
                     int scomp, int ncomp, DTOS const& dtos, Proj const& proj)
{
#ifdef AMREX_USE_MPI
    BL_PROFILE("FillBoundary_finish(cmd)");
    Comm_finish(mf, scomp, ncomp, cmd, std::move(handler), dtos, proj);
#else
    amrex::ignore_unused(handler,mf,cmd,scomp,ncomp,dtos,proj);
#endif
}

template <typename DTOS>
Vector<std::pair<Box,Box>>
get_src_dst_boxes (DTOS const& dtos, Box const& dstbox, Box const& domain)
{
    Vector<std::pair<Box,Box>> r;
    IntVect mapped_smallend(dtos(amrex::lbound(dstbox)));
    IntVect mapped_bigend  (dtos(amrex::ubound(dstbox)));
    if (!domain.contains(mapped_smallend) || !domain.contains(mapped_bigend)) {
        return r;
    }

    auto sign = dtos.sign(amrex::lbound(dstbox));
    auto perm = dtos.permutation(amrex::lbound(dstbox));
    auto dtype = dstbox.type();
    IntVect stype{AMREX_D_DECL(dtype[perm[0]],
                               dtype[perm[1]],
                               dtype[perm[2]])};
    Array<Array<std::pair<int,int>,2>,AMREX_SPACEDIM> ends;
    Array<Array<std::pair<int,int>,2>,AMREX_SPACEDIM> dst_ends;
    Array<int,AMREX_SPACEDIM> nboxes;
    for (int ddim = 0; ddim < AMREX_SPACEDIM; ++ddim) {
        int sdim = perm[ddim];
        auto mm = std::minmax(mapped_smallend[sdim],mapped_bigend[sdim]);
        if (((sign[ddim] > 0) && (mapped_smallend[sdim] <= mapped_bigend[sdim])) ||
            ((sign[ddim] < 0) && (mapped_bigend[sdim] <= mapped_smallend[sdim])))
        {
            nboxes[sdim] = 1;
            ends[sdim][0] = mm;
            dst_ends[ddim][0] = std::make_pair(dstbox.smallEnd(ddim),
                                               dstbox.bigEnd(ddim));
        } else {
            nboxes[sdim] = 2;
            ends[sdim][0].first  = domain.smallEnd(sdim);
            ends[sdim][0].second = mm.first;
            ends[sdim][1].first  = mm.second;
            ends[sdim][1].second = domain.bigEnd(sdim);
            int n0 = ends[sdim][0].second - ends[sdim][0].first;
            int n1 = ends[sdim][1].second - ends[sdim][1].first;
            if (mm.first == mapped_smallend[sdim]) {
                dst_ends[ddim][0] = std::make_pair(dstbox.smallEnd(ddim),
                                                   dstbox.smallEnd(ddim)+n0);
                dst_ends[ddim][1] = std::make_pair(dstbox.bigEnd(ddim)-n1,
                                                   dstbox.bigEnd(ddim));
            } else {
                dst_ends[ddim][0] = std::make_pair(dstbox.bigEnd(ddim)-n0,
                                                   dstbox.bigEnd(ddim));
                dst_ends[ddim][1] = std::make_pair(dstbox.smallEnd(ddim),
                                                   dstbox.smallEnd(ddim)+n1);
            }
        }
    }

    r.reserve(AMREX_D_TERM(nboxes[0],*nboxes[1],*nboxes[2]));

#if (AMREX_SPACEDIM == 3)
    for (int kbox = 0; kbox < nboxes[2]; ++kbox) {
#endif
#if (AMREX_SPACEDIM >=2 )
    for (int jbox = 0; jbox < nboxes[1]; ++jbox) {
#endif
    for (int ibox = 0; ibox < nboxes[0]; ++ibox)
    {
        IntVect siv(AMREX_D_DECL(ibox,jbox,kbox));
        IntVect div(AMREX_D_DECL(siv[perm[0]],siv[perm[1]],siv[perm[2]]));
        r.emplace_back(Box(IntVect(AMREX_D_DECL(ends[0][ibox].first,
                                                ends[1][jbox].first,
                                                ends[2][kbox].first)),
                           IntVect(AMREX_D_DECL(ends[0][ibox].second,
                                                ends[1][jbox].second,
                                                ends[2][kbox].second)),
                           stype),
                       Box(IntVect(AMREX_D_DECL(dst_ends[0][div[0]].first,
                                                dst_ends[1][div[1]].first,
                                                dst_ends[2][div[2]].first)),
                           IntVect(AMREX_D_DECL(dst_ends[0][div[0]].second,
                                                dst_ends[1][div[1]].second,
                                                dst_ends[2][div[2]].second)),
                           dtype));
    AMREX_D_TERM(},},})

    return r;
}

template <typename DTOS>
Box get_dst_subbox (DTOS const& dtos, std::pair<Box,Box> const& sdboxes,
                    Box const& srcsubbox)
{
    Box const& srcbox = sdboxes.first;
    Box const& dstbox = sdboxes.second;
    if (srcbox == srcsubbox) {
        return dstbox;
    } else {
        auto sign = dtos.sign(amrex::lbound(dstbox));
        auto perm = dtos.permutation(amrex::lbound(dstbox));
        Box dstsubbox = dstbox;
        for (int ddim = 0; ddim < AMREX_SPACEDIM; ++ddim) {
            int sdim = perm[ddim];
            if (sign[ddim] > 0) {
                dstsubbox.growLo(ddim, srcbox.smallEnd(sdim)-srcsubbox.smallEnd(sdim));
                dstsubbox.growHi(ddim, srcsubbox.bigEnd(sdim)-srcbox.bigEnd(sdim));
            } else {
                dstsubbox.growLo(ddim, srcsubbox.bigEnd(sdim)-srcbox.bigEnd(sdim));
                dstsubbox.growHi(ddim, srcbox.smallEnd(sdim)-srcsubbox.smallEnd(sdim));
            }
        }
        return dstsubbox;
    }
}

namespace detail {
    void split_boxes (BoxList& bl, Box const& domain);
}

template <typename FAB, typename DTOS>
std::enable_if_t<IsBaseFab<FAB>() && IsCallableR<Dim3,DTOS,Dim3>(),
                 FabArrayBase::CommMetaData>
makeFillBoundaryMetaData (FabArray<FAB>& mf, IntVect const& nghost,
                          Geometry const& geom, DTOS const& dtos)
{
    FabArrayBase::CommMetaData cmd;
    cmd.m_LocTags = std::make_unique<FabArrayBase::CopyComTagsContainer>();
    cmd.m_SndTags = std::make_unique<FabArrayBase::MapOfCopyComTagContainers>();
    cmd.m_RcvTags = std::make_unique<FabArrayBase::MapOfCopyComTagContainers>();

    // Normal FillBoundary part
    mf.define_fb_metadata(cmd, nghost, false, geom.periodicity(), false);

    BoxArray const& ba = mf.boxArray();
    DistributionMapping const& dm = mf.DistributionMap();
    Box dombox = amrex::convert(geom.Domain(), ba.ixType());
    Box pdombox = amrex::convert(geom.growPeriodicDomain(nghost), ba.ixType());

    const int myproc = ParallelDescriptor::MyProc();
    const auto nboxes = static_cast<int>(ba.size());
    std::vector<std::pair<int,Box> > isects;

    for (int i = 0; i < nboxes; ++i) {
        Box const& gbx = amrex::grow(ba[i], nghost);
        BoxList bl = amrex::boxDiff(gbx, pdombox);
        if (bl.isEmpty()) { continue; }

        detail::split_boxes(bl, dombox);

        const int dst_owner = dm[i];
        for (auto const& dst_box : bl) {
            auto const& src_dst_boxes = get_src_dst_boxes(dtos, dst_box, dombox);
            for (auto const& sd_box_pair : src_dst_boxes) {
                ba.intersections(sd_box_pair.first, isects);
                for (auto const& is : isects) {
                    int const     k = is.first;
                    Box const src_b = is.second;
                    int const src_owner = dm[k];
                    if (dst_owner == myproc || src_owner == myproc) {
                        Box const& dst_b = get_dst_subbox(dtos, sd_box_pair, src_b);
                        if (src_owner == dst_owner) {
                            cmd.m_LocTags->emplace_back(dst_b, src_b, i, k);
                        } else {
                            auto& tags = (dst_owner == myproc) ?
                                (*cmd.m_RcvTags)[src_owner] :
                                (*cmd.m_SndTags)[dst_owner];
                            tags.emplace_back(dst_b, src_b, i, k);
                        }
                    }
                }
            }
        }
    }

    return cmd;
}

struct SphThetaPhiRIndexMapping
{
    SphThetaPhiRIndexMapping (Box const& a_domain)
        : nx(a_domain.length(0)),
          ny(a_domain.length(1)),
          nz(a_domain.length(2))
    {
        AMREX_ASSERT(a_domain.smallEnd() == 0);
    }

    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    Dim3 operator() (Dim3 const& ijk) const noexcept
    {
        const int i = ijk.x;
        const int j = ijk.y;
        const int k = ijk.z;
        bool ilo = i < 0;
        bool ihi = i >= nx;
        bool imd = i >= 0 && i < nx;
        bool jlo = j < 0;
        bool jhi = j >= ny;
        bool jmd = j >= 0 && j < ny;
        bool klo = k < 0;
        bool kmd = k >= 0 && k < nz;
        // We do not need to do anything at the theta-lo/r-lo edge,
        // theta-hi/r-lo edge, and r > r-hi.
        if (ilo && jmd && kmd)
        {
            return {-1-i, (j+ny/2)%ny, k};
        }
        else if (ihi && jmd && kmd)
        {
            return {2*nx-1-i, (j+ny/2)%ny, k};
        }
        else if (imd && jlo && kmd)
        {
            return {i, j+ny, k};
        }
        else if (imd && jhi && kmd)
        {
            return {i, j-ny, k};
        }
        else if (imd && jmd & klo)
        {
            return {nx-1-i, (j+ny/2)%ny, -1-k};
        }
        else if (ilo && jlo && kmd)
        {
            return {-1-i, (j+ny/2)%ny, k};
        }
        else if (ihi && jlo && kmd)
        {
            return {2*nx-1-i, (j+ny/2)%ny, k};
        }
        else if (ilo && jhi && kmd)
        {
            return {-1-i, (j+ny/2)%ny, k};
        }
        else if (ihi && jhi && kmd)
        {
            return {2*nx-1-i, (j+ny/2)%ny, k};
        }
        else if (imd && jlo && klo)
        {
            return {nx-1-i, (j+ny/2)%ny, -1-k};
        }
        else if (imd && jhi && klo)
        {
            return {nx-1-i, (j+ny/2)%ny, -1-k};
        }
        else
        {
            return ijk;
        }
    }

    [[nodiscard]] IntVect sign (Dim3 const& ijk) const noexcept
    {
        if (ijk.z < 0) {
            return IntVect{AMREX_D_DECL(-1, 1,-1)};
        } else if (ijk.z >=0 && ijk.z < nz &&
                   (ijk.x < 0 || ijk.x >= nx)) {
            return IntVect{AMREX_D_DECL(-1, 1, 1)};
        } else {
            return IntVect{AMREX_D_DECL( 1, 1, 1)};
        }
    }

    [[nodiscard]] IntVect permutation (Dim3 const& /*ijk*/) const noexcept // NOLINT(readability-convert-member-functions-to-static)
    {
        return IntVect(AMREX_D_DECL(0,1,2));
    }

private:
    int nx, ny, nz;
};

struct SphThetaPhiRComponentMapping
{
    SphThetaPhiRComponentMapping (Box const& a_domain, int a_start_index)
        : nx(a_domain.length(0)),
          ny(a_domain.length(1)),
          nz(a_domain.length(2)),
          scomp(a_start_index) {}

    template <typename T>
    [[nodiscard]] AMREX_GPU_HOST_DEVICE
    T operator()(Array4<const T> const& a, Dim3 const& ijk, int n) const noexcept
    {
        const int i = ijk.x;
        const int j = ijk.y;
        const int k = ijk.z;
        auto r = a(i,j,k,n);
        if (n == scomp) {
            if ((i >= 0 && i <  nx) &&
                (j <  0 || j >= ny) &&
                (k >= 0 && k <  nz)) {
                return r;
            } else {
                // We do not need to worry about the theta-lo/r-lo edge,
                // theta-hi/r-lo edge, and r > r-hi.
                return -r;
            }
        } else if (n == scomp+2) {
            if (k < 0) {
                return -r;
            } else {
                return r;
            }
        } else {
            return r;
        }
    }
private:
    int nx, ny, nz;
    int scomp;
};

extern template MultiBlockCommMetaData ParallelCopy(FabArray<FArrayBox>& dest, const Box& destbox,
                                                    const FabArray<FArrayBox>& src, int destcomp,
                                                    int srccomp, int numcomp, const IntVect& ngrow,
                                                    MultiBlockIndexMapping const&, Identity const&);
}

#endif
