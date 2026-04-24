#ifndef AMREX_NEIGHBORPARTICLESCPUIMPL_H_
#define AMREX_NEIGHBORPARTICLESCPUIMPL_H_
#include <AMReX_Config.H>

#include <array>
#include <cstring>

//! \cond CPU

namespace amrex {

namespace detail {

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void applyNeighborPeriodicShift (T& value, int dir, const IntVect& periodic_shift,
                                 const Periodicity& periodicity,
                                 const RealBox& prob_domain) noexcept
{
    if (!periodicity.isPeriodic(dir)) { return; }
    if (periodic_shift[dir] < 0) {
        value += static_cast<T>(prob_domain.length(dir));
    } else if (periodic_shift[dir] > 0) {
        value -= static_cast<T>(prob_domain.length(dir));
    }
}

template <typename PTD>
AMREX_FORCE_INLINE
void shiftNeighborParticlePositions (PTD const& ptd, int index, const IntVect& periodic_shift,
                                     const Periodicity& periodicity,
                                     const RealBox& prob_domain) noexcept
{
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        applyNeighborPeriodicShift(ptd.pos(dir, index), dir, periodic_shift, periodicity, prob_domain);
    }
}

template <typename PTD>
AMREX_FORCE_INLINE
void packNeighborParticleData (char*& dst, PTD const& src, int src_index,
                               int num_real_comps, int num_int_comps,
                               const Vector<int>& ghost_real_comp,
                               const Vector<int>& ghost_int_comp,
                               const IntVect& periodic_shift,
                               const Periodicity& periodicity,
                               const RealBox& prob_domain) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;
    constexpr int real_comm_comp_start =
        ParticleType::is_soa_particle ? 0 : AMREX_SPACEDIM + ParticleType::NReal;
    constexpr int int_comm_comp_start = 2 + ParticleType::NInt;

    if constexpr (!ParticleType::is_soa_particle) {
        auto const p = src[src_index];
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (!ghost_real_comp[dir]) { continue; }
            auto value = p.pos(dir);
            applyNeighborPeriodicShift(value, dir, periodic_shift, periodicity, prob_domain);
            std::memcpy(dst, &value, sizeof(typename ParticleType::RealType));
            dst += sizeof(typename ParticleType::RealType);
        }
        for (int ii = 0; ii < ParticleType::NReal; ++ii) {
            if (!ghost_real_comp[AMREX_SPACEDIM + ii]) { continue; }
            auto value = p.rdata(ii);
            std::memcpy(dst, &value, sizeof(typename ParticleType::RealType));
            dst += sizeof(typename ParticleType::RealType);
        }
    }

    for (int ii = 0; ii < num_real_comps; ++ii) {
        if (!ghost_real_comp[real_comm_comp_start + ii]) { continue; }

        ParticleReal value;
        if (ii < std::decay_t<PTD>::NAR) {
            value = src.rdata(ii)[src_index];
        } else {
            value = src.m_runtime_rdata[ii - std::decay_t<PTD>::NAR][src_index];
        }
        if constexpr (ParticleType::is_soa_particle) {
            if (ii < AMREX_SPACEDIM) {
                applyNeighborPeriodicShift(value, ii, periodic_shift, periodicity, prob_domain);
            }
        }
        std::memcpy(dst, &value, sizeof(ParticleReal));
        dst += sizeof(ParticleReal);
    }

    std::array<int,2> idcpu_parts{};
    auto const idcpu = src.idcpu(src_index);
    std::memcpy(idcpu_parts.data(), &idcpu, sizeof(idcpu));
    for (int ii = 0; ii < 2; ++ii) {
        if (!ghost_int_comp[ii]) { continue; }
        std::memcpy(dst, &idcpu_parts[ii], sizeof(int));
        dst += sizeof(int);
    }

    if constexpr (!ParticleType::is_soa_particle) {
        auto const p = src[src_index];
        for (int ii = 0; ii < ParticleType::NInt; ++ii) {
            if (!ghost_int_comp[2 + ii]) { continue; }
            auto value = p.idata(ii);
            std::memcpy(dst, &value, sizeof(int));
            dst += sizeof(int);
        }
    }

    for (int ii = 0; ii < num_int_comps; ++ii) {
        if (!ghost_int_comp[int_comm_comp_start + ii]) { continue; }

        int value;
        if (ii < std::decay_t<PTD>::NAI) {
            value = src.idata(ii)[src_index];
        } else {
            value = src.m_runtime_idata[ii - std::decay_t<PTD>::NAI][src_index];
        }
        std::memcpy(dst, &value, sizeof(int));
        dst += sizeof(int);
    }
}

template <typename PTD>
AMREX_FORCE_INLINE
void unpackNeighborParticleData (PTD const& dst_ptd, int dst_index, const char*& src,
                                 int num_real_comps, int num_int_comps,
                                 const Vector<int>& ghost_real_comp,
                                 const Vector<int>& ghost_int_comp) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;
    constexpr int real_comm_comp_start =
        ParticleType::is_soa_particle ? 0 : AMREX_SPACEDIM + ParticleType::NReal;
    constexpr int int_comm_comp_start = 2 + ParticleType::NInt;

    if constexpr (!ParticleType::is_soa_particle) {
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            if (!ghost_real_comp[dir]) { continue; }
            std::memcpy(&(dst_ptd.pos(dir, dst_index)), src, sizeof(typename ParticleType::RealType));
            src += sizeof(typename ParticleType::RealType);
        }
        for (int ii = 0; ii < ParticleType::NReal; ++ii) {
            if (!ghost_real_comp[AMREX_SPACEDIM + ii]) { continue; }
            std::memcpy(&(dst_ptd[dst_index].rdata(ii)), src, sizeof(typename ParticleType::RealType));
            src += sizeof(typename ParticleType::RealType);
        }
    }

    for (int ii = 0; ii < num_real_comps; ++ii) {
        if (!ghost_real_comp[real_comm_comp_start + ii]) { continue; }

        if (ii < std::decay_t<PTD>::NAR) {
            std::memcpy(dst_ptd.rdata(ii) + dst_index, src, sizeof(ParticleReal));
        } else {
            std::memcpy(dst_ptd.m_runtime_rdata[ii - std::decay_t<PTD>::NAR] + dst_index, src,
                        sizeof(ParticleReal));
        }
        src += sizeof(ParticleReal);
    }

    std::array<int,2> idcpu_parts{};
    std::memcpy(idcpu_parts.data(), &(dst_ptd.idcpu(dst_index)), sizeof(uint64_t));
    for (int ii = 0; ii < 2; ++ii) {
        if (!ghost_int_comp[ii]) { continue; }
        std::memcpy(&idcpu_parts[ii], src, sizeof(int));
        src += sizeof(int);
    }
    std::memcpy(&(dst_ptd.idcpu(dst_index)), idcpu_parts.data(), sizeof(uint64_t));

    if constexpr (!ParticleType::is_soa_particle) {
        for (int ii = 0; ii < ParticleType::NInt; ++ii) {
            if (!ghost_int_comp[2 + ii]) { continue; }
            std::memcpy(&(dst_ptd[dst_index].idata(ii)), src, sizeof(int));
            src += sizeof(int);
        }
    }

    for (int ii = 0; ii < num_int_comps; ++ii) {
        if (!ghost_int_comp[int_comm_comp_start + ii]) { continue; }

        if (ii < std::decay_t<PTD>::NAI) {
            std::memcpy(dst_ptd.idata(ii) + dst_index, src, sizeof(int));
        } else {
            std::memcpy(dst_ptd.m_runtime_idata[ii - std::decay_t<PTD>::NAI] + dst_index, src,
                        sizeof(int));
        }
        src += sizeof(int);
    }
}

template <typename PTD>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ParticleReal getSummedRealComp (PTD const& ptd, int index, int comp) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;

    if constexpr (!ParticleType::is_soa_particle) {
        if (comp < ParticleType::NReal) {
            return ptd[index].rdata(comp);
        }
        comp -= ParticleType::NReal;
    }

    if (comp < std::decay_t<PTD>::NAR) {
        return ptd.rdata(comp)[index];
    }

    return ptd.m_runtime_rdata[comp - std::decay_t<PTD>::NAR][index];
}

template <typename PTD>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int getSummedIntComp (PTD const& ptd, int index, int comp) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;

    if constexpr (!ParticleType::is_soa_particle) {
        if (comp < ParticleType::NInt) {
            return ptd[index].idata(comp);
        }
        comp -= ParticleType::NInt;
    }

    if (comp < std::decay_t<PTD>::NAI) {
        return ptd.idata(comp)[index];
    }

    return ptd.m_runtime_idata[comp - std::decay_t<PTD>::NAI][index];
}

template <typename PTD>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void addSummedRealComp (PTD const& ptd, int index, int comp, ParticleReal value) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;

    if constexpr (!ParticleType::is_soa_particle) {
        if (comp < ParticleType::NReal) {
            ptd[index].rdata(comp) += value;
            return;
        }
        comp -= ParticleType::NReal;
    }

    if (comp < std::decay_t<PTD>::NAR) {
        ptd.rdata(comp)[index] += value;
        return;
    }

    ptd.m_runtime_rdata[comp - std::decay_t<PTD>::NAR][index] += value;
}

template <typename PTD>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void addSummedIntComp (PTD const& ptd, int index, int comp, int value) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;

    if constexpr (!ParticleType::is_soa_particle) {
        if (comp < ParticleType::NInt) {
            ptd[index].idata(comp) += value;
            return;
        }
        comp -= ParticleType::NInt;
    }

    if (comp < std::decay_t<PTD>::NAI) {
        ptd.idata(comp)[index] += value;
        return;
    }

    ptd.m_runtime_idata[comp - std::decay_t<PTD>::NAI][index] += value;
}

}

template <typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
NeighborParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::fillNeighborsCPU () {
    BL_PROFILE("NeighborParticleContainer::fillNeighborsCPU");
    if (!areMasksValid()) {
        BuildMasks();
        GetNeighborCommTags();
    }
    cacheNeighborInfo();
    updateNeighborsCPU(false);
}

template <typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
NeighborParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::sumNeighborsCPU (int real_start_comp, int real_num_comp,
                   int int_start_comp,  int int_num_comp)
{
    BL_PROFILE("NeighborParticleContainer::sumNeighborsCPU");

    if ( ! enableInverse() )
    {
        amrex::Abort("Need to enable inverse to true to use sumNeighbors. \n");
    }

    const int MyProc = ParallelContext::MyProcSub();

    std::map<int, Vector<char> > isend_data;

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
        {
            PairIndex src_index(pti.index(), pti.LocalTileIndex());
            const auto& tags = inverse_tags[lev][src_index];
            const auto& neighb_tile = neighbors[lev][src_index];
            const auto neighbs = neighb_tile.getConstParticleTileData();
            const Long num_neighbs = static_cast<Long>(neighb_tile.size());
            AMREX_ASSERT(tags.size() == num_neighbs);

            for (int i = 0; i < static_cast<int>(num_neighbs); ++i)
            {
                const auto& tag = tags[i];
                const int dst_grid = tag.src_grid;
                const int global_rank = this->ParticleDistributionMap(lev)[dst_grid];
                const int dst_proc = ParallelContext::global_to_local_rank(global_rank);
                const int dst_tile = tag.src_tile;
                const int dst_index = tag.src_index;
                const int dst_level = tag.src_level;

                if (dst_proc == MyProc)
                {
                    auto pair = std::make_pair(dst_grid, dst_tile);
                    auto& dst_ptile = this->GetParticles(dst_level)[pair];
                    auto dst_ptd = dst_ptile.getParticleTileData();

                    for (int comp = real_start_comp; comp < real_start_comp + real_num_comp; ++comp)
                    {
                        detail::addSummedRealComp(
                            dst_ptd, dst_index, comp, detail::getSummedRealComp(neighbs, i, comp));
                    }

                    for (int comp = int_start_comp; comp < int_start_comp + int_num_comp; ++comp)
                    {
                        detail::addSummedIntComp(
                            dst_ptd, dst_index, comp, detail::getSummedIntComp(neighbs, i, comp));
                    }
                }

                else
                {
                    auto& sdata = isend_data[dst_proc];
                    auto old_size = sdata.size();
                    auto new_size = old_size + real_num_comp*sizeof(ParticleReal) + int_num_comp*sizeof(int)
                        + 4*sizeof(int);
                    sdata.resize(new_size);
                    char* dst = &sdata[old_size];
                    std::memcpy(dst, &dst_grid, sizeof(int)); dst += sizeof(int);
                    std::memcpy(dst, &dst_tile, sizeof(int)); dst += sizeof(int);
                    std::memcpy(dst, &dst_index, sizeof(int)); dst += sizeof(int);
                    std::memcpy(dst, &dst_level, sizeof(int)); dst += sizeof(int);
                    for (int comp = real_start_comp; comp < real_start_comp + real_num_comp; ++comp)
                    {
                        ParticleReal data = detail::getSummedRealComp(neighbs, i, comp);
                        std::memcpy(dst, &data, sizeof(ParticleReal));
                        dst += sizeof(ParticleReal);
                    }
                    for (int comp = int_start_comp; comp < int_start_comp + int_num_comp; ++comp)
                    {
                        int data = detail::getSummedIntComp(neighbs, i, comp);
                        std::memcpy(dst, &data, sizeof(int));
                        dst += sizeof(int);
                    }
                }
            }
        }
    }

    sumNeighborsMPI(isend_data, real_start_comp, real_num_comp, int_start_comp, int_num_comp);
}

template <typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
NeighborParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
sumNeighborsMPI (std::map<int, Vector<char> >& not_ours,
                 int real_start_comp, int real_num_comp,
                 int int_start_comp, int int_num_comp)
{
    BL_PROFILE("NeighborParticleContainer::sumNeighborsMPI");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();

    Vector<Long> isnds(NProcs, 0);
    Vector<Long> ircvs(NProcs, 0);
    for (int i = 0; i < NProcs; ++i) {
        ircvs[i] = 0;
    }

    {
        // each proc figures out how many bytes it will send, and how
        // many it will receive

        Long num_isnds = 0;
        for (const auto& kv : not_ours)
        {
            num_isnds      += kv.second.size();
            isnds[kv.first] = kv.second.size();
        }

        ParallelAllReduce::Max(num_isnds, ParallelContext::CommunicatorSub());

        if (num_isnds == 0) { return; }

        const int num_ircvs = neighbor_procs.size();
        Vector<MPI_Status>  stats(num_ircvs);
        Vector<MPI_Request> rreqs(num_ircvs);

        const int SeqNum = ParallelDescriptor::SeqNum();

        // Post receives
        for (int i = 0; i < num_ircvs; ++i)
        {
            const int Who = neighbor_procs[i];
            const Long Cnt = 1;

            AMREX_ASSERT(Who >= 0 && Who < NProcs);

            rreqs[i] = ParallelDescriptor::Arecv(&ircvs[Who], Cnt, Who, SeqNum,
                                                 ParallelContext::CommunicatorSub()).req();
        }

        // Send.
        for (int i = 0; i < num_ircvs; ++i) {
        const int Who = neighbor_procs[i];
        const Long Cnt = 1;

        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        ParallelDescriptor::Send(&isnds[Who], Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
        }

        if (num_ircvs > 0) { ParallelDescriptor::Waitall(rreqs, stats); }
    }

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (ircvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += ircvs[i];
        }
    }

    const auto nrcvs = int(RcvProc.size());
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();

    // Allocate data for rcvs as one big chunk.
    Vector<char> recvdata(TotRcvBytes);

    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = ircvs[Who];

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum,
                                             ParallelContext::CommunicatorSub()).req();
    }

    // Send.
    for (const auto& kv : not_ours) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
    }

    // unpack the received data and put them into the proper neighbor buffers
    if (nrcvs > 0)
    {
        ParallelDescriptor::Waitall(rreqs, stats);

        const size_t data_size =
            real_num_comp*sizeof(ParticleReal) + int_num_comp*sizeof(int) + 4 * sizeof(int);

        if (recvdata.size() % data_size != 0) {
            amrex::Print() << recvdata.size() << " " << data_size << "\n";
            if (this->m_verbose) {
                amrex::AllPrint() << "NeighborParticles::sumNeighbors: sizes = "
                                  << recvdata.size() << ", " << data_size << "\n";
            }
            amrex::Abort("NeighborParticles::sumNeighbors: How did this happen?");
        }

        auto npart = int(recvdata.size() / data_size);

        char* buffer = recvdata.data();
        for (int j = 0; j < npart; ++j)
        {
            int grid, tile, index, lev;
            std::memcpy(&grid,  buffer, sizeof(int)); buffer += sizeof(int);
            std::memcpy(&tile,  buffer, sizeof(int)); buffer += sizeof(int);
            std::memcpy(&index,  buffer, sizeof(int)); buffer += sizeof(int);
            std::memcpy(&lev, buffer, sizeof(int)); buffer += sizeof(int);

            auto pair = std::make_pair(grid, tile);
            auto& ptile = this->GetParticles(lev)[pair];
            auto dst_ptd = ptile.getParticleTileData();

            for (int comp = real_start_comp; comp < real_start_comp + real_num_comp; ++comp)
            {
                ParticleReal data;
                std::memcpy(&data, buffer, sizeof(ParticleReal));
                detail::addSummedRealComp(dst_ptd, index, comp, data);
                buffer += sizeof(ParticleReal);
            }

            for (int comp = int_start_comp; comp < int_start_comp + int_num_comp; ++comp)
            {
                int data;
                std::memcpy(&data, buffer, sizeof(int));
                detail::addSummedIntComp(dst_ptd, index, comp, data);
                buffer += sizeof(int);
            }
        }
    }
#else
    amrex::ignore_unused(not_ours, real_start_comp, real_num_comp, int_start_comp, int_num_comp);
#endif
}

template <typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
NeighborParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::updateNeighborsCPU (bool reuse_rcv_counts) {

    BL_PROFILE_VAR("NeighborParticleContainer::updateNeighborsCPU", update);

    const int MyProc = ParallelContext::MyProcSub();

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        const Periodicity& periodicity = this->Geom(lev).periodicity();
        const RealBox& prob_domain = this->Geom(lev).ProbDomain();

        int num_threads = OpenMP::get_max_threads();

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex src_index(pti.index(), pti.LocalTileIndex());
            auto const src = pti.GetParticleTile().getConstParticleTileData();
            for (int j = 0; j < num_threads; ++j) {
                auto& tags = buffer_tag_cache[lev][src_index][j];
                int num_tags = tags.size();
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
                for (int i = 0; i < num_tags; ++i) {
                    const NeighborCopyTag& tag = tags[i];
                    const int global_who = this->ParticleDistributionMap(tag.level)[tag.grid];
                    const int who = ParallelContext::global_to_local_rank(global_who);
                    if (who == MyProc) {
                        PairIndex dst_index(tag.grid, tag.tile);
                        auto dst = neighbors[tag.level][dst_index].getParticleTileData();
                        copyParticle(dst, src, tag.src_index, tag.dst_index);
                        if (periodicity.isAnyPeriodic()) {
                            detail::shiftNeighborParticlePositions(dst, tag.dst_index, tag.periodic_shift,
                                                                  periodicity, prob_domain);
                        }

                        if ( enableInverse() )
                        {
                            auto& itags = inverse_tags[tag.level][dst_index];
                            AMREX_ASSERT(tag.dst_index < itags.size());
                            itags[tag.dst_index].src_grid = src_index.first;
                            itags[tag.dst_index].src_tile = src_index.second;
                            itags[tag.dst_index].src_index = tag.src_index;
                            itags[tag.dst_index].src_level = lev;
                        }
                    } else {
                        char* dst_ptr = &send_data[who][tag.dst_index];
                        detail::packNeighborParticleData(dst_ptr, src, tag.src_index,
                                                         this->NumRealComps(), this->NumIntComps(),
                                                         ghost_real_comp, ghost_int_comp,
                                                         tag.periodic_shift, periodicity, prob_domain);
                        if ( enableInverse() )
                        {
                            std::memcpy(dst_ptr,&(src_index.first),sizeof(int)); dst_ptr += sizeof(int);
                            std::memcpy(dst_ptr,&(src_index.second),sizeof(int)); dst_ptr += sizeof(int);
                            std::memcpy(dst_ptr,&(tag.src_index),sizeof(int)); dst_ptr += sizeof(int);
                            std::memcpy(dst_ptr,&(lev),sizeof(int)); dst_ptr += sizeof(int);
                        }
                    }
                }
            }
        }

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
            const int grid = mfi.index();
            const int tile = mfi.LocalTileIndex();
            PairIndex dst_index(grid, tile);
            neighbors[lev][dst_index].resize(local_neighbor_sizes[lev][dst_index]);
            if ( enableInverse() ) {
               inverse_tags[lev][dst_index].resize(local_neighbor_sizes[lev][dst_index]);
            }
        }
    }
    BL_PROFILE_VAR_STOP(update);

    fillNeighborsMPI(reuse_rcv_counts);

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int src_grid = mfi.index();
            int src_tile = mfi.LocalTileIndex();
            auto index = std::make_pair(src_grid, src_tile);
            auto& ptile = this->GetParticles(lev)[index];
            ptile.setNumNeighbors(neighbors[lev][index].size());
            amrex::copyParticles(ptile, neighbors[lev][index], 0,
                                 ptile.numRealParticles(), ptile.numNeighborParticles());
        }
    }

}

template <typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
NeighborParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>
::clearNeighborsCPU ()
{
    BL_PROFILE("NeighborParticleContainer::clearNeighborsCPU");

    resizeContainers(this->numLevels());
    for (int lev = 0; lev < this->numLevels(); ++lev) {
        neighbors[lev].clear();
        if ( enableInverse() ) { inverse_tags[lev].clear(); }
        buffer_tag_cache[lev].clear();

        for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int src_grid = mfi.index();
            int src_tile = mfi.LocalTileIndex();
            auto index = std::make_pair(src_grid, src_tile);
            auto& ptile = this->GetParticles(lev)[index];
            ptile.setNumNeighbors(0);
        }
    }

    send_data.clear();
}

template <typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
NeighborParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
getRcvCountsMPI () {

    BL_PROFILE("NeighborParticleContainer::getRcvCountsMPI");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();

    // each proc figures out how many bytes it will send, and how
    // many it will receive
    Vector<Long> snds(NProcs, 0);
    rcvs.resize(NProcs);
    for (int i = 0; i < NProcs; ++i) {
        rcvs[i] = 0;
    }

    num_snds = 0;
    for (const auto& kv : send_data) {
        num_snds      += kv.second.size();
        snds[kv.first] = kv.second.size();
    }

    ParallelAllReduce::Max(num_snds, ParallelContext::CommunicatorSub());

    if (num_snds == 0) { return; }

    const int num_rcvs = neighbor_procs.size();
    Vector<MPI_Status>  stats(num_rcvs);
    Vector<MPI_Request> rreqs(num_rcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();

    // Post receives
    for (int i = 0; i < num_rcvs; ++i) {
        const int Who = neighbor_procs[i];
        const Long Cnt = 1;

        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&rcvs[Who], Cnt, Who, SeqNum,
                                             ParallelContext::CommunicatorSub()).req();
    }

    // Send.
    for (int i = 0; i < num_rcvs; ++i) {
        const int Who = neighbor_procs[i];
        const Long Cnt = 1;

        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        ParallelDescriptor::Send(&snds[Who], Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
    }

    if (num_rcvs > 0) { ParallelDescriptor::Waitall(rreqs, stats); }

#endif // AMREX_USE_MPI
}

template <typename T_ParticleType, int NArrayReal, int NArrayInt,
          template<class> class Allocator, class CellAssignor>
void
NeighborParticleContainer_impl<T_ParticleType, NArrayReal, NArrayInt, Allocator, CellAssignor>::
fillNeighborsMPI (bool reuse_rcv_counts) {

    BL_PROFILE("NeighborParticleContainer::fillNeighborsMPI");

#ifdef AMREX_USE_MPI
    const int NProcs = ParallelContext::NProcsSub();

    // each proc figures out how many bytes it will send, and how
    // many it will receive
    if (!reuse_rcv_counts) { getRcvCountsMPI(); }
    if (num_snds == 0) { return; }

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += rcvs[i];
        }
    }

    const auto nrcvs = int(RcvProc.size());
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();

    // Allocate data for rcvs as one big chunk.
    Vector<char> recvdata(TotRcvBytes);

    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = rcvs[Who];

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());
        AMREX_ASSERT(Who >= 0 && Who < NProcs);

        rreqs[i] = ParallelDescriptor::Arecv(&recvdata[offset], Cnt, Who, SeqNum,
                                             ParallelContext::CommunicatorSub()).req();
    }

    // Send.
    for (const auto& kv : send_data) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum,
                                 ParallelContext::CommunicatorSub());
    }

    // unpack the received data and put them into the proper neighbor buffers
    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);
        for (int i = 0; i < nrcvs; ++i) {
            const std::size_t offset = rOffset[i];
            char* buffer = &recvdata[offset];
            int num_tiles, lev, gid, tid, size, np;
            std::memcpy(&num_tiles, buffer, sizeof(int)); buffer += sizeof(int);
            for (int j = 0; j < num_tiles; ++j) {
                std::memcpy(&lev,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&gid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&tid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&size, buffer, sizeof(int)); buffer += sizeof(int);

                if (size == 0) { continue; }

                np = size / cdata_size;

                AMREX_ASSERT(size % cdata_size == 0);

                PairIndex dst_index(gid, tid);
                size_t old_size = neighbors[lev][dst_index].size();
                size_t new_size = neighbors[lev][dst_index].size() + np;
                if ( enableInverse() )
                {
                    AMREX_ASSERT(neighbors[lev][dst_index].size() ==
                                 size_t(inverse_tags[lev][dst_index].size()));
                    inverse_tags[lev][dst_index].resize(new_size);
                }
                neighbors[lev][dst_index].resize(new_size);

                auto dst_ptd = neighbors[lev][dst_index].getParticleTileData();
                char const* src = buffer;
                for (int n = 0; n < np; ++n) {
                    detail::unpackNeighborParticleData(dst_ptd, int(old_size+n), src,
                                                       this->NumRealComps(), this->NumIntComps(),
                                                       ghost_real_comp, ghost_int_comp);

                    if ( enableInverse() )
                    {
                        auto& tag = inverse_tags[lev][dst_index][old_size+n];
                        std::memcpy(&(tag.src_grid),src,sizeof(int));
                        src += sizeof(int);

                        std::memcpy(&(tag.src_tile),src,sizeof(int));
                        src += sizeof(int);

                        std::memcpy(&(tag.src_index),src,sizeof(int));
                        src += sizeof(int);

                        std::memcpy(&(tag.src_level),src,sizeof(int));
                        src += sizeof(int);
                    }
                }
                buffer += size;
            }
        }
    }
#else
    amrex::ignore_unused(reuse_rcv_counts);
#endif
}

}

//! \endcond

#endif
