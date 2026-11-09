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
        value = static_cast<T>(value + prob_domain.length(dir));
    } else if (periodic_shift[dir] > 0) {
        value = static_cast<T>(value - prob_domain.length(dir));
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

/// Same-rank ghost refresh: copy particle si of src to slot di of dst, but only
/// the components in the plan (the ones flagged for ghost communication), and
/// from the slot group only the occupied slots. An occupied slot gets its two
/// keys and its payload. A slot whose two keys are both -1 is empty: its keys
/// are set to -1 in dst if they were not already, and its payload is left alone,
/// since nothing reads the payload of an empty slot before writing it.
template <typename DTD, typename STD>
AMREX_FORCE_INLINE
void copyGhostParticle (DTD const& dst, STD const& src, int si, int di,
                        GhostCopyPlan const& plan) noexcept
{
    using ParticleType = typename std::decay_t<DTD>::ParticleType;
    constexpr int NAR = std::decay_t<DTD>::NAR;
    constexpr int NAI = std::decay_t<DTD>::NAI;

    if constexpr (ParticleType::is_soa_particle) {
        dst.m_idcpu[di] = src.m_idcpu[si];
    } else {
        dst.m_aos[di] = src.m_aos[si];
    }
    if constexpr (NAR > 0) {
        for (int j : plan.ct_real) { dst.m_rdata[j][di] = src.m_rdata[j][si]; }
    }
    for (int j : plan.rt_real) { dst.m_runtime_rdata[j][di] = src.m_runtime_rdata[j][si]; }
    if constexpr (NAI > 0) {
        for (int j : plan.ct_int) { dst.m_idata[j][di] = src.m_idata[j][si]; }
    }
    for (int j : plan.rt_int) { dst.m_runtime_idata[j][di] = src.m_runtime_idata[j][si]; }

    int const ns = plan.nslots;
    for (int s = 0; s < ns; ++s) {
        int const pid = src.m_runtime_idata[plan.key_id + s][si];
        int const cpu = src.m_runtime_idata[plan.key_cpu + s][si];
        int* const dpid = dst.m_runtime_idata[plan.key_id + s] + di;
        int* const dcpu = dst.m_runtime_idata[plan.key_cpu + s] + di;
        if (pid == -1 && cpu == -1) {
            if (*dpid != -1 || *dcpu != -1) {
                *dpid = -1;
                *dcpu = -1;
            }
        } else {
            *dpid = pid;
            *dcpu = cpu;
            for (int b = 0; b < plan.nblocks; ++b) {
                int const jr = plan.real0 + b*ns + s;
                dst.m_runtime_rdata[jr][di] = src.m_runtime_rdata[jr][si];
            }
        }
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
                               const RealBox& prob_domain,
                               GhostCopyPlan const* wire = nullptr) noexcept
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

    // Masked MPI path: the payload of an empty slot (both keys -1) is not
    // written, its bytes in the record stay as they are. The receiver does not
    // read them (see unpackNeighborParticleData).
    unsigned char occ[MaxWireSlots];
    bool const wm = (wire != nullptr) && wire->wire;
    if (wm) {
        for (int s = 0; s < wire->nslots; ++s) {
            int const pid = src.m_runtime_idata[wire->key_id  + s][src_index];
            int const cpu = src.m_runtime_idata[wire->key_cpu + s][src_index];
            occ[s] = static_cast<unsigned char>(pid != -1 || cpu != -1);
        }
    }

    for (int ii = 0; ii < num_real_comps; ++ii) {
        if (!ghost_real_comp[real_comm_comp_start + ii]) { continue; }

        if (wm) {
            int const s = wire->wire_real_slot[ii];
            if (s >= 0 && !occ[s]) {
                dst += sizeof(ParticleReal);
                continue;
            }
        }

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
                                 const Vector<int>& ghost_int_comp,
                                 GhostCopyPlan const* wire = nullptr) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;
    constexpr int real_comm_comp_start =
        ParticleType::is_soa_particle ? 0 : AMREX_SPACEDIM + ParticleType::NReal;
    constexpr int int_comm_comp_start = 2 + ParticleType::NInt;

    // Masked MPI path: the keys of the slots sit at known offsets of the record,
    // so the occupancy is known before the reals are read; the payload of an
    // empty slot is skipped (its bytes were not written by the sender).
    char const* const record = src;
    unsigned char occ[MaxWireSlots];
    bool const wm = (wire != nullptr) && wire->wire;
    if (wm) {
        for (int s = 0; s < wire->nslots; ++s) {
            int pid, cpu;
            std::memcpy(&pid, record + wire->wire_key_id_off[s],  sizeof(int));
            std::memcpy(&cpu, record + wire->wire_key_cpu_off[s], sizeof(int));
            occ[s] = static_cast<unsigned char>(pid != -1 || cpu != -1);
        }
    }

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

        if (wm) {
            int const s = wire->wire_real_slot[ii];
            if (s >= 0 && !occ[s]) {
                src += sizeof(ParticleReal);
                continue;
            }
        }

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

/// Batched packNeighborParticleData for the n remote ghosts of one source tile:
/// the same records, written one component at a time for all of them (rec[t] is
/// the record of ghost t, si[t] its source index, shift[t] its periodic shift).
/// The payload of the empty slots is not written, as in the one-ghost version.
/// Only for the plan of the masked MPI path with plan.batch.
template <typename PTD>
void packNeighborBatch (PTD const& src, int n, char* const* rec, int const* si, IntVect const* shift,
                        int num_real_comps, int num_int_comps,
                        const Periodicity& periodicity, const RealBox& prob_domain,
                        GhostCopyPlan const& pl, GhostBatchScratch& sc) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;
    using RT = typename ParticleType::RealType;
    constexpr int NAR = std::decay_t<PTD>::NAR;
    constexpr int NAI = std::decay_t<PTD>::NAI;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        int const po = pl.wire_pos_off[dir];
        if (po < 0) { continue; }
        for (int t = 0; t < n; ++t) {
            RT value = src[si[t]].pos(dir);
            applyNeighborPeriodicShift(value, dir, shift[t], periodicity, prob_domain);
            std::memcpy(rec[t] + po, &value, sizeof(RT));
        }
    }

    for (int ii = 0; ii < num_real_comps; ++ii) {
        int const off = pl.wire_real_off[ii];
        if (off < 0 || pl.wire_real_slot[ii] >= 0) { continue; }
        ParticleReal const* arr = (ii < NAR) ? src.rdata(ii) : src.m_runtime_rdata[ii - NAR];
        for (int t = 0; t < n; ++t) {
            std::memcpy(rec[t] + off, arr + si[t], sizeof(ParticleReal));
        }
    }

    // ghosts with an occupied slot, per slot
    int const ns = pl.nslots;
    sc.cnt.assign(ns, 0);
    sc.lists.resize(std::size_t(ns) * std::size_t(n));
    for (int s = 0; s < ns; ++s) {
        int const* kid  = src.m_runtime_idata[pl.key_id  + s];
        int const* kcpu = src.m_runtime_idata[pl.key_cpu + s];
        int* const lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
        int c = 0;
        for (int t = 0; t < n; ++t) {
            int const i = si[t];
            lst[c] = t;
            c += (kid[i] != -1 || kcpu[i] != -1) ? 1 : 0;
        }
        sc.cnt[s] = c;
    }
    for (int b = 0; b < pl.nblocks; ++b) {
        for (int s = 0; s < ns; ++s) {
            int const jr = pl.real0 + b*ns + s;
            int const off = pl.wire_real_off[jr + NAR];
            ParticleReal const* arr = src.m_runtime_rdata[jr];
            int const* lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
            int const c = sc.cnt[s];
            for (int k = 0; k < c; ++k) {
                int const t = lst[k];
                std::memcpy(rec[t] + off, arr + si[t], sizeof(ParticleReal));
            }
        }
    }

    if (pl.wire_idcpu_off[0] >= 0 || pl.wire_idcpu_off[1] >= 0) {
        for (int t = 0; t < n; ++t) {
            std::array<int,2> idcpu_parts{};
            auto const idcpu = src.idcpu(si[t]);
            std::memcpy(idcpu_parts.data(), &idcpu, sizeof(idcpu));
            for (int ii = 0; ii < 2; ++ii) {
                if (pl.wire_idcpu_off[ii] >= 0) {
                    std::memcpy(rec[t] + pl.wire_idcpu_off[ii], &idcpu_parts[ii], sizeof(int));
                }
            }
        }
    }

    for (int ii = 0; ii < num_int_comps; ++ii) {
        int const off = pl.wire_int_off[ii];
        if (off < 0) { continue; }
        int const* arr = (ii < NAI) ? src.idata(ii) : src.m_runtime_idata[ii - NAI];
        for (int t = 0; t < n; ++t) {
            std::memcpy(rec[t] + off, arr + si[t], sizeof(int));
        }
    }
}

/// Batched unpackNeighborParticleData for the n records base + k*rs of one
/// message section, into the ghosts d0 + k of dst, one component at a time. The
/// payload of an empty slot (both keys -1 in the record) is not read.
template <typename PTD>
void unpackNeighborBatch (PTD const& dst, int d0, int n, const char* base, std::size_t rs,
                          int num_real_comps, int num_int_comps,
                          GhostCopyPlan const& pl, GhostBatchScratch& sc) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;
    using RT = typename ParticleType::RealType;
    constexpr int NAR = std::decay_t<PTD>::NAR;
    constexpr int NAI = std::decay_t<PTD>::NAI;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        int const po = pl.wire_pos_off[dir];
        if (po < 0) { continue; }
        for (int k = 0; k < n; ++k) {
            std::memcpy(&(dst.pos(dir, d0 + k)), base + std::size_t(k)*rs + po, sizeof(RT));
        }
    }

    for (int ii = 0; ii < num_real_comps; ++ii) {
        int const off = pl.wire_real_off[ii];
        if (off < 0 || pl.wire_real_slot[ii] >= 0) { continue; }
        ParticleReal* const arr = ((ii < NAR) ? dst.rdata(ii) : dst.m_runtime_rdata[ii - NAR]) + d0;
        for (int k = 0; k < n; ++k) {
            std::memcpy(arr + k, base + std::size_t(k)*rs + off, sizeof(ParticleReal));
        }
    }

    int const ns = pl.nslots;
    sc.cnt.assign(ns, 0);
    sc.lists.resize(std::size_t(ns) * std::size_t(n));
    for (int s = 0; s < ns; ++s) {
        int const oid  = pl.wire_key_id_off[s];
        int const ocpu = pl.wire_key_cpu_off[s];
        int* const lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
        int c = 0;
        for (int k = 0; k < n; ++k) {
            char const* r = base + std::size_t(k)*rs;
            int pid, cpu;
            std::memcpy(&pid, r + oid,  sizeof(int));
            std::memcpy(&cpu, r + ocpu, sizeof(int));
            lst[c] = k;
            c += (pid != -1 || cpu != -1) ? 1 : 0;
        }
        sc.cnt[s] = c;
    }
    for (int b = 0; b < pl.nblocks; ++b) {
        for (int s = 0; s < ns; ++s) {
            int const jr = pl.real0 + b*ns + s;
            int const off = pl.wire_real_off[jr + NAR];
            ParticleReal* const arr = dst.m_runtime_rdata[jr] + d0;
            int const* lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
            int const c = sc.cnt[s];
            for (int q = 0; q < c; ++q) {
                int const k = lst[q];
                std::memcpy(arr + k, base + std::size_t(k)*rs + off, sizeof(ParticleReal));
            }
        }
    }

    if (pl.wire_idcpu_off[0] >= 0 || pl.wire_idcpu_off[1] >= 0) {
        for (int k = 0; k < n; ++k) {
            std::array<int,2> idcpu_parts{};
            std::memcpy(idcpu_parts.data(), &(dst.idcpu(d0 + k)), sizeof(uint64_t));
            for (int ii = 0; ii < 2; ++ii) {
                if (pl.wire_idcpu_off[ii] >= 0) {
                    std::memcpy(&idcpu_parts[ii], base + std::size_t(k)*rs + pl.wire_idcpu_off[ii], sizeof(int));
                }
            }
            std::memcpy(&(dst.idcpu(d0 + k)), idcpu_parts.data(), sizeof(uint64_t));
        }
    }

    for (int ii = 0; ii < num_int_comps; ++ii) {
        int const off = pl.wire_int_off[ii];
        if (off < 0) { continue; }
        int* const arr = ((ii < NAI) ? dst.idata(ii) : dst.m_runtime_idata[ii - NAI]) + d0;
        for (int k = 0; k < n; ++k) {
            std::memcpy(arr + k, base + std::size_t(k)*rs + off, sizeof(int));
        }
    }
}

/// Compact version of packNeighborBatch (plan with compact set): the head of ghost
/// t (rec[t], chead_size bytes) holds the components that are not the payload of a
/// slot and the offset of the payload of the ghost in the payload area of its
/// message; the payload of the occupied slots (cpay_bytes each, in slot order) is
/// appended to that area, at *sinfo[t]->used, which advances. The values written
/// are the ones of packNeighborBatch; only their places differ.
template <typename PTD>
void packNeighborBatchCompact (PTD const& src, int n, char* const* rec, int const* si, IntVect const* shift,
                               GhostSendInfo* const* sinfo,
                               int num_real_comps, int num_int_comps,
                               const Periodicity& periodicity, const RealBox& prob_domain,
                               GhostCopyPlan const& pl, GhostBatchScratch& sc) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;
    using RT = typename ParticleType::RealType;
    constexpr int NAR = std::decay_t<PTD>::NAR;
    constexpr int NAI = std::decay_t<PTD>::NAI;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        int const po = pl.cwire_pos_off[dir];
        if (po < 0) { continue; }
        for (int t = 0; t < n; ++t) {
            RT value = src[si[t]].pos(dir);
            applyNeighborPeriodicShift(value, dir, shift[t], periodicity, prob_domain);
            std::memcpy(rec[t] + po, &value, sizeof(RT));
        }
    }

    for (int ii = 0; ii < num_real_comps; ++ii) {
        int const off = pl.cwire_real_off[ii];
        if (off < 0) { continue; }
        ParticleReal const* arr = (ii < NAR) ? src.rdata(ii) : src.m_runtime_rdata[ii - NAR];
        for (int t = 0; t < n; ++t) {
            std::memcpy(rec[t] + off, arr + si[t], sizeof(ParticleReal));
        }
    }

    // ghosts with an occupied slot, per slot, and the place of the slot among the occupied slots of the ghost
    int const ns = pl.nslots;
    sc.cnt.assign(ns, 0);
    sc.ord.assign(n, 0);
    sc.lists.resize(std::size_t(ns) * std::size_t(n));
    sc.kord.resize(std::size_t(ns) * std::size_t(n));
    for (int s = 0; s < ns; ++s) {
        int const* kid  = src.m_runtime_idata[pl.key_id  + s];
        int const* kcpu = src.m_runtime_idata[pl.key_cpu + s];
        int* const lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
        int* const kor = sc.kord.data()  + std::size_t(s) * std::size_t(n);
        int c = 0;
        for (int t = 0; t < n; ++t) {
            int const i = si[t];
            int const o = (kid[i] != -1 || kcpu[i] != -1) ? 1 : 0;
            lst[c] = t;
            kor[c] = sc.ord[t];
            c += o;
            sc.ord[t] += o;
        }
        sc.cnt[s] = c;
    }

    // place the payload of every ghost in the payload area of its message
    sc.poff.resize(n);
    for (int t = 0; t < n; ++t) {
        GhostSendInfo* const inf = sinfo[t];
        int const po = int(inf->used);
        sc.poff[t] = po;
        inf->used += std::size_t(sc.ord[t]) * std::size_t(pl.cpay_bytes);
        std::memcpy(rec[t] + pl.cpoff_off, &po, sizeof(int));
    }

    // payload of the occupied slots: the nblocks reals of one slot are contiguous
    ParticleReal const* parr[MaxPayloadBlocks];
    int const nbk = pl.nblocks;
    for (int s = 0; s < ns; ++s) {
        int const c = sc.cnt[s];
        if (c == 0) { continue; }
        for (int b = 0; b < nbk; ++b) { parr[b] = src.m_runtime_rdata[pl.real0 + b*ns + s]; }
        int const* lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
        int const* kor = sc.kord.data()  + std::size_t(s) * std::size_t(n);
        for (int q = 0; q < c; ++q) {
            int const t = lst[q];
            int const i = si[t];
            char* d = sinfo[t]->pay + sc.poff[t] + std::size_t(kor[q]) * std::size_t(pl.cpay_bytes);
            for (int b = 0; b < nbk; ++b) {
                std::memcpy(d + std::size_t(b) * sizeof(ParticleReal), parr[b] + i, sizeof(ParticleReal));
            }
        }
    }

    if (pl.cwire_idcpu_off[0] >= 0 || pl.cwire_idcpu_off[1] >= 0) {
        for (int t = 0; t < n; ++t) {
            std::array<int,2> idcpu_parts{};
            auto const idcpu = src.idcpu(si[t]);
            std::memcpy(idcpu_parts.data(), &idcpu, sizeof(idcpu));
            for (int ii = 0; ii < 2; ++ii) {
                if (pl.cwire_idcpu_off[ii] >= 0) {
                    std::memcpy(rec[t] + pl.cwire_idcpu_off[ii], &idcpu_parts[ii], sizeof(int));
                }
            }
        }
    }

    for (int ii = 0; ii < num_int_comps; ++ii) {
        int const off = pl.cwire_int_off[ii];
        if (off < 0) { continue; }
        int const* arr = (ii < NAI) ? src.idata(ii) : src.m_runtime_idata[ii - NAI];
        for (int t = 0; t < n; ++t) {
            std::memcpy(rec[t] + off, arr + si[t], sizeof(int));
        }
    }
}

/// Compact version of unpackNeighborBatch: the n heads base + k*hs of one message
/// section go into the ghosts d0 + k of dst, and the payload of their occupied
/// slots is read from the payload area pay of the message.
template <typename PTD>
void unpackNeighborBatchCompact (PTD const& dst, int d0, int n, const char* base, std::size_t hs,
                                 const char* pay,
                                 int num_real_comps, int num_int_comps,
                                 GhostCopyPlan const& pl, GhostBatchScratch& sc) noexcept
{
    using ParticleType = typename std::decay_t<PTD>::ParticleType;
    using RT = typename ParticleType::RealType;
    constexpr int NAR = std::decay_t<PTD>::NAR;
    constexpr int NAI = std::decay_t<PTD>::NAI;

    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        int const po = pl.cwire_pos_off[dir];
        if (po < 0) { continue; }
        for (int k = 0; k < n; ++k) {
            std::memcpy(&(dst.pos(dir, d0 + k)), base + std::size_t(k)*hs + po, sizeof(RT));
        }
    }

    for (int ii = 0; ii < num_real_comps; ++ii) {
        int const off = pl.cwire_real_off[ii];
        if (off < 0) { continue; }
        ParticleReal* const arr = ((ii < NAR) ? dst.rdata(ii) : dst.m_runtime_rdata[ii - NAR]) + d0;
        for (int k = 0; k < n; ++k) {
            std::memcpy(arr + k, base + std::size_t(k)*hs + off, sizeof(ParticleReal));
        }
    }

    int const ns = pl.nslots;
    sc.cnt.assign(ns, 0);
    sc.ord.assign(n, 0);
    sc.lists.resize(std::size_t(ns) * std::size_t(n));
    sc.kord.resize(std::size_t(ns) * std::size_t(n));
    for (int s = 0; s < ns; ++s) {
        int const oid  = pl.cwire_key_id_off[s];
        int const ocpu = pl.cwire_key_cpu_off[s];
        int* const lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
        int* const kor = sc.kord.data()  + std::size_t(s) * std::size_t(n);
        int c = 0;
        for (int k = 0; k < n; ++k) {
            char const* r = base + std::size_t(k)*hs;
            int pid, cpu;
            std::memcpy(&pid, r + oid,  sizeof(int));
            std::memcpy(&cpu, r + ocpu, sizeof(int));
            int const o = (pid != -1 || cpu != -1) ? 1 : 0;
            lst[c] = k;
            kor[c] = sc.ord[k];
            c += o;
            sc.ord[k] += o;
        }
        sc.cnt[s] = c;
    }

    sc.poff.resize(n);
    for (int k = 0; k < n; ++k) {
        std::memcpy(&sc.poff[k], base + std::size_t(k)*hs + pl.cpoff_off, sizeof(int));
    }

    ParticleReal* parr[MaxPayloadBlocks];
    int const nbk = pl.nblocks;
    for (int s = 0; s < ns; ++s) {
        int const c = sc.cnt[s];
        if (c == 0) { continue; }
        for (int b = 0; b < nbk; ++b) { parr[b] = dst.m_runtime_rdata[pl.real0 + b*ns + s] + d0; }
        int const* lst = sc.lists.data() + std::size_t(s) * std::size_t(n);
        int const* kor = sc.kord.data()  + std::size_t(s) * std::size_t(n);
        for (int q = 0; q < c; ++q) {
            int const k = lst[q];
            char const* d = pay + sc.poff[k] + std::size_t(kor[q]) * std::size_t(pl.cpay_bytes);
            for (int b = 0; b < nbk; ++b) {
                std::memcpy(parr[b] + k, d + std::size_t(b) * sizeof(ParticleReal), sizeof(ParticleReal));
            }
        }
    }

    if (pl.cwire_idcpu_off[0] >= 0 || pl.cwire_idcpu_off[1] >= 0) {
        for (int k = 0; k < n; ++k) {
            std::array<int,2> idcpu_parts{};
            std::memcpy(idcpu_parts.data(), &(dst.idcpu(d0 + k)), sizeof(uint64_t));
            for (int ii = 0; ii < 2; ++ii) {
                if (pl.cwire_idcpu_off[ii] >= 0) {
                    std::memcpy(&idcpu_parts[ii], base + std::size_t(k)*hs + pl.cwire_idcpu_off[ii], sizeof(int));
                }
            }
            std::memcpy(&(dst.idcpu(d0 + k)), idcpu_parts.data(), sizeof(uint64_t));
        }
    }

    for (int ii = 0; ii < num_int_comps; ++ii) {
        int const off = pl.cwire_int_off[ii];
        if (off < 0) { continue; }
        int* const arr = ((ii < NAI) ? dst.idata(ii) : dst.m_runtime_idata[ii - NAI]) + d0;
        for (int k = 0; k < n; ++k) {
            std::memcpy(arr + k, base + std::size_t(k)*hs + off, sizeof(int));
        }
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
::fillNeighborsCPU (amrex::Real search_radius) {
    BL_PROFILE("NeighborParticleContainer::fillNeighborsCPU");
    if (!areMasksValid(search_radius)) {
        BuildMasks(search_radius);
        GetNeighborCommTags(search_radius);
    }
    cacheNeighborInfo(search_radius);
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
                const int dst_level = tag.src_level;
                const int global_rank = this->ParticleDistributionMap(dst_level)[dst_grid];
                const int dst_proc = ParallelContext::global_to_local_rank(global_rank);
                const int dst_tile = tag.src_tile;
                const int dst_index = tag.src_index;

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

    // Copy only the flagged components (and the occupied contact slots) for
    // the ghosts that stay on this rank, as the MPI path does for the others.
    const bool masked = ghost_copy_masked;
    // With the persistent buffers the component lists are made at a
    // rebuild and kept for the updates that follow it.
    detail::GhostCopyPlan plan_local;
    if (masked) {
        if (ghost_persistent_recv) {
            if (!reuse_rcv_counts || !ghost_plan_valid) {
                ghost_plan_cache = makeGhostCopyPlan();
                ghost_plan_valid = true;
            }
        } else {
            plan_local = makeGhostCopyPlan();
        }
    }
    const detail::GhostCopyPlan& plan = (masked && ghost_persistent_recv) ? ghost_plan_cache : plan_local;

    const bool inverse = enableInverse();
    // ghosts written straight into the tiles (updates only: the tile already has its ghost count)
    const bool direct = ghost_direct && reuse_rcv_counts && !inverse;
    // records packed and unpacked one component at a time (one thread: the batch is collected serially)
    const bool batched = plan.batch && !inverse && (OpenMP::get_max_threads() == 1);
    // compact records: the layout of the message buffers is the one made at the last rebuild
    // (cacheNeighborInfo), which has to be the one this update packs
    const bool compact = ghost_compact_layout;
    if (compact && !(batched && plan.compact)) {
        amrex::Abort("NeighborParticleContainer::updateNeighborsCPU: compact ghost records need the batched layout of the rebuild");
    }
    detail::GhostBatchScratch& sc = ghost_scratch;
    if (!reuse_rcv_counts && !ghost_switches_told && amrex::Verbose() > 0 &&
        (ghost_persistent_recv || ghost_direct || ghost_batched || ghost_isend || ghost_compact)) {
        ghost_switches_told = true;
        amrex::Print() << "Ghost refresh switches: persistent receive buffer "
                       << (ghost_persistent_recv ? "on" : "off")
                       << ", ghosts written straight into the tiles "
                       << (ghost_direct ? (inverse ? "off (inverse tags)" : "on") : "off")
                       << ", records packed per component "
                       << (batched ? "on" : (ghost_batched ? "off (not available for this layout)" : "off"))
                       << ", sends " << (ghost_isend ? "nonblocking" : "blocking")
                       << ", compact records "
                       << (compact ? "on" : (ghost_compact ? "off (not available for this layout)" : "off"))
                       << " (record " << plan.rec_size << " B, communication size " << cdata_size << " B";
        if (compact) { amrex::Print() << ", head " << plan.chead_size << " B and " << plan.cpay_bytes << " B per occupied slot"; }
        amrex::Print() << ")\n";
    }

    if (compact) {
        // the payload area of the message to each rank starts after its heads and is empty again
        ghost_send_info.assign(ParallelContext::NProcsSub(), detail::GhostSendInfo());
        for (auto& kv : send_data) {
            auto& inf = ghost_send_info[kv.first];
            inf.pay  = kv.second.data() + ghost_head_bytes.at(kv.first);
            inf.used = 0;
        }
    }

    if (direct) {
        // the tile keeps the ghost count of the last rebuild; the staging tile still holds it
        for (int lev = 0; lev < this->numLevels(); ++lev) {
            for (MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi) {
                auto index = std::make_pair(mfi.index(), mfi.LocalTileIndex());
                auto& ptile = this->GetParticles(lev)[index];
                const int nn = int(neighbors[lev][index].size());
                if (ptile.numNeighborParticles() != nn) { ptile.setNumNeighbors(nn); }
            }
        }
    }

    for (int lev = 0; lev < this->numLevels(); ++lev) {
        const Periodicity& periodicity = this->Geom(lev).periodicity();
        const RealBox& prob_domain = this->Geom(lev).ProbDomain();

        int num_threads = OpenMP::get_max_threads();

        for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
            PairIndex src_index(pti.index(), pti.LocalTileIndex());
            auto const src = pti.GetParticleTile().getConstParticleTileData();
            // pack the collected remote ghosts of this tile, one component at a time
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
                        // staging tile, or the ghost region of the particle tile itself
                        auto& dtile = direct ? this->GetParticles(tag.level)[dst_index]
                                             : neighbors[tag.level][dst_index];
                        const int di = direct ? dtile.numRealParticles() + tag.dst_index : tag.dst_index;
                        auto dst = dtile.getParticleTileData();
                        if (masked) {
                            detail::copyGhostParticle(dst, src, tag.src_index, di, plan);
                        } else {
                            copyParticle(dst, src, tag.src_index, di);
                        }
                        if (periodicity.isAnyPeriodic()) {
                            detail::shiftNeighborParticlePositions(dst, di, tag.periodic_shift,
                                                                  periodicity, prob_domain);
                        }

                        if ( inverse )
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
                        if (batched) {
                            sc.rec.push_back(dst_ptr);
                            sc.si.push_back(tag.src_index);
                            sc.shift.push_back(tag.periodic_shift);
                            if (compact) { sc.sinfo.push_back(&ghost_send_info[who]); }
                        } else {
                        detail::packNeighborParticleData(dst_ptr, src, tag.src_index,
                                                         this->NumRealComps(), this->NumIntComps(),
                                                         ghost_real_comp, ghost_int_comp,
                                                         tag.periodic_shift, periodicity, prob_domain,
                                                         &plan);
                        if ( inverse )
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
            if (batched && !sc.si.empty()) {
                const int nb = int(sc.si.size());
                if (compact) {
                    detail::packNeighborBatchCompact(src, nb, sc.rec.data(), sc.si.data(), sc.shift.data(),
                                                     sc.sinfo.data(),
                                                     this->NumRealComps(), this->NumIntComps(),
                                                     periodicity, prob_domain, plan, sc);
                } else {
                    detail::packNeighborBatch(src, nb, sc.rec.data(), sc.si.data(), sc.shift.data(),
                                              this->NumRealComps(), this->NumIntComps(),
                                              periodicity, prob_domain, plan, sc);
                }
                sc.rec.clear();
                sc.si.clear();
                sc.shift.clear();
                sc.sinfo.clear();
            }
        }

        if (!direct) {
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
    }
    BL_PROFILE_VAR_STOP(update);

    fillNeighborsMPI(reuse_rcv_counts, plan, direct, batched);

    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            int src_grid = mfi.index();
            int src_tile = mfi.LocalTileIndex();
            auto index = std::make_pair(src_grid, src_tile);
            auto& ptile = this->GetParticles(lev)[index];
            if (direct) {
                // the ghosts are in the tile already
                continue;
            }
            ptile.setNumNeighbors(neighbors[lev][index].size());
            if (masked) {
                auto const sdata = neighbors[lev][index].getConstParticleTileData();
                auto const ddata = ptile.getParticleTileData();
                const int nrp = ptile.numRealParticles();
                const int nn = ptile.numNeighborParticles();
#ifdef AMREX_USE_OMP
#pragma omp parallel for
#endif
                for (int i = 0; i < nn; ++i) {
                    detail::copyGhostParticle(ddata, sdata, i, nrp + i, plan);
                }
            } else {
                amrex::copyParticles(ptile, neighbors[lev][index], 0,
                                     ptile.numRealParticles(), ptile.numNeighborParticles());
            }
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
fillNeighborsMPI (bool reuse_rcv_counts, detail::GhostCopyPlan const& plan, bool direct, bool batched) {

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

    // Compact records: a message is its heads followed by the used part of its payload area
    // (the receive sizes of the rebuild are the capacities, so a shorter message is fine).
    const bool compact = ghost_compact_layout;
    if (compact && !batched) {
        amrex::Abort("NeighborParticleContainer::fillNeighborsMPI: compact ghost records need the batched layout of the rebuild");
    }

    // Allocate data for rcvs as one big chunk (kept between updates with the persistent buffer).
    Vector<char> recv_local;
    Vector<char>& recvdata = ghost_persistent_recv ? ghost_recv_buf : recv_local;
    if (ghost_persistent_recv) {
        if (recvdata.size() < static_cast<Long>(TotRcvBytes)) { recvdata.resize(TotRcvBytes); }
    } else {
        recvdata.resize(TotRcvBytes);
    }

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

    // Send (nonblocking with ghost_isend: the requests are waited for below).
    Vector<MPI_Request> sreqs;
    for (const auto& kv : send_data) {
        const auto Who = kv.first;
        const std::size_t Cnt = compact ? (ghost_head_bytes.at(Who) + ghost_send_info[Who].used)
                                        : kv.second.size();

        AMREX_ASSERT(Cnt > 0);
        AMREX_ASSERT(Who >= 0 && Who < NProcs);
        AMREX_ASSERT(Cnt < std::numeric_limits<int>::max());

        if (ghost_isend) {
            sreqs.push_back(ParallelDescriptor::Asend(kv.second.data(), Cnt, Who, SeqNum,
                                                      ParallelContext::CommunicatorSub()).req());
        } else {
            ParallelDescriptor::Send(kv.second.data(), Cnt, Who, SeqNum,
                                     ParallelContext::CommunicatorSub());
        }
    }
    if (!sreqs.empty()) {
        Vector<MPI_Status> sstats(sreqs.size());
        ParallelDescriptor::Waitall(sreqs, sstats);
    }

    // unpack the received data and put them into the proper neighbor buffers
    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);
        // batched unpack of the np records of a message section
        const std::size_t rs = compact ? std::size_t(plan.chead_size) : cdata_size;   // bytes of a record or a head
        auto unpack_batched = [&] (auto const& dptd, int d0, int cnt, const char* buf, const char* pay) {
            if (compact) {
                detail::unpackNeighborBatchCompact(dptd, d0, cnt, buf, rs, pay,
                                                   this->NumRealComps(), this->NumIntComps(),
                                                   plan, ghost_scratch);
            } else {
                detail::unpackNeighborBatch(dptd, d0, cnt, buf, cdata_size,
                                            this->NumRealComps(), this->NumIntComps(),
                                            plan, ghost_scratch);
            }
        };
        // direct mode: index of the next remote ghost of each tile (after the ghosts of this rank)
        Vector<std::map<PairIndex, int> > fill_pos;
        if (direct) { fill_pos = local_neighbor_sizes; }
        for (int i = 0; i < nrcvs; ++i) {
            const std::size_t offset = rOffset[i];
            char* buffer = &recvdata[offset];
            int num_tiles, lev, gid, tid, size, np;
            std::memcpy(&num_tiles, buffer, sizeof(int)); buffer += sizeof(int);
            const char* pay_base = nullptr;
            if (compact) {
                // the payload area starts after the heads of the last tile of the message
                const char* p = buffer;
                for (int j = 0; j < num_tiles; ++j) {
                    int sz;
                    std::memcpy(&sz, p + 3*sizeof(int), sizeof(int));
                    p += 4*sizeof(int) + std::size_t(sz);
                }
                pay_base = p;
            }
            for (int j = 0; j < num_tiles; ++j) {
                std::memcpy(&lev,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&gid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&tid,  buffer, sizeof(int)); buffer += sizeof(int);
                std::memcpy(&size, buffer, sizeof(int)); buffer += sizeof(int);

                if (size == 0) { continue; }

                np = int(std::size_t(size) / rs);

                AMREX_ASSERT(std::size_t(size) % rs == 0);

                PairIndex dst_index(gid, tid);

                if (direct) {
                    // straight into the ghost region of the particle tile
                    auto& ptile = this->GetParticles(lev)[dst_index];
                    int& pos = fill_pos[lev][dst_index];
                    const int d0 = ptile.numRealParticles() + pos;
                    pos += np;
                    auto dst_tile_ptd = ptile.getParticleTileData();
                    if (batched) {
                        unpack_batched(dst_tile_ptd, d0, np, buffer, pay_base);
                    } else {
                        char const* dsrc = buffer;
                        for (int n = 0; n < np; ++n) {
                            detail::unpackNeighborParticleData(dst_tile_ptd, d0 + n, dsrc,
                                                               this->NumRealComps(), this->NumIntComps(),
                                                               ghost_real_comp, ghost_int_comp, &plan);
                        }
                    }
                    buffer += size;
                    continue;
                }

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
                if (batched) {
                    unpack_batched(dst_ptd, int(old_size), np, buffer, pay_base);
                    buffer += size;
                    continue;
                }
                char const* src = buffer;
                for (int n = 0; n < np; ++n) {
                    detail::unpackNeighborParticleData(dst_ptd, int(old_size+n), src,
                                                       this->NumRealComps(), this->NumIntComps(),
                                                       ghost_real_comp, ghost_int_comp, &plan);

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
        if (direct) {
            // every tile has received exactly the ghosts it held after the rebuild
            for (int lev = 0; lev < this->numLevels(); ++lev) {
                for (const auto& kv : fill_pos[lev]) {
                    if (size_t(kv.second) != neighbors[lev][kv.first].size()) {
                        amrex::Abort("NeighborParticleContainer::fillNeighborsMPI: ghost count of a tile changed since the rebuild");
                    }
                }
            }
        }
    }
#else
    amrex::ignore_unused(reuse_rcv_counts, plan, direct, batched);
#endif
}

}

//! \endcond

#endif
