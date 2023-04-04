#ifndef AMREX_NEIGHBORPARTICLESGPUIMPL_H_
#define AMREX_NEIGHBORPARTICLESGPUIMPL_H_
#include <AMReX_Config.H>

namespace detail
{
    inline Vector<Box> getBoundaryBoxes(const Box& box, int ncells)
    {
        AMREX_ASSERT_WITH_MESSAGE(box.size() > 2*IntVect(AMREX_D_DECL(ncells, ncells, ncells)),
                                  "Too many cells requested in getBoundaryBoxes");

        AMREX_ASSERT_WITH_MESSAGE(box.ixType().cellCentered(),
                                  "Box must be cell-centered");

        Vector<Box> bl;
        for (int i = 0; i < AMREX_SPACEDIM; ++i) {
            BoxList face_boxes;
            Box hi_face_box = adjCellHi(box, i, ncells);
            Box lo_face_box = adjCellLo(box, i, ncells);
            face_boxes.push_back(hi_face_box); bl.push_back(hi_face_box);
            face_boxes.push_back(lo_face_box); bl.push_back(lo_face_box);
            for (auto face_box : face_boxes) {
                for (int j = 0; j < AMREX_SPACEDIM; ++j) {
                    if (i == j) continue;
                    BoxList edge_boxes;
                    Box hi_edge_box = adjCellHi(face_box, j, ncells);
                    Box lo_edge_box = adjCellLo(face_box, j, ncells);
                    edge_boxes.push_back(hi_edge_box); bl.push_back(hi_edge_box);
                    edge_boxes.push_back(lo_edge_box); bl.push_back(lo_edge_box);
                    for (auto edge_box : edge_boxes) {
                        for (int k = 0; k < AMREX_SPACEDIM; ++k) {
                            if ((j == k) || (i == k)) continue;
                            Box hi_corner_box = adjCellHi(edge_box, k, ncells);
                            Box lo_corner_box = adjCellLo(edge_box, k, ncells);
                            bl.push_back(hi_corner_box);
                            bl.push_back(lo_corner_box);
                        }
                    }
                }
            }
        }

        RemoveDuplicates(bl);
        return bl;
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
buildNeighborMask ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::buildNeighborMask");
    m_neighbor_mask_initialized = true;
    const int lev = 0;
    const Geometry& geom = this->Geom(lev);
    const BoxArray& ba = this->ParticleBoxArray(lev);
    const DistributionMapping& dmap = this->ParticleDistributionMap(lev);

    if (ba.size() == 1 && (! geom.isAnyPeriodic()) ) return;

    if (m_neighbor_mask_ptr == nullptr ||
        ! BoxArray::SameRefs(m_neighbor_mask_ptr->boxArray(), ba) ||
        ! DistributionMapping::SameRefs(m_neighbor_mask_ptr->DistributionMap(), dmap))
    {
        const Periodicity& periodicity = geom.periodicity();
        const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();

        for (MFIter mfi(ba, dmap); mfi.isValid(); ++mfi)
        {
            int grid = mfi.index();

            std::set<NeighborTask> neighbor_grids;
            for (auto pit=pshifts.cbegin(); pit!=pshifts.cend(); ++pit)
            {
                const Box box = ba[mfi] + *pit;

                const bool first_only = false;
                auto isecs = ba.intersections(box, first_only, m_num_neighbor_cells);

                for (auto& isec : isecs)
                {
                    int nbor_grid = isec.first;
                    const Box isec_box = isec.second - *pit;
                    if ( (grid == nbor_grid) && (*pit == IntVect(AMREX_D_DECL(0, 0, 0)))) continue;
                    neighbor_grids.insert(NeighborTask(nbor_grid, isec_box, *pit));
                    const int global_rank = dmap[nbor_grid];
                    neighbor_procs.push_back(ParallelContext::global_to_local_rank(global_rank));
                }
            }

            Gpu::HostVector<Box>          h_isec_boxes;
            Gpu::HostVector<NeighborCode> h_code_arr;
            for (auto nbor_grid : neighbor_grids)
            {
                NeighborCode code;
                code.grid_id        = nbor_grid.grid_id;
                code.periodic_shift = nbor_grid.periodic_shift;
                h_code_arr.push_back(code);
                h_isec_boxes.push_back(nbor_grid.box);
            }

            m_code_array[grid].resize(h_code_arr.size());
            Gpu::copyAsync(Gpu::hostToDevice, h_code_arr.begin(), h_code_arr.end(),
                      m_code_array[grid].begin());
            m_isec_boxes[grid].resize(h_isec_boxes.size());
            Gpu::copyAsync(Gpu::hostToDevice, h_isec_boxes.begin(), h_isec_boxes.end(),
                      m_isec_boxes[grid].begin());

            Gpu::streamSynchronize();
        }

        RemoveDuplicates(neighbor_procs);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
buildNeighborCopyOp (bool use_boundary_neighbor)
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::buildNeighborCopyOp()");

    AMREX_ASSERT(!hasNeighbors() || use_boundary_neighbor);

    const int lev = 0;
    const auto& geom = this->Geom(lev);
    const auto dxi = this->Geom(lev).InvCellSizeArray();
    const auto plo = this->Geom(lev).ProbLoArray();
    const auto domain = this->Geom(lev).Domain();
    auto& plev  = this->GetParticles(lev);
    auto& ba = this->ParticleBoxArray(lev);

    if (ba.size() == 1 && (! geom.isAnyPeriodic()) ) return;

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& src_tile = plev[index];
        auto& aos   = src_tile.GetArrayOfStructs();
        const size_t np_real = aos.numParticles();

        size_t np = np_real;
        if (use_boundary_neighbor) {
            np = m_boundary_particle_ids[lev][index].size();
        }
        else {
            m_boundary_particle_ids.resize(1);
            m_boundary_particle_ids[lev][index];
        }

        const auto*  p_bndry_pid = m_boundary_particle_ids[lev][index].dataPtr();

        Gpu::DeviceVector<int> counts(np + 1, 0);
        Gpu::DeviceVector<int> offsets(np + 1);
        auto p_counts = counts.dataPtr();
        auto p_offsets = offsets.dataPtr();

        ParticleType* p_ptr = &(aos[0]);
        auto p_code_array   = m_code_array[gid].dataPtr();
        auto p_isec_boxes   = m_isec_boxes[gid].dataPtr();
        const int nisec_box = m_isec_boxes[gid].size();
        // auto p_code_offsets = m_code_offsets[gid].dataPtr();

        AMREX_FOR_1D ( np, i,
        {
            // note that cannot use ternary statement here because p_bndry is not
            // properly allocated when use_boundary_neighbor=false
            int pid = i;
            if (use_boundary_neighbor) {
                pid = p_bndry_pid[i];
            }
            IntVect iv = getParticleCell(p_ptr[pid], plo, dxi, domain);
            for (int j=0; j<nisec_box; ++j) {
                if (p_isec_boxes[j].contains(iv)) {
                    ++p_counts[i];
                }
            }
        });

        amrex::Gpu::exclusive_scan(counts.begin(), counts.end(), offsets.begin());

        int num_copies;
        Gpu::dtoh_memcpy_async(&num_copies, offsets.data()+np, sizeof(int));
        Gpu::streamSynchronize();

        neighbor_copy_op.resize(gid, lev, num_copies);

        auto p_boxes = neighbor_copy_op.m_boxes[lev][gid].dataPtr();
        auto p_levs = neighbor_copy_op.m_levels[lev][gid].dataPtr();
        auto p_src_indices = neighbor_copy_op.m_src_indices[lev][gid].dataPtr();
        auto p_periodic_shift = neighbor_copy_op.m_periodic_shift[lev][gid].dataPtr();

        Gpu::streamSynchronize();
        AMREX_FOR_1D ( np, i,
        {
            int pid = i;
            if (use_boundary_neighbor) {
                pid = p_bndry_pid[i];
            }
            IntVect iv = getParticleCell(p_ptr[pid], plo, dxi, domain);
            int      k = p_offsets[i];
            for (int j=0; j<nisec_box; ++j) {
                if (p_isec_boxes[j].contains(iv)) {
                    p_boxes[k]          = p_code_array[j].grid_id;
                    p_levs[k]           = 0;
                    p_periodic_shift[k] = p_code_array[j].periodic_shift;
                    p_src_indices[k]    = pid;
                    ++k;
                }
            }
            AMREX_ALWAYS_ASSERT(k == p_offsets[i+1]);
        });

        Gpu::streamSynchronize();
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
fillNeighborsGPU ()
{
    BL_PROFILE("NeighborParticleContainer::fillNeighbors");

    AMREX_ASSERT(numParticlesOutOfRange(*this, 0) == 0);

    buildNeighborMask();
    this->defineBufferMap();

    neighbor_copy_op.clear();
    neighbor_copy_plan.clear();
    buildNeighborCopyOp();
    neighbor_copy_plan.build(*this, neighbor_copy_op, ghost_int_comp,
                             ghost_real_comp, true);
    updateNeighborsGPU(false);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
updateNeighborsGPU (bool boundary_neighbors_only)
{
    BL_PROFILE("NeighborParticleContainer::updateNeighborsGPU");

    if (boundary_neighbors_only) {
        neighbor_copy_op.clear();
        neighbor_copy_plan.clear();
        buildNeighborCopyOp(true);
        neighbor_copy_plan.build(*this, neighbor_copy_op, ghost_int_comp,
                                 ghost_real_comp, true);
    }

    clearNeighbors();
    packBuffer(*this, neighbor_copy_op, neighbor_copy_plan, snd_buffer);
    if (ParallelDescriptor::UseGpuAwareMpi())
    {
        neighbor_copy_plan.buildMPIFinish(this->BufferMap());
        communicateParticlesStart(*this, neighbor_copy_plan, snd_buffer, rcv_buffer);
        unpackBuffer(*this, neighbor_copy_plan, snd_buffer, NeighborUnpackPolicy());
        communicateParticlesFinish(neighbor_copy_plan);
        unpackRemotes(*this, neighbor_copy_plan, rcv_buffer, NeighborUnpackPolicy());
    }
    else
    {
        Gpu::streamSynchronize();
        if (snd_buffer.arena()->isPinned()) {
            neighbor_copy_plan.buildMPIFinish(this->BufferMap());
            Gpu::streamSynchronize();
            communicateParticlesStart(*this, neighbor_copy_plan, snd_buffer, pinned_rcv_buffer);
        } else {
            pinned_snd_buffer.resize(snd_buffer.size());
            Gpu::dtoh_memcpy_async(pinned_snd_buffer.dataPtr(), snd_buffer.dataPtr(), snd_buffer.size());
            neighbor_copy_plan.buildMPIFinish(this->BufferMap());
            Gpu::streamSynchronize();
            communicateParticlesStart(*this, neighbor_copy_plan, pinned_snd_buffer, pinned_rcv_buffer);
        }

        rcv_buffer.resize(pinned_rcv_buffer.size());
        unpackBuffer(*this, neighbor_copy_plan, snd_buffer, NeighborUnpackPolicy());
        communicateParticlesFinish(neighbor_copy_plan);
        Gpu::htod_memcpy_async(rcv_buffer.dataPtr(), pinned_rcv_buffer.dataPtr(), pinned_rcv_buffer.size());
        unpackRemotes(*this, neighbor_copy_plan, rcv_buffer, NeighborUnpackPolicy());
    }

    Gpu::streamSynchronize();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
NeighborParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
clearNeighborsGPU()
{
    BL_PROFILE("NeighborParticleContainer::clearNeighborsGPU");

    this->reserveData();
    this->resizeData();
    for (int lev = 0; lev < this->numLevels(); ++lev)
    {
        for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            auto& ptile = this->DefineAndReturnParticleTile(lev, mfi);
            ptile.setNumNeighbors(0);
        }
    }
}

#endif
