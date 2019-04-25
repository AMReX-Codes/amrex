#ifndef AMREX_NEIGHBORPARTICLESGPUIMPL_H_
#define AMREX_NEIGHBORPARTICLESGPUIMPL_H_

namespace detail
{
    inline Vector<Box> getBoundaryBoxes(const Box& box, const int ncells)
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
                            if ((j == k) or (i == k)) continue;
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

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
buildNeighborMask()
{    
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::buildNeighborMask");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);
    
    m_neighbor_mask_initialized = true;

    const int lev = 0;
    const Geometry& geom = this->Geom(lev);
    const BoxArray& ba = this->ParticleBoxArray(lev);
    const DistributionMapping& dmap = this->ParticleDistributionMap(lev);

    if (ba.size() == 1 and (not geom.isAnyPeriodic()) ) return;

    if (m_neighbor_mask_ptr == nullptr ||
        ! BoxArray::SameRefs(m_neighbor_mask_ptr->boxArray(), ba) ||
        ! DistributionMapping::SameRefs(m_neighbor_mask_ptr->DistributionMap(), dmap))
    {
        m_neighbor_mask_ptr.reset(new iMultiFab(ba, dmap, 1, 0));
        m_neighbor_mask_ptr->setVal(-1);
        
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
                    if ( (grid == nbor_grid) and (*pit == IntVect(AMREX_D_DECL(0, 0, 0)))) continue; 
                    neighbor_grids.insert(NeighborTask(nbor_grid, isec_box, *pit));
                }
            }
            
            BoxList isec_bl;
            std::vector<int> isec_grids;
            std::vector<IntVect> isec_pshifts;
            for (auto nbor_grid : neighbor_grids)
            {
                isec_grids.push_back(nbor_grid.grid_id);
                isec_bl.push_back(nbor_grid.box);
                isec_pshifts.push_back(nbor_grid.periodic_shift);
            }
            BoxArray isec_ba(isec_bl);
            
            Vector<Box> bl = detail::getBoundaryBoxes(
                amrex::grow(ba[mfi], -m_num_neighbor_cells), m_num_neighbor_cells);
            
            m_grid_map[grid].resize(bl.size());
            for (int i = 0; i < static_cast<int>(bl.size()); ++i)
            {
                const Box& box = bl[i];
            
                const int nGrow = 0;
                const bool first_only = false;
                auto isecs = isec_ba.intersections(box, first_only, nGrow);
                
                if (! isecs.empty() ) (*m_neighbor_mask_ptr)[mfi].setVal(i, box);

                for (auto& isec : isecs)
                {
                    NeighborCode code;
                    code.grid_id = isec_grids[isec.first];
                    code.periodic_shift = isec_pshifts[isec.first];
                    m_grid_map[grid][i].push_back(code);
                }
            }
        }
    }
}

template <int NStructReal, int NStructInt>
void 
NeighborParticleContainer<NStructReal, NStructInt>::
sortParticlesByNeighborDestGPU()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::sortParticlesByNeighborDestGPU");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    const int lev = 0;
    const auto& geom = this->Geom(lev);
    const auto dxi = this->Geom(lev).InvCellSizeArray();
    const auto plo = this->Geom(lev).ProbLoArray();
    const auto domain = this->Geom(lev).Domain();
    auto& plev  = this->GetParticles(lev);
    auto& ba = this->ParticleBoxArray(lev);
    auto& dmap = this->ParticleDistributionMap(lev);

    if (ba.size() == 1 and (not geom.isAnyPeriodic()) ) return;

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();        
        auto index = std::make_pair(gid, tid);

        const Box& bx = mfi.tilebox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        auto& src_tile = plev[index];
        auto& aos   = src_tile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        Gpu::ManagedDeviceVector<int> neighbor_codes(np);
        int* pcodes = neighbor_codes.dataPtr();

        BaseFab<int>* mask_ptr = m_neighbor_mask_ptr->fabPtr(mfi);
        
        ParticleType* p_ptr = &(aos[0]);
        AMREX_FOR_1D ( np, i,
        {
	    IntVect iv = IntVect(
                AMREX_D_DECL(floor((p_ptr[i].pos(0)-plo[0])*dxi[0]),
                             floor((p_ptr[i].pos(1)-plo[1])*dxi[1]),
                             floor((p_ptr[i].pos(2)-plo[2])*dxi[2]))
                );
            
            iv += domain.smallEnd();
            
	    int code = (*mask_ptr)(iv);
            
            pcodes[i] = code;
        });

        Gpu::ManagedDeviceVector<int> index_sequence_tmp(np);
        thrust::sequence(thrust::device, index_sequence_tmp.begin(), index_sequence_tmp.end());

        ParticleVector aos_r(np);
        thrust::sort_by_key(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
                            neighbor_codes.begin(),
                            neighbor_codes.end(),
                            index_sequence_tmp.begin());

        //
        // Reorder the particle data 
        //
        {
	    // reorder structs
	    aos_r.resize(np);
            thrust::gather(thrust::device,
                           index_sequence_tmp.begin(), index_sequence_tmp.end(),
                           aos().begin(), aos_r.begin());
            aos().swap(aos_r);
        }

        int num_codes = m_grid_map[gid].size();
        Gpu::ManagedDeviceVector<int> neighbor_code_begin(num_codes + 1);
        Gpu::ManagedDeviceVector<int> neighbor_code_end  (num_codes + 1);
        
        thrust::counting_iterator<int> search_begin(-1);
        thrust::lower_bound(thrust::device,
                            neighbor_codes.begin(),
                            neighbor_codes.end(),
                            search_begin,
                            search_begin + num_codes + 1,
                            neighbor_code_begin.begin());
        
        thrust::upper_bound(thrust::device,
                            neighbor_codes.begin(),
                            neighbor_codes.end(),
                            search_begin,
                            search_begin + num_codes + 1,
                            neighbor_code_end.begin());
        
        m_start[gid].resize(num_codes + 1);
        m_stop[gid].resize(num_codes + 1);
        
        Cuda::thrust_copy(neighbor_code_begin.begin(),
                          neighbor_code_begin.end(),
                          m_start[gid].begin());
        
        Cuda::thrust_copy(neighbor_code_end.begin(),
                          neighbor_code_end.end(),
                          m_stop[gid].begin());
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
fillNeighborsGPU ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::fillNeighbors");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    buildNeighborMask();
    sortParticlesByNeighborDestGPU();
    updateNeighborsGPU();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
updateNeighborsGPU ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::updateNeighborsGPU");

    AMREX_ASSERT(this->numLevels() == 1);

    clearNeighbors();

    const int lev = 0;
    auto& plev  = this->GetParticles(lev);
    auto& ba = this->ParticleBoxArray(lev);
    auto& dmap = this->ParticleDistributionMap(lev);
    
    if (ba.size() == 1 and (not this->Geom(lev).isAnyPeriodic()) ) return;

    std::map<int, SendBuffer> not_ours;
    std::map<int, size_t> grid_counts;
    
    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int src_grid = mfi.index();
        int num_codes = m_grid_map[src_grid].size();
        for (int i = 1; i < num_codes + 1; ++i)
        {
            for (auto dst_code : m_grid_map[src_grid][i-1])
            {            
                int dst_grid = dst_code.grid_id;
                const int dest_proc = dmap[dst_grid];
                const int num_to_add = m_stop[src_grid][i] - m_start[src_grid][i];
                if ( (dest_proc != ParallelDescriptor::MyProc() ) and (num_to_add > 0) )
                {
                    grid_counts[dest_proc] += 1;
                }
            }
        }
    }

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int src_grid = mfi.index();
        int src_tile = mfi.LocalTileIndex();
        AMREX_ASSERT(src_tile == 0);
        auto index = std::make_pair(src_grid, src_tile);

        const Box& bx = mfi.tilebox();
        const auto lo = amrex::lbound(bx);
        const auto hi = amrex::ubound(bx);

        auto& src_ptile = plev[index];

        int num_codes = m_grid_map[src_grid].size();
        for (int i = 1; i < num_codes + 1; ++i)
        {
            const size_t num_to_add = m_stop[src_grid][i] - m_start[src_grid][i];
            if (num_to_add == 0) continue;
            for (auto dst_code : m_grid_map[src_grid][i-1]) {
                int dst_grid = dst_code.grid_id;
                const int dst_tile = 0;                
                auto pair_index = std::make_pair(dst_grid, dst_tile);
                const int dest_proc = dmap[dst_grid];
                if (dest_proc == ParallelDescriptor::MyProc())  // this is a local copy
                {
                    auto& dst_ptile = plev[pair_index];
                    const int nRParticles = dst_ptile.numRealParticles();
                    const int nNParticles = dst_ptile.getNumNeighbors();
                    const int new_num_neighbors = nNParticles + num_to_add;
                    dst_ptile.setNumNeighbors(new_num_neighbors);
                    
                    // copy structs
                    {
                        auto& src = src_ptile.GetArrayOfStructs();
                        auto& dst = dst_ptile.GetArrayOfStructs();
                        thrust::copy(thrust::device,
                                     src.begin() + m_start[src_grid][i], src.begin() + m_stop[src_grid][i], 
                                     dst().begin() + nRParticles + nNParticles);

                        // apply periodic shift so that a naive distance calculation returns
                        // the shortest distance between real and neighbor particles
                        if (this->Geom(lev).isAnyPeriodic()) {
                            thrust::transform(thrust::device, 
                                              dst().begin() + nRParticles + nNParticles,
                                              dst().begin() + nRParticles + nNParticles + num_to_add,
                                              dst().begin() + nRParticles + nNParticles,
                                              functors::periodicTransformation(this->Geom(lev).data(),
                                                                               dst_code.periodic_shift));
                        }
                    }
                }
                else { // this is the non-local case
                    char* dst;
                    const size_t old_size = not_ours[dest_proc].size();
                    const size_t new_size
                        = old_size + num_to_add*sizeof(ParticleType) + sizeof(size_t) + 2*sizeof(int);
                
                    if (old_size == 0)
                    {
                        not_ours[dest_proc].resize(new_size + sizeof(size_t));
                        cudaMemcpyAsync(thrust::raw_pointer_cast(not_ours[dest_proc].data()),
                                        &grid_counts[dest_proc], sizeof(size_t), cudaMemcpyHostToHost);
                        dst = thrust::raw_pointer_cast(
                            not_ours[dest_proc].data() + old_size + sizeof(size_t));
                    } else
                    {
                        not_ours[dest_proc].resize(new_size);
                        dst = thrust::raw_pointer_cast(not_ours[dest_proc].data() + old_size);
                    }
                    
                    cudaMemcpyAsync(thrust::raw_pointer_cast(dst), 
                                    &num_to_add, sizeof(size_t), cudaMemcpyHostToHost);
                    dst += sizeof(size_t);
                    
                    cudaMemcpyAsync(thrust::raw_pointer_cast(dst), &dst_grid, sizeof(int), 
                                    cudaMemcpyHostToHost);
                    dst += sizeof(int);
                    
                    cudaMemcpyAsync(thrust::raw_pointer_cast(dst), 
                                    &dest_proc, sizeof(int), cudaMemcpyHostToHost);
                    dst += sizeof(int);
                    
                    // pack structs
                    {
                        auto& aos = src_ptile.GetArrayOfStructs();
                        cudaMemcpyAsync(thrust::raw_pointer_cast(dst), 
                                        thrust::raw_pointer_cast(aos().data() + m_start[src_grid][i]),
                                        num_to_add*sizeof(ParticleType), cudaMemcpyDeviceToHost);

                        // apply periodic shift so that a naive distance calculation returns
                        // the shortest distance between real and neighbor particles
                        if (this->Geom(lev).isAnyPeriodic()) {
                            thrust::transform(thrust::device, 
                                              thrust::raw_pointer_cast(reinterpret_cast<ParticleType*>(dst)),
                                              thrust::raw_pointer_cast(reinterpret_cast<ParticleType*>(dst) + num_to_add),
                                              thrust::raw_pointer_cast(reinterpret_cast<ParticleType*>(dst)),
                                              functors::periodicTransformation(this->Geom(lev).data(),
                                                                               dst_code.periodic_shift));
                        }
                        dst += num_to_add*sizeof(ParticleType);
                    }
                }
            }
        }
    }
    
    if (ParallelDescriptor::NProcs() == 1) {
        BL_ASSERT(not_ours.empty());
    }
    else {
        fillNeighborsMPIGPU(not_ours);
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
fillNeighborsMPIGPU (std::map<int, SendBuffer>& not_ours)
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::fillNeighborsMPIGPU()");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    const int lev = 0;
    
    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;

    for (const auto& kv : not_ours)
    {
        const size_t nbytes = kv.second.size();
        Snds[kv.first] = nbytes;
        NumSnds += nbytes;
    }
    
    ParallelDescriptor::ReduceLongMax(NumSnds);

    if (NumSnds == 0) return;

    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::BeforeCall());
    
    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<long>::type(),
                                 ParallelDescriptor::Communicator()) );
    
    BL_ASSERT(Rcvs[ParallelDescriptor::MyProc()] == 0);
    
    BL_COMM_PROFILE(BLProfiler::Alltoall, sizeof(long),
                    ParallelDescriptor::MyProc(), BLProfiler::AfterCall());

    Vector<int> RcvProc;
    Vector<std::size_t> rOffset; // Offset (in bytes) in the receive buffer
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes);
            TotRcvBytes += Rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();
    
    // Allocate data for rcvs as one big chunk.    
    char* rcv_buffer;
    if (ParallelDescriptor::UseGpuAwareMpi()) 
    {
        rcv_buffer = static_cast<char*>(amrex::The_Device_Arena()->alloc(TotRcvBytes));
    }
    else 
    {
        rcv_buffer = static_cast<char*>(amrex::The_Pinned_Arena()->alloc(TotRcvBytes));
    }
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = Rcvs[Who];
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv(rcv_buffer + offset,
                                             Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (const auto& kv : not_ours) {
        const auto Who = kv.first;
        const auto Cnt = kv.second.size();

        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send(thrust::raw_pointer_cast(kv.second.data()),
                                 Cnt, Who, SeqNum);
    }

    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);

        for (int i = 0; i < nrcvs; ++i) {
            const int offset = rOffset[i];
            char* buffer = thrust::raw_pointer_cast(rcv_buffer + offset);
            size_t num_grids, num_particles;
            int gid, pid;
            cudaMemcpy(&num_grids, buffer, sizeof(size_t), cudaMemcpyHostToHost);
            buffer += sizeof(size_t);

            for (int g = 0; g < num_grids; ++g) {
                cudaMemcpyAsync(&num_particles, buffer, sizeof(size_t), cudaMemcpyHostToHost);
                buffer += sizeof(size_t);
                cudaMemcpyAsync(&gid, buffer, sizeof(int), cudaMemcpyHostToHost);
                buffer += sizeof(int);
                cudaMemcpyAsync(&pid, buffer, sizeof(int), cudaMemcpyHostToHost);
                buffer += sizeof(int);

                Gpu::Device::streamSynchronize();

                if (num_particles == 0) continue;

                AMREX_ALWAYS_ASSERT(pid == ParallelDescriptor::MyProc());
                {
                    const int tid = 0;
                    auto pair_index = std::make_pair(gid, tid);
                    auto& ptile = this->GetParticles(lev)[pair_index];
                    const int nRParticles = ptile.numRealParticles();
                    const int nNParticles = ptile.getNumNeighbors();
                    const int new_num_neighbors = nNParticles + num_particles;
                    ptile.setNumNeighbors(new_num_neighbors);

                    //copy structs
                    auto& aos = ptile.GetArrayOfStructs();
                    cudaMemcpyAsync(static_cast<ParticleType*>(aos().data()) + nRParticles + nNParticles,
                                    buffer, num_particles*sizeof(ParticleType),
                                    cudaMemcpyHostToDevice);
                    buffer += num_particles*sizeof(ParticleType);
                }
            }
        }
    }

    if (ParallelDescriptor::UseGpuAwareMpi())
    {
        amrex::The_Device_Arena()->free(rcv_buffer);
    } else {
        amrex::The_Pinned_Arena()->free(rcv_buffer);
    }

#endif // MPI    
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
clearNeighborsGPU()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::clearNeighborsGPU");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    const int lev = 0;

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int src_grid = mfi.index();
        int src_tile = mfi.LocalTileIndex();
        AMREX_ASSERT(src_tile == 0);
        auto index = std::make_pair(src_grid, src_tile);
        auto& ptile = this->GetParticles(lev)[index];
        ptile.setNumNeighbors(0);
    }
}

template <int NStructReal, int NStructInt>
template <class CheckPair>
void 
NeighborParticleContainer<NStructReal, NStructInt>::
buildNeighborListGPU(CheckPair check_pair)
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::buildNeighborListGPU");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    const int lev = 0;
    const Geometry& geom = this->Geom(lev);
    auto& plev  = this->GetParticles(lev);

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();        
        auto index = std::make_pair(gid, tid);

        auto& ptile = plev[index];
        auto& aos   = ptile.GetArrayOfStructs();

        Box bx = mfi.tilebox();
        bx.grow(m_num_neighbor_cells);

        m_neighbor_list[index].build(aos(), bx, geom, check_pair);
    }
}

template <int NStructReal, int NStructInt>
void 
NeighborParticleContainer<NStructReal, NStructInt>::
printNeighborListGPU()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::printNeighborListGPU");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    const int lev = 0;
    const Geometry& geom = this->Geom(lev);
    auto& plev  = this->GetParticles(lev);

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        m_neighbor_list[index].print();
    }
}

#endif
