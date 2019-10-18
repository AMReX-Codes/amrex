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
buildNeighborMask ()
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
	    int num_codes = 0;
       
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
	    m_code_offsets[grid].push_back(0);
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
		    m_code_array[grid].push_back(code);
                    neighbor_procs.push_back(dmap[code.grid_id]);
		    ++num_codes;
                }
   	        m_code_offsets[grid].push_back(m_code_array[grid].size());
            }        
	}
        RemoveDuplicates(neighbor_procs);
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
buildNeighborCopyOp ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::buildNeighborCopyOp()");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);
    AMREX_ASSERT(hasNeighbors() == false);

    const int lev = 0;
    const auto& geom = this->Geom(lev);
    const auto dxi = this->Geom(lev).InvCellSizeArray();
    const auto plo = this->Geom(lev).ProbLoArray();
    const auto domain = this->Geom(lev).Domain();
    auto& plev  = this->GetParticles(lev);
    auto& ba = this->ParticleBoxArray(lev);

    if (ba.size() == 1 and (not geom.isAnyPeriodic()) ) return;

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& src_tile = plev[index];
        auto& aos   = src_tile.GetArrayOfStructs();
        const size_t np = aos.numParticles();

        Array4<const int> const& mask_arr = m_neighbor_mask_ptr->array(mfi);

	Gpu::DeviceVector<int> counts(np + 1, 0);
	Gpu::DeviceVector<int> offsets(np + 1);
	auto p_counts = counts.dataPtr();
	auto p_offsets = offsets.dataPtr();

        ParticleType* p_ptr = &(aos[0]);
	auto p_code_array = m_code_array[gid].dataPtr();
	auto p_code_offsets = m_code_offsets[gid].dataPtr();	
	AMREX_FOR_1D ( np, i,
        {
	    IntVect iv = IntVect(
                AMREX_D_DECL(floor((p_ptr[i].pos(0)-plo[0])*dxi[0]),
                             floor((p_ptr[i].pos(1)-plo[1])*dxi[1]),
                             floor((p_ptr[i].pos(2)-plo[2])*dxi[2]))
                );
            
            iv += domain.smallEnd();
            
	    int code = mask_arr(iv);
	    if (code >= 0)
	    {
	      int start = p_code_offsets[code];
	      int stop =  p_code_offsets[code+1];
	      Gpu::Atomic::Add(&p_counts[i], stop-start);
	    }
	});

        amrex::Gpu::exclusive_scan(counts.begin(), counts.end() + 1, offsets.begin());

	int num_copies;
        Gpu::dtoh_memcpy(&num_copies, thrust::raw_pointer_cast(offsets.data())+np, sizeof(int));

	neighbor_copy_op.resize(gid, lev, num_copies);

	auto p_boxes = neighbor_copy_op.m_boxes[lev][gid].dataPtr();
        auto p_levs = neighbor_copy_op.m_levels[lev][gid].dataPtr();
	auto p_src_indices = neighbor_copy_op.m_src_indices[lev][gid].dataPtr();
	auto p_periodic_shift = neighbor_copy_op.m_periodic_shift[lev][gid].dataPtr();

	AMREX_FOR_1D ( np, i,
        {
	    IntVect iv = IntVect(
                AMREX_D_DECL(floor((p_ptr[i].pos(0)-plo[0])*dxi[0]),
                             floor((p_ptr[i].pos(1)-plo[1])*dxi[1]),
                             floor((p_ptr[i].pos(2)-plo[2])*dxi[2]))
                );
            
            iv += domain.smallEnd();
            	   
	    int code = mask_arr(iv);
	    if (code >= 0)
	    {
	        int k = 0;
	        int start = p_code_offsets[code];
	        int stop =  p_code_offsets[code+1];
	        for (int j = start; j < stop; ++j)
	        {
		  int index = p_offsets[i]+k;
		  p_boxes[index] = p_code_array[j].grid_id;
                  p_levs[index] = 0;
		  p_periodic_shift[index] = p_code_array[j].periodic_shift;
		  p_src_indices[index] = i;
		  ++k;
		}
	    }
	});		
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
fillNeighborsGPU ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::fillNeighbors");

    AMREX_ASSERT(this->numLevels() == 1);

    buildNeighborMask();
    this->defineBufferMap();

    neighbor_copy_op.clear();
    neighbor_copy_plan.clear();
    buildNeighborCopyOp();
    neighbor_copy_plan.build(*this, neighbor_copy_op, true);
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
        pinned_snd_buffer.resize(snd_buffer.size());
        cudaMemcpyAsync(pinned_snd_buffer.dataPtr(), snd_buffer.dataPtr(), snd_buffer.size()*sizeof(SuperParticleType), cudaMemcpyDeviceToHost);        
        neighbor_copy_plan.buildMPIFinish(this->BufferMap());
        communicateParticlesStart(*this, neighbor_copy_plan, pinned_snd_buffer, pinned_rcv_buffer);
        rcv_buffer.resize(pinned_rcv_buffer.size());
        unpackBuffer(*this, neighbor_copy_plan, snd_buffer, NeighborUnpackPolicy());
        communicateParticlesFinish(neighbor_copy_plan);
        cudaMemcpyAsync(rcv_buffer.dataPtr(), pinned_rcv_buffer.dataPtr(), pinned_rcv_buffer.size()*sizeof(SuperParticleType), cudaMemcpyHostToDevice);
        unpackRemotes(*this, neighbor_copy_plan, rcv_buffer, NeighborUnpackPolicy());
    }
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
