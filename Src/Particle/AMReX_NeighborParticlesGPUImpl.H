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
		    ++num_codes;
                }
   	        m_code_offsets[grid].push_back(m_code_array[grid].size());
            }        
	}    
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
computeHaloParticles ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::computeHaloParticles()");

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

	Gpu::DeviceVector<int> counts(np, 0);
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

	neighbor_copy_op.resize(gid, num_copies);

	auto p_boxes = neighbor_copy_op.m_boxes[gid].dataPtr();
	auto p_src_indices = neighbor_copy_op.m_src_indices[gid].dataPtr();
	auto p_periodic_shift = neighbor_copy_op.m_periodic_shift[gid].dataPtr();

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
clearCopyPlan ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::clearCopyPlan");
    clearNeighbors();
    neighbor_copy_op.clear();
    neighbor_copy_plan.clear();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
computeCopyPlan ()
{
    neighbor_copy_plan.build(*this, neighbor_copy_op);
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
packBuffer ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::packBuffer");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    const int lev = 0;
    const auto& geom = this->Geom(lev);
    auto& plev  = this->GetParticles(lev);
    auto& ba = this->ParticleBoxArray(lev);
    auto& dmap = this->ParticleDistributionMap(lev);
    const auto plo = this->Geom(lev).ProbLoArray();
    const auto phi = this->Geom(lev).ProbHiArray();
    const auto is_per = this->Geom(lev).isPeriodicArray();

    int num_boxes = ba.size();
    if (num_boxes == 1 and (not geom.isAnyPeriodic()) ) return;
    
    int total_buffer_size = neighbor_copy_plan.m_box_offsets[num_boxes];

    snd_buffer.resize(total_buffer_size);    

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& src_tile = plev[index];
        auto& aos   = src_tile.GetArrayOfStructs();
        ParticleType* p_ptr = &(aos[0]);

	int num_copies = neighbor_copy_op.m_boxes[gid].size();

	auto p_boxes = neighbor_copy_op.m_boxes[gid].dataPtr();
	auto p_src_indices = neighbor_copy_op.m_src_indices[gid].dataPtr();
	auto p_dst_indices = neighbor_copy_plan.m_dst_indices[gid].dataPtr();
	auto p_periodic_shift = neighbor_copy_op.m_periodic_shift[gid].dataPtr();
	auto p_box_offsets = neighbor_copy_plan.m_box_offsets.dataPtr();
        auto p_box_perm = this->BufferMap().gridToBucketPtr();
	auto p_snd_buffer = snd_buffer.dataPtr();

	AMREX_FOR_1D ( num_copies, i,
        {
	  int box_offset = p_box_offsets[p_box_perm[p_boxes[i]]];
	  int dst_index = box_offset + p_dst_indices[i];
	  int src_index = p_src_indices[i];
	  p_snd_buffer[dst_index] = p_ptr[src_index];

	  ParticleType& p = p_snd_buffer[dst_index];
	  const IntVect& pshift = p_periodic_shift[i];
	  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
	  {
	    if (not is_per[idim]) continue;
	    if (pshift[idim] > 0) 
	      p.pos(idim) += phi[idim] - plo[idim]; 
	    else if (pshift[idim] < 0)
	      p.pos(idim) -= phi[idim] - plo[idim];
	  }
	});
    }
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
unpackBuffer ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::unpackBuffer");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    const int lev = 0;
    const auto& geom = this->Geom(lev);
    auto& plev  = this->GetParticles(lev);
    auto& ba = this->ParticleBoxArray(lev);
    auto& dmap = this->ParticleDistributionMap(lev);

    int num_boxes = ba.size();
    if (num_boxes == 1 and (not geom.isAnyPeriodic()) ) return;
    
    // local unpack
    auto p_dst_box_counts = neighbor_copy_plan.m_box_counts.dataPtr();
    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& tile = plev[index];
        auto& aos   = tile.GetArrayOfStructs();
	auto np = tile.numTotalParticles();

	int num_copies = neighbor_copy_plan.m_box_counts[this->BufferMap().gridToBucket(gid)];

	auto p_box_offsets = neighbor_copy_plan.m_box_offsets.dataPtr();
	auto p_box_perm = this->BufferMap().gridToBucketPtr();
	auto p_snd_buffer = snd_buffer.dataPtr();

        int nn = tile.getNumNeighbors();
        tile.setNumNeighbors(nn + num_copies);	

        ParticleType* p_ptr = aos().dataPtr();

        AMREX_FOR_1D ( num_copies, i,
        {
	  int box_offset = p_box_offsets[p_box_perm[gid]];
	  int src_index = box_offset + i;
	  int dst_index = np + i;
	  p_ptr[dst_index] = p_snd_buffer[src_index];
        });
    }    
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
fillNeighborsGPU ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::fillNeighbors");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    clearNeighbors();
    clearCopyPlan();
    buildNeighborMask();
    this->defineBufferMap();
    computeHaloParticles();
    computeCopyPlan();
    getBoxRcvCounts();
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
    packBuffer();
    unpackBuffer();
    communicateParticles();
}



template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
getBoxRcvCounts ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::getBoxRcvCounts()");
    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    const int MyProc = ParallelDescriptor::MyProc();
    
    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // number of particles per snd / recieve

    std::map<int, Vector<int> > snd_data;
    long NumSnds = 0;
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        auto box_buffer_indices = this->BufferMap().allBucketsOnProc(i);
        int nboxes = box_buffer_indices.size();
        for (auto bucket : box_buffer_indices)
	{
            int dst = this->BufferMap().bucketToGrid(bucket);
            int npart = neighbor_copy_plan.m_box_counts[bucket];
            snd_data[i].push_back(npart);
            snd_data[i].push_back(dst);
	}
	long nbytes = 2*nboxes*sizeof(int);
	Snds[i] = nbytes;
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
    Vector<std::size_t> rOffset;
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes/sizeof(int));
            TotRcvBytes += Rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();

    Gpu::ManagedDeviceVector<int> rcv_data(TotRcvBytes/sizeof(int));

    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = Rcvs[Who];
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv((char*) thrust::raw_pointer_cast(&rcv_data[offset]),
                                             Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        const auto Who = i;
        const auto Cnt = Snds[i];

        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        ParallelDescriptor::Send((char*) thrust::raw_pointer_cast(snd_data[i].data()),
                                 Cnt, Who, SeqNum);
    }

    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);
	
	neighbor_copy_plan.m_rcv_box_offsets.push_back(0);
	for (int i = 0; i < rcv_data.size(); i +=2)
	{
	  neighbor_copy_plan.m_rcv_box_counts.push_back(rcv_data[i]);
	  AMREX_ASSERT(MyProc == this->ParticleDistributionMap(0)[neighbor_copy_plan.m_rcv_data[i+1]]);
	  neighbor_copy_plan.m_rcv_box_ids.push_back(rcv_data[i+1]);
	  neighbor_copy_plan.m_rcv_box_offsets.push_back(neighbor_copy_plan.m_rcv_box_offsets.back() + neighbor_copy_plan.m_rcv_box_counts.back());
	}
    }
    
#endif // MPI
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
communicateParticles ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::communicateParticles()");
    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

#ifdef BL_USE_MPI
    const int NProcs = ParallelDescriptor::NProcs();
    const int MyProc = ParallelDescriptor::MyProc();
    
    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        auto box_buffer_indices = this->BufferMap().allBucketsOnProc(i);
        int nboxes = box_buffer_indices.size();
        long nbytes = 0;
        for (auto bucket : box_buffer_indices)
        {
            nbytes += neighbor_copy_plan.m_box_counts[bucket]*sizeof(ParticleType); 
        }
	Snds[i] = nbytes;
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
    Vector<std::size_t> rOffset;
    
    std::size_t TotRcvBytes = 0;
    for (int i = 0; i < NProcs; ++i) {
        if (Rcvs[i] > 0) {
            RcvProc.push_back(i);
            rOffset.push_back(TotRcvBytes/sizeof(ParticleType));
            TotRcvBytes += Rcvs[i];
        }
    }
    
    const int nrcvs = RcvProc.size();
    Vector<MPI_Status>  stats(nrcvs);
    Vector<MPI_Request> rreqs(nrcvs);

    const int SeqNum = ParallelDescriptor::SeqNum();

    rcv_buffer.resize(TotRcvBytes/sizeof(ParticleType));
    
    // Post receives.
    for (int i = 0; i < nrcvs; ++i) {
        const auto Who    = RcvProc[i];
        const auto offset = rOffset[i];
        const auto Cnt    = Rcvs[Who];
        
        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        BL_ASSERT(Who >= 0 && Who < NProcs);
        
        rreqs[i] = ParallelDescriptor::Arecv((char*) thrust::raw_pointer_cast(&rcv_buffer[offset]),
                                             Cnt, Who, SeqNum).req();
    }
    
    // Send.
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
        const auto Who = i;
        const auto Cnt = Snds[i];

        BL_ASSERT(Cnt > 0);
        BL_ASSERT(Who >= 0 && Who < NProcs);
        BL_ASSERT(Cnt < std::numeric_limits<int>::max());
        
        int bucket = this->BufferMap().firstBucketOnProc(i);
        const auto snd_offset = neighbor_copy_plan.m_box_offsets[bucket];

        ParallelDescriptor::Send(
	    (char*) thrust::raw_pointer_cast(&snd_buffer[snd_offset]),
	    Cnt, Who, SeqNum);
    }

    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);

	const int lev = 0;
	auto& plev  = this->GetParticles(lev);

	auto p_rcv_box_counts = neighbor_copy_plan.m_rcv_box_counts.dataPtr();
	auto p_rcv_box_offsets = neighbor_copy_plan.m_rcv_box_offsets.dataPtr();
	auto p_rcv_box_ids = neighbor_copy_plan.m_rcv_box_ids.dataPtr();
	auto p_rcv_buffer = rcv_buffer.dataPtr();

	for (int i = 0; i < neighbor_copy_plan.m_rcv_box_counts.size(); ++i)
	{
	  int copy_size = neighbor_copy_plan.m_rcv_box_counts[i];
	  int gid = neighbor_copy_plan.m_rcv_box_ids[i];
	  int tid = 0;
	  int offset = neighbor_copy_plan.m_rcv_box_offsets[i];

	  auto index = std::make_pair(gid, tid);
	  
	  auto& tile = plev[index];
	  auto& aos   = tile.GetArrayOfStructs();
	  auto np = tile.numTotalParticles();

	  AMREX_ALWAYS_ASSERT(MyProc == this->ParticleDistributionMap(0)[gid]);

	  int nn = tile.getNumNeighbors();
	  tile.setNumNeighbors(nn + neighbor_copy_plan.m_rcv_box_counts[i]);	

	  ParticleType* p_ptr = &(aos[0]);
	  AMREX_FOR_1D ( copy_size, ip,
	  {
	    int src_index = offset + ip;
	    int dst_index = np + ip;
	    p_ptr[dst_index] = p_rcv_buffer[src_index];
	  });	 
	}	
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
