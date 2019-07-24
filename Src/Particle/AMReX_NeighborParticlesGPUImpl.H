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
        
	// for (int i = 0; i < m_code_offsets[grid].size()-1; ++i)
	// {
	//   int start = m_code_offsets[grid][i];
	//   int stop  = m_code_offsets[grid][i+1];
	//   amrex::Print() << "Code " << i << " corresponds to: ";
	//   for (int j = start; j < stop; ++j)
	//     amrex::Print() << m_code_array[grid][j].grid_id << " ";
	//   amrex::Print() << "\n";
	// }
	// amrex::Print() << "\n";
	}
    
	BuildMasks();
	GetNeighborCommTags();

	int num_boxes = ba.size();
	bucket_to_gid.resize(num_boxes);
        gid_to_bucket.resize(num_boxes);

        Gpu::DeviceVector<int> box_proc_ids(num_boxes);
	Gpu::thrust_copy(dmap.ProcessorMap().begin(), dmap.ProcessorMap().end(), box_proc_ids.begin());

	thrust::sequence(thrust::device, bucket_to_gid.begin(), bucket_to_gid.end());
	
	thrust::sort_by_key(thrust::cuda::par(Cuda::The_ThrustCachedAllocator()),
			    box_proc_ids.begin(),
			    box_proc_ids.end(),
			    bucket_to_gid.begin());

        auto b_to_gid_ptr = bucket_to_gid.dataPtr();
        auto gid_to_b_ptr = gid_to_bucket.dataPtr();
        AMREX_FOR_1D ( num_boxes, i,
	{
            gid_to_b_ptr[b_to_gid_ptr[i]] = i;
        });	 

	proc_box_counts.resize(ParallelDescriptor::NProcs(), 0);
	proc_box_offsets.resize(ParallelDescriptor::NProcs());
	for (auto& val : dmap.ProcessorMap() )
        {
	    proc_box_counts[val]++;
	} 

        amrex::Gpu::exclusive_scan(proc_box_counts.begin(), proc_box_counts.end(),
                                   proc_box_offsets.begin());
 
	proc_box_offsets.push_back(num_boxes);
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
    auto& dmap = this->ParticleDistributionMap(lev);

    int num_boxes = ba.size();
    if (num_boxes == 1 and (not geom.isAnyPeriodic()) ) return;

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

	copy_plan_boxes[gid].resize(num_copies);
	copy_plan_src_indices[gid].resize(num_copies);
	copy_plan_dst_indices[gid].resize(num_copies);
	copy_plan_periodic_shift[gid].resize(num_copies);

	auto p_copy_plan_boxes = copy_plan_boxes[gid].dataPtr();
	auto p_copy_plan_src_indices = copy_plan_src_indices[gid].dataPtr();
	auto p_copy_plan_periodic_shift = copy_plan_periodic_shift[gid].dataPtr();

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
		  p_copy_plan_boxes[index] = p_code_array[j].grid_id;
		  p_copy_plan_periodic_shift[index] = p_code_array[j].periodic_shift;
		  p_copy_plan_src_indices[index] = i;
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

    box_counts.resize(0);
    box_offsets.resize(0);

    copy_plan_boxes.clear();
    copy_plan_src_indices.clear();
    copy_plan_dst_indices.clear();
    copy_plan_periodic_shift.clear();
}

template <int NStructReal, int NStructInt>
void
NeighborParticleContainer<NStructReal, NStructInt>::
computeCopyPlan ()
{
    BL_PROFILE("NeighborParticleContainer<NStructReal, NStructInt>::computeCopyPlan");

    AMREX_ALWAYS_ASSERT(this->numLevels() == 1);

    AMREX_ASSERT(hasNeighbors() == false);

    const int lev = 0;
    const auto& geom = this->Geom(lev);
    const auto dxi = this->Geom(lev).InvCellSizeArray();
    const auto plo = this->Geom(lev).ProbLoArray();
    const auto domain = this->Geom(lev).Domain();
    auto& plev  = this->GetParticles(lev);
    auto& ba = this->ParticleBoxArray(lev);
    auto& dmap = this->ParticleDistributionMap(lev);

    int num_boxes = ba.size();
    if (num_boxes == 1 and (not geom.isAnyPeriodic()) ) return;

    box_counts.resize(0);
    box_counts.resize(num_boxes, 0);
    box_offsets.resize(num_boxes+1);
    auto p_dst_box_counts = box_counts.dataPtr();
    auto p_box_perm = gid_to_bucket.dataPtr();

    constexpr unsigned int max_unsigned_int = std::numeric_limits<unsigned int>::max();

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
	int num_copies = copy_plan_boxes[gid].size();

	auto p_copy_plan_boxes = copy_plan_boxes[gid].dataPtr();
	auto p_copy_plan_dst_indices = copy_plan_dst_indices[gid].dataPtr();

	AMREX_FOR_1D ( num_copies, i,
        {
	    int dst_box = p_copy_plan_boxes[i];
	    int index = Gpu::Atomic::Inc(&p_dst_box_counts[p_box_perm[dst_box]], max_unsigned_int);
	    p_copy_plan_dst_indices[i] = index;
	});
    }

    amrex::Gpu::exclusive_scan(box_counts.begin(), box_counts.end()+1, box_offsets.begin());
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
    
    int total_buffer_size = box_offsets[num_boxes];

    snd_buffer.resize(total_buffer_size);    

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& src_tile = plev[index];
        auto& aos   = src_tile.GetArrayOfStructs();
        ParticleType* p_ptr = &(aos[0]);

	int num_copies = copy_plan_boxes[gid].size();

	auto p_copy_plan_boxes = copy_plan_boxes[gid].dataPtr();
	auto p_copy_plan_src_indices = copy_plan_src_indices[gid].dataPtr();
	auto p_copy_plan_dst_indices = copy_plan_dst_indices[gid].dataPtr();
	auto p_copy_plan_periodic_shift = copy_plan_periodic_shift[gid].dataPtr();
	auto p_box_offsets = box_offsets.dataPtr();
        auto p_box_perm = gid_to_bucket.dataPtr();
	auto p_snd_buffer = snd_buffer.dataPtr();

	AMREX_FOR_1D ( num_copies, i,
        {
	  int box_offset = p_box_offsets[p_box_perm[p_copy_plan_boxes[i]]];
	  int dst_index = box_offset + p_copy_plan_dst_indices[i];
	  int src_index = p_copy_plan_src_indices[i];
	  p_snd_buffer[dst_index] = p_ptr[src_index];

	  ParticleType& p = p_snd_buffer[dst_index];
	  const IntVect& pshift = p_copy_plan_periodic_shift[i];
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
    // for (int i = 0; i < snd_buffer.size(); ++i)
    // {
    //   amrex::Print() << snd_buffer[i] << "\n";
    // }
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
    
    auto p_dst_box_counts = box_counts.dataPtr();

    for(MFIter mfi = this->MakeMFIter(lev); mfi.isValid(); ++mfi)
    {
        int gid = mfi.index();
        int tid = mfi.LocalTileIndex();
        auto index = std::make_pair(gid, tid);

        auto& tile = plev[index];
        auto& aos   = tile.GetArrayOfStructs();
	auto np = tile.numTotalParticles();
	
	int num_copies = box_counts[gid_to_bucket[gid]];

	auto p_box_offsets = box_offsets.dataPtr();
	auto p_box_perm = gid_to_bucket.dataPtr();
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
	long nboxes = proc_box_offsets[i+1] - proc_box_offsets[i];
	for (int j = 0; j < nboxes; ++j)
	{
	  int dst = bucket_to_gid[proc_box_offsets[i]+j];
	  int npart = box_counts[proc_box_offsets[i]+j];
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
	
	rcv_box_offsets.push_back(0);
	for (int i = 0; i < rcv_data.size(); i +=2)
	{
	  rcv_box_counts.push_back(rcv_data[i]);
	  AMREX_ASSERT(MyProc == this->ParticleDistributionMap(0)[rcv_data[i+1]]);
	  rcv_box_ids.push_back(rcv_data[i+1]);
	  rcv_box_offsets.push_back(rcv_box_offsets.back() + rcv_box_counts.back());
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
    
    std::map<int, std::vector<int> > boxes_rcv;
    for (int i = 0; i < local_neighbors.size(); ++i)
    {
        int proc = local_neighbors[i].proc_id;
	if (proc == MyProc) continue;
	int grid = local_neighbors[i].grid_id;
	boxes_rcv[proc].push_back(grid);
    }

    // We may now have particles that are rightfully owned by another CPU.
    Vector<long> Snds(NProcs, 0), Rcvs(NProcs, 0);  // bytes!

    long NumSnds = 0;
    for (int i = 0; i < NProcs; ++i)
    {
        if (i == MyProc) continue;
	long nboxes = proc_box_offsets[i+1] - proc_box_offsets[i];
        long nbytes = 0;
	for (int j = 0; j < nboxes; ++j)
        {
            nbytes += box_counts[proc_box_offsets[i]+j]*sizeof(ParticleType); 
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
        
        const auto snd_offset = box_offsets[proc_box_offsets[i]];

        ParallelDescriptor::Send(
	    (char*) thrust::raw_pointer_cast(&snd_buffer[snd_offset]),
	    Cnt, Who, SeqNum);
    }

    if (nrcvs > 0) {
        ParallelDescriptor::Waitall(rreqs, stats);

	const int lev = 0;
	auto& plev  = this->GetParticles(lev);

	auto p_rcv_box_counts = rcv_box_counts.dataPtr();
	auto p_rcv_box_offsets = rcv_box_offsets.dataPtr();
	auto p_rcv_box_ids = rcv_box_ids.dataPtr();
	auto p_rcv_buffer = rcv_buffer.dataPtr();

	for (int i = 0; i < rcv_box_counts.size(); ++i)
	{
	  int copy_size = rcv_box_counts[i];
	  int gid = rcv_box_ids[i];
	  int tid = 0;
	  int offset = rcv_box_offsets[i];

	  auto index = std::make_pair(gid, tid);
	  
	  auto& tile = plev[index];
	  auto& aos   = tile.GetArrayOfStructs();
	  auto np = tile.numTotalParticles();

	  AMREX_ALWAYS_ASSERT(MyProc == this->ParticleDistributionMap(0)[gid]);

	  int nn = tile.getNumNeighbors();
	  tile.setNumNeighbors(nn + rcv_box_counts[i]);	

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
