
#include <AMReX_MFIter.H>
#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>

namespace amrex {

int MFIter::nextDynamicIndex = std::numeric_limits<int>::min();

MFIter::MFIter (const FabArrayBase& fabarray_, 
		unsigned char       flags_)
    :
    fabArray(fabarray_),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    Initialize();
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		bool                do_tiling_)
    :
    fabArray(fabarray_),
    tile_size((do_tiling_) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(do_tiling_ ? Tiling : 0),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    Initialize();
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		const IntVect&      tilesize_, 
		unsigned char       flags_)
    :
    fabArray(fabarray_),
    tile_size(tilesize_),
    flags(flags_ | Tiling),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    Initialize();
}

MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm, unsigned char flags_)
    :
    m_fa(new FabArray<FArrayBox>(ba, dm, 1, 0,
                                 MFInfo().SetAlloc(false),
                                 FArrayBoxFactory())),
    fabArray(*m_fa),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    Initialize();
}

MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm, bool do_tiling_)
    :
    m_fa(new FabArray<FArrayBox>(ba, dm, 1, 0,
                                 MFInfo().SetAlloc(false),
                                 FArrayBoxFactory())),
    fabArray(*m_fa),
    tile_size((do_tiling_) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(do_tiling_ ? Tiling : 0),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    Initialize();
}


MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm,
                const IntVect& tilesize_, unsigned char flags_)
    :
    m_fa(new FabArray<FArrayBox>(ba, dm, 1, 0,
                                 MFInfo().SetAlloc(false),
                                 FArrayBoxFactory())),
    fabArray(*m_fa),
    tile_size(tilesize_),
    flags(flags_ | Tiling),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    Initialize();
}


MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm, const MFItInfo& info)
    :
    m_fa(new FabArray<FArrayBox>(ba, dm, 1, 0,
                                 MFInfo().SetAlloc(false),
                                 FArrayBoxFactory())),
    fabArray(*m_fa),
    tile_size(info.tilesize),
    flags(info.do_tiling ? Tiling : 0),
#ifdef _OPENMP
    dynamic(info.dynamic && (omp_get_num_threads() > 1)),
#else
    dynamic(false),
#endif
    device_sync(info.device_sync),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef _OPENMP
    if (dynamic) {
#pragma omp barrier
#pragma omp single
        nextDynamicIndex = omp_get_num_threads();
        // yes omp single has an implicit barrier and we need it because nextDynamicIndex is static.
    }
#endif

    Initialize();
}

MFIter::MFIter (const FabArrayBase& fabarray_, const MFItInfo& info)
    :
    fabArray(fabarray_),
    tile_size(info.tilesize),
    flags(info.do_tiling ? Tiling : 0),
#ifdef _OPENMP
    dynamic(info.dynamic && (omp_get_num_threads() > 1)),
#else
    dynamic(false),
#endif
    device_sync(info.device_sync),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef _OPENMP
    if (dynamic) {
#pragma omp barrier
#pragma omp single
        nextDynamicIndex = omp_get_num_threads();
        // yes omp single has an implicit barrier and we need it because nextDynamicIndex is static.
    }
#endif

    Initialize();
}


MFIter::~MFIter ()
{
#ifdef BL_USE_TEAM
    if ( ! (flags & NoTeamBarrier) )
	ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#ifdef AMREX_USE_GPU
    if (device_sync) Gpu::Device::synchronize();
#endif

#ifdef AMREX_USE_GPU
    reduce();
#endif

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        for (int i = 0; i < real_reduce_list.size(); ++i)
            for (int j = 0; j < real_reduce_list[i].size(); ++j)
                amrex::The_MFIter_Arena()->free(real_device_reduce_list[i][j]);
    }
#endif

#ifdef AMREX_USE_GPU
    AMREX_GPU_ERROR_CHECK();
    Gpu::Device::resetStreamIndex();
#endif
}

void 
MFIter::Initialize ()
{
    if (flags & SkipInit) {
	return;
    }
    else if (flags & AllBoxes)  // a very special case
    {
	index_map    = &(fabArray.IndexArray());
	currentIndex = 0;
	beginIndex   = 0;
	endIndex     = index_map->size();
    }
    else
    {
	const FabArrayBase::TileArray* pta = fabArray.getTileArray(tile_size);
	
	index_map            = &(pta->indexMap);
	local_index_map      = &(pta->localIndexMap);
	tile_array           = &(pta->tileArray);
	local_tile_index_map = &(pta->localTileIndexMap);
	num_local_tiles      = &(pta->numLocalTiles);

	{
	    int rit = 0;
	    int nworkers = 1;
#ifdef BL_USE_TEAM
	    if (ParallelDescriptor::TeamSize() > 1) {
		if ( tile_size == IntVect::TheZeroVector() ) {
		    // In this case the TileArray contains only boxes owned by this worker.
		    // So there is no sharing going on.
		    rit = 0;
		    nworkers = 1;
		} else {
		    rit = ParallelDescriptor::MyRankInTeam();
		    nworkers = ParallelDescriptor::TeamSize();
		}
	    }
#endif

	    int ntot = index_map->size();
	    
	    if (nworkers == 1)
	    {
		beginIndex = 0;
		endIndex = ntot;
	    }
	    else
	    {
		int nr   = ntot / nworkers;
		int nlft = ntot - nr * nworkers;
		if (rit < nlft) {  // get nr+1 items
		    beginIndex = rit * (nr + 1);
		    endIndex = beginIndex + nr + 1;
		} else {           // get nr items
		    beginIndex = rit * nr + nlft;
		    endIndex = beginIndex + nr;
		}
	    }
	}
	
#ifdef _OPENMP
	int nthreads = omp_get_num_threads();
	if (nthreads > 1)
	{
            if (dynamic)
            {
                beginIndex = omp_get_thread_num();
            }
            else
            {
                int tid = omp_get_thread_num();
                int ntot = endIndex - beginIndex;
                int nr   = ntot / nthreads;
                int nlft = ntot - nr * nthreads;
                if (tid < nlft) {  // get nr+1 items
                    beginIndex += tid * (nr + 1);
                    endIndex = beginIndex + nr + 1;
                } else {           // get nr items
                    beginIndex += tid * nr + nlft;
                    endIndex = beginIndex + nr;
                }
            }
	}
#endif

	currentIndex = beginIndex;

#ifdef AMREX_USE_GPU
	Gpu::Device::setStreamIndex(currentIndex);
#endif

	typ = fabArray.boxArray().ixType();
    }
}

Box 
MFIter::tilebox () const noexcept
{ 
    BL_ASSERT(tile_array != 0);
    Box bx((*tile_array)[currentIndex]);
    if (! typ.cellCentered())
    {
	bx.convert(typ);
	const Box& vbx = validbox();
	const IntVect& Big = vbx.bigEnd();
	for (int d=0; d<AMREX_SPACEDIM; ++d) {
	    if (typ.nodeCentered(d)) { // validbox should also be nodal in d-direction.
		if (bx.bigEnd(d) < Big[d]) {
		    bx.growHi(d,-1);
		}
	    }
	}
    }
    return bx;
}

Box
MFIter::tilebox (const IntVect& nodal) const noexcept
{
    BL_ASSERT(tile_array != 0);
    Box bx((*tile_array)[currentIndex]);
    const IndexType new_typ {nodal};
    if (! new_typ.cellCentered())
    {
	bx.setType(new_typ);
	const Box& valid_cc_box = amrex::enclosedCells(validbox());
	const IntVect& Big = valid_cc_box.bigEnd();
	for (int d=0; d<AMREX_SPACEDIM; ++d) {
	    if (new_typ.nodeCentered(d)) { // validbox should also be nodal in d-direction.
		if (bx.bigEnd(d) == Big[d]) {
		    bx.growHi(d,1);
		}
	    }
	}
    }
    return bx;
}

Box
MFIter::tilebox (const IntVect& nodal, const IntVect& ngrow) const noexcept
{
    Box bx = tilebox(nodal);
    const Box& vbx = validbox();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, ngrow[d]);
	}
	if (bx.bigEnd(d) >= vbx.bigEnd(d)) {
	    bx.growHi(d, ngrow[d]);
	}
    }
    return bx;
}

Box
MFIter::nodaltilebox (int dir) const noexcept
{ 
    BL_ASSERT(dir < AMREX_SPACEDIM);
    BL_ASSERT(tile_array != 0);
    Box bx((*tile_array)[currentIndex]);
    bx.convert(typ);
    const Box& vbx = validbox();
    const IntVect& Big = vbx.bigEnd();
    int d0, d1;
    if (dir < 0) {
	d0 = 0;
	d1 = AMREX_SPACEDIM-1;
    } else {
	d0 = d1 = dir;
    }
    for (int d=d0; d<=d1; ++d) {
	if (typ.cellCentered(d)) { // validbox should also be cell-centered in d-direction.
	    bx.surroundingNodes(d);
	    if (bx.bigEnd(d) <= Big[d]) {
		bx.growHi(d,-1);
	    }
	}
    }
    return bx;
}

// Note that a small negative ng is supported.
Box 
MFIter::growntilebox (int a_ng) const noexcept
{
    Box bx = tilebox();
    IntVect ngv{a_ng};
    if (a_ng < -100) ngv = fabArray.nGrowVect();
    const Box& vbx = validbox();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, ngv[d]);
	}
	if (bx.bigEnd(d) == vbx.bigEnd(d)) {
	    bx.growHi(d, ngv[d]);
	}
    }
    return bx;
}

Box
MFIter::growntilebox (const IntVect& ng) const noexcept
{
    Box bx = tilebox();
    const Box& vbx = validbox();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, ng[d]);
	}
	if (bx.bigEnd(d) == vbx.bigEnd(d)) {
	    bx.growHi(d, ng[d]);
	}
    }
    return bx;
}

Box
MFIter::grownnodaltilebox (int dir, int a_ng) const noexcept
{
    IntVect ngv(a_ng);
    if (a_ng < -100) ngv = fabArray.nGrowVect();
    return grownnodaltilebox(dir, ngv);
}

Box
MFIter::grownnodaltilebox (int dir, IntVect const& a_ng) const noexcept
{
    BL_ASSERT(dir < AMREX_SPACEDIM);
    Box bx = nodaltilebox(dir);
    const Box& vbx = validbox();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, a_ng[d]);
	}
	if (bx.bigEnd(d) >= vbx.bigEnd(d)) {
	    bx.growHi(d, a_ng[d]);
	}
    }
    return bx;
}

void
MFIter::operator++ () noexcept
{
#ifdef _OPENMP
    int numOmpThreads = omp_get_num_threads();
#else
    int numOmpThreads = 1;
#endif

#ifdef _OPENMP
    if (dynamic)
    {
#pragma omp atomic capture
        currentIndex = nextDynamicIndex++;
    }
    else
#endif
    {
#ifdef AMREX_USE_GPU
        bool use_gpu = (numOmpThreads == 1) && Gpu::inLaunchRegion();
        if (use_gpu) {
            if (!real_reduce_list.empty()) {
                for (int i = 0; i < real_reduce_list[currentIndex].size(); ++i) {
                    Gpu::Device::dtoh_memcpy_async(&real_reduce_list[currentIndex][i],
                                                   real_device_reduce_list[currentIndex][i],
                                                   sizeof(Real));
                }
            }
        }
#endif

        ++currentIndex;

#ifdef AMREX_USE_GPU
        if (use_gpu) {
            Gpu::Device::setStreamIndex(currentIndex);
            AMREX_GPU_ERROR_CHECK();
#ifdef DEBUG
//            Gpu::Device::synchronize();
#endif
        }
#endif
    }
}

#ifdef AMREX_USE_GPU
Real*
MFIter::add_reduce_value(Real* val, MFReducer r)
{

    if (Gpu::inLaunchRegion()) {

        reducer = r;

        // For the reduce lists, the outer vector is length
        // (endIndex - beginIndex) and has a contribution from
        // every tile. The inner vector is all of the individual
        // variables to reduce over. While the former is a known
        // quantity, the latter is not, so we'll push_back for
        // each new quantity to reduce. Since the elements will
        // be added in the same order in every MFIter iteration,
        // we can access them in the same order in each entry
        // of the reduce_list.

        if (real_reduce_list.empty()) {
            real_reduce_list.resize(length());
            real_device_reduce_list.resize(length());
        }

        // Store the current value of the data.

        Real reduce_val = *val;
        real_reduce_list[currentIndex].push_back(reduce_val);

        // Create a device copy of the data to update within
        // the kernel.

        Real* dval = static_cast<Real*>(amrex::The_MFIter_Arena()->alloc(sizeof(Real)));
        real_device_reduce_list[currentIndex].push_back(dval);

        // Queue up a host to device copy of the input data,
        // so that we start from the correct value.

        const int list_idx = real_reduce_list[currentIndex].size() - 1;

        Gpu::Device::htod_memcpy_async(real_device_reduce_list[currentIndex][list_idx],
                                       &real_reduce_list[currentIndex][list_idx],
                                       sizeof(Real));

        // If we haven't already, store the address to the variable
        // we will update at the end.

        if (real_reduce_val.size() < real_reduce_list[currentIndex].size())
            real_reduce_val.push_back(val);

        return dval;

    }
    else {

        return val;

    }

}
#endif

#ifdef AMREX_USE_GPU
// Reduce over the values in the list.
void
MFIter::reduce()
{

    // Do nothing if we're not currently executing on the device.

    if (Gpu::notInLaunchRegion()) return;

    // Do nothing if we don't have enough values to reduce on.

    if (real_reduce_list.empty()) return;

    // Assume that the number of reductions we want is fixed
    // in each vector, and just grab the number from the first
    // entry.

    const int num_reductions = real_reduce_list[0].size();

    BL_ASSERT(real_reduce_list.size() == length());
    
    for (int j = 0; j < num_reductions; ++j) {

        Real result;

        if (reducer == MFReducer::SUM) result = 0.0e0;
        if (reducer == MFReducer::MIN) result = 1.e200;
        if (reducer == MFReducer::MAX) result = -1.e200;

        for (int i = 0; i < real_reduce_list.size(); ++i) {

            // Double check our assumption from above.

            BL_ASSERT(real_reduce_list[i].size() == num_reductions);

            if (reducer == MFReducer::SUM) {
                result += real_reduce_list[i][j];
            }
            else if (reducer == MFReducer::MIN) {
                result = std::min(result, real_reduce_list[i][j]);
            }
            else if (reducer == MFReducer::MAX) {
                result = std::max(result, real_reduce_list[i][j]);
            }

        }

        *(real_reduce_val[j]) = result;

    }

}
#endif

MFGhostIter::MFGhostIter (const FabArrayBase& fabarray)
    :
    MFIter(fabarray, (unsigned char)(SkipInit|Tiling))
{
    Initialize();
}

void
MFGhostIter::Initialize ()
{
    int rit = 0;
    int nworkers = 1;
#ifdef BL_USE_TEAM
    if (ParallelDescriptor::TeamSize() > 1) {
	rit = ParallelDescriptor::MyRankInTeam();
	nworkers = ParallelDescriptor::TeamSize();
    }
#endif

    int tid = 0;
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = omp_get_num_threads();
    if (nthreads > 1)
	tid = omp_get_thread_num();
#endif

    int npes = nworkers*nthreads;
    int pid = rit*nthreads+tid;

    BoxList alltiles;
    Vector<int> allindex;
    Vector<int> alllocalindex;

    for (int i=0; i < fabArray.IndexArray().size(); ++i) {
	int K = fabArray.IndexArray()[i];
	const Box& vbx = fabArray.box(K);
	const Box& fbx = fabArray.fabbox(K);

	const BoxList& diff = amrex::boxDiff(fbx, vbx);
	
	for (BoxList::const_iterator bli = diff.begin(); bli != diff.end(); ++bli) {
	    BoxList tiles(*bli, FabArrayBase::mfghostiter_tile_size);
	    int nt = tiles.size();
	    for (int it=0; it<nt; ++it) {
		allindex.push_back(K);
		alllocalindex.push_back(i);
	    }
	    alltiles.catenate(tiles);
	}
    }

    int n_tot_tiles = alltiles.size();
    int navg = n_tot_tiles / npes;
    int nleft = n_tot_tiles - navg*npes;
    int ntiles = navg;
    if (pid < nleft) ntiles++;

    // how many tiles should we skip?
    int nskip = pid*navg + std::min(pid,nleft);
    BoxList::const_iterator bli = alltiles.begin();
    for (int i=0; i<nskip; ++i) ++bli;

    lta.indexMap.reserve(ntiles);
    lta.localIndexMap.reserve(ntiles);
    lta.tileArray.reserve(ntiles);

    for (int i=0; i<ntiles; ++i) {
	lta.indexMap.push_back(allindex[i+nskip]);
	lta.localIndexMap.push_back(alllocalindex[i+nskip]);
	lta.tileArray.push_back(*bli++);
    }

    currentIndex = beginIndex = 0;
    endIndex = lta.indexMap.size();

    lta.nuse = 0;
    index_map       = &(lta.indexMap);
    local_index_map = &(lta.localIndexMap);
    tile_array      = &(lta.tileArray);
}

}
