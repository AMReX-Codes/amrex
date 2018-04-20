
#include <AMReX_MFIter.H>
#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#ifdef AMREX_USE_DEVICE
#include <AMReX_Device.H>
#endif

namespace amrex {

#ifdef AMREX_USE_CUDA
int MFIter_init::m_cnt = 0;

namespace
{
    Arena* the_mfiter_arena = 0;
}

MFIter_init::MFIter_init ()
{
    if (m_cnt++ == 0)
    {
        BL_ASSERT(the_mfiter_arena == 0);

        the_mfiter_arena = new CArena;

	the_mfiter_arena->SetDeviceMemory();
    }
}

MFIter_init::~MFIter_init ()
{
    if (--m_cnt == 0)
        delete the_mfiter_arena;
}

Arena*
The_MFIter_Arena ()
{
    BL_ASSERT(the_mfiter_arena != 0);

    return the_mfiter_arena;
}
#endif

int MFIter::nextDynamicIndex = std::numeric_limits<int>::min();

MFIter::MFIter (const FabArrayBase& fabarray_, 
		unsigned char       flags_)
    :
    fabArray(fabarray_),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    dynamic(false),
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
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    Initialize();
}


MFIter::MFIter (const FabArrayBase& fabarray_, const MFItInfo& info)
    :
    fabArray(fabarray_),
    tile_size(info.tilesize),
    flags(info.do_tiling ? Tiling : 0),
    dynamic(info.dynamic),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
    if (dynamic) {
#ifdef _OPENMP
#pragma omp single
        nextDynamicIndex = omp_get_num_threads();
        // yes omp single has an implicit barrier and we need it because nextDynamicIndex is static.
#else
        dynamic = false;  // dynamic doesn't make sense if OMP is not used.
#endif
    }

    Initialize();
}


MFIter::~MFIter ()
{
#if BL_USE_TEAM
    if ( ! (flags & NoTeamBarrier) )
	ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#ifdef AMREX_USE_DEVICE
    Device::synchronize();
#endif

#ifdef AMREX_USE_CUDA
    reduce();
#endif

#ifdef AMREX_USE_CUDA
    for (int i = 0; i < real_reduce_list.size(); ++i)
        amrex::The_MFIter_Arena()->free(real_device_reduce_list[i]);
#endif

#ifdef AMREX_USE_DEVICE
    Device::check_for_errors();
    Device::set_stream_index(-1);
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

#ifdef AMREX_USE_DEVICE
	Device::set_stream_index(currentIndex);
#endif

	typ = fabArray.boxArray().ixType();
    }
}

Box 
MFIter::tilebox () const
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
    registerBox(bx);
    return bx;
}

Box
MFIter::tilebox (const IntVect& nodal) const
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
    registerBox(bx);
    return bx;
}

Box
MFIter::tilebox (const IntVect& nodal, const IntVect& ngrow) const
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
MFIter::nodaltilebox (int dir) const 
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
    registerBox(bx);
    return bx;
}

// Note that a small negative ng is supported.
Box 
MFIter::growntilebox (int ng) const 
{
    Box bx = tilebox();
    if (ng < -100) ng = fabArray.nGrow();
    const Box& vbx = validbox();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, ng);
	}
	if (bx.bigEnd(d) == vbx.bigEnd(d)) {
	    bx.growHi(d, ng);
	}
    }
    registerBox(bx);
    return bx;
}

Box
MFIter::grownnodaltilebox (int dir, int ng) const
{
    BL_ASSERT(dir < AMREX_SPACEDIM);
    Box bx = nodaltilebox(dir);
    if (ng < -100) ng = fabArray.nGrow();
    const Box& vbx = validbox();
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, ng);
	}
	if (bx.bigEnd(d) >= vbx.bigEnd(d)) {
	    bx.growHi(d, ng);
	}
    }
    registerBox(bx);
    return bx;
}

#ifndef _OPENMP
void
MFIter::operator++ () {

#ifdef AMREX_USE_CUDA
    if (real_reduce_list.size() == currentIndex + 1) {
        Device::device_dtoh_memcpy_async(&real_reduce_list[currentIndex],
                                         real_device_reduce_list[currentIndex],
                                         sizeof(Real));
    }
#endif

    ++currentIndex;

#ifdef AMREX_USE_DEVICE
    Device::set_stream_index(currentIndex);
    Device::check_for_errors();
#ifdef DEBUG
    Device::synchronize();
#endif
#endif

}
#endif

#ifdef AMREX_USE_CUDA
Real*
MFIter::add_reduce_value(Real* val, MFReducer r)
{

    real_reduce_val = val;

    reducer = r;

    Real reduce_val = *val;
    real_reduce_list.push_back(reduce_val);

    Real* dval = static_cast<Real*>(amrex::The_MFIter_Arena()->alloc(sizeof(Real)));
    real_device_reduce_list.push_back(dval);

    Device::device_htod_memcpy_async(real_device_reduce_list[currentIndex],
                                     &real_reduce_list[currentIndex],
                                     sizeof(Real));

    return dval;

}
#endif

#ifdef AMREX_USE_CUDA
// Reduce over the values in the list.
void
MFIter::reduce()
{

    // Do nothing if we don't have enough values to reduce on.

    if (real_reduce_list.empty()) return;

    if (real_reduce_list.size() < length()) return;

    Real result;

    if (reducer == MFReducer::SUM) {
        result = 0.0;
        for (int i = 0; i < real_reduce_list.size(); ++i) {
            result += real_reduce_list[i];
        }
    }
    else if (reducer == MFReducer::MIN) {
        result = 1.e200;
        for (int i = 0; i < real_reduce_list.size(); ++i) {
            result = std::min(result, real_reduce_list[i]);
        }
    }
    else if (reducer == MFReducer::MAX) {
        result = -1.e200;
        for (int i = 0; i < real_reduce_list.size(); ++i) {
            result = std::max(result, real_reduce_list[i]);
        }
    }

    *real_reduce_val = result;

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
