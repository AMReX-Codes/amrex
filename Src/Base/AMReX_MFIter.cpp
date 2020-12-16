
#include <AMReX_MFIter.H>
#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_OpenMP.H>

namespace amrex {

int MFIter::nextDynamicIndex = std::numeric_limits<int>::min();
int MFIter::depth = 0;
int MFIter::allow_multiple_mfiters = 0;

int
MFIter::allowMultipleMFIters (int allow)
{
    std::swap(allow, allow_multiple_mfiters);
    return allow;
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		unsigned char       flags_)
    :
    fabArray(fabarray_),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    streams(Gpu::numGpuStreams()),
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
    streams(Gpu::numGpuStreams()),
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
    streams(Gpu::numGpuStreams()),
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
    m_fa(new FabArrayBase(ba,dm,1,0)),
    fabArray(*m_fa),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef _OPENMP
#pragma omp single
#endif
    {
        m_fa->addThisBD();
    }
    Initialize();
}

MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm, bool do_tiling_)
    :
    m_fa(new FabArrayBase(ba,dm,1,0)),
    fabArray(*m_fa),
    tile_size((do_tiling_) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(do_tiling_ ? Tiling : 0),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef _OPENMP
#pragma omp single
#endif
    {
        m_fa->addThisBD();
    }
    Initialize();
}


MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm,
                const IntVect& tilesize_, unsigned char flags_)
    :
    m_fa(new FabArrayBase(ba,dm,1,0)),
    fabArray(*m_fa),
    tile_size(tilesize_),
    flags(flags_ | Tiling),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(true),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef _OPENMP
#pragma omp single
#endif
    {
        m_fa->addThisBD();
    }
    Initialize();
}


MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm, const MFItInfo& info)
    :
    m_fa(new FabArrayBase(ba, dm, 1, 0)),
    fabArray(*m_fa),
    tile_size(info.tilesize),
    flags(info.do_tiling ? Tiling : 0),
    streams(info.num_streams),
    dynamic(info.dynamic && (OpenMP::get_num_threads() > 1)),
    device_sync(info.device_sync),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef _OPENMP
#pragma omp single
#endif
    {
        m_fa->addThisBD();
    }
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
    streams(info.num_streams),
    dynamic(info.dynamic && (OpenMP::get_num_threads() > 1)),
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
#ifdef _OPENMP
#pragma omp master
#endif
    {
        depth = 0;
    }

#ifdef BL_USE_TEAM
    if ( ! (flags & NoTeamBarrier) )
	ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#ifdef AMREX_USE_GPU
    if (device_sync) Gpu::synchronize();
#endif

#ifdef AMREX_USE_GPU
    AMREX_GPU_ERROR_CHECK();
    Gpu::Device::resetStreamIndex();
    Gpu::resetNumCallbacks();
    if (!OpenMP::in_parallel() && Gpu::inFuseRegion()) {
        Gpu::LaunchFusedKernels();
    }
#endif

    if (m_fa) {
#ifdef _OPENMP
#pragma omp barrier
#pragma omp single
#endif
        m_fa->clearThisBD();
    }
}

void
MFIter::Initialize ()
{
#ifdef _OPENMP
#pragma omp master
#endif
    {
        ++depth;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(depth == 1 || MFIter::allow_multiple_mfiters,
            "Nested or multiple active MFIters is not supported by default.  This can be changed by calling MFIter::allowMultipleMFIters(true)".);
    }

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
	Gpu::Device::setStreamIndex((streams > 0) ? currentIndex%streams : -1);
        Gpu::resetNumCallbacks();
        if (!OpenMP::in_parallel()) {
            if (index_map->size() >= Gpu::getFuseNumKernelsThreshold()) {
                gpu_fsg.reset(new Gpu::FuseSafeGuard(true));
            }
        }
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
    const Box& vccbx = amrex::enclosedCells(validbox());
    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vccbx.smallEnd(d)) {
	    bx.growLo(d, ngrow[d]);
	}
	if (bx.bigEnd(d) >= vccbx.bigEnd(d)) {
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
    if (dir < 0) return tilebox(IntVect::TheNodeVector(), a_ng);
    return tilebox(IntVect::TheDimensionVector(dir), a_ng);
}

void
MFIter::operator++ () noexcept
{
#ifdef _OPENMP
    if (dynamic)
    {
#pragma omp atomic capture
        currentIndex = nextDynamicIndex++;
    }
    else
#endif
    {
        ++currentIndex;

#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            Gpu::Device::setStreamIndex((streams > 0) ? currentIndex%streams : -1);
            AMREX_GPU_ERROR_CHECK();
#ifdef AMREX_DEBUG
//            Gpu::synchronize();
#endif
        }
#endif
    }
}

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

    int tid = OpenMP::get_thread_num();
    int nthreads = OpenMP::get_num_threads();

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
