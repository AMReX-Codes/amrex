
#include <AMReX_MFIter.H>
#include <AMReX_FabArray.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_OpenMP.H>

#include <algorithm>

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

int
MFIter::currentDepth ()
{
    int r;
#ifdef AMREX_USE_OMP
#pragma omp atomic read
#endif
    r = MFIter::depth;
    return r;
}

MFIter::MFIter (const FabArrayBase& fabarray_,
                unsigned char       flags_)
    :
    fabArray(&fabarray_),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(!Gpu::inNoSyncRegion()),
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
    fabArray(&fabarray_),
    tile_size((do_tiling_) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(do_tiling_ ? Tiling : 0),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(!Gpu::inNoSyncRegion()),
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
    fabArray(&fabarray_),
    tile_size(tilesize_),
    flags(flags_ | Tiling),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(!Gpu::inNoSyncRegion()),
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
    m_fa(std::make_unique<FabArrayBase>(ba,dm,1,0)),
    fabArray(m_fa.get()),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(!Gpu::inNoSyncRegion()),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef AMREX_USE_OMP
#pragma omp single
#endif
    {
        m_fa->addThisBD();
    }
    Initialize();
}

MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm, bool do_tiling_)
    :
    m_fa(std::make_unique<FabArrayBase>(ba,dm,1,0)),
    fabArray(m_fa.get()),
    tile_size((do_tiling_) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(do_tiling_ ? Tiling : 0),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(!Gpu::inNoSyncRegion()),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef AMREX_USE_OMP
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
    m_fa(std::make_unique<FabArrayBase>(ba,dm,1,0)),
    fabArray(m_fa.get()),
    tile_size(tilesize_),
    flags(flags_ | Tiling),
    streams(Gpu::numGpuStreams()),
    dynamic(false),
    device_sync(!Gpu::inNoSyncRegion()),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef AMREX_USE_OMP
#pragma omp single
#endif
    {
        m_fa->addThisBD();
    }
    Initialize();
}


MFIter::MFIter (const BoxArray& ba, const DistributionMapping& dm, const MFItInfo& info)
    :
    m_fa(std::make_unique<FabArrayBase>(ba, dm, 1, 0)),
    fabArray(m_fa.get()),
    tile_size(info.tilesize),
    flags(info.do_tiling ? Tiling : 0),
    streams(std::max(1,std::min(Gpu::numGpuStreams(),info.num_streams))),
    dynamic(info.dynamic && (OpenMP::get_num_threads() > 1)),
    device_sync(info.device_sync),
    device_sync_pre(info.device_sync_pre),
    device_sync_post(info.device_sync_post),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef AMREX_USE_OMP
#pragma omp single
#endif
    {
        m_fa->addThisBD();
    }
#ifdef AMREX_USE_OMP
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
    fabArray(&fabarray_),
    tile_size(info.tilesize),
    flags(info.do_tiling ? Tiling : 0),
    streams(std::max(1,std::min(Gpu::numGpuStreams(),info.num_streams))),
    dynamic(info.dynamic && (OpenMP::get_num_threads() > 1)),
    device_sync(info.device_sync),
    device_sync_pre(info.device_sync_pre),
    device_sync_post(info.device_sync_post),
    index_map(nullptr),
    local_index_map(nullptr),
    tile_array(nullptr),
    local_tile_index_map(nullptr),
    num_local_tiles(nullptr)
{
#ifdef AMREX_USE_OMP
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
    Finalize();
}

void
MFIter::Finalize ()
{
    // avoid double finalize
    if (finalized) { return; }
    finalized = true;

    // mark as invalid
    currentIndex = endIndex;

#ifdef BL_USE_TEAM
    if ( ! (flags & NoTeamBarrier) )
        ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#ifdef AMREX_USE_GPU
    // Orders the work queued after the loop after the loop's work (on the
    // host if the loop synchronizes), and makes the stream that was current
    // before the loop current again.
    if (m_fan.isOpen()) {
        Gpu::theStreamSequencer().join(m_fan);
    } else if (!m_fan_tried && device_sync && device_sync_post && beginIndex < endIndex) {
        // An iterator with local boxes that was never iterated (e.g., one
        // used only to index a FabArray, or a FillPatchIterator used only for
        // its data) still returns with the GPU work queued so far completed.
        Gpu::streamSynchronize();
    }
    m_fan_tried.flag = true; // isValid() does not open the fan from now on
    AMREX_GPU_ERROR_CHECK();
#endif

    if (m_fa) {
#ifdef AMREX_USE_OMP
#pragma omp barrier
#pragma omp single
#endif
        m_fa->clearThisBD();
    }
    if (m_fa) {
        m_fa.reset(nullptr);
    }

    if (depth_counted) {  // only on the master thread
        depth_counted.flag = false;
#ifdef AMREX_USE_OMP
#pragma omp atomic update
#endif
        --depth;
    }
}

void
MFIter::Initialize ()
{
#ifdef AMREX_USE_OMP
#pragma omp master
#endif
    {
        int d;
#ifdef AMREX_USE_OMP
#pragma omp atomic capture
#endif
        d = ++depth;
        depth_counted.flag = true;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(d == 1 || MFIter::allow_multiple_mfiters,
            "Nested or multiple active MFIters is not supported by default.  This can be changed by calling MFIter::allowMultipleMFIters(true)");
    }
    if (flags & AllBoxes)  // a very special case
    {
        index_map    = &(fabArray->IndexArray());
        currentIndex = 0;
        beginIndex   = 0;
        endIndex     = static_cast<int>(index_map->size());
    }
    else
    {
        const FabArrayBase::TileArray* pta = fabArray->getTileArray(tile_size);

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

            int ntot = static_cast<int>(index_map->size());

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

#ifdef AMREX_USE_OMP
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

        typ = fabArray->boxArray().ixType();
    }
}

#ifdef AMREX_USE_GPU
void
MFIter::openFan () const
{
    m_fan_tried.flag = true;

    using Edge = Gpu::StreamEdge;
    const bool launch = Gpu::inLaunchRegion();
    const bool all_boxes = (flags & AllBoxes) != 0;
    const int nitems = static_cast<int>(index_map->size());
#ifdef AMREX_USE_OMP
    const bool in_team = omp_in_parallel() && omp_get_num_threads() > 1;
#else
    const bool in_team = false;
#endif

    Gpu::StreamSequencer::ForkInfo info{.nitems = nitems, .shared = in_team};
    info.ordered = Gpu::inStreamOrderedRegion();
    if (info.ordered) {
        // The loop is kept in program order. If it does not synchronize the
        // host, its streams wait on the GPU for the work queued before the
        // loop, and the work queued after it waits for them. Outside the
        // launch region the body queues no kernels, so the loop stays on the
        // current stream.
        info.max_width = launch ? streams : 1;
        info.fork = (device_sync && device_sync_pre) ? Edge::Host : Edge::Gpu;
        info.join = (device_sync && device_sync_post) ? Edge::Host : Edge::Gpu;
    } else {
        // Outside a Gpu::StreamOrderedRegion, a loop that synchronizes waits
        // on the host before (if it uses more than one stream) and after. A
        // loop that does not synchronize is not ordered.
        info.max_width = streams;
        info.fork = (launch && device_sync && device_sync_pre && !all_boxes &&
                     streams > 1 && nitems > 1) ? Edge::Host : Edge::None;
        info.join = (device_sync && device_sync_post) ? Edge::Host : Edge::None;
    }
    Gpu::theStreamSequencer().fork(m_fan, info);
    // Outside a Gpu::StreamOrderedRegion, an AllBoxes loop starts on the
    // current stream.
    if (info.ordered || !all_boxes) {
        Gpu::theStreamSequencer().lane(m_fan, currentIndex);
    }
}
#endif

Box
MFIter::tilebox () const noexcept
{
    BL_ASSERT(tile_array != nullptr);
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
    BL_ASSERT(tile_array != nullptr);
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
    if (dir < 0) { return tilebox(IntVect::TheNodeVector()); }
    return tilebox(IntVect::TheDimensionVector(dir));
}

// Note that a small negative ng is supported.
Box
MFIter::growntilebox (int a_ng) const noexcept
{
    Box bx = tilebox();
    IntVect ngv{a_ng};
    if (a_ng < -100) { ngv = fabArray->nGrowVect(); }
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
    if (a_ng < -100) { ngv = fabArray->nGrowVect(); }
    return grownnodaltilebox(dir, ngv);
}

Box
MFIter::grownnodaltilebox (int dir, IntVect const& a_ng) const noexcept
{
    BL_ASSERT(dir < AMREX_SPACEDIM);
    if (dir < 0) { return tilebox(IntVect::TheNodeVector(), a_ng); }
    return tilebox(IntVect::TheDimensionVector(dir), a_ng);
}

void
MFIter::operator++ () noexcept
{
#ifdef AMREX_USE_OMP
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
        // A loop that advances without having called isValid() (e.g., a
        // derived iterator with its own validity check) starts here.
        if (!m_fan_tried) { openFan(); }
        // (With more than one lane outside the launch region, which only
        // happens outside a Gpu::StreamOrderedRegion, the loop stays on its
        // first stream.)
        if (m_fan.isOpen() && currentIndex < endIndex &&
            (Gpu::inLaunchRegion() || m_fan.width() == 1)) {
            Gpu::theStreamSequencer().lane(m_fan, currentIndex);
        }
#endif
    }
}

}
