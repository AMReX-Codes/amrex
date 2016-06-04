#include <winstd.H>

#include <iterator>
#include <numeric>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#include <Utility.H>
#include <FabArray.H>
#include <ParmParse.H>

#ifdef BL_MEM_PROFILING
#include <MemProfiler.H>
#endif

//
// Set default values in Initialize()!!!
//
bool    FabArrayBase::do_async_sends;
int     FabArrayBase::MaxComp;
#if BL_SPACEDIM == 1
IntVect FabArrayBase::mfiter_tile_size(1024000);
#elif BL_SPACEDIM == 2
IntVect FabArrayBase::mfiter_tile_size(1024000,1024000);
#else
IntVect FabArrayBase::mfiter_tile_size(1024000,8,8);
#endif
IntVect FabArrayBase::comm_tile_size(D_DECL(1024000, 8, 8));
IntVect FabArrayBase::mfghostiter_tile_size(D_DECL(1024000, 8, 8));

int FabArrayBase::nFabArrays(0);

FabArrayBase::TACache              FabArrayBase::m_TheTileArrayCache;
FabArrayBase::FBCache              FabArrayBase::m_TheFBCache;
FabArrayBase::CPCCache             FabArrayBase::m_TheCopyCache;
FabArrayBase::FPCCache             FabArrayBase::m_TheFillPatchCache;

FabArrayBase::CacheStats           FabArrayBase::m_TAC_stats("TileArrayCache");
FabArrayBase::CacheStats           FabArrayBase::m_FBC_stats("SICache");
FabArrayBase::CacheStats           FabArrayBase::m_CPC_stats("CopyCache");
FabArrayBase::CacheStats           FabArrayBase::m_FPC_stats("FillPatchCache");

std::map<FabArrayBase::BDKey, int> FabArrayBase::m_BD_count;

FabArrayBase::FabArrayStats        FabArrayBase::m_FA_stats;

namespace
{
    bool initialized = false;
    //
    // Set default values in Initialize()!!!
    //
    int fb_cache_max_size;
    int copy_cache_max_size;
}

void
FabArrayBase::Initialize ()
{
    if (initialized) return;
    //
    // Set default values here!!!
    //
    FabArrayBase::do_async_sends    = true;
    FabArrayBase::MaxComp           = 25;

    copy_cache_max_size = 25;
    fb_cache_max_size   = 25;

    ParmParse pp("fabarray");

    Array<int> tilesize(BL_SPACEDIM);

    if (pp.queryarr("mfiter_tile_size", tilesize, 0, BL_SPACEDIM))
    {
	for (int i=0; i<BL_SPACEDIM; i++) FabArrayBase::mfiter_tile_size[i] = tilesize[i];
    }

    if (pp.queryarr("mfghostiter_tile_size", tilesize, 0, BL_SPACEDIM))
    {
	for (int i=0; i<BL_SPACEDIM; i++) FabArrayBase::mfghostiter_tile_size[i] = tilesize[i];
    }

    if (pp.queryarr("comm_tile_size", tilesize, 0, BL_SPACEDIM))
    {
        for (int i=0; i<BL_SPACEDIM; i++) FabArrayBase::comm_tile_size[i] = tilesize[i];
    }

    pp.query("maxcomp",             FabArrayBase::MaxComp);
    pp.query("do_async_sends",      FabArrayBase::do_async_sends);
    pp.query("fb_cache_max_size",   fb_cache_max_size);
    pp.query("copy_cache_max_size", copy_cache_max_size);
    //
    // Don't let the caches get too small. This simplifies some logic later.
    //
    if (fb_cache_max_size < 1)
        fb_cache_max_size = 1;
    if (copy_cache_max_size < 1)
        copy_cache_max_size = 1;
    if (MaxComp < 1)
        MaxComp = 1;

    FabArrayBase::nFabArrays = 0;

    BoxLib::ExecOnFinalize(FabArrayBase::Finalize);

#ifdef BL_MEM_PROFILING
    MemProfiler::add(m_TAC_stats.name, [] () -> MemProfiler::MemInfo {
	    return {m_TAC_stats.bytes, m_TAC_stats.bytes_hwm};
	});
    MemProfiler::add(m_FBC_stats.name, [] () -> MemProfiler::MemInfo {
	    return {m_FBC_stats.bytes, m_FBC_stats.bytes_hwm};
	});
    MemProfiler::add(m_CPC_stats.name, [] () -> MemProfiler::MemInfo {
	    return {m_CPC_stats.bytes, m_CPC_stats.bytes_hwm};
	});
    MemProfiler::add(m_FPC_stats.name, [] () -> MemProfiler::MemInfo {
	    return {m_FPC_stats.bytes, m_FPC_stats.bytes_hwm};
	});
#endif

    initialized = true;
}

FabArrayBase::FabArrayBase ()
{
    Initialize();
    faID = nFabArrays++;
}

FabArrayBase::~FabArrayBase () {}

Box
FabArrayBase::fabbox (int K) const
{
    return BoxLib::grow(boxarray[K], n_grow);
}

long
FabArrayBase::bytesOfMapOfCopyComTagContainers (const FabArrayBase::MapOfCopyComTagContainers& m)
{
    long r = sizeof(MapOfCopyComTagContainers);
    for (MapOfCopyComTagContainers::const_iterator it = m.begin(); it != m.end(); ++it) {
	r += sizeof(it->first) + BoxLib::bytesOf(it->second)
	    + BoxLib::gcc_map_node_extra_bytes;
    }
    return r;
}

long
FabArrayBase::CPC::bytes () const
{
    long cnt = sizeof(FabArrayBase::CPC);

    if (m_LocTags)
	cnt += BoxLib::bytesOf(*m_LocTags);

    if (m_SndTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_SndTags);

    if (m_RcvTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_RcvTags);

    if (m_SndVols)
	cnt += BoxLib::bytesOf(*m_SndVols);

    if (m_RcvVols)
	cnt += BoxLib::bytesOf(*m_RcvVols);

    return cnt;
}

long
FabArrayBase::bytesOfCPCCache ()
{
    long r;
    if (m_TheCopyCache.empty()) {
	r = 0L;
    } else {
	r = sizeof(m_TheCopyCache);
	for (CPCCacheIter it = m_TheCopyCache.begin(), End = m_TheCopyCache.end();
	     it != End; ++it)
	{
	    r += sizeof(it->first) + it->second.bytes() + BoxLib::gcc_map_node_extra_bytes;
	}
    }
    return r;
}

long
FabArrayBase::SI::bytes () const
{
    int cnt = sizeof(FabArrayBase::SI);

    if (m_LocTags)
	cnt += BoxLib::bytesOf(*m_LocTags);

    if (m_SndTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_SndTags);

    if (m_RcvTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_RcvTags);

    if (m_SndVols)
	cnt += BoxLib::bytesOf(*m_SndVols);

    if (m_RcvVols)
	cnt += BoxLib::bytesOf(*m_RcvVols);

    return cnt;
}

long
FabArrayBase::bytesOfFBCache ()
{
    long r;
    if (m_TheFBCache.empty()) {
	r = 0L;
    } else {
	r = sizeof(m_TheFBCache);
	for (FBCacheIter it = m_TheFBCache.begin(), End = m_TheFBCache.end();
	     it != End; ++it)
	{
	    r += sizeof(it->first) + it->second.bytes() + BoxLib::gcc_map_node_extra_bytes;
	}
    }
    return r;
}

long
FabArrayBase::TileArray::bytes () const
{
    return sizeof(*this) 
	+ (BoxLib::bytesOf(this->indexMap)      - sizeof(this->indexMap))
	+ (BoxLib::bytesOf(this->localIndexMap) - sizeof(this->localIndexMap))
	+ (BoxLib::bytesOf(this->tileArray)     - sizeof(this->tileArray));
}

long
FabArrayBase::bytesOfTACache ()
{
    long r;
    if (m_TheTileArrayCache.empty()) {
	r = 0L;
    } else {
	r = sizeof(m_TheTileArrayCache);
	for (TACache::const_iterator it = m_TheTileArrayCache.begin();
	     it != m_TheTileArrayCache.end(); ++it)
	{
	    r += sizeof(it->first) + sizeof(it->second)
		+ BoxLib::gcc_map_node_extra_bytes;
	    for (TAMap::const_iterator it2 = it->second.begin();
		 it2 != it->second.end(); ++it2)
	    {
		r += sizeof(it2->first) + (it2->second).bytes()
		    + BoxLib::gcc_map_node_extra_bytes;
	    }
	}
    }
    return r;
}

//
// Stuff used for copy() caching.
//

FabArrayBase::CPC::CPC (const BoxArray&            dstba,
                        const BoxArray&            srcba,
                        const DistributionMapping& dstdm,
                        const DistributionMapping& srcdm,
			int                        dstng,
			int                        srcng)
    :
    m_dstba(dstba),
    m_srcba(srcba),
    m_dstdm(dstdm),
    m_srcdm(srcdm),
    m_dstng(dstng),
    m_srcng(srcng),
    m_nuse(0),
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0) {}

FabArrayBase::CPC::~CPC ()
{
    delete m_LocTags;
    delete m_SndTags;
    delete m_RcvTags;
    delete m_SndVols;
    delete m_RcvVols;
}

bool
FabArrayBase::CPC::operator== (const CPC& rhs) const
{
    return
        m_dstba == rhs.m_dstba && m_srcba == rhs.m_srcba && m_dstdm == rhs.m_dstdm && m_srcdm == rhs.m_srcdm && m_dstng == rhs.m_dstng && m_srcng == rhs.m_srcng;
}

FabArrayBase::CPCCacheIter
FabArrayBase::TheCPC (const CPC&          cpc,
                      const FabArrayBase& dst,
                      const FabArrayBase& src)
{
    BL_PROFILE("FabArrayBase::TheCPC()");

    BL_ASSERT(cpc.m_dstba.size() > 0 && cpc.m_srcba.size() > 0);
    //
    // We want to choose our keys wisely to minimize search time.
    // We'd like to distinguish between copies of the same length
    // but with different edgeness of boxes.  We also want to
    // differentiate dst.copy(src) from src.copy(dst).
    //
    CPCCache&      TheCopyCache = FabArrayBase::m_TheCopyCache;
    const IntVect& Typ          = cpc.m_dstba[0].type();
    const int      Scale        = D_TERM(Typ[0],+3*Typ[1],+5*Typ[2]) + 11;

    int Key = cpc.m_dstba.size() + cpc.m_srcba.size() + Scale;
    Key    += cpc.m_dstba[0].numPts() + cpc.m_dstba[cpc.m_dstba.size()-1].numPts();
    Key    += cpc.m_dstdm[0] + cpc.m_dstdm[cpc.m_dstdm.size()-1];

    std::pair<CPCCacheIter,CPCCacheIter> er_it = TheCopyCache.equal_range(Key);

    for (CPCCacheIter it = er_it.first; it != er_it.second; ++it)
    {
        if (it->second == cpc)
        {
	    ++it->second.m_nuse;
	    m_CPC_stats.recordUse();
            return it;
        }
    }

    if (TheCopyCache.size() >= copy_cache_max_size)
    {
        //
        // Don't let the size of the cache get too big.
        // Get rid of entries with the biggest largest key that haven't been reused.
        // Otherwise just remove the entry with the largest key.
        //
        CPCCache::iterator End      = TheCopyCache.end();
        CPCCache::iterator last_it  = End;
        CPCCache::iterator erase_it = End;

        for (CPCCache::iterator it = TheCopyCache.begin(); it != End; ++it)
        {
            last_it = it;

            if (it->second.m_nuse <= 1)
                erase_it = it;
        }

        if (erase_it != End)
        {
#ifdef BL_MEM_PROFILING
	    m_CPC_stats.bytes -= erase_it->second.bytes();
#endif
	    m_CPC_stats.recordErase(erase_it->second.m_nuse);
            TheCopyCache.erase(erase_it);
        }
        else if (last_it != End)
        {
#ifdef BL_MEM_PROFILING
	    m_CPC_stats.bytes -= last_it->second.bytes();
#endif
	    m_CPC_stats.recordErase(last_it->second.m_nuse);
            TheCopyCache.erase(last_it);
        }
    }
    //
    // Got to insert one & then build it.
    //
    CPCCacheIter cache_it = TheCopyCache.insert(CPCCache::value_type(Key,cpc));
    CPC&         TheCPC   = cache_it->second;
    const int    MyProc   = ParallelDescriptor::MyProc();
    //
    // Here's where we allocate memory for the cache innards.
    // We do this so we don't have to build objects of these types
    // each time we search the cache.  Otherwise we'd be constructing
    // and destroying said objects quite frequently.
    //
    TheCPC.m_LocTags = new CopyComTag::CopyComTagsContainer;
    TheCPC.m_SndTags = new CopyComTag::MapOfCopyComTagContainers;
    TheCPC.m_RcvTags = new CopyComTag::MapOfCopyComTagContainers;
    TheCPC.m_SndVols = new std::map<int,int>;
    TheCPC.m_RcvVols = new std::map<int,int>;

    TheCPC.m_nuse = 1;

    m_CPC_stats.recordBuild();
    m_CPC_stats.recordUse();

    if (dst.IndexMap().empty() && src.IndexMap().empty()) {
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
#ifdef BL_MEM_PROFILING
	m_CPC_stats.bytes += TheCPC.bytes();
	m_CPC_stats.bytes_hwm = std::max(m_CPC_stats.bytes_hwm, m_CPC_stats.bytes);
#endif
        return cache_it;
    }

    const BoxArray& ba_src = TheCPC.m_srcba;
    const DistributionMapping& dm_src = TheCPC.m_srcdm;
    const Array<int>& imap_src = src.IndexMap();
    const int nlocal_src = imap_src.size();
    const int ng_src = TheCPC.m_srcng;

    const BoxArray& ba_dst = TheCPC.m_dstba;
    const DistributionMapping& dm_dst = TheCPC.m_dstdm;
    const Array<int>& imap_dst = dst.IndexMap();
    const int nlocal_dst = dst.IndexMap().size();
    const int ng_dst = TheCPC.m_dstng;

    std::vector< std::pair<int,Box> > isects;

    CopyComTag::MapOfCopyComTagContainers send_tags; // temp copy

    for (int i = 0; i < nlocal_src; ++i)
    {
	const int   k_src = imap_src[i];
	const Box& bx_src = BoxLib::grow(ba_src[k_src], ng_src);

	ba_dst.intersections(bx_src, isects, ng_dst);

	for (int j = 0, M = isects.size(); j < M; ++j)
        {
	    const int k_dst     = isects[j].first;
	    const Box& bx       = isects[j].second;
	    const int dst_owner = dm_dst[k_dst];

	    if (ParallelDescriptor::sameTeam(dst_owner)) {
		continue; // local copy will be dealt with later
	    } else if (MyProc == dm_src[k_src]) {
		send_tags[dst_owner].push_back(CopyComTag(bx, k_dst, k_src));
	    }
	}
    }

    CopyComTag::MapOfCopyComTagContainers recv_tags; // temp copy

    BaseFab<int> localtouch, remotetouch;
    bool check_local = false, check_remote = false;
#ifdef _OPENMP
    if (omp_get_max_threads() > 1) {
        check_local = true;
        check_remote = true;
    }
#endif    

    if (ParallelDescriptor::TeamSize() > 1) {
	check_local = true;
    }

    for (int i = 0; i < nlocal_dst; ++i)
    {
	const int   k_dst = imap_dst[i];
	const Box& bx_dst = BoxLib::grow(ba_dst[k_dst], ng_dst);

	if (check_local) {
	    localtouch.resize(bx_dst);
	    localtouch.setVal(0);
	}

	if (check_remote) {
	    remotetouch.resize(bx_dst);
	    remotetouch.setVal(0);
	}

	ba_src.intersections(bx_dst, isects, ng_src);

	for (int j = 0, M = isects.size(); j < M; ++j)
        {
	    const int k_src     = isects[j].first;
	    const Box& bx       = isects[j].second;
	    const int src_owner = dm_src[k_src];

	    if (ParallelDescriptor::sameTeam(src_owner)) { // local copy
		const BoxList tilelist(bx, FabArrayBase::comm_tile_size);
		for (BoxList::const_iterator
			 it_tile  = tilelist.begin(),
			 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
		{
		    TheCPC.m_LocTags->push_back(CopyComTag(*it_tile, k_dst, k_src));
		}
		if (check_local) {
		    localtouch.plus(1, bx);
		}
	    } else if (MyProc == dm_dst[k_dst]) {
		recv_tags[src_owner].push_back(CopyComTag(bx, k_dst, k_src));
		if (check_remote) {
		    remotetouch.plus(1, bx);
		}
	    }
	}

	if (check_local) {  
	    // safe if a cell is touched no more than once 
	    // keep checking thread safety if it is safe so far
            check_local = TheCPC.m_threadsafe_loc = localtouch.max() <= 1;
        }

	if (check_remote) {
            check_remote = TheCPC.m_threadsafe_rcv = remotetouch.max() <= 1;
        }
    }

//    ba_src.clear_hash_bin();
//    ba_dst.clear_hash_bin();

    for (int ipass = 0; ipass < 2; ++ipass) // pass 0: send; pass 1: recv
    {
	CopyComTag::MapOfCopyComTagContainers & Tags
	    = (ipass == 0) ? *TheCPC.m_SndTags : *TheCPC.m_RcvTags;
	CopyComTag::MapOfCopyComTagContainers & tmpTags
	    = (ipass == 0) ?         send_tags :         recv_tags;
	std::map<int,int> & Vols
	    = (ipass == 0) ? *TheCPC.m_SndVols : *TheCPC.m_RcvVols;

	for (CopyComTag::MapOfCopyComTagContainers::iterator 
		 it  = tmpTags.begin(), 
		 End = tmpTags.end();   it != End; ++it)
	{
	    const int key = it->first;
	    std::vector<CopyComTag>& cctv = it->second;

	    // We need to fix the order so that the send and recv processes match.
	    std::sort(cctv.begin(), cctv.end());

	    std::vector<CopyComTag> new_cctv;
	    new_cctv.reserve(cctv.size());

	    for (std::vector<CopyComTag>::const_iterator 
		     it2  = cctv.begin(),
		     End2 = cctv.end();   it2 != End2; ++it2)
	    {
		const Box& bx = it2->box;

		Vols[key] += bx.numPts();

		const BoxList tilelist(bx, FabArrayBase::comm_tile_size);
		for (BoxList::const_iterator 
			 it_tile  = tilelist.begin(), 
			 End_tile = tilelist.end();    it_tile != End_tile; ++it_tile)
                {
		    new_cctv.push_back(CopyComTag(*it_tile, it2->fabIndex, it2->srcIndex));
		}
	    }

	    Tags[key].swap(new_cctv);
	}
    }    

#ifdef BL_MEM_PROFILING
    m_CPC_stats.bytes += TheCPC.bytes();
    m_CPC_stats.bytes_hwm = std::max(m_CPC_stats.bytes_hwm, m_CPC_stats.bytes);
#endif

    return cache_it;
}

void
FabArrayBase::CPC::FlushCache ()
{
    long stats[3] = {0,0,0}; // size, reused, bytes

    stats[0] = m_TheCopyCache.size();

    for (CPCCacheIter it = m_TheCopyCache.begin(), End = m_TheCopyCache.end();
         it != End;
         ++it)
    {
        if (it->second.m_nuse >= 2)
            stats[1]++;
	m_CPC_stats.recordErase(it->second.m_nuse);
    }

    stats[2] = bytesOfCPCCache();

    if (BoxLib::verbose)
    {
#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceLongMax(&stats[0], 3, ParallelDescriptor::IOProcessorNumber());
        if (stats[0] > 0 && ParallelDescriptor::IOProcessor())
        {
            std::cout << "CPC::m_TheCopyCache: max size: "
                      << stats[0]
                      << ", max # reused: "
                      << stats[1]
                      << ", max bytes used: "
                      << stats[2]
                      << std::endl;
        }
#ifdef BL_LAZY
	});
#endif
    }

    m_TheCopyCache.clear();
#ifdef BL_MEM_PROFILING
    m_CPC_stats.bytes = 0L;
#endif
}

FabArrayBase::SI::SI ()
    :
    m_ngrow(-1),
    m_nuse(0),
    m_cross(false),
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0) {}

FabArrayBase::SI::SI (const BoxArray&            ba,
                      const DistributionMapping& dm,
                      int                        ngrow,
                      bool                       cross)
    :
    m_ba(ba),
    m_dm(dm),
    m_ngrow(ngrow),
    m_nuse(0),
    m_cross(cross),
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0)
{
    BL_ASSERT(ngrow >= 0);
}

FabArrayBase::SI::~SI ()
{
    delete m_LocTags;
    delete m_SndTags;
    delete m_RcvTags;
    delete m_SndVols;
    delete m_RcvVols;
}

bool
FabArrayBase::SI::operator== (const SI& rhs) const
{
    return
        m_ngrow == rhs.m_ngrow && m_cross == rhs.m_cross && m_ba == rhs.m_ba && m_dm == rhs.m_dm;
}

FabArrayBase::FBCacheIter
FabArrayBase::TheFB (bool                cross,
                     const FabArrayBase& mf)
{
    BL_PROFILE("FabArray::TheFB");

    BL_ASSERT(mf.size() > 0);

    const FabArrayBase::SI si(mf.boxArray(), mf.DistributionMap(), mf.nGrow(), cross);

    const IntVect& Typ   = mf.boxArray()[0].type();
    const int      Scale = D_TERM(Typ[0],+3*Typ[1],+5*Typ[2]) + 11;
    const int      Key   = mf.size() + mf.boxArray()[0].numPts() + mf.nGrow() + Scale + cross;

    std::pair<FBCacheIter,FBCacheIter> er_it = m_TheFBCache.equal_range(Key);

    for (FBCacheIter it = er_it.first; it != er_it.second; ++it)
    {
        if (it->second == si)
        {
	    ++it->second.m_nuse;
	    m_FBC_stats.recordUse();
            return it;
        }
    }

    if (m_TheFBCache.size() >= fb_cache_max_size)
    {
        //
        // Don't let the size of the cache get too big.
        // Get rid of entries with the biggest largest key that haven't been reused.
        // Otherwise just remove the entry with the largest key.
        //
        FBCacheIter End      = m_TheFBCache.end();
        FBCacheIter last_it  = End;
        FBCacheIter erase_it = End;

        for (FBCacheIter it = m_TheFBCache.begin(); it != End; ++it)
        {
            last_it = it;

            if (it->second.m_nuse <= 1)
                erase_it = it;
        }

        if (erase_it != End)
        {
#ifdef BL_MEM_PROFILING
	    m_FBC_stats.bytes -= erase_it->second.bytes();
#endif
	    m_FBC_stats.recordErase(erase_it->second.m_nuse);
            m_TheFBCache.erase(erase_it);
        }
        else if (last_it != End)
        {
#ifdef BL_MEM_PROFILING
	    m_FBC_stats.bytes -= last_it->second.bytes();
#endif
	    m_FBC_stats.recordErase(last_it->second.m_nuse);
	    m_TheFBCache.erase(last_it);
        }
    }
    //
    // Got to insert one & then build it.
    //
    FBCacheIter                cache_it = m_TheFBCache.insert(FBCache::value_type(Key,si));
    SI&                        TheFB    = cache_it->second;
    const int                  MyProc   = ParallelDescriptor::MyProc();
    const BoxArray&            ba       = mf.boxArray();
    const DistributionMapping& dm       = mf.DistributionMap();
    const Array<int>&          imap     = mf.IndexMap();
    //
    // Here's where we allocate memory for the cache innards.
    // We do this so we don't have to build objects of these types
    // each time we search the cache.  Otherwise we'd be constructing
    // and destroying said objects quite frequently.
    //
    TheFB.m_LocTags = new CopyComTag::CopyComTagsContainer;
    TheFB.m_SndTags = new CopyComTag::MapOfCopyComTagContainers;
    TheFB.m_RcvTags = new CopyComTag::MapOfCopyComTagContainers;
    TheFB.m_SndVols = new std::map<int,int>;
    TheFB.m_RcvVols = new std::map<int,int>;

    TheFB.m_nuse = 1;

    m_FBC_stats.recordBuild();
    m_FBC_stats.recordUse();

    if (imap.empty()) {
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
#ifdef BL_MEM_PROFILING
	m_FBC_stats.bytes += TheFB.bytes();
	m_FBC_stats.bytes_hwm = std::max(m_FBC_stats.bytes_hwm, m_FBC_stats.bytes);
#endif
        return cache_it;
    }

    // For local copy, all workers in the same team will have the identical copy of tags
    // so that they can share work.  But for remote communication, they are all different.

    const int nlocal = imap.size();
    const int ng = si.m_ngrow;
    std::vector< std::pair<int,Box> > isects;

    CopyComTag::MapOfCopyComTagContainers send_tags; // temp copy

    for (int i = 0; i < nlocal; ++i)
    {
	const int ksnd = imap[i];
	const Box& vbx = ba[ksnd];

	ba.intersections(vbx, isects, ng);

	for (int j = 0, M = isects.size(); j < M; ++j)
	{
	    const int krcv      = isects[j].first;
	    const Box& bx       = isects[j].second;
	    const int dst_owner = dm[krcv];

	    if (krcv == ksnd) continue;  // same box

	    if (ParallelDescriptor::sameTeam(dst_owner)) {
		continue;  // local copy will be dealt with later
	    } else if (MyProc == dm[ksnd]) {
		const BoxList& bl = BoxLib::boxDiff(bx, ba[krcv]);
		for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
		    send_tags[dst_owner].push_back(CopyComTag(*lit, krcv, ksnd));
	    }
	}
    }

    CopyComTag::MapOfCopyComTagContainers recv_tags; // temp copy

    BaseFab<int> localtouch, remotetouch;
    bool check_local = false, check_remote = false;
#ifdef _OPENMP
    if (omp_get_max_threads() > 1) {
        check_local = true;
        check_remote = true;
    }
#endif

    if (ParallelDescriptor::TeamSize() > 1) {
	check_local = true;
    }

    if (ba.ixType().cellCentered()) {
	TheFB.m_threadsafe_loc = true;
	TheFB.m_threadsafe_rcv = true;
        check_local = false;
        check_remote = false;
    }

    for (int i = 0; i < nlocal; ++i)
    {
	const int   krcv = imap[i];
	const Box& vbx   = ba[krcv];
	const Box& bxrcv = BoxLib::grow(vbx, ng);

	if (check_local) {
	    localtouch.resize(bxrcv);
	    localtouch.setVal(0);
	}

	if (check_remote) {
	    remotetouch.resize(bxrcv);
	    remotetouch.setVal(0);
	}

	ba.intersections(bxrcv, isects);

	for (int j = 0, M = isects.size(); j < M; ++j)
	{
	    const int ksnd      = isects[j].first;
	    const Box& bx       = isects[j].second;
	    const int src_owner = dm[ksnd];

	    if (krcv == ksnd) continue;  // same box

	    const BoxList& bl = BoxLib::boxDiff(bx, vbx);
	    for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
	    {
		const Box& blbx = *lit;

		if (ParallelDescriptor::sameTeam(src_owner)) { // local copy
		    const BoxList tilelist(blbx, FabArrayBase::comm_tile_size);
		    for (BoxList::const_iterator
			     it_tile  = tilelist.begin(),
			     End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
		    {
			TheFB.m_LocTags->push_back(CopyComTag(*it_tile, krcv, ksnd));
		    }
		    if (check_local) {
			localtouch.plus(1, blbx);
		    }
		} else if (MyProc == dm[krcv]) {
		    recv_tags[src_owner].push_back(CopyComTag(blbx, krcv, ksnd));
		    if (check_remote) {
			remotetouch.plus(1, blbx);
		    }
		}
	    }
	}

	if (check_local) {  
	    // safe if a cell is touched no more than once 
	    // keep checking thread safety if it is safe so far
            check_local = TheFB.m_threadsafe_loc = localtouch.max() <= 1;
        }

	if (check_remote) {
            check_remote = TheFB.m_threadsafe_rcv = remotetouch.max() <= 1;
        }
    }

//    ba.clear_hash_bin();

    for (int ipass = 0; ipass < 2; ++ipass) // pass 0: send; pass 1: recv
    {
	CopyComTag::MapOfCopyComTagContainers & Tags
	    = (ipass == 0) ? *TheFB.m_SndTags : *TheFB.m_RcvTags;
	CopyComTag::MapOfCopyComTagContainers & tmpTags
	    = (ipass == 0) ?        send_tags :        recv_tags;
	std::map<int,int> & Vols
	    = (ipass == 0) ? *TheFB.m_SndVols : *TheFB.m_RcvVols;

	for (CopyComTag::MapOfCopyComTagContainers::iterator 
		 it  = tmpTags.begin(), 
		 End = tmpTags.end();   it != End; ++it)
	{
	    const int key = it->first;
	    std::vector<CopyComTag>& cctv = it->second;

	    // We need to fix the order so that the send and recv processes match.
	    std::sort(cctv.begin(), cctv.end());

	    std::vector<CopyComTag> new_cctv;
	    new_cctv.reserve(cctv.size());

	    for (std::vector<CopyComTag>::const_iterator 
		     it2  = cctv.begin(),
		     End2 = cctv.end();   it2 != End2; ++it2)
	    {
		const Box& bx = it2->box;

		std::vector<Box> boxes;
		int vol = 0;

		if (si.m_cross) {
		    const Box& dstfabbx = ba[it2->fabIndex];
		    for (int dir = 0; dir < BL_SPACEDIM; dir++)
	            {
			Box lo = dstfabbx;
			lo.setSmall(dir, dstfabbx.smallEnd(dir) - ng);
			lo.setBig  (dir, dstfabbx.smallEnd(dir) - 1);
			lo &= bx;
			if (lo.ok()) {
			    boxes.push_back(lo);
			    vol += lo.numPts();
			}

			Box hi = dstfabbx;
			hi.setSmall(dir, dstfabbx.bigEnd(dir) + 1);
			hi.setBig  (dir, dstfabbx.bigEnd(dir) + ng);
			hi &= bx;
			if (hi.ok()) {
			    boxes.push_back(hi);
			    vol += hi.numPts();
			}
		    }
		} else {
		    boxes.push_back(bx);
		    vol += bx.numPts();
		}

		if (vol > 0) 
		{
		    Vols[key] += vol;

		    for (std::vector<Box>::const_iterator 
			     it_bx  = boxes.begin(),
			     End_bx = boxes.end();    it_bx != End_bx; ++it_bx)
		    {
			const BoxList tilelist(*it_bx, FabArrayBase::comm_tile_size);
			for (BoxList::const_iterator 
				 it_tile  = tilelist.begin(), 
				 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
			{
			    new_cctv.push_back(CopyComTag(*it_tile, it2->fabIndex, it2->srcIndex));
			}
		    }
		}
	    }

	    if (!new_cctv.empty())
		Tags[key].swap(new_cctv);
	}
    }

#ifdef BL_MEM_PROFILING
	m_FBC_stats.bytes += TheFB.bytes();
	m_FBC_stats.bytes_hwm = std::max(m_FBC_stats.bytes_hwm, m_FBC_stats.bytes);
#endif

    return cache_it;
}

FabArrayBase::FPC::FPC (const FabArrayBase& srcfa,
			const FabArrayBase& dstfa,
			Box                 dstdomain,
			int                 dstng,
			const BoxConverter& coarsener)
    : m_srcbdk   (srcfa.getBDKey()),
      m_dstbdk   (dstfa.getBDKey()),
      m_dstdomain(dstdomain),
      m_dstng    (dstng),
      m_coarsener(coarsener.clone()),
      m_nuse     (0)
{ 
    BL_PROFILE("FPC::FPC()");

    const BoxArray& srcba = srcfa.boxArray();
    const BoxArray& dstba = dstfa.boxArray();
    BL_ASSERT(srcba.ixType() == dstba.ixType());

    const IndexType& boxtype = dstba.ixType();
    BL_ASSERT(boxtype == dstdomain.ixType());
     
    BL_ASSERT(dstng <= dstfa.nGrow());

    const DistributionMapping& dstdm = dstfa.DistributionMap();
    
    const int myproc = ParallelDescriptor::MyProc();

    BoxList bl(boxtype);
    Array<int> iprocs;

    for (int i = 0, N = dstba.size(); i < N; ++i)
    {
	Box bx = dstba[i];
	bx.grow(m_dstng);
	bx &= m_dstdomain;

	BoxList leftover = srcba.complement(bx);

	bool ismybox = (dstdm[i] == myproc);
	for (BoxList::const_iterator bli = leftover.begin(); bli != leftover.end(); ++bli)
	{
	    bl.push_back(m_coarsener->doit(*bli));
	    if (ismybox) {
		dst_boxes.push_back(*bli);
		dst_idxs.push_back(i);
	    }
	    iprocs.push_back(dstdm[i]);
	}
    }

    if (!iprocs.empty()) {
	ba_crse_patch.define(bl);
	iprocs.push_back(myproc);
	dm_crse_patch.define(iprocs);
    }
}

FabArrayBase::FPC::~FPC ()
{
    delete m_coarsener;
}

long
FabArrayBase::FPC::bytes () const
{
    long cnt = sizeof(FabArrayBase::FPC);
    cnt += sizeof(Box) * (ba_crse_patch.capacity() + dst_boxes.capacity());
    cnt += sizeof(int) * (dm_crse_patch.capacity() + dst_idxs.capacity());
    return cnt;
}

const FabArrayBase::FPC&
FabArrayBase::TheFPC (const FabArrayBase& srcfa,
		      const FabArrayBase& dstfa,
		      Box                 dstdomain,
		      int                 dstng,
		      const BoxConverter& coarsener)
{
    BL_PROFILE("FabArrayBase::TheFPC()");

    const BDKey& srckey = srcfa.getBDKey();
    const BDKey& dstkey = dstfa.getBDKey();

    std::pair<FPCCacheIter,FPCCacheIter> er_it = m_TheFillPatchCache.equal_range(dstkey);

    for (FPCCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	if (it->second->m_srcbdk    == srckey    &&
	    it->second->m_dstdomain == dstdomain &&
	    it->second->m_dstng     == dstng     &&
	    it->second->m_dstdomain.ixType() == dstdomain.ixType() &&
	    it->second->m_coarsener->doit(it->second->m_dstdomain) == coarsener.doit(dstdomain))
	{
	    ++it->second->m_nuse;
	    m_FPC_stats.recordUse();
	    return *(it->second);
	}
    }

    // Have to build a new one
    FPC* new_fpc = new FPC(srcfa, dstfa, dstdomain, dstng, coarsener);

#ifdef BL_MEM_PROFILING
    m_FPC_stats.bytes += new_fpc->bytes();
    m_FPC_stats.bytes_hwm = std::max(m_FPC_stats.bytes_hwm, m_FPC_stats.bytes);
#endif
    
    new_fpc->m_nuse = 1;
    m_FPC_stats.recordBuild();
    m_FPC_stats.recordUse();

    m_TheFillPatchCache.insert(er_it.second, FPCCache::value_type(dstkey,new_fpc));
    if (srckey != dstkey)
	m_TheFillPatchCache.insert(          FPCCache::value_type(srckey,new_fpc));

    return *new_fpc;
}

void
FabArrayBase::flushFPC ()
{
    BL_ASSERT(getBDKey() == m_bdkey);

    std::vector<FPCCacheIter> others;

    std::pair<FPCCacheIter,FPCCacheIter> er_it = m_TheFillPatchCache.equal_range(m_bdkey);

    for (FPCCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	const BDKey& srckey = it->second->m_srcbdk;
	const BDKey& dstkey = it->second->m_dstbdk;

	BL_ASSERT((srckey==dstkey && srckey==m_bdkey) || 
		  (m_bdkey==srckey) || (m_bdkey==dstkey));

	if (srckey != dstkey) {
	    const BDKey& otherkey = (m_bdkey == srckey) ? dstkey : srckey;
	    std::pair<FPCCacheIter,FPCCacheIter> o_er_it = m_TheFillPatchCache.equal_range(otherkey);

	    for (FPCCacheIter oit = o_er_it.first; oit != o_er_it.second; ++oit)
	    {
		if (it->second == oit->second)
		    others.push_back(oit);
	    }
	} 

#ifdef BL_MEM_PROFILING
	m_FPC_stats.bytes -= it->second->bytes();
#endif
	m_FPC_stats.recordErase(it->second->m_nuse);
	delete it->second;
    }
    
    m_TheFillPatchCache.erase(er_it.first, er_it.second);

    for (std::vector<FPCCacheIter>::iterator it = others.begin(),
	     End = others.end(); it != End; ++it)
    {
	m_TheFillPatchCache.erase(*it);
    }
}

void
FabArrayBase::Finalize ()
{
    FabArrayBase::FlushSICache();
    FabArrayBase::CPC::FlushCache();

    FabArrayBase::flushTileArrayCache();

    if (ParallelDescriptor::IOProcessor() && BoxLib::verbose) {
	m_FA_stats.print();
	m_TAC_stats.print();
	m_FBC_stats.print();
	m_CPC_stats.print();
	m_FPC_stats.print();
    }

    initialized = false;
}

void
FabArrayBase::FlushSICache ()
{
    long stats[3] = {0,0,0}; // size, reused, bytes

    stats[0] = m_TheFBCache.size();

    for (FBCacheIter it = m_TheFBCache.begin(), End = m_TheFBCache.end();
         it != End;
         ++it)
    {
        if (it->second.m_nuse >= 2)
            stats[1]++;
	m_FBC_stats.recordErase(it->second.m_nuse);
    }

    stats[2] = bytesOfFBCache();

    if (BoxLib::verbose)
    {
#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceLongMax(&stats[0], 3, ParallelDescriptor::IOProcessorNumber());
        if (stats[0] > 0 && ParallelDescriptor::IOProcessor())
        {
            std::cout << "SI::TheFBCache: max size: "
                      << stats[0]
                      << ", max # reused: "
                      << stats[1]
                      << ", max bytes used: "
                      << stats[2]
                      << std::endl;
        }
#ifdef BL_LAZY
	});
#endif
    }

    m_TheFBCache.clear();
#ifdef BL_MEM_PROFILING
    m_FBC_stats.bytes = 0L;
#endif
}

int
FabArrayBase::SICacheSize ()
{
    return m_TheFBCache.size();
}

const FabArrayBase::TileArray* 
FabArrayBase::getTileArray (const IntVect& tilesize) const
{
    TileArray* p;

#ifdef _OPENMP
#pragma omp critical(gettilearray)
#endif
    {
	BL_ASSERT(getBDKey() == m_bdkey);
	p = &FabArrayBase::m_TheTileArrayCache[m_bdkey][tilesize];
	if (p->nuse == -1) {
	    buildTileArray(tilesize, *p);
	    p->nuse = 0;
	    m_TAC_stats.recordBuild();
#ifdef BL_MEM_PROFILING
	    m_TAC_stats.bytes += p->bytes();
	    m_TAC_stats.bytes_hwm = std::max(m_TAC_stats.bytes_hwm,
					     m_TAC_stats.bytes);
#endif
	}
#ifdef _OPENMP
#pragma omp master
#endif
	{
	    ++(p->nuse);
	    m_TAC_stats.recordUse();
        }
    }

    return p;
}

void
FabArrayBase::buildTileArray (const IntVect& tileSize, TileArray& ta) const
{
    // Note that we store Tiles always as cell-centered boxes, even if the boxarray is nodal.
    const int N = indexMap.size();

    if (tileSize == IntVect::TheZeroVector())
    {
	for (int i = 0; i < N; ++i)
	{
	    if (isOwner(i))
	    {
		const int K = indexMap[i]; 
		const Box& bx = boxarray.getCellCenteredBox(K);
		ta.indexMap.push_back(K);
		ta.localIndexMap.push_back(i);
		ta.tileArray.push_back(bx);
	    }
	}
    }
    else
    {
#if defined(BL_USE_TEAM) && !defined(__INTEL_COMPILER)
	std::vector<int> local_idxs(N);
	std::iota(std::begin(local_idxs), std::end(local_idxs), 0);
#else
	std::vector<int> local_idxs;
	for (int i = 0; i < N; ++i)
	    local_idxs.push_back(i);
#endif

#if defined(BL_USE_TEAM)
	const int nworkers = ParallelDescriptor::TeamSize();
	if (nworkers > 1) {
	    // reorder it so that each worker will be more likely to work on their own fabs
	    std::stable_sort(local_idxs.begin(), local_idxs.end(), [this](int i, int j) 
			     { return  this->distributionMap[this->indexMap[i]] 
				     < this->distributionMap[this->indexMap[j]]; });
	}
#endif	

	for (std::vector<int>::const_iterator it = local_idxs.begin(); it != local_idxs.end(); ++it)
	{
	    const int i = *it;         // local index 
	    const int K = indexMap[i]; // global index
	    const Box& bx = boxarray.getCellCenteredBox(K);
	    
	    IntVect nt_in_fab, tsize, nleft;
	    int ntiles = 1;
	    for (int d=0; d<BL_SPACEDIM; d++) {
		int ncells = bx.length(d);
		nt_in_fab[d] = std::max(ncells/tileSize[d], 1);
		tsize    [d] = ncells/nt_in_fab[d];
		nleft    [d] = ncells - nt_in_fab[d]*tsize[d];
		ntiles *= nt_in_fab[d];
	    }
	    
	    IntVect small, big, ijk;  // note that the initial values are all zero.
	    ijk[0] = -1;
	    for (int t = 0; t < ntiles; ++t) {
		ta.indexMap.push_back(K);
		ta.localIndexMap.push_back(i);
		
		for (int d=0; d<BL_SPACEDIM; d++) {
		    if (ijk[d]<nt_in_fab[d]-1) {
			ijk[d]++;
			break;
		    } else {
			ijk[d] = 0;
		    }
		}
		
		for (int d=0; d<BL_SPACEDIM; d++) {
		    if (ijk[d] < nleft[d]) {
			small[d] = ijk[d]*(tsize[d]+1);
			big[d] = small[d] + tsize[d];
		    } else {
			small[d] = ijk[d]*tsize[d] + nleft[d];
			big[d] = small[d] + tsize[d] - 1;
		    }
		}
		
		Box tbx(small, big, IndexType::TheCellType());
		tbx.shift(bx.smallEnd());
		
		ta.tileArray.push_back(tbx);
	    }
	}
    }
}

void
FabArrayBase::flushTileArray (const IntVect& tileSize) const
{
    BL_ASSERT(getBDKey() == m_bdkey);

    TACache& tao = m_TheTileArrayCache;
    TACache::iterator tao_it = tao.find(m_bdkey);
    if(tao_it != tao.end()) 
    {
	if (tileSize == IntVect::TheZeroVector()) 
	{
	    for (TAMap::const_iterator tai_it = tao_it->second.begin();
		 tai_it != tao_it->second.end(); ++tai_it)
	    {
#ifdef BL_MEM_PROFILING
		m_TAC_stats.bytes -= tai_it->second.bytes();
#endif		
		m_TAC_stats.recordErase(tai_it->second.nuse);
	    }
	    tao.erase(tao_it);
	} 
	else 
	{
	    TAMap& tai = tao_it->second;
	    TAMap::iterator tai_it = tai.find(tileSize);
	    if (tai_it != tai.end()) {
#ifdef BL_MEM_PROFILING
		m_TAC_stats.bytes -= tai_it->second.bytes();
#endif		
		m_TAC_stats.recordErase(tai_it->second.nuse);
		tai.erase(tai_it);
	    }
	}
    }
}

void
FabArrayBase::flushTileArrayCache ()
{
    for (TACache::const_iterator tao_it = m_TheTileArrayCache.begin();
	 tao_it != m_TheTileArrayCache.end(); ++tao_it)
    {
	for (TAMap::const_iterator tai_it = tao_it->second.begin();
	     tai_it != tao_it->second.end(); ++tai_it)
	{
	    m_TAC_stats.recordErase(tai_it->second.nuse);
	}
    }
    m_TheTileArrayCache.clear();
#ifdef BL_MEM_PROFILING
    m_TAC_stats.bytes = 0L;
#endif
}

void
FabArrayBase::clearThisBD (bool no_assertion)
{
    if (! boxarray.empty() ) 
    {
	BL_ASSERT(no_assertion || getBDKey() == m_bdkey);

	std::map<BDKey, int>::iterator cnt_it = m_BD_count.find(m_bdkey);
	if (cnt_it != m_BD_count.end()) 
	{
	    --(cnt_it->second);
	    if (cnt_it->second == 0) 
	    {
		m_BD_count.erase(cnt_it);
		
		// Since this is the last one built with these BoxArray 
		// and DistributionMapping, erase it from caches.
		flushTileArray();
		flushFPC();
	    }
	}
    }
}

void
FabArrayBase::addThisBD ()
{
    m_bdkey = getBDKey();
    int cnt = ++(m_BD_count[m_bdkey]);
    if (cnt == 1) { // new one
	m_FA_stats.recordMaxNumBoxArrays(m_BD_count.size());
    } else {
	m_FA_stats.recordMaxNumBAUse(cnt);
    }
}

void
FabArrayBase::updateBDKey ()
{
    if (getBDKey() != m_bdkey) {
	clearThisBD(true);
	addThisBD();
    }
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		unsigned char       flags_)
    :
    fabArray(fabarray_),
    tile_size((flags_ & Tiling) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(flags_),
    index_map(0),
    local_index_map(0),
    tile_array(0)
{
    Initialize();
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		bool                do_tiling_)
    :
    fabArray(fabarray_),
    tile_size((do_tiling_) ? FabArrayBase::mfiter_tile_size : IntVect::TheZeroVector()),
    flags(do_tiling_ ? Tiling : 0),
    index_map(0),
    local_index_map(0),
    tile_array(0)
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
    index_map(0),
    local_index_map(0),
    tile_array(0)
{
    Initialize();
}

MFIter::~MFIter ()
{
#if BL_USE_TEAM
    if ( ! (flags & NoTeamBarrier) )
	ParallelDescriptor::MyTeam().MemoryBarrier();
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
	index_map    = &(fabArray.IndexMap());
	currentIndex = 0;
	beginIndex   = 0;
	endIndex     = index_map->size();
    }
    else
    {
	const FabArrayBase::TileArray* pta = fabArray.getTileArray(tile_size);
	
	index_map       = &(pta->indexMap);
	local_index_map = &(pta->localIndexMap);
	tile_array      = &(pta->tileArray);

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
#endif

	currentIndex = beginIndex;

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
	const IntVect& Big = validbox().bigEnd();
	for (int d=0; d<BL_SPACEDIM; ++d) {
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
MFIter::nodaltilebox (int dir) const 
{ 
    BL_ASSERT(tile_array != 0);
    Box bx((*tile_array)[currentIndex]);
    bx.convert(typ);
    const IntVect& Big = validbox().bigEnd();
    int d0, d1;
    if (dir < 0) {
	d0 = 0;
	d1 = BL_SPACEDIM-1;
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

Box 
MFIter::growntilebox (int ng) const 
{
    Box bx = tilebox();
    if (ng < -100) ng = fabArray.nGrow();
    const Box& vbx = validbox();
    for (int d=0; d<BL_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, ng);
	}
	if (bx.bigEnd(d) == vbx.bigEnd(d)) {
	    bx.growHi(d, ng);
	}
    }
    return bx;
}

Box
MFIter::grownnodaltilebox (int dir, int ng) const
{
    Box bx = nodaltilebox(dir);
    if (ng < -100) ng = fabArray.nGrow();
    const Box& vbx = validbox();
    for (int d=0; d<BL_SPACEDIM; ++d) {
	if (bx.smallEnd(d) == vbx.smallEnd(d)) {
	    bx.growLo(d, ng);
	}
	if (bx.bigEnd(d) >= vbx.bigEnd(d)) {
	    bx.growHi(d, ng);
	}
    }
    return bx;
}

// For whatever reason, gcc sometimes doesn't like these two functions in FabArray.H,
// unless '-fno-inline' is used.  Moving it here seems to solve the problem.
Box
#ifdef __GNUG__ 
__attribute__ ((noinline))
#endif
MFIter::validbox () const 
{ 
    return fabArray.box((*index_map)[currentIndex]);
}
Box
MFIter::fabbox () const 
{ 
    return fabArray.fabbox((*index_map)[currentIndex]);
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
    Array<int> allindex;
    Array<int> alllocalindex;

    for (int i=0; i < fabArray.IndexMap().size(); ++i) {
	int K = fabArray.IndexMap()[i];
	const Box& vbx = fabArray.box(K);
	const Box& fbx = fabArray.fabbox(K);

	const BoxList& diff = BoxLib::boxDiff(fbx, vbx);
	
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
