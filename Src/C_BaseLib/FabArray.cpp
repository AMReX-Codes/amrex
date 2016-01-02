#include <winstd.H>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#include <FabArray.H>
#include <ParmParse.H>
//
// Set default values in Initialize()!!!
//
bool    FabArrayBase::Verbose;
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
IntVect FabArrayBase::mfiter_huge_box_size(D_DECL(1024000,1024000,1024000));

int FabArrayBase::nFabArrays(0);

FabArrayBase::TACache              FabArrayBase::m_TheTileArrayCache;
FabArrayBase::FBCache              FabArrayBase::m_TheFBCache;
FabArrayBase::CPCCache             FabArrayBase::m_TheCopyCache;

FabArrayBase::CacheStats           FabArrayBase::m_TAC_stats("Tile Array Cache");
FabArrayBase::CacheStats           FabArrayBase::m_FBC_stats("Fill Boundary Cache");
FabArrayBase::CacheStats           FabArrayBase::m_CPC_stats("Copy Cache");

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
    FabArrayBase::Verbose           = true;
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

    Array<int> commtilesize(BL_SPACEDIM);
    if (pp.queryarr("comm_tile_size", commtilesize, 0, BL_SPACEDIM))
    {
        for (int i=0; i<BL_SPACEDIM; i++) FabArrayBase::comm_tile_size[i] = commtilesize[i];
    }

    pp.query("verbose",             FabArrayBase::Verbose);
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

//
// Stuff used for copy() caching.
//

FabArrayBase::CPC::CPC ()
    :
    m_nuse(0),
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0) {}

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
        m_dstba == rhs.m_dstba && m_srcba == rhs.m_srcba && m_dstdm == rhs.m_dstdm && m_srcdm == rhs.m_srcdm;
}

int
FabArrayBase::CPC::bytes () const
{
    int cnt = sizeof(FabArrayBase::CPC);

    if (m_LocTags)
    {
        cnt += sizeof(CopyComTagsContainer) + m_LocTags->size()*sizeof(CopyComTag);
    }

    if (m_SndTags)
    {
        cnt += sizeof(MapOfCopyComTagContainers);

        cnt += m_SndTags->size()*sizeof(MapOfCopyComTagContainers::value_type);

        for (MapOfCopyComTagContainers::const_iterator it = m_SndTags->begin(),
                 m_End = m_SndTags->end();
             it != m_End;
             ++it)
        {
            cnt += it->second.size()*sizeof(CopyComTag);
        }
    }

    if (m_RcvTags)
    {
        cnt += sizeof(MapOfCopyComTagContainers);

        cnt += m_RcvTags->size()*sizeof(MapOfCopyComTagContainers::value_type);

        for (MapOfCopyComTagContainers::const_iterator it = m_RcvTags->begin(),
                 m_End = m_RcvTags->end();
             it != m_End;
             ++it)
        {
            cnt += it->second.size()*sizeof(CopyComTag);
        }
    }

    if (m_SndVols)
    {
        cnt += sizeof(std::map<int,int>) + m_SndVols->size()*sizeof(std::map<int,int>::value_type);
    }

    if (m_RcvVols)
    {
        cnt += sizeof(std::map<int,int>) + m_RcvVols->size()*sizeof(std::map<int,int>::value_type);
    }

    return cnt;
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
	    m_CPC_stats.recordErase(erase_it->second.m_nuse);
            TheCopyCache.erase(erase_it);
        }
        else if (last_it != End)
        {
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

    if (dst.IndexMap().empty() && src.IndexMap().empty())
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
        return cache_it;

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

	    if (dst_owner == MyProc) {
		const BoxList tilelist(bx, FabArrayBase::comm_tile_size);
		for (BoxList::const_iterator
			 it_tile  = tilelist.begin(),
			 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
		{
		    TheCPC.m_LocTags->push_back(CopyComTag(*it_tile, k_dst, k_src));
		}
	    } else {
		send_tags[dst_owner].push_back(CopyComTag(bx, k_dst, k_src));
	    }
	}
    }

    CopyComTag::MapOfCopyComTagContainers recv_tags; // temp copy

    for (int i = 0; i < nlocal_dst; ++i)
    {
	const int   k_dst = imap_dst[i];
	const Box& bx_dst = BoxLib::grow(ba_dst[k_dst], ng_dst);

	ba_src.intersections(bx_dst, isects, ng_src);

	for (int j = 0, M = isects.size(); j < M; ++j)
        {
	    const int k_src     = isects[j].first;
	    const Box& bx       = isects[j].second;
	    const int src_owner = dm_src[k_src];

	    if (src_owner == MyProc) continue; // local copy has been dealt with already

	    recv_tags[src_owner].push_back(CopyComTag(bx, k_dst, k_src));
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

    //
    // set thread safety
    //
#ifdef _OPENMP
    TheCPC.m_threadsafe_loc = LocThreadSafety(TheCPC.m_LocTags);
    TheCPC.m_threadsafe_rcv = RcvThreadSafety(TheCPC.m_RcvTags);
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
        stats[2] += it->second.bytes();
        if (it->second.m_nuse >= 2)
            stats[1]++;
	m_CPC_stats.recordErase(it->second.m_nuse);
    }

    if (FabArrayBase::Verbose)
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

int
FabArrayBase::SI::bytes () const
{
    int cnt = sizeof(FabArrayBase::SI);

    if (m_LocTags)
    {
        cnt += sizeof(CopyComTagsContainer) + m_LocTags->size()*sizeof(CopyComTag);
    }

    if (m_SndTags)
    {
        cnt += sizeof(MapOfCopyComTagContainers);

        cnt += m_SndTags->size()*sizeof(MapOfCopyComTagContainers::value_type);

        for (MapOfCopyComTagContainers::const_iterator it = m_SndTags->begin(),
                 m_End = m_SndTags->end();
             it != m_End;
             ++it)
        {
            cnt += it->second.size()*sizeof(CopyComTag);
        }
    }

    if (m_RcvTags)
    {
        cnt += sizeof(MapOfCopyComTagContainers);

        cnt += m_RcvTags->size()*sizeof(MapOfCopyComTagContainers::value_type);

        for (MapOfCopyComTagContainers::const_iterator it = m_RcvTags->begin(),
                 m_End = m_RcvTags->end();
             it != m_End;
             ++it)
        {
            cnt += it->second.size()*sizeof(CopyComTag);
        }
    }

    if (m_SndVols)
    {
        cnt += sizeof(std::map<int,int>) + m_SndVols->size()*sizeof(std::map<int,int>::value_type);
    }

    if (m_RcvVols)
    {
        cnt += sizeof(std::map<int,int>) + m_RcvVols->size()*sizeof(std::map<int,int>::value_type);
    }

    return cnt;
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
	    m_FBC_stats.recordErase(erase_it->second.m_nuse);
            m_TheFBCache.erase(erase_it);
        }
        else if (last_it != End)
        {
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

    if (imap.empty())
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
        return cache_it;

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

	    if (dst_owner == MyProc) {
		const BoxList tilelist(bx, FabArrayBase::comm_tile_size);
		for (BoxList::const_iterator
			 it_tile  = tilelist.begin(),
			 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
		{
		    TheFB.m_LocTags->push_back(CopyComTag(*it_tile, krcv, ksnd));
		}
	    } else {
		send_tags[dst_owner].push_back(CopyComTag(bx, krcv, ksnd));
	    }
	}
    }

    CopyComTag::MapOfCopyComTagContainers recv_tags; // temp copy

    for (int i = 0; i < nlocal; ++i)
    {
	const int krcv = imap[i];
	const Box& vbx = ba[krcv];

	ba.intersections(BoxLib::grow(vbx,ng), isects);

	for (int j = 0, M = isects.size(); j < M; ++j)
	{
	    const int ksnd      = isects[j].first;
	    const Box& bx       = isects[j].second;
	    const int src_owner = dm[ksnd];

	    if (krcv == ksnd) continue;  // same box

	    if (src_owner == MyProc) continue;  // local copy has been dealt with already

	    recv_tags[src_owner].push_back(CopyComTag(bx, krcv, ksnd));
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

	    Tags[key].swap(new_cctv);
	}
    }

    //
    // set thread safety
    //
#ifdef _OPENMP
    if (ba.ixType().cellCentered()) {
	TheFB.m_threadsafe_loc = true;
	TheFB.m_threadsafe_rcv = true;
    } else {
	TheFB.m_threadsafe_loc = false;
	TheFB.m_threadsafe_rcv = false;
    }
#endif

    return cache_it;
}

void
FabArrayBase::Finalize ()
{
    FabArrayBase::FlushSICache();
    FabArrayBase::CPC::FlushCache();

    FabArrayBase::flushTileArrayCache();

    if (ParallelDescriptor::IOProcessor()) {
	m_FA_stats.print();
	m_TAC_stats.print();
	m_FBC_stats.print();
	m_CPC_stats.print();
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
        stats[2] += it->second.bytes();
        if (it->second.m_nuse >= 2)
            stats[1]++;
	m_FBC_stats.recordErase(it->second.m_nuse);
    }

    if (FabArrayBase::Verbose)
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
}

int
FabArrayBase::SICacheSize ()
{
    return m_TheFBCache.size();
}

bool
FabArrayBase::LocThreadSafety(const CopyComTagsContainer* LocTags)
{
#ifdef _OPENMP
    bool tsall = true;
    int N_loc = (*LocTags).size();
    if (N_loc > 0) {
#pragma omp parallel reduction(&&:tsall)
	{
	    bool tsthis = true;
#pragma omp for schedule(static,1)
	    for (int i=0; i<N_loc-1; ++i) {
		if (tsthis) {
		    const CopyComTag& tagi = (*LocTags)[i];
		    for (int j=i+1; j<N_loc; ++j) {
			const CopyComTag& tagj = (*LocTags)[j];
			if ( tagi.fabIndex == tagj.fabIndex &&
			     tagi.box.intersects(tagj.box) ) {
			    tsthis = false;
			    break;
			}
		    }
		}
	    }
	    tsall = tsall && tsthis;
	}
    }
    return tsall;
#else
    return true;
#endif
}

bool 
FabArrayBase::RcvThreadSafety(const MapOfCopyComTagContainers* RcvTags)
{
#ifdef _OPENMP
    bool tsall = true;
    const int N_rcvs = RcvTags->size();
    if (N_rcvs > 0) {
	Array<const CopyComTagsContainer*> recv_cctc;
	recv_cctc.reserve(N_rcvs);
	
	for (MapOfCopyComTagContainers::const_iterator m_it = RcvTags->begin(),
		 m_End = RcvTags->end();
	     m_it != m_End;
	     ++m_it)
	{
	    recv_cctc.push_back(&(m_it->second));
	}
	
#pragma omp parallel reduction(&&:tsall)
	{
	    bool tsthis = true;
#pragma omp for schedule(static,1)
	    for (int i=0; i<N_rcvs-1; ++i) {
		if (tsthis) {
		    const CopyComTagsContainer& cctci = *recv_cctc[i];
		    for (CopyComTagsContainer::const_iterator iti = cctci.begin();
			 iti != cctci.end(); ++iti)
		    {
			for (int j=i+1; j<N_rcvs; ++j) {
			    const CopyComTagsContainer& cctcj = *recv_cctc[j];
			    for (CopyComTagsContainer::const_iterator itj = cctcj.begin();
				 itj != cctcj.end(); ++itj)
			    {
				if ( iti->fabIndex == itj->fabIndex &&
				     (iti->box).intersects(itj->box) ) {
				    tsthis = false;
				    goto labelRTS;
				}
			    }			    
			}
		    }
		}
	    labelRTS: ;
	    }
	    tsall = tsall && tsthis;
	}
    }
    return tsall;
#else
    return true;
#endif
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

    for (int i = 0; i < indexMap.size(); ++i)
    {
	const int K = indexMap[i]; 
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
		m_TAC_stats.recordErase(tai_it->second.nuse);
	    }
	    tao.erase(tao_it);
	} 
	else 
	{
	    TAMap& tai = tao_it->second;
	    TAMap::iterator tai_it = tai.find(tileSize);
	    if (tai_it != tai.end()) {
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
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		unsigned char       flags_)
    :
    fabArray(fabarray_),
    pta(fabarray_.getTileArray((flags_ & Tiling) 
			       ? FabArrayBase::mfiter_tile_size
			       : FabArrayBase::mfiter_huge_box_size)),
    flags(flags_)
{
    Initialize();
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		bool                do_tiling_)
    :
    fabArray(fabarray_),
    pta(fabarray_.getTileArray(do_tiling_ 
			       ? FabArrayBase::mfiter_tile_size
			       : FabArrayBase::mfiter_huge_box_size)),
    flags(do_tiling_ ? Tiling : 0)
{
    Initialize();
}

MFIter::MFIter (const FabArrayBase& fabarray_, 
		const IntVect&      tilesize_, 
		unsigned char       flags_)
    :
    fabArray(fabarray_),
    pta(fabarray_.getTileArray(tilesize_)),
    flags(flags_ | Tiling)
{
    Initialize();
}

MFIter::MFIter (const FabArrayBase&            fabarray_, 
		const FabArrayBase::TileArray* pta_,
		unsigned char                  flags_)
    :
    fabArray(fabarray_),
    pta(pta_),
    flags(flags_ | Tiling)
{
    Initialize();
}

MFIter::~MFIter ()
{
    ;
}

void 
MFIter::Initialize ()
{
    if (pta == 0) return;

    int rit = 0;
    int nworkers = 1;
    
#ifdef _OPENMP
    int nosharing = flags & NoSharing;
    if (omp_in_parallel() && !nosharing) {
	rit = omp_get_thread_num();
	nworkers = omp_get_num_threads();
    }
#endif

    int ntot = pta->indexMap.size();

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
    currentIndex = beginIndex;

    typ = fabArray.boxArray().ixType();
}

Box 
MFIter::tilebox () const
{ 
    Box bx(pta->tileArray[currentIndex]);
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
    Box bx(pta->tileArray[currentIndex]);
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

