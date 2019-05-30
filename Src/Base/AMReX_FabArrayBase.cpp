
#include <AMReX_FabArrayBase.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Utility.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>

#include <AMReX_BArena.H>
#include <AMReX_CArena.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#endif

namespace amrex {

//
// Set default values in Initialize()!!!
//
int     FabArrayBase::MaxComp;

#if defined(AMREX_USE_GPU) && defined(AMREX_USE_GPU_PRAGMA)

#if AMREX_SPACEDIM == 1
IntVect FabArrayBase::mfiter_tile_size(1024000);
#elif AMREX_SPACEDIM == 2
IntVect FabArrayBase::mfiter_tile_size(1024000,1024000);
#else
IntVect FabArrayBase::mfiter_tile_size(1024000,1024000,1024000);
#endif

#else

#if AMREX_SPACEDIM == 1
IntVect FabArrayBase::mfiter_tile_size(1024000);
#elif AMREX_SPACEDIM == 2
IntVect FabArrayBase::mfiter_tile_size(1024000,1024000);
#else
IntVect FabArrayBase::mfiter_tile_size(1024000,8,8);
#endif

#endif

#ifdef AMREX_USE_GPU
IntVect FabArrayBase::comm_tile_size(AMREX_D_DECL(1024000, 1024000, 1024000));
IntVect FabArrayBase::mfghostiter_tile_size(AMREX_D_DECL(1024000, 1024000, 1024000));
#else
IntVect FabArrayBase::comm_tile_size(AMREX_D_DECL(1024000, 8, 8));
IntVect FabArrayBase::mfghostiter_tile_size(AMREX_D_DECL(1024000, 8, 8));
#endif

FabArrayBase::TACache              FabArrayBase::m_TheTileArrayCache;
FabArrayBase::FBCache              FabArrayBase::m_TheFBCache;
FabArrayBase::CPCache              FabArrayBase::m_TheCPCache;
FabArrayBase::FPinfoCache          FabArrayBase::m_TheFillPatchCache;
FabArrayBase::CFinfoCache          FabArrayBase::m_TheCrseFineCache;

FabArrayBase::CacheStats           FabArrayBase::m_TAC_stats("TileArrayCache");
FabArrayBase::CacheStats           FabArrayBase::m_FBC_stats("FBCache");
FabArrayBase::CacheStats           FabArrayBase::m_CPC_stats("CopyCache");
FabArrayBase::CacheStats           FabArrayBase::m_FPinfo_stats("FillPatchCache");
FabArrayBase::CacheStats           FabArrayBase::m_CFinfo_stats("CrseFineCache");

std::map<FabArrayBase::BDKey, int> FabArrayBase::m_BD_count;

FabArrayBase::FabArrayStats        FabArrayBase::m_FA_stats;

namespace
{
    Arena* the_fa_arena = nullptr;
    bool initialized = false;
}

void
FabArrayBase::Initialize ()
{
    if (initialized) return;
    initialized = true;

    //
    // Set default values here!!!
    //
    FabArrayBase::MaxComp           = 25;

    ParmParse pp("fabarray");

    Vector<int> tilesize(AMREX_SPACEDIM);

    if (pp.queryarr("mfiter_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
	for (int i=0; i<AMREX_SPACEDIM; i++) FabArrayBase::mfiter_tile_size[i] = tilesize[i];
    }

    if (pp.queryarr("mfghostiter_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
	for (int i=0; i<AMREX_SPACEDIM; i++) FabArrayBase::mfghostiter_tile_size[i] = tilesize[i];
    }

    if (pp.queryarr("comm_tile_size", tilesize, 0, AMREX_SPACEDIM))
    {
        for (int i=0; i<AMREX_SPACEDIM; i++) FabArrayBase::comm_tile_size[i] = tilesize[i];
    }

    pp.query("maxcomp",             FabArrayBase::MaxComp);

    if (MaxComp < 1) {
        MaxComp = 1;
    }

    if (ParallelDescriptor::UseGpuAwareMpi()) {
        the_fa_arena = The_Device_Arena();
    } else {
        the_fa_arena = The_Pinned_Arena();
    }

    amrex::ExecOnFinalize(FabArrayBase::Finalize);

#ifdef BL_MEM_PROFILING
    MemProfiler::add(m_TAC_stats.name, std::function<MemProfiler::MemInfo()>
		     ([] () -> MemProfiler::MemInfo {
			 return {m_TAC_stats.bytes, m_TAC_stats.bytes_hwm};
		     }));
    MemProfiler::add(m_FBC_stats.name, std::function<MemProfiler::MemInfo()>
		     ([] () -> MemProfiler::MemInfo {
			 return {m_FBC_stats.bytes, m_FBC_stats.bytes_hwm};
		     }));
    MemProfiler::add(m_CPC_stats.name, std::function<MemProfiler::MemInfo()>
		     ([] () -> MemProfiler::MemInfo {
			 return {m_CPC_stats.bytes, m_CPC_stats.bytes_hwm};
		     }));
    MemProfiler::add(m_FPinfo_stats.name, std::function<MemProfiler::MemInfo()>
		     ([] () -> MemProfiler::MemInfo {
			 return {m_FPinfo_stats.bytes, m_FPinfo_stats.bytes_hwm};
		     }));
    MemProfiler::add(m_CFinfo_stats.name, std::function<MemProfiler::MemInfo()>
		     ([] () -> MemProfiler::MemInfo {
			 return {m_CFinfo_stats.bytes, m_CFinfo_stats.bytes_hwm};
		     }));
#endif
}

Arena*
The_FA_Arena ()
{
    return the_fa_arena;
}

FabArrayBase::FabArrayBase ()
{
}

FabArrayBase::~FabArrayBase () {}

void
FabArrayBase::define (const BoxArray&            bxs,
                      const DistributionMapping& dm,
                      int                        nvar,
                      int                        ngrow)
{
    define(bxs, dm, nvar, IntVect(ngrow));
}

void
FabArrayBase::define (const BoxArray&            bxs,
                      const DistributionMapping& dm,
                      int                        nvar,
                      const IntVect&             ngrow)
{
    BL_ASSERT(ngrow.allGE(IntVect::TheZeroVector()));
    BL_ASSERT(boxarray.size() == 0);
    indexArray.clear();
    ownership.clear();
    n_grow = ngrow;
    n_comp = nvar;
    
    boxarray = bxs;
    
    BL_ASSERT(dm.ProcessorMap().size() == bxs.size());
    distributionMap = dm;

    indexArray = distributionMap.getIndexArray();
    ownership = distributionMap.getOwnerShip();    
}

void
FabArrayBase::clear ()
{
    boxarray.clear();
    distributionMap = DistributionMapping();
    indexArray.clear();
    ownership.clear();
    m_bdkey = BDKey();
}

Box
FabArrayBase::fabbox (int K) const noexcept
{
    return amrex::grow(boxarray[K], n_grow);
}

long
FabArrayBase::bytesOfMapOfCopyComTagContainers (const FabArrayBase::MapOfCopyComTagContainers& m)
{
    long r = sizeof(MapOfCopyComTagContainers);
    for (MapOfCopyComTagContainers::const_iterator it = m.begin(); it != m.end(); ++it) {
	r += sizeof(it->first) + amrex::bytesOf(it->second)
	    + amrex::gcc_map_node_extra_bytes;
    }
    return r;
}

long
FabArrayBase::CPC::bytes () const
{
    long cnt = sizeof(FabArrayBase::CPC);

    if (m_LocTags)
	cnt += amrex::bytesOf(*m_LocTags);

    if (m_SndTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_SndTags);

    if (m_RcvTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_RcvTags);

    return cnt;
}

long
FabArrayBase::FB::bytes () const
{
    int cnt = sizeof(FabArrayBase::FB);

    if (m_LocTags)
	cnt += amrex::bytesOf(*m_LocTags);

    if (m_SndTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_SndTags);

    if (m_RcvTags)
	cnt += FabArrayBase::bytesOfMapOfCopyComTagContainers(*m_RcvTags);

    return cnt;
}

long
FabArrayBase::TileArray::bytes () const
{
    return sizeof(*this) 
	+ (amrex::bytesOf(this->numLocalTiles)     - sizeof(this->numLocalTiles))
	+ (amrex::bytesOf(this->indexMap)          - sizeof(this->indexMap))
	+ (amrex::bytesOf(this->localIndexMap)     - sizeof(this->localIndexMap))
	+ (amrex::bytesOf(this->localTileIndexMap) - sizeof(this->localTileIndexMap))
	+ (amrex::bytesOf(this->tileArray)         - sizeof(this->tileArray));
}

//
// Stuff used for copy() caching.
//

FabArrayBase::CPC::CPC (const FabArrayBase& dstfa, const IntVect& dstng,
			const FabArrayBase& srcfa, const IntVect& srcng,
			const Periodicity& period)
    : m_srcbdk(srcfa.getBDKey()), 
      m_dstbdk(dstfa.getBDKey()), 
      m_srcng(srcng), 
      m_dstng(dstng), 
      m_period(period),
      m_srcba(srcfa.boxArray()), 
      m_dstba(dstfa.boxArray()),
      m_threadsafe_loc(false), m_threadsafe_rcv(false),
      m_LocTags(0), m_SndTags(0), m_RcvTags(0), m_nuse(0)
{
    this->define(m_dstba, dstfa.DistributionMap(), dstfa.IndexArray(), 
		 m_srcba, srcfa.DistributionMap(), srcfa.IndexArray());
}

FabArrayBase::CPC::CPC (const BoxArray& dstba, const DistributionMapping& dstdm, 
			const Vector<int>& dstidx, const IntVect& dstng,
			const BoxArray& srcba, const DistributionMapping& srcdm, 
			const Vector<int>& srcidx, const IntVect& srcng,
			const Periodicity& period, int myproc)
    : m_srcbdk(), 
      m_dstbdk(), 
      m_srcng(srcng), 
      m_dstng(dstng), 
      m_period(period),
      m_srcba(srcba), 
      m_dstba(dstba),
      m_threadsafe_loc(false), m_threadsafe_rcv(false),
      m_LocTags(0), m_SndTags(0), m_RcvTags(0), m_nuse(0)
{
    this->define(dstba, dstdm, dstidx, srcba, srcdm, srcidx, myproc);
}

FabArrayBase::CPC::~CPC ()
{
    delete m_LocTags;
    delete m_SndTags;
    delete m_RcvTags;
}

void
FabArrayBase::CPC::define (const BoxArray& ba_dst, const DistributionMapping& dm_dst,
			   const Vector<int>& imap_dst,
			   const BoxArray& ba_src, const DistributionMapping& dm_src,
			   const Vector<int>& imap_src,
			   int MyProc)
{
    BL_PROFILE("FabArrayBase::CPC::define()");

    BL_ASSERT(ba_dst.size() > 0 && ba_src.size() > 0);
    BL_ASSERT(ba_dst.ixType() == ba_src.ixType());
    
    m_LocTags = new CopyComTag::CopyComTagsContainer;
    m_SndTags = new CopyComTag::MapOfCopyComTagContainers;
    m_RcvTags = new CopyComTag::MapOfCopyComTagContainers;

    if (!(imap_dst.empty() && imap_src.empty())) 
    {
	const int nlocal_src = imap_src.size();
	const IntVect& ng_src = m_srcng;
	const int nlocal_dst = imap_dst.size();
	const IntVect& ng_dst = m_dstng;

	std::vector< std::pair<int,Box> > isects;

	const std::vector<IntVect>& pshifts = m_period.shiftIntVect();

	auto& send_tags = *m_SndTags;
	
	for (int i = 0; i < nlocal_src; ++i)
	{
	    const int   k_src = imap_src[i];
	    const Box& bx_src = amrex::grow(ba_src[k_src], ng_src);

	    for (std::vector<IntVect>::const_iterator pit=pshifts.begin(); pit!=pshifts.end(); ++pit)
	    {
		ba_dst.intersections(bx_src+(*pit), isects, false, ng_dst);
	    
		for (int j = 0, M = isects.size(); j < M; ++j)
		{
		    const int k_dst     = isects[j].first;
		    const Box& bx       = isects[j].second;
		    const int dst_owner = dm_dst[k_dst];
		
		    if (ParallelDescriptor::sameTeam(dst_owner)) {
			continue; // local copy will be dealt with later
		    } else if (MyProc == dm_src[k_src]) {
			send_tags[dst_owner].push_back(CopyComTag(bx, bx-(*pit), k_dst, k_src));
		    }
		}
	    }
	}

	auto& recv_tags = *m_RcvTags;

	BaseFab<int,CpuDataAllocator<int> > localtouch, remotetouch;
	bool check_local = false, check_remote = false;
#if defined(_OPENMP)
	if (omp_get_max_threads() > 1) {
	    check_local = true;
	    check_remote = true;
	}
#elif defined(AMREX_USE_GPU)
        check_local = true;
        check_remote = true;
#endif    
	
	if (ParallelDescriptor::TeamSize() > 1) {
	    check_local = true;
	}
	
	for (int i = 0; i < nlocal_dst; ++i)
	{
	    const int   k_dst = imap_dst[i];
	    const Box& bx_dst = amrex::grow(ba_dst[k_dst], ng_dst);
	    
	    if (check_local) {
		localtouch.resize(bx_dst);
		localtouch.setVal(0);
	    }
	    
	    if (check_remote) {
		remotetouch.resize(bx_dst);
		remotetouch.setVal(0);
	    }
	    
	    for (std::vector<IntVect>::const_iterator pit=pshifts.begin(); pit!=pshifts.end(); ++pit)
	    {
		ba_src.intersections(bx_dst+(*pit), isects, false, ng_src);
	    
		for (int j = 0, M = isects.size(); j < M; ++j)
		{
		    const int k_src     = isects[j].first;
		    const Box& bx       = isects[j].second - *pit;
		    const int src_owner = dm_src[k_src];
		
		    if (ParallelDescriptor::sameTeam(src_owner, MyProc)) { // local copy
			const BoxList tilelist(bx, FabArrayBase::comm_tile_size);
			for (BoxList::const_iterator
				 it_tile  = tilelist.begin(),
				 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
			{
			    m_LocTags->push_back(CopyComTag(*it_tile, (*it_tile)+(*pit), k_dst, k_src));
			}
			if (check_local) {
			    localtouch.plus(1, bx);
			}
		    } else if (MyProc == dm_dst[k_dst]) {
			recv_tags[src_owner].push_back(CopyComTag(bx, bx+(*pit), k_dst, k_src));
			if (check_remote) {
			    remotetouch.plus(1, bx);
			}
		    }
		}
	    }
	    
	    if (check_local) {  
		// safe if a cell is touched no more than once 
		// keep checking thread safety if it is safe so far
		check_local = m_threadsafe_loc = localtouch.max() <= 1;
	    }
	    
	    if (check_remote) {
		check_remote = m_threadsafe_rcv = remotetouch.max() <= 1;
	    }
	}
	
	for (int ipass = 0; ipass < 2; ++ipass) // pass 0: send; pass 1: recv
	{
	    CopyComTag::MapOfCopyComTagContainers & Tags = (ipass == 0) ? *m_SndTags : *m_RcvTags;
            for (auto& kv : Tags)
	    {
		std::vector<CopyComTag>& cctv = kv.second;		
		// We need to fix the order so that the send and recv processes match.
		std::sort(cctv.begin(), cctv.end());
	    }
	}    
    }
}

FabArrayBase::CPC::CPC (const BoxArray& ba, const IntVect& ng,
                        const DistributionMapping& dstdm, const DistributionMapping& srcdm)
    : m_srcbdk(), 
      m_dstbdk(), 
      m_srcng(ng), 
      m_dstng(ng), 
      m_period(),
      m_srcba(ba), 
      m_dstba(ba),
      m_threadsafe_loc(true), m_threadsafe_rcv(true),
      m_LocTags(0), m_SndTags(0), m_RcvTags(0), m_nuse(0)
{
    BL_ASSERT(ba.size() > 0);

    m_LocTags = new CopyComTag::CopyComTagsContainer;
    m_SndTags = new CopyComTag::MapOfCopyComTagContainers;
    m_RcvTags = new CopyComTag::MapOfCopyComTagContainers;

    const int myproc = ParallelDescriptor::MyProc();

    for (int i = 0, N = ba.size(); i < N; ++i)
    {
        const int src_owner = srcdm[i];
        const int dst_owner = dstdm[i];
        if (src_owner == myproc || dst_owner == myproc)
        {
            const Box& bx = amrex::grow(ba[i], ng);
            const BoxList tilelist(bx, FabArrayBase::comm_tile_size);
            if (src_owner == myproc && dst_owner == myproc)
            {
                for (const Box& tbx : tilelist)
                {
                    m_LocTags->push_back(CopyComTag(tbx, tbx, i, i));
                }
            }
            else
            {
                auto& Tags = (src_owner == myproc) ? (*m_SndTags)[dst_owner] : (*m_RcvTags)[src_owner];
                Tags.push_back(CopyComTag(bx, bx, i, i));
            }
        }
    }
}

void
FabArrayBase::flushCPC (bool no_assertion) const
{
    BL_ASSERT(no_assertion || getBDKey() == m_bdkey);

    std::vector<CPCacheIter> others;

    std::pair<CPCacheIter,CPCacheIter> er_it = m_TheCPCache.equal_range(m_bdkey);

    for (CPCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	const BDKey& srckey = it->second->m_srcbdk;
	const BDKey& dstkey = it->second->m_dstbdk;

	BL_ASSERT((srckey==dstkey && srckey==m_bdkey) || 
		  (m_bdkey==srckey) || (m_bdkey==dstkey));

	if (srckey != dstkey) {
	    const BDKey& otherkey = (m_bdkey == srckey) ? dstkey : srckey;
	    std::pair<CPCacheIter,CPCacheIter> o_er_it = m_TheCPCache.equal_range(otherkey);

	    for (CPCacheIter oit = o_er_it.first; oit != o_er_it.second; ++oit)
	    {
		if (it->second == oit->second)
		    others.push_back(oit);
	    }
	}

#ifdef BL_MEM_PROFILING
	m_CPC_stats.bytes -= it->second->bytes();
#endif
	m_CPC_stats.recordErase(it->second->m_nuse);
	delete it->second;
    }

    m_TheCPCache.erase(er_it.first, er_it.second);

    for (std::vector<CPCacheIter>::iterator it = others.begin(),
	     End = others.end(); it != End; ++it)
    {
	m_TheCPCache.erase(*it);
    }    
}

void
FabArrayBase::flushCPCache ()
{
    for (CPCacheIter it = m_TheCPCache.begin(); it != m_TheCPCache.end(); ++it)
    {
	if (it->first == it->second->m_srcbdk) {
	    m_CPC_stats.recordErase(it->second->m_nuse);
	    delete it->second;
	}
    }
    m_TheCPCache.clear();
#ifdef BL_MEM_PROFILING
    m_CPC_stats.bytes = 0L;
#endif
}

const FabArrayBase::CPC&
FabArrayBase::getCPC (const IntVect& dstng, const FabArrayBase& src, const IntVect& srcng, const Periodicity& period) const
{
    BL_PROFILE("FabArrayBase::getCPC()");

    BL_ASSERT(getBDKey() == m_bdkey);
    BL_ASSERT(src.getBDKey() == src.m_bdkey);
    BL_ASSERT(boxArray().ixType() == src.boxArray().ixType());

    const BDKey& srckey = src.getBDKey();
    const BDKey& dstkey =     getBDKey();

    std::pair<CPCacheIter,CPCacheIter> er_it = m_TheCPCache.equal_range(dstkey);

    for (CPCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	if (it->second->m_srcng  == srcng &&
	    it->second->m_dstng  == dstng &&
	    it->second->m_srcbdk == srckey &&
	    it->second->m_dstbdk == dstkey &&
	    it->second->m_period == period &&
	    it->second->m_srcba  == src.boxArray() &&
	    it->second->m_dstba  == boxArray())
	{
	    ++(it->second->m_nuse);
	    m_CPC_stats.recordUse();
	    return *(it->second);
	}
    }
    
    // Have to build a new one
    CPC* new_cpc = new CPC(*this, dstng, src, srcng, period);

#ifdef BL_MEM_PROFILING
    m_CPC_stats.bytes += new_cpc->bytes();
    m_CPC_stats.bytes_hwm = std::max(m_CPC_stats.bytes_hwm, m_CPC_stats.bytes);
#endif    

    new_cpc->m_nuse = 1;
    m_CPC_stats.recordBuild();
    m_CPC_stats.recordUse();

    m_TheCPCache.insert(er_it.second, CPCache::value_type(dstkey,new_cpc));
    if (srckey != dstkey)
	m_TheCPCache.insert(          CPCache::value_type(srckey,new_cpc));

    return *new_cpc;
}

//
// Some stuff for fill boundary
//

FabArrayBase::FB::FB (const FabArrayBase& fa, const IntVect& nghost,
                      bool cross, const Periodicity& period, 
                      bool enforce_periodicity_only)
    : m_typ(fa.boxArray().ixType()), m_crse_ratio(fa.boxArray().crseRatio()),
      m_ngrow(nghost), m_cross(cross),
      m_epo(enforce_periodicity_only), m_period(period),
      m_threadsafe_loc(false), m_threadsafe_rcv(false),
      m_LocTags(new CopyComTag::CopyComTagsContainer),
      m_SndTags(new CopyComTag::MapOfCopyComTagContainers),
      m_RcvTags(new CopyComTag::MapOfCopyComTagContainers),
      m_nuse(0)
{
    BL_PROFILE("FabArrayBase::FB::FB()");

    if (!fa.IndexArray().empty()) {
	if (enforce_periodicity_only) {
	    BL_ASSERT(m_cross==false);
	    define_epo(fa);
	} else {
	    define_fb(fa);
	}
    }
}

void
FabArrayBase::FB::define_fb(const FabArrayBase& fa)
{
    const int                  MyProc   = ParallelDescriptor::MyProc();
    const BoxArray&            ba       = fa.boxArray();
    const DistributionMapping& dm       = fa.DistributionMap();
    const Vector<int>&         imap     = fa.IndexArray();

    // For local copy, all workers in the same team will have the identical copy of tags
    // so that they can share work.  But for remote communication, they are all different.
    
    const int nlocal = imap.size();
    const IntVect& ng = m_ngrow;
    std::vector< std::pair<int,Box> > isects;
    
    const std::vector<IntVect>& pshifts = m_period.shiftIntVect();
    
    auto& send_tags = *m_SndTags;
    
    for (int i = 0; i < nlocal; ++i)
    {
	const int ksnd = imap[i];
	const Box& vbx = ba[ksnd];
	
	for (auto pit=pshifts.cbegin(); pit!=pshifts.cend(); ++pit)
	{
	    ba.intersections(vbx+(*pit), isects, false, ng);

	    for (int j = 0, M = isects.size(); j < M; ++j)
	    {
		const int krcv      = isects[j].first;
		const Box& bx       = isects[j].second;
		const int dst_owner = dm[krcv];
		
		if (ParallelDescriptor::sameTeam(dst_owner)) {
		    continue;  // local copy will be dealt with later
		} else if (MyProc == dm[ksnd]) {
		    const BoxList& bl = amrex::boxDiff(bx, ba[krcv]);
		    for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
			send_tags[dst_owner].push_back(CopyComTag(*lit, (*lit)-(*pit), krcv, ksnd));
		}
	    }
	}
    }

    auto& recv_tags = *m_RcvTags;

    BaseFab<int,CpuDataAllocator<int> > localtouch, remotetouch;
    bool check_local = false, check_remote = false;
#if defined(_OPENMP)
    if (omp_get_max_threads() > 1) {
	check_local = true;
	check_remote = true;
    }
#elif defined(AMREX_USE_GPU)
    check_local = true;
    check_remote = true;
#endif

    if (ParallelDescriptor::TeamSize() > 1) {
	check_local = true;
    }

    for (int i = 0; i < nlocal; ++i)
    {
	const int   krcv = imap[i];
	const Box& vbx   = ba[krcv];
	const Box& bxrcv = amrex::grow(vbx, ng);
	
	if (check_local) {
	    localtouch.resize(bxrcv);
	    localtouch.setVal(0);
	}
	
	if (check_remote) {
	    remotetouch.resize(bxrcv);
	    remotetouch.setVal(0);
	}
	
	for (auto pit=pshifts.cbegin(); pit!=pshifts.cend(); ++pit)
	{
	    ba.intersections(bxrcv+(*pit), isects);

	    for (int j = 0, M = isects.size(); j < M; ++j)
	    {
		const int ksnd      = isects[j].first;
		const Box& dst_bx   = isects[j].second - *pit;
		const int src_owner = dm[ksnd];
		
		const BoxList& bl = amrex::boxDiff(dst_bx, vbx);
		for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
		{
		    const Box& blbx = *lit;
			
		    if (ParallelDescriptor::sameTeam(src_owner)) { // local copy
			const BoxList tilelist(blbx, FabArrayBase::comm_tile_size);
			for (BoxList::const_iterator
				 it_tile  = tilelist.begin(),
				 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
			{
			    m_LocTags->push_back(CopyComTag(*it_tile, (*it_tile)+(*pit), krcv, ksnd));
			}
			if (check_local) {
			    localtouch.plus(1, blbx);
			}
		    } else if (MyProc == dm[krcv]) {
			recv_tags[src_owner].push_back(CopyComTag(blbx, blbx+(*pit), krcv, ksnd));
			if (check_remote) {
			    remotetouch.plus(1, blbx);
			}
		    }
		}
	    }
	}

	if (check_local) {  
	    // safe if a cell is touched no more than once 
	    // keep checking thread safety if it is safe so far
	    check_local = m_threadsafe_loc = localtouch.max() <= 1;
	}

	if (check_remote) {
	    check_remote = m_threadsafe_rcv = remotetouch.max() <= 1;
	}
    }

    for (int ipass = 0; ipass < 2; ++ipass) // pass 0: send; pass 1: recv
    {
	CopyComTag::MapOfCopyComTagContainers & Tags = (ipass == 0) ? *m_SndTags : *m_RcvTags;

        Vector<int> to_be_deleted;
	    
        for (auto& kv : Tags)
	{
            std::vector<CopyComTag>& cctv = kv.second;
		
	    // We need to fix the order so that the send and recv processes match.
	    std::sort(cctv.begin(), cctv.end());
		
            std::vector<CopyComTag> cctv_tags_cross;
            cctv_tags_cross.reserve(cctv.size());

            for (auto const& tag : cctv)
            {
		const Box& bx = tag.dbox;
		const IntVect& d2s = tag.sbox.smallEnd() - tag.dbox.smallEnd();

		std::vector<Box> boxes;
		if (m_cross) {
		    const Box& dstvbx = ba[tag.dstIndex];
		    for (int dir = 0; dir < AMREX_SPACEDIM; dir++)
		    {
			Box lo = dstvbx;
			lo.setSmall(dir, dstvbx.smallEnd(dir) - ng[dir]);
			lo.setBig  (dir, dstvbx.smallEnd(dir) - 1);
			lo &= bx;
			if (lo.ok()) {
			    boxes.push_back(lo);
			}
			    
			Box hi = dstvbx;
			hi.setSmall(dir, dstvbx.bigEnd(dir) + 1);
			hi.setBig  (dir, dstvbx.bigEnd(dir) + ng[dir]);
			hi &= bx;
			if (hi.ok()) {
			    boxes.push_back(hi);
			}
		    }
		} else {
		    boxes.push_back(bx);
		}
		
		if (!boxes.empty()) 
		{
                    for (auto const& cross_box : boxes)
                    {
                        if (m_cross)
                        {
                            cctv_tags_cross.push_back(CopyComTag(cross_box, cross_box+d2s, 
                                                                 tag.dstIndex, tag.srcIndex));
                        }
		    }
		}
	    }
		
            if (!cctv_tags_cross.empty()) {
                cctv.swap(cctv_tags_cross);
            }
	}

        for (int key : to_be_deleted) {
            Tags.erase(key);
        }
    }
}

void
FabArrayBase::FB::define_epo (const FabArrayBase& fa)
{
    const int                  MyProc   = ParallelDescriptor::MyProc();
    const BoxArray&            ba       = fa.boxArray();
    const DistributionMapping& dm       = fa.DistributionMap();
    const Vector<int>&         imap     = fa.IndexArray();

    // For local copy, all workers in the same team will have the identical copy of tags
    // so that they can share work.  But for remote communication, they are all different.
    
    const int nlocal = imap.size();
    const IntVect& ng = m_ngrow;
    const IndexType& typ = ba.ixType();
    std::vector< std::pair<int,Box> > isects;
    
    const std::vector<IntVect>& pshifts = m_period.shiftIntVect();
    
    auto& send_tags = *m_SndTags;

    Box pdomain = m_period.Domain();
    pdomain.convert(typ);
    
    for (int i = 0; i < nlocal; ++i)
    {
	const int ksnd = imap[i];
	Box bxsnd = amrex::grow(ba[ksnd],ng);
	bxsnd &= pdomain; // source must be inside the periodic domain.

	if (!bxsnd.ok()) continue;

	for (auto pit=pshifts.cbegin(); pit!=pshifts.cend(); ++pit)
	{
	    if (*pit != IntVect::TheZeroVector())
	    {
		ba.intersections(bxsnd+(*pit), isects, false, ng);
		
		for (int j = 0, M = isects.size(); j < M; ++j)
		{
		    const int krcv      = isects[j].first;
		    const Box& bx       = isects[j].second;
		    const int dst_owner = dm[krcv];
		    
		    if (ParallelDescriptor::sameTeam(dst_owner)) {
			continue;  // local copy will be dealt with later
		    } else if (MyProc == dm[ksnd]) {
			const BoxList& bl = amrex::boxDiff(bx, pdomain);
			for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit) {
			    send_tags[dst_owner].push_back(CopyComTag(*lit, (*lit)-(*pit), krcv, ksnd));
			}
		    }
		}
	    }
	}
    }

    auto& recv_tags = *m_RcvTags;

    BaseFab<int,CpuDataAllocator<int> > localtouch, remotetouch;
    bool check_local = false, check_remote = false;
#if defined(_OPENMP)
    if (omp_get_max_threads() > 1) {
	check_local = true;
	check_remote = true;
    }
#elif defined(AMREX_USE_GPU)
    check_local = true;
    check_remote = true;
#endif

    if (ParallelDescriptor::TeamSize() > 1) {
	check_local = true;
    }

    for (int i = 0; i < nlocal; ++i)
    {
	const int   krcv = imap[i];
	const Box& vbx   = ba[krcv];
	const Box& bxrcv = amrex::grow(vbx, ng);
	
	if (pdomain.contains(bxrcv)) continue;

	if (check_local) {
	    localtouch.resize(bxrcv);
	    localtouch.setVal(0);
	}
	
	if (check_remote) {
	    remotetouch.resize(bxrcv);
	    remotetouch.setVal(0);
	}
	
	for (std::vector<IntVect>::const_iterator pit=pshifts.begin(); pit!=pshifts.end(); ++pit)
	{
	    if (*pit != IntVect::TheZeroVector())
	    {
		ba.intersections(bxrcv+(*pit), isects, false, ng);

		for (int j = 0, M = isects.size(); j < M; ++j)
		{
		    const int ksnd      = isects[j].first;
		    const Box& dst_bx   = isects[j].second - *pit;
		    const int src_owner = dm[ksnd];
		    
		    const BoxList& bl = amrex::boxDiff(dst_bx, pdomain);

		    for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
		    {
			Box sbx = (*lit) + (*pit);
			sbx &= pdomain; // source must be inside the periodic domain.
			
			if (sbx.ok()) {
			    Box dbx = sbx - (*pit);
			    if (ParallelDescriptor::sameTeam(src_owner)) { // local copy
				const BoxList tilelist(dbx, FabArrayBase::comm_tile_size);
				for (BoxList::const_iterator
					 it_tile  = tilelist.begin(),
					 End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
				{
				    m_LocTags->push_back(CopyComTag(*it_tile, (*it_tile)+(*pit), krcv, ksnd));
				}
				if (check_local) {
				    localtouch.plus(1, dbx);
				}
			    } else if (MyProc == dm[krcv]) {
				recv_tags[src_owner].push_back(CopyComTag(dbx, sbx, krcv, ksnd));
				if (check_remote) {
				    remotetouch.plus(1, dbx);
				}
			    }
			}
		    }
		}
	    }
	}

	if (check_local) {  
	    // safe if a cell is touched no more than once 
	    // keep checking thread safety if it is safe so far
	    check_local = m_threadsafe_loc = localtouch.max() <= 1;
	}

	if (check_remote) {
	    check_remote = m_threadsafe_rcv = remotetouch.max() <= 1;
	}
    }

    for (int ipass = 0; ipass < 2; ++ipass) // pass 0: send; pass 1: recv
    {
	CopyComTag::MapOfCopyComTagContainers & Tags = (ipass == 0) ? *m_SndTags : *m_RcvTags;
        for (auto& kv : Tags)
	{
	    std::vector<CopyComTag>& cctv = kv.second;
	    // We need to fix the order so that the send and recv processes match.
	    std::sort(cctv.begin(), cctv.end());
	}
    }
}

FabArrayBase::FB::~FB ()
{
    delete m_LocTags;
    delete m_SndTags;
    delete m_RcvTags;
}

void
FabArrayBase::flushFB (bool no_assertion) const
{
    BL_ASSERT(no_assertion || getBDKey() == m_bdkey);
    std::pair<FBCacheIter,FBCacheIter> er_it = m_TheFBCache.equal_range(m_bdkey);
    for (FBCacheIter it = er_it.first; it != er_it.second; ++it)
    {
#ifdef BL_MEM_PROFILING
	m_FBC_stats.bytes -= it->second->bytes();
#endif
	m_FBC_stats.recordErase(it->second->m_nuse);
	delete it->second;
    }
    m_TheFBCache.erase(er_it.first, er_it.second);
}

void
FabArrayBase::flushFBCache ()
{
    for (FBCacheIter it = m_TheFBCache.begin(); it != m_TheFBCache.end(); ++it)
    {
	m_FBC_stats.recordErase(it->second->m_nuse);
	delete it->second;
    }
    m_TheFBCache.clear();
#ifdef BL_MEM_PROFILING
    m_FBC_stats.bytes = 0L;
#endif
}

const FabArrayBase::FB&
FabArrayBase::getFB (const IntVect& nghost, const Periodicity& period,
                     bool cross, bool enforce_periodicity_only) const
{
    BL_PROFILE("FabArrayBase::getFB()");

    BL_ASSERT(getBDKey() == m_bdkey);
    std::pair<FBCacheIter,FBCacheIter> er_it = m_TheFBCache.equal_range(m_bdkey);
    for (FBCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	if (it->second->m_typ        == boxArray().ixType()      &&
            it->second->m_crse_ratio == boxArray().crseRatio()   &&
	    it->second->m_ngrow      == nghost                   &&
	    it->second->m_cross      == cross                    &&
	    it->second->m_epo        == enforce_periodicity_only &&
	    it->second->m_period     == period              )
	{
	    ++(it->second->m_nuse);
	    m_FBC_stats.recordUse();
	    return *(it->second);
	}
    }

    // Have to build a new one
    FB* new_fb = new FB(*this, nghost, cross, period, enforce_periodicity_only);

#ifdef BL_PROFILE
    m_FBC_stats.bytes += new_fb->bytes();
    m_FBC_stats.bytes_hwm = std::max(m_FBC_stats.bytes_hwm, m_FBC_stats.bytes);
#endif

    new_fb->m_nuse = 1;
    m_FBC_stats.recordBuild();
    m_FBC_stats.recordUse();

    m_TheFBCache.insert(er_it.second, FBCache::value_type(m_bdkey,new_fb));

    return *new_fb;
}

FabArrayBase::FPinfo::FPinfo (const FabArrayBase& srcfa,
			      const FabArrayBase& dstfa,
			      const Box&          dstdomain,
			      const IntVect&      dstng,
			      const BoxConverter& coarsener,
                              const Box&          cdomain,
                              const EB2::IndexSpace* index_space)
    : m_srcbdk   (srcfa.getBDKey()),
      m_dstbdk   (dstfa.getBDKey()),
      m_dstdomain(dstdomain),
      m_dstng    (dstng),
      m_coarsener(coarsener.clone()),
      m_nuse     (0)
{ 
    BL_PROFILE("FPinfo::FPinfo()");
    const BoxArray& srcba = srcfa.boxArray();
    const BoxArray& dstba = dstfa.boxArray();
    BL_ASSERT(srcba.ixType() == dstba.ixType());

    const IndexType& boxtype = dstba.ixType();
    BL_ASSERT(boxtype == dstdomain.ixType());
     
    BL_ASSERT(dstng.allLE(dstfa.nGrowVect()));

    const DistributionMapping& dstdm = dstfa.DistributionMap();
    
    const int myproc = ParallelDescriptor::MyProc();

    BoxList bl(boxtype);
    Vector<int> iprocs;

    for (int i = 0, N = dstba.size(); i < N; ++i)
    {
        Box bx = dstba[i];
        bx.grow(m_dstng);
        bx &= m_dstdomain;

        BoxList leftover = srcba.complementIn(bx);

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
        dm_crse_patch.define(std::move(iprocs));
#ifdef AMREX_USE_EB
        if (index_space)
        {
            fact_crse_patch = makeEBFabFactory(index_space, Geometry(cdomain),
                                               ba_crse_patch,
                                               dm_crse_patch,
                                               {0,0,0}, EBSupport::basic);
        }
        else
#endif
        {
            fact_crse_patch.reset(new FArrayBoxFactory());
        }
    }
}

FabArrayBase::FPinfo::~FPinfo ()
{
    delete m_coarsener;
}

long
FabArrayBase::FPinfo::bytes () const
{
    long cnt = sizeof(FabArrayBase::FPinfo);
    cnt += sizeof(Box) * (ba_crse_patch.capacity() + dst_boxes.capacity());
    cnt += sizeof(int) * (dm_crse_patch.capacity() + dst_idxs.capacity());
    return cnt;
}

const FabArrayBase::FPinfo&
FabArrayBase::TheFPinfo (const FabArrayBase& srcfa,
                         const FabArrayBase& dstfa,
                         const Box&          dstdomain,
                         const IntVect&      dstng,
                         const BoxConverter& coarsener,
                         const Box&          cdomain,
                         const EB2::IndexSpace* index_space)
{
    BL_PROFILE("FabArrayBase::TheFPinfo()");

    const BDKey& srckey = srcfa.getBDKey();
    const BDKey& dstkey = dstfa.getBDKey();

    std::pair<FPinfoCacheIter,FPinfoCacheIter> er_it = m_TheFillPatchCache.equal_range(dstkey);

    for (FPinfoCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	if (it->second->m_srcbdk    == srckey    &&
	    it->second->m_dstbdk    == dstkey    &&
	    it->second->m_dstdomain == dstdomain &&
	    it->second->m_dstng     == dstng     &&
	    it->second->m_dstdomain.ixType() == dstdomain.ixType() &&
	    it->second->m_coarsener->doit(it->second->m_dstdomain) == coarsener.doit(dstdomain))
	{
	    ++(it->second->m_nuse);
	    m_FPinfo_stats.recordUse();
	    return *(it->second);
	}
    }

    // Have to build a new one
    FPinfo* new_fpc = new FPinfo(srcfa, dstfa, dstdomain, dstng, coarsener, cdomain, index_space);

#ifdef BL_MEM_PROFILING
    m_FPinfo_stats.bytes += new_fpc->bytes();
    m_FPinfo_stats.bytes_hwm = std::max(m_FPinfo_stats.bytes_hwm, m_FPinfo_stats.bytes);
#endif
    
    new_fpc->m_nuse = 1;
    m_FPinfo_stats.recordBuild();
    m_FPinfo_stats.recordUse();

    m_TheFillPatchCache.insert(er_it.second, FPinfoCache::value_type(dstkey,new_fpc));
    if (srckey != dstkey)
	m_TheFillPatchCache.insert(          FPinfoCache::value_type(srckey,new_fpc));

    return *new_fpc;
}

void
FabArrayBase::flushFPinfo (bool no_assertion)
{
    BL_ASSERT(no_assertion || getBDKey() == m_bdkey);

    std::vector<FPinfoCacheIter> others;

    std::pair<FPinfoCacheIter,FPinfoCacheIter> er_it = m_TheFillPatchCache.equal_range(m_bdkey);

    for (FPinfoCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	const BDKey& srckey = it->second->m_srcbdk;
	const BDKey& dstkey = it->second->m_dstbdk;

	BL_ASSERT((srckey==dstkey && srckey==m_bdkey) || 
		  (m_bdkey==srckey) || (m_bdkey==dstkey));

	if (srckey != dstkey) {
	    const BDKey& otherkey = (m_bdkey == srckey) ? dstkey : srckey;
	    std::pair<FPinfoCacheIter,FPinfoCacheIter> o_er_it = m_TheFillPatchCache.equal_range(otherkey);

	    for (FPinfoCacheIter oit = o_er_it.first; oit != o_er_it.second; ++oit)
	    {
		if (it->second == oit->second)
		    others.push_back(oit);
	    }
	} 

#ifdef BL_MEM_PROFILING
	m_FPinfo_stats.bytes -= it->second->bytes();
#endif
	m_FPinfo_stats.recordErase(it->second->m_nuse);
	delete it->second;
    }
    
    m_TheFillPatchCache.erase(er_it.first, er_it.second);

    for (std::vector<FPinfoCacheIter>::iterator it = others.begin(),
	     End = others.end(); it != End; ++it)
    {
	m_TheFillPatchCache.erase(*it);
    }
}

FabArrayBase::CFinfo::CFinfo (const FabArrayBase& finefa,
                              const Geometry&     finegm,
                              const IntVect&      ng,
                              bool                include_periodic,
                              bool                include_physbndry)
    : m_fine_bdk (finefa.getBDKey()),
      m_ng       (ng),
      m_include_periodic(include_periodic),
      m_include_physbndry(include_physbndry),
      m_nuse     (0)
{
    BL_PROFILE("CFinfo::CFinfo()");
    
    m_fine_domain = Domain(finegm, ng, include_periodic, include_physbndry);

    const BoxArray& fba = amrex::convert(finefa.boxArray(), IndexType::TheCellType());
    const DistributionMapping& fdm = finefa.DistributionMap();

    BoxList bl(fba.ixType());
    Vector<int> iprocs;
    const int myproc = ParallelDescriptor::MyProc();

    for (int i = 0, N = fba.size(); i < N; ++i)
    {
        Box bx = fba[i];
        bx.grow(m_ng);
        bx &= m_fine_domain;

        const BoxList& noncovered = fba.complementIn(bx);
        for (const Box& b : noncovered) {
            bl.push_back(b);
            iprocs.push_back(fdm[i]);
            if (fdm[i] == myproc) {
                fine_grid_idx.push_back(i);
            }
        }
    }

    if (!iprocs.empty())
    {
        ba_cfb.define(bl);
        dm_cfb.define(std::move(iprocs));
    }
}

Box
FabArrayBase::CFinfo::Domain (const Geometry& geom, const IntVect& ng,
                              bool include_periodic, bool include_physbndry)
{
#if !defined(BL_NO_FORT)
    Box bx = geom.Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (Geometry::isPeriodic(idim)) {
            if (include_periodic) {
                bx.grow(idim, ng[idim]);
            }
        } else {
            if (include_physbndry) {
                bx.grow(idim, ng[idim]);
            }
        }
    }
    return bx;
#else
    return Box();
#endif
}

long
FabArrayBase::CFinfo::bytes () const
{
    long cnt = sizeof(FabArrayBase::CFinfo);
    cnt += sizeof(Box) * ba_cfb.capacity();
    cnt += sizeof(int) * (dm_cfb.capacity() + fine_grid_idx.capacity());
    return cnt;
}

const FabArrayBase::CFinfo&
FabArrayBase::TheCFinfo (const FabArrayBase& finefa,
                         const Geometry&     finegm,
                         const IntVect&      ng,
                         bool                include_periodic,
                         bool                include_physbndry)
{
    BL_PROFILE("FabArrayBase::TheCFinfo()");

    const BDKey& key = finefa.getBDKey();
    auto er_it = m_TheCrseFineCache.equal_range(key);
    for (auto it = er_it.first; it != er_it.second; ++it)
    {
        if (it->second->m_fine_bdk    == key                        &&
            it->second->m_fine_domain == CFinfo::Domain(finegm, ng,
                                                        include_periodic,
                                                        include_physbndry) &&
            it->second->m_ng          == ng)
        {
            ++(it->second->m_nuse);
            m_CFinfo_stats.recordUse();
            return *(it->second);
        }
    }

    // Have to build a new one
    CFinfo* new_cfinfo = new CFinfo(finefa, finegm, ng, include_periodic, include_physbndry);

#ifdef BL_MEM_PROFILING
    m_CFinfo_stats.bytes += new_cfinfo->bytes();
    m_CFinfo_stats.bytes_hwm = std::max(m_CFinfo_stats.bytes_hwm, m_CFinfo_stats.bytes);
#endif

    new_cfinfo->m_nuse = 1;
    m_CFinfo_stats.recordBuild();
    m_CFinfo_stats.recordUse();

    m_TheCrseFineCache.insert(er_it.second, CFinfoCache::value_type(key,new_cfinfo));

    return *new_cfinfo;
}

void
FabArrayBase::flushCFinfo (bool no_assertion)
{
    BL_ASSERT(no_assertion || getBDKey() == m_bdkey);
    auto er_it = m_TheCrseFineCache.equal_range(m_bdkey);
    for (auto it = er_it.first; it != er_it.second; ++it)
    {
#ifdef BL_MEM_PROFILING
        m_CFinfo_stats.bytes -= it->second->bytes();
#endif
        m_CFinfo_stats.recordErase(it->second->m_nuse);
        delete it->second;
    }
    m_TheCrseFineCache.erase(er_it.first, er_it.second);
}

void
FabArrayBase::Finalize ()
{
    FabArrayBase::flushFBCache();
    FabArrayBase::flushCPCache();
    FabArrayBase::flushTileArrayCache();

    if (ParallelDescriptor::IOProcessor() && amrex::system::verbose > 1) {
	m_FA_stats.print();
	m_TAC_stats.print();
	m_FBC_stats.print();
	m_CPC_stats.print();
	m_FPinfo_stats.print();
	m_CFinfo_stats.print();
    }

    m_TAC_stats = CacheStats("TileArrayCache");
    m_FBC_stats = CacheStats("FBCache");
    m_CPC_stats = CacheStats("CopyCache");
    m_FPinfo_stats = CacheStats("FillPatchCache");
    m_CFinfo_stats = CacheStats("CrseFineCache");

    m_BD_count.clear();
    
    m_FA_stats = FabArrayStats();

    the_fa_arena = nullptr;

    initialized = false;
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

        const IntVect& crse_ratio = boxArray().crseRatio();
	p = &FabArrayBase::m_TheTileArrayCache[m_bdkey][std::pair<IntVect,IntVect>(tilesize,crse_ratio)];
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
    const int N = indexArray.size();

    if (tileSize == IntVect::TheZeroVector())
    {
	for (int i = 0; i < N; ++i)
	{
	    if (isOwner(i))
	    {
		const int K = indexArray[i]; 
		const Box& bx = boxarray.getCellCenteredBox(K);
		ta.indexMap.push_back(K);
		ta.localIndexMap.push_back(i);
		ta.localTileIndexMap.push_back(0);
		ta.numLocalTiles.push_back(1);
		ta.tileArray.push_back(bx);
	    }
	}
    }
    else
    {
	std::vector<int> local_idxs(N);
	std::iota(std::begin(local_idxs), std::end(local_idxs), 0);

#if defined(BL_USE_TEAM)
	const int nworkers = ParallelDescriptor::TeamSize();
	if (nworkers > 1) {
	    // reorder it so that each worker will be more likely to work on their own fabs
	    std::stable_sort(local_idxs.begin(), local_idxs.end(), [this](int i, int j) 
			     { return  this->distributionMap[this->indexArray[i]] 
				     < this->distributionMap[this->indexArray[j]]; });
	}
#endif	

	for (std::vector<int>::const_iterator it = local_idxs.begin(); it != local_idxs.end(); ++it)
	{
	    const int i = *it;         // local index 
	    const int K = indexArray[i]; // global index
	    const Box& bx = boxarray.getCellCenteredBox(K);

            //
            //  This must be consistent with ParticleContainer::getTileIndex function!!!
            //
	    
	    IntVect nt_in_fab, tsize, nleft;
	    int ntiles = 1;
	    for (int d=0; d<AMREX_SPACEDIM; d++) {
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
		ta.localTileIndexMap.push_back(t);
		ta.numLocalTiles.push_back(ntiles);

		for (int d=0; d<AMREX_SPACEDIM; d++) {
		    if (ijk[d]<nt_in_fab[d]-1) {
			ijk[d]++;
			break;
		    } else {
			ijk[d] = 0;
		    }
		}
		
		for (int d=0; d<AMREX_SPACEDIM; d++) {
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
FabArrayBase::flushTileArray (const IntVect& tileSize, bool no_assertion) const
{
    BL_ASSERT(no_assertion || getBDKey() == m_bdkey);

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
            const IntVect& crse_ratio = boxArray().crseRatio();
	    TAMap::iterator tai_it = tai.find(std::pair<IntVect,IntVect>(tileSize,crse_ratio));
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
    BL_ASSERT(boxarray.empty() || no_assertion || getBDKey() == m_bdkey);

    std::map<BDKey, int>::iterator cnt_it = m_BD_count.find(m_bdkey);
    if (cnt_it != m_BD_count.end()) 
    {
        --(cnt_it->second);
        if (cnt_it->second == 0) 
        {
            m_BD_count.erase(cnt_it);
            
            // Since this is the last one built with these BoxArray 
            // and DistributionMapping, erase it from caches.
            flushTileArray(IntVect::TheZeroVector(), no_assertion);
            flushFPinfo(no_assertion);
            flushCFinfo(no_assertion);
            flushFB(no_assertion);
            flushCPC(no_assertion);
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


void
FabArrayBase::WaitForAsyncSends (int                 N_snds,
                                 Vector<MPI_Request>& send_reqs,
                                 Vector<char*>&       send_data,
                                 Vector<MPI_Status>&  stats)
{
#ifdef BL_USE_MPI
    BL_ASSERT(N_snds > 0);

    stats.resize(N_snds);

    BL_ASSERT(send_reqs.size() == N_snds);
    BL_ASSERT(send_data.size() == N_snds);

    ParallelDescriptor::Waitall(send_reqs, stats);
#endif /*BL_USE_MPI*/
}


#ifdef BL_USE_MPI
bool
FabArrayBase::CheckRcvStats(Vector<MPI_Status>& recv_stats,
			    const Vector<int>& recv_size,
			    MPI_Datatype datatype, int tag)
{
    bool r = true;
    for (int i = 0, n = recv_size.size(); i < n; ++i) {
	if (recv_size[i] > 0) {
	    int count;

	    MPI_Get_count(&recv_stats[i], datatype, &count);

	    if (count != recv_size[i]) {
		r = false;
                if (amrex::Verbose()) {
                    amrex::AllPrint() << "ERROR: Proc. " << ParallelContext::MyProcSub()
                                      << " received " << count << " counts of data from Proc. "
                                      << recv_stats[i].MPI_SOURCE
                                      << " with tag " << recv_stats[i].MPI_TAG
                                      << " error " << recv_stats[i].MPI_ERROR
                                      << ", but the expected counts is " << recv_size[i]
                                      << " with tag " << tag << "\n";
                }
	    }
	}
    }
    return r;
}
#endif

std::ostream&
operator<< (std::ostream& os, const FabArrayBase::BDKey& id)
{
    os << "(" << id.m_ba_id << ", " << id.m_dm_id << ")";
    return os;
}

}
