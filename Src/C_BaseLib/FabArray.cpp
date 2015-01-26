#include <winstd.H>

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

    BoxLib::ExecOnFinalize(FabArrayBase::Finalize);

    initialized = true;
}

FabArrayBase::FabArrayBase ()
{
    Initialize();
}

FabArrayBase::~FabArrayBase () {}

const Box
FabArrayBase::fabbox (int K) const
{
    return BoxLib::grow(boxarray[K], n_grow);
}

//
// Stuff used for copy() caching.
//

FabArrayBase::CPC::CPC ()
    :
    m_reused(false),
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
                        const DistributionMapping& srcdm)
    :
    m_dstba(dstba),
    m_srcba(srcba),
    m_dstdm(dstdm),
    m_srcdm(srcdm),
    m_reused(false),
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

//
// The copy() cache.
//
FabArrayBase::CPCCache FabArrayBase::m_TheCopyCache;

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
            it->second.m_reused = true;

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

            if (!it->second.m_reused)
                erase_it = it;
        }

        if (erase_it != End)
        {
            TheCopyCache.erase(erase_it);
        }
        else if (last_it != End)
        {
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

    if (dst.IndexMap().empty() && src.IndexMap().empty())
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
        return cache_it;

    std::vector< std::pair<int,Box> > isects;

    for (int i = 0, N = TheCPC.m_dstba.size(); i < N; i++)
    {
        TheCPC.m_srcba.intersections(TheCPC.m_dstba[i],isects);

        const int dst_owner = TheCPC.m_dstdm[i];

        for (int j = 0, M = isects.size(); j < M; j++)
        {
            const Box& bx        = isects[j].second;
            const int  k         = isects[j].first;
            const int  src_owner = TheCPC.m_srcdm[k];

            if (dst_owner != MyProc && src_owner != MyProc) continue;

            CopyComTag tag;

            tag.box      = bx;
            tag.fabIndex = i;
            tag.srcIndex = k;

            if (dst_owner == MyProc)
            {
                if (src_owner == MyProc)
                {
                    TheCPC.m_LocTags->push_back(tag);
                }
                else
                {
                    FabArrayBase::SetRecvTag(*TheCPC.m_RcvTags,src_owner,tag,*TheCPC.m_RcvVols,bx);
                }
            }
            else if (src_owner == MyProc)
            {
                FabArrayBase::SetSendTag(*TheCPC.m_SndTags,dst_owner,tag,*TheCPC.m_SndVols,bx);
            }
        }
    }
    //
    // Squeeze out any unused memory ...
    //
    CopyComTagsContainer tmp(*TheCPC.m_LocTags); 

    TheCPC.m_LocTags->swap(tmp);

    for (MapOfCopyComTagContainers::iterator it = TheCPC.m_SndTags->begin(), End = TheCPC.m_SndTags->end(); it != End; ++it)
    {
        CopyComTagsContainer tmp(it->second);

        it->second.swap(tmp);
    }

    for (MapOfCopyComTagContainers::iterator it = TheCPC.m_RcvTags->begin(), End = TheCPC.m_RcvTags->end(); it != End; ++it)
    {
        CopyComTagsContainer tmp(it->second);

        it->second.swap(tmp);
    }

    TheCPC.m_srcba.clear_hash_bin();

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
FabArrayBase::EraseFromTheCPC (CPCCacheIter& cache_it)
{
    FabArrayBase::m_TheCopyCache.erase(cache_it);    
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
        if (it->second.m_reused)
            stats[1]++;
    }

    if (FabArrayBase::Verbose)
    {
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
    }

    m_TheCopyCache.clear();
}

FabArrayBase::SI::SI ()
    :
    m_ngrow(-1),
    m_cross(false),
    m_reused(false),
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
    m_cross(cross),
    m_reused(false),
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

FabArrayBase::FBCache FabArrayBase::m_TheFBCache;

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
            it->second.m_reused = true;

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

            if (!it->second.m_reused)
                erase_it = it;
        }

        if (erase_it != End)
        {
            m_TheFBCache.erase(erase_it);
        }
        else if (last_it != End)
        {
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

    if (mf.IndexMap().empty())
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
        return cache_it;

    std::vector<Box>                  boxes;
    std::vector< std::pair<int,Box> > isects;

    boxes.resize(si.m_cross ? 2*BL_SPACEDIM : 1);

    for (int i = 0, N = ba.size(); i < N; i++)
    {
        const Box& vbx = ba[i];

        if (si.m_cross)
        {
            for (int dir = 0; dir < BL_SPACEDIM; dir++)
            {
                Box lo = vbx;
                lo.setSmall(dir, vbx.smallEnd(dir) - si.m_ngrow);
                lo.setBig  (dir, vbx.smallEnd(dir) - 1);
                boxes[2*dir+0] = lo;

                Box hi = vbx;
                hi.setSmall(dir, vbx.bigEnd(dir) + 1);
                hi.setBig  (dir, vbx.bigEnd(dir) + si.m_ngrow);
                boxes[2*dir+1] = hi;
            }
        }
        else
        {
            boxes[0] = BoxLib::grow(vbx,si.m_ngrow);
        }

        const int dst_owner = dm[i];

        for (std::vector<Box>::const_iterator it = boxes.begin(),
                 End = boxes.end();
             it != End;
             ++it)
        {
            ba.intersections(*it,isects);

            for (int j = 0, M = isects.size(); j < M; j++)
            {
                const int  k         = isects[j].first;
                const Box& bx        = isects[j].second;
                const int  src_owner = dm[k];

                if ( (k == i) || (dst_owner != MyProc && src_owner != MyProc) ) continue;

                CopyComTag tag;

                tag.box      = bx;
                tag.fabIndex = i;
                tag.srcIndex = k;

                if (dst_owner == MyProc)
                {
                    if (src_owner == MyProc)
                    {
                        TheFB.m_LocTags->push_back(tag);
                    }
                    else
                    {
                        FabArrayBase::SetRecvTag(*TheFB.m_RcvTags,src_owner,tag,*TheFB.m_RcvVols,bx);
                    }
                }
                else if (src_owner == MyProc)
                {
                    FabArrayBase::SetSendTag(*TheFB.m_SndTags,dst_owner,tag,*TheFB.m_SndVols,bx);
                }
            }
        }
    }
    //
    // Squeeze out any unused memory ...
    //
    CopyComTagsContainer tmp(*TheFB.m_LocTags); 

    TheFB.m_LocTags->swap(tmp);

    for (MapOfCopyComTagContainers::iterator it = TheFB.m_SndTags->begin(), End = TheFB.m_SndTags->end(); it != End; ++it)
    {
        CopyComTagsContainer tmp(it->second);

        it->second.swap(tmp);
    }

    for (MapOfCopyComTagContainers::iterator it = TheFB.m_RcvTags->begin(), End = TheFB.m_RcvTags->end(); it != End; ++it)
    {
        CopyComTagsContainer tmp(it->second);

        it->second.swap(tmp);
    }

    ba.clear_hash_bin();

    //
    // set thread safety
    //
#ifdef _OPENMP
    if (ba[0].cellCentered()) {
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
        if (it->second.m_reused)
            stats[1]++;
    }

    if (FabArrayBase::Verbose)
    {
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


MFIter::MFIter (const FabArrayBase& fabarray, int sharing)
    :
    fabArray(fabarray),
    currentIndex(0),
    tileSize(D_DECL(1024000,1024000,1024000))
{
    Initialize(sharing,0);
}

MFIter::MFIter (const FabArrayBase& fabarray, bool do_tiling, int chunksize)
    :
    fabArray(fabarray),
    currentIndex(0),
    tileSize(D_DECL(1024000,1024000,1024000))
{
    if (do_tiling) 
	tileSize = FabArrayBase::mfiter_tile_size;
    Initialize(1,chunksize);
}

MFIter::MFIter (const FabArrayBase& fabarray, const IntVect& tilesize, int chunksize)
    :
    fabArray(fabarray),
    currentIndex(0),
    tileSize(tilesize)
{
    Initialize(1,chunksize);
}

void 
MFIter::Initialize (int sharing, int chunksize) 
{
    int tid = 0;
    int nthreads = 1;
    
#ifdef _OPENMP
    if (omp_in_parallel() && sharing) {
	tid = omp_get_thread_num();
	nthreads = omp_get_num_threads();
    }
#endif
    
    // First let's figure out how many tiles we are going to have
    int n_tot_tiles=0;
    Array<IntVect> nt_in_fab;
    for (int i=0; i<fabArray.IndexMap().size(); i++) {
	int K = fabArray.IndexMap()[i]; 
	const Box& bx = fabArray.box(K);
	
	int ntiles = 1;
	IntVect nt;
	for (int d=0; d<BL_SPACEDIM; d++) {
	    nt[d] = std::max(bx.length(d)/tileSize[d], 1);
	    ntiles *= nt[d];
	}
	
	nt_in_fab.push_back(nt);
	n_tot_tiles += ntiles;
    }

    int tlo, thi, sz;
    if (chunksize <= 0) {
	// figure out the tile no range, tlo and thi for this thread
	if (n_tot_tiles < nthreads) { // there are more threads than tiles
	    if (tid < n_tot_tiles) {
		tlo = thi = tid;
	    } else {
		return;
	    }
	} else {
	    int tiles_per_thread = n_tot_tiles/nthreads;
	    int nleft = n_tot_tiles - tiles_per_thread*nthreads;
	    if (tid < nleft) {
		tlo = tid*(tiles_per_thread+1);
		thi = tlo + tiles_per_thread;
	    } else {
	    tlo = tid*tiles_per_thread + nleft;
	    thi = tlo + tiles_per_thread - 1;
	    }
	}

	sz = thi-tlo+1;
    }
    else {
	tlo = -1;
	thi = n_tot_tiles;
	int nchunks = (n_tot_tiles+chunksize-1) / chunksize;
	int chunks_per_thread = nchunks/nthreads;
	int nleft = nchunks-chunks_per_thread*nthreads;
	if (tid < nleft) chunks_per_thread++;
	sz = chunks_per_thread*chunksize;
    }

    indexMap.reserve(sz);
    localIndexMap.reserve(sz);
    tileArray.reserve(sz);
	
    int it = 0;
    for (int i=0; i<fabArray.IndexMap().size(); i++) {
	
	int ntiles = 1;
	for (int d=0; d<BL_SPACEDIM; d++) {
	    ntiles *= nt_in_fab[i][d];
	}
	
	if (it+ntiles-1 < tlo) {
	    it += ntiles;
	    continue;
	}
	
	int K = fabArray.IndexMap()[i]; 
	const Box& bx = fabArray.box(K);
	
	IntVect tsize, nleft;
	for (int d=0; d<BL_SPACEDIM; d++) {
	    int ncells = bx.length(d);
	    tsize[d] = ncells/nt_in_fab[i][d];
	    nleft[d] = ncells - nt_in_fab[i][d]*tsize[d];
	}
	
	IntVect small, big, ijk;  // note that the initial values are all zero.
	ijk[0] = -1;
	for (int t=0; t<ntiles; t++) {
	    for (int d=0; d<BL_SPACEDIM; d++) {
		if (ijk[d]<nt_in_fab[i][d]-1) {
		    ijk[d]++;
		    break;
		} else {
		    ijk[d] = 0;
		}
	    }
	    
	    if (chunksize > 0 && (it/chunksize)%nthreads == tid
		|| chunksize <=0 && it >= tlo) {
		for (int d=0; d<BL_SPACEDIM; d++) {
		    if (ijk[d] < nleft[d]) {
			small[d] = ijk[d]*(tsize[d]+1);
			big[d] = small[d] + tsize[d];
		    } else {
			small[d] = ijk[d]*tsize[d] + nleft[d];
			big[d] = small[d] + tsize[d] - 1;
		    }
		}
		
		indexMap.push_back(K);
		localIndexMap.push_back(i);
		
		Box tbx(small, big, bx.ixType());
		tbx.shift(bx.smallEnd());
		tileArray.push_back(tbx);
	    }
	    
	    it++;
	    if (it > thi) return;
	}
    }
}

Box
MFIter::nodaltilebox (int dir) const 
{ 
    Box bx(tileArray[currentIndex]);
    if (bx.type(dir) == IndexType::CELL) {
	bx.surroundingNodes(dir);
	if (bx.bigEnd(dir) != validbox().bigEnd(dir)+1) {
	    bx.growHi(dir,-1);
	}
    }
    return bx;
}

Box 
MFIter::growntilebox (int ng) const 
{
    if (ng < 0) ng = fabArray.nGrow();
    if (ng == 0) {
	return tileArray[currentIndex];
    } else {
	Box bx(tileArray[currentIndex]);
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
}
