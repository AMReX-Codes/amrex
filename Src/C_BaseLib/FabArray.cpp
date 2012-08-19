#include <winstd.H>

#include <FabArray.H>
#include <ParmParse.H>
//
// Set default values in Initialize()!!!
//
bool FabArrayBase::verbose;
bool FabArrayBase::do_async_sends;
bool FabArrayBase::do_not_use_cache;

namespace
{
    bool initialized = false;
    //
    // Set default values in Initialize()!!!
    //
    bool use_copy_cache;
    int  copy_cache_max_size;
    bool use_fb_cache;
    int  fb_cache_max_size;
}

void
FabArrayBase::Initialize ()
{
    if (initialized) return;
    //
    // Set default values here!!!
    //
    FabArrayBase::verbose          = false;
    FabArrayBase::do_async_sends   = false;
    FabArrayBase::do_not_use_cache = false;

    use_copy_cache      = true;
    copy_cache_max_size = 50;   // -1 ==> no maximum size
    use_fb_cache        = true;
    fb_cache_max_size   = 50;   // -1 ==> no maximum size

    ParmParse pp("fabarray");

    pp.query("verbose",          FabArrayBase::verbose);
    pp.query("do_async_sends",   FabArrayBase::do_async_sends);
    pp.query("do_not_use_cache", FabArrayBase::do_not_use_cache);

    pp.query("use_copy_cache",      use_copy_cache);
    pp.query("copy_cache_max_size", copy_cache_max_size);
    pp.query("use_fb_cache",        use_fb_cache);
    pp.query("fb_cache_max_size",   fb_cache_max_size);

    if (fb_cache_max_size <= 0 && fb_cache_max_size != -1)
        use_fb_cache = false;

    if (copy_cache_max_size <= 0 && copy_cache_max_size != -1)
        use_copy_cache = false;

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
// Used to cache some CommData stuff in CollectData().
//

FabArrayBase::CommDataCache::CommDataCache ()
    :
    m_valid(false)
{}

void
FabArrayBase::CommDataCache::operator= (const Array<ParallelDescriptor::CommData>& rhs)
{
    m_commdata = rhs;
    m_valid    = true;
}

//
// Stuff used for copy() caching.
//

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
    return m_dstba == rhs.m_dstba && m_srcba == rhs.m_srcba && m_dstdm == rhs.m_dstdm && m_srcdm == rhs.m_srcdm;
}

typedef std::multimap<int,FabArrayBase::CPC> CPCCache;

typedef CPCCache::iterator CPCCacheIter;

static CPCCache TheCopyCache;

FabArrayBase::CPC&
FabArrayBase::CPC::TheCPC (const CPC& cpc)
{
    const int Key = cpc.m_dstba.size() + cpc.m_srcba.size();

    if (use_copy_cache)
    {
        std::pair<CPCCacheIter,CPCCacheIter> er_it = TheCopyCache.equal_range(Key);

        for (CPCCacheIter it = er_it.first; it != er_it.second; ++it)
        {
            if (it->second == cpc)
            {
                it->second.m_reused = true;

                return it->second;
            }
        }

        if (TheCopyCache.size() >= copy_cache_max_size && copy_cache_max_size != -1)
        {
            //
            // Don't let the size of the cache get too big.
            //
            for (CPCCacheIter it = TheCopyCache.begin(); it != TheCopyCache.end(); )
            {
                if (!it->second.m_reused)
                {
                    TheCopyCache.erase(it++);

                    if (TheCopyCache.size() < copy_cache_max_size)
                        //
                        // Only delete enough entries to stay under limit.
                        //
                        break;
                }
                else
                {
                    ++it;
                }
            }

            if (TheCopyCache.size() >= copy_cache_max_size)
            {
                //
                // Get rid of entry with the smallest key.
                //
                TheCopyCache.erase(TheCopyCache.begin());
            }
        }
    }
    else
    {
        TheCopyCache.clear();
    }
    //
    // Got to build it.
    //
    const int    MyProc = ParallelDescriptor::MyProc();
    CPCCacheIter it     = TheCopyCache.insert(std::make_pair(Key,cpc));
    CPC&         TheCPC = it->second;
    //
    // Here is where we allocate space for the stuff used in the cache.
    //
    TheCPC.m_LocTags = new CPC::CopyComTagsContainer;
    TheCPC.m_SndTags = new CPC::MapOfCopyComTagContainers;
    TheCPC.m_RcvTags = new CPC::MapOfCopyComTagContainers;
    TheCPC.m_SndVols = new std::map<int,int>;
    TheCPC.m_RcvVols = new std::map<int,int>;

    CopyComTag                        tag;
    std::vector< std::pair<int,Box> > isects;

    for (int i = 0, N = TheCPC.m_dstba.size(); i < N; i++)
    {
        TheCPC.m_srcba.intersections(TheCPC.m_dstba[i],isects);

        const int d_owner = TheCPC.m_dstdm[i];

        for (int j = 0, M = isects.size(); j < M; j++)
        {
            const Box& bx      = isects[j].second;
            const int  k       = isects[j].first;
            const int  s_owner = TheCPC.m_srcdm[k];

            if (d_owner != MyProc && s_owner != MyProc) continue;

            tag.box = bx;

            const int vol = bx.numPts();

            if (d_owner == MyProc)
            {
                tag.fabIndex = i;

                if (s_owner == MyProc)
                {
                    tag.srcIndex = k;

                    TheCPC.m_LocTags->push_back(tag);
                }
                else
                {
                    (*TheCPC.m_RcvTags)[s_owner].push_back(tag);

                    if (TheCPC.m_RcvVols->count(s_owner) > 0)
                    {
                        (*TheCPC.m_RcvVols)[s_owner] += vol;
                    }
                    else
                    {
                        (*TheCPC.m_RcvVols)[s_owner] = vol;
                    }
                }
            }
            else if (s_owner == MyProc)
            {
                tag.fabIndex = k;

                (*TheCPC.m_SndTags)[d_owner].push_back(tag);

                if (TheCPC.m_SndVols->count(d_owner) > 0)
                {
                    (*TheCPC.m_SndVols)[d_owner] += vol;
                }
                else
                {
                    (*TheCPC.m_SndVols)[d_owner] = vol;
                }
            }
        }
    }

    return TheCPC;
}

void
FabArrayBase::CPC::FlushCache ()
{
    if (ParallelDescriptor::IOProcessor() && !TheCopyCache.empty() && FabArrayBase::verbose)
    {
        int reused = 0;

        for (CPCCacheIter it = TheCopyCache.begin(), End = TheCopyCache.end(); it != End; ++it)
            if (it->second.m_reused)
                reused++;

        std::cout << "CPC::TheCopyCache.size() = " << TheCopyCache.size() << ", # reused = " << reused << '\n';
    }
    TheCopyCache.clear();
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
    return m_ngrow == rhs.m_ngrow && m_cross == rhs.m_cross && m_ba == rhs.m_ba && m_dm == rhs.m_dm;
}

typedef std::multimap<int,FabArrayBase::SI> SIMMap;

typedef SIMMap::iterator SIMMapIter;

static SIMMap SICache;

void
FabArrayBase::Finalize ()
{
    SICache.clear();

    TheCopyCache.clear();

    initialized = false;
}

void
FabArrayBase::FlushSICache ()
{
    if (ParallelDescriptor::IOProcessor() && !SICache.empty() && FabArrayBase::verbose)
    {
        int reused = 0;

        for (SIMMapIter it = SICache.begin(), End = SICache.end(); it != End; ++it)
            if (it->second.m_reused)
                reused++;

        std::cout << "FabArrayBase::SICache.size() = " << SICache.size() << ", # reused = " << reused << '\n';
    }

    SICache.clear();
}

int
FabArrayBase::SICacheSize ()
{
    return SICache.size();
}

FabArrayBase::SI&
FabArrayBase::TheFBsirec (bool                cross,
                          const FabArrayBase& mf)
{
    const FabArrayBase::SI si(mf.boxArray(), mf.DistributionMap(), mf.nGrow(), cross);

    const int Key = mf.size() + mf.nGrow() + cross;

    if (use_fb_cache)
    {
        std::pair<SIMMapIter,SIMMapIter> er_it = SICache.equal_range(Key);
    
        for (SIMMapIter it = er_it.first; it != er_it.second; ++it)
        {
            if (it->second == si)
            {
                it->second.m_reused = true;

                return it->second;
            }
        }

        if (SICache.size() >= fb_cache_max_size && fb_cache_max_size != -1)
        {
            //
            // Don't let the size of the cache get too big.
            //
            for (SIMMapIter it = SICache.begin(); it != SICache.end(); )
            {
                if (!it->second.m_reused)
                {
                    SICache.erase(it++);

                    if (SICache.size() < fb_cache_max_size)
                        //
                        // Only delete enough entries to stay under limit.
                        //
                        break;
                }
                else
                {
                    ++it;
                }
            }

            if (SICache.size() >= fb_cache_max_size)
            {
                //
                // Get rid of entry with the smallest key.
                //
                SICache.erase(SICache.begin());
            }
        }
    }
    else
    {
        SICache.clear();
    }
    //
    // Got to build one.
    //
    SIMMapIter                 it     = SICache.insert(std::make_pair(Key,si));
    const BoxArray&            ba     = mf.boxArray();
    const DistributionMapping& dm     = mf.DistributionMap();
    const int                  MyProc = ParallelDescriptor::MyProc();
    SI&                        TheSI  = it->second;
    //
    // Here is where we allocate space for the stuff used in the cache.
    //
    TheSI.m_LocTags = new SI::CopyComTagsContainer;
    TheSI.m_SndTags = new SI::MapOfCopyComTagContainers;
    TheSI.m_RcvTags = new SI::MapOfCopyComTagContainers;
    TheSI.m_SndVols = new std::map<int,int>;
    TheSI.m_RcvVols = new std::map<int,int>;

    CopyComTag                        tag;
    std::vector<Box>                  boxes;
    std::vector< std::pair<int,Box> > isects;

    boxes.resize(si.m_cross ? 2*BL_SPACEDIM : 1);

    for (int i = 0, N = ba.size(); i < N; i++)
    {
        if (si.m_cross)
        {
            const Box& vbx = ba[i];

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
            boxes[0] = BoxLib::grow(ba[i],si.m_ngrow);
        }

        const int d_owner = dm[i];

        for (std::vector<Box>::const_iterator it = boxes.begin(),
                 End = boxes.end();
             it != End;
             ++it)
        {
            ba.intersections(*it,isects);

            for (int j = 0, M = isects.size(); j < M; j++)
            {
                const Box& bx      = isects[j].second;
                const int  k       = isects[j].first;
                const int  s_owner = dm[k];

                if (d_owner != MyProc && s_owner != MyProc) continue;

                if (k == i) continue;

                const int vol = bx.numPts();

                tag.box = bx;

                if (d_owner == MyProc)
                {
                    tag.fabIndex = i;

                    if (s_owner == MyProc)
                    {
                        tag.srcIndex = k;

                        TheSI.m_LocTags->push_back(tag);
                    }
                    else
                    {
                        (*TheSI.m_RcvTags)[s_owner].push_back(tag);

                        if (TheSI.m_RcvVols->count(s_owner) > 0)
                        {
                            (*TheSI.m_RcvVols)[s_owner] += vol;
                        }
                        else
                        {
                            (*TheSI.m_RcvVols)[s_owner] = vol;
                        }
                    }
                }
                else if (s_owner == MyProc)
                {
                    tag.fabIndex = k;

                    (*TheSI.m_SndTags)[d_owner].push_back(tag);

                    if (TheSI.m_SndVols->count(d_owner) > 0)
                    {
                        (*TheSI.m_SndVols)[d_owner] += vol;
                    }
                    else
                    {
                        (*TheSI.m_SndVols)[d_owner] = vol;
                    }
                }
            }
        }
    }

    return TheSI;
}
