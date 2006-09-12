#include <winstd.H>

#include <map>

#include <FabArray.H>

FabArrayBase::FabArrayBase ()
{}

FabArrayBase::FabArrayBase (const BoxArray& bxs, int nvar, int ngrow)
    :
    boxarray(bxs),
    distributionMap(boxarray, ParallelDescriptor::NProcsCFD()),
    n_grow(ngrow),
    n_comp(nvar)
{}

FabArrayBase::FabArrayBase (const BoxArray&            bxs,
                            int                        nvar,
                            int                        ngrow,
                            const DistributionMapping& map)
    :
    boxarray(bxs),
    distributionMap(map),
    n_grow(ngrow),
    n_comp(nvar)
{}

FabArrayBase::~FabArrayBase ()
{}

Box
FabArrayBase::fabbox (int K) const
{
    //
    // Do not use fabparray[K] because it may not be valid in parallel.
    //
    return BoxLib::grow(boxarray[K], n_grow);
}

bool MFIter::g_debugging = false;

MFIter::MFIter (const FabArrayBase& fabarray)
    :
    fabArray(fabarray),
    currentIndex(0),
    m_debugging(g_debugging)
{
    //
    // Increment the currentIndex to start at the first valid index
    // for this ParallelDescriptor::MyProc.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    while (fabArray.DistributionMap()[currentIndex] != MyProc)
    {
        ++currentIndex;
    }
}

void
MFIter::operator++ ()
{
    const int MyProc = ParallelDescriptor::MyProc();
    //
    // Go to the next index on this processor.
    //
    do
    {
        ++currentIndex;
    }
    while (fabArray.DistributionMap()[currentIndex] != MyProc);
}

void
MFIter::setDebugging (bool debugging)
{
    g_debugging = debugging;
}

bool
MFIter::isValid ()
{
    BL_ASSERT(currentIndex >= 0);

    bool rc = currentIndex < fabArray.size();

    if (m_debugging)
    {
        if (rc) return true;
        ParallelDescriptor::Barrier();
        return false;
    }

    return rc;
}

const Box&
MFIter::validbox () const
{
    return fabArray.box(currentIndex);
}

Box
MFIter::fabbox () const
{
    return fabArray.fabbox(currentIndex);
}

FillBoxId::FillBoxId ()
    :
    m_fillBoxId(-1),
    m_fabIndex(-1)
{}

FillBoxId::FillBoxId (int        newid,
		      const Box& fillbox)
    :
    m_fillBox(fillbox),
    m_fillBoxId(newid),
    m_fabIndex(-1)
{}

//
// Used to cache some CommData stuff in CollectData().
//

CommDataCache::CommDataCache ()
    :
    m_valid(false)
{}

void
CommDataCache::operator= (const Array<CommData>& rhs)
{
    m_commdata = rhs;
    m_valid    = true;
}

//
// Stuff used for copy() caching.
//

CopyComTag::CopyComTag () {}

CopyComTag::CopyComTag (const CopyComTag& cct)
    :
    box(cct.box),
    fabIndex(cct.fabIndex),
    srcIndex(cct.srcIndex)
{}

CPC::CPC ()
    :
    m_reused(false)
{}

CPC::CPC (const BoxArray&            dstba,
          const BoxArray&            srcba,
          const DistributionMapping& dstdm,
          const DistributionMapping& srcdm)
    :
    m_dstba(dstba),
    m_srcba(srcba),
    m_dstdm(dstdm),
    m_srcdm(srcdm),
    m_reused(false)
{}

CPC::CPC (const CPC& rhs)
    :
    m_dstba(rhs.m_dstba),
    m_srcba(rhs.m_srcba),
    m_dstdm(rhs.m_dstdm),
    m_srcdm(rhs.m_srcdm),
    m_reused(rhs.m_reused)
{}

CPC::~CPC () {}

bool
CPC::operator== (const CPC& rhs) const
{
    return m_dstba == rhs.m_dstba &&
           m_srcba == rhs.m_srcba &&
           m_dstdm == rhs.m_dstdm &&
           m_srcdm == rhs.m_srcdm;
}

bool
CPC::operator!= (const CPC& rhs) const
{
    return !operator==(rhs);
}

typedef std::multimap<int,CPC> CPCCache;

typedef CPCCache::iterator CPCCacheIter;

static CPCCache TheCopyCache;

CPC&
CPC::TheCPC (const CPC& cpc, bool& got_from_cache)
{
    static bool first          = false;
    static bool use_copy_cache = true;

    if (!first)
    {
        first = true;
        ParmParse pp("fabarray");
        pp.query("use_copy_cache", use_copy_cache);
    }

    got_from_cache = false;

    const int key = cpc.m_dstba.size() + cpc.m_srcba.size();

    if (use_copy_cache)
    {
        std::pair<CPCCacheIter,CPCCacheIter> er_it = TheCopyCache.equal_range(key);

        for (CPCCacheIter it = er_it.first; it != er_it.second; ++it)
        {
            if (it->second == cpc)
            {
                it->second.m_reused = true;
                got_from_cache = true;
                return it->second;
            }
        }
        //
        // Don't let the size of the cache get too big.
        //
        const int N = 50;

        if (TheCopyCache.size() >= N)
        {
            for (CPCCacheIter it = TheCopyCache.begin(); it != TheCopyCache.end(); )
            {
                if (!it->second.m_reused)
                {
                    TheCopyCache.erase(it++);
                }
                else
                {
                    ++it;
                }
            }

            if (TheCopyCache.size() >= N)
                //
                // Get rid of first entry which is the one with the smallest key.
                //
                TheCopyCache.erase(TheCopyCache.begin());
        }
    }
    else
    {
        TheCopyCache.clear();
    }

    CPCCacheIter it = TheCopyCache.insert(std::make_pair(key,cpc));

    return it->second;
}

void
CPC::FlushCache ()
{
    TheCopyCache.clear();
}
