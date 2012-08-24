#include <winstd.H>

#include <FabArray.H>
#include <ParmParse.H>
//
// Set default values in Initialize()!!!
//
bool FabArrayBase::verbose;
bool FabArrayBase::do_async_sends;

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
    FabArrayBase::verbose          = true;
    FabArrayBase::do_async_sends   = false;

    copy_cache_max_size = 75;
    fb_cache_max_size   = 75;

    ParmParse pp("fabarray");

    pp.query("verbose",             FabArrayBase::verbose);
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

int
FabArrayBase::CPC::bytes () const
{
    //
    // Get a estimate on number of bytes used by a CPC.
    // This doesn't count any "overhead" in the STL containers.
    //
    int cnt = 0;

    cnt += m_LocTags->size()*sizeof(CopyComTag);

    for (MapOfCopyComTagContainers::const_iterator it = m_SndTags->begin(),
             m_End = m_SndTags->end();
         it != m_End;
         ++it)
    {
        cnt += it->second.size()*sizeof(CopyComTag);
    }

    for (MapOfCopyComTagContainers::const_iterator it = m_RcvTags->begin(),
             m_End = m_RcvTags->end();
         it != m_End;
         ++it)
    {
        cnt += it->second.size()*sizeof(CopyComTag);
    }

    cnt += 2*m_SndVols->size()*sizeof(int);
    cnt += 2*m_RcvVols->size()*sizeof(int);

    return cnt;
}

void
FabArrayBase::CopyComTag::GrokAsyncSends (const MapOfCopyComTagContainers& m_SndTags,
                                          Array<MPI_Request>&              send_reqs,
                                          Array<double*>&                  send_data,
                                          Array<MPI_Status>&               stats)
{
#ifdef BL_USE_MPI
    BL_ASSERT(FabArrayBase::do_async_sends && !m_SndTags.empty());

    const int N_snds = m_SndTags.size();

    stats.resize(N_snds);

    BL_ASSERT(send_reqs.size() == N_snds);
    BL_ASSERT(send_data.size() == N_snds);

    BL_MPI_REQUIRE( MPI_Waitall(N_snds, send_reqs.dataPtr(), stats.dataPtr()) );

    for (int i = 0; i < N_snds; i++)
        BoxLib::The_Arena()->free(send_data[i]);
#endif /*BL_USE_MPI*/
}

void
FabArrayBase::CopyComTag::SetRecvTag (MapOfCopyComTagContainers& m_RcvTags,
                                      int                        src_owner,
                                      CopyComTag&                tag,
                                      std::map<int,int>&         m_RcvVols,
                                      int                        vol)
{
    m_RcvTags[src_owner].push_back(tag);

    std::map<int,int>::iterator vol_it = m_RcvVols.find(src_owner);

    if (vol_it != m_RcvVols.end())
    {
        vol_it->second += vol;
    }
    else
    {
        m_RcvVols[src_owner] = vol;
    }
}

void
FabArrayBase::CopyComTag::SetSendTag (MapOfCopyComTagContainers& m_SndTags,
                                      int                        dst_owner,
                                      CopyComTag&                tag,
                                      std::map<int,int>&         m_SndVols,
                                      int                        vol)
{
    m_SndTags[dst_owner].push_back(tag);

    std::map<int,int>::iterator vol_it = m_SndVols.find(dst_owner);

    if (vol_it != m_SndVols.end())
    {
        vol_it->second += vol;
    }
    else
    {
        m_SndVols[dst_owner] = vol;
    }
}

void
FabArrayBase::CopyComTag::PostRcvs (const MapOfCopyComTagContainers& m_RcvTags,
                                    const std::map<int,int>&         m_RcvVols,
                                    double*&                         the_recv_data,
                                    Array<double*>&                  recv_data,
                                    Array<int>&                      recv_from,
                                    Array<MPI_Request>&              recv_reqs,
                                    int                              ncomp,
                                    int                              SeqNum)
{
    int TotalRcvsVolume = 0;

    for (std::map<int,int>::const_iterator it = m_RcvVols.begin(),
             End = m_RcvVols.end();
         it != End;
         ++it)
    {
        TotalRcvsVolume += it->second;
    }

    TotalRcvsVolume *= ncomp;

    BL_ASSERT((TotalRcvsVolume*sizeof(double)) < std::numeric_limits<int>::max());

    the_recv_data = static_cast<double*>(BoxLib::The_Arena()->alloc(TotalRcvsVolume*sizeof(double)));

    int Offset = 0;

    for (MapOfCopyComTagContainers::const_iterator m_it = m_RcvTags.begin(),
             m_End = m_RcvTags.end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = m_RcvVols.find(m_it->first);

        BL_ASSERT(vol_it != m_RcvVols.end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        recv_data.push_back(&the_recv_data[Offset]);
        recv_from.push_back(m_it->first);
        recv_reqs.push_back(ParallelDescriptor::Arecv(recv_data.back(),N,m_it->first,SeqNum).req());

        Offset += N;
    }
}

//
// The copy() cache.
//
FabArrayBase::CPCCache FabArrayBase::m_TheCopyCache;

FabArrayBase::CPCCacheIter
FabArrayBase::TheCPC (const CPC& cpc)
{
    const int Key = cpc.m_dstba.size() + cpc.m_srcba.size();

    std::pair<CPCCacheIter,CPCCacheIter> er_it = m_TheCopyCache.equal_range(Key);

    for (CPCCacheIter it = er_it.first; it != er_it.second; ++it)
    {
        if (it->second == cpc)
        {
            it->second.m_reused = true;

            return it;
        }
    }

    if (m_TheCopyCache.size() >= copy_cache_max_size)
    {
        //
        // Don't let the size of the cache get too big.
        //
        for (CPCCacheIter it = m_TheCopyCache.begin(); it != m_TheCopyCache.end(); )
        {
            if (!it->second.m_reused)
            {
                m_TheCopyCache.erase(it++);

                if (m_TheCopyCache.size() < copy_cache_max_size)
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

        if (m_TheCopyCache.size() >= copy_cache_max_size && !m_TheCopyCache.empty())
            //
            // Get rid of first entry which is the one with the smallest key.
            //
            m_TheCopyCache.erase(m_TheCopyCache.begin());
    }
    //
    // Got to insert one & then build it.
    //
    CPCCacheIter cache_it = m_TheCopyCache.insert(std::make_pair(Key,cpc));
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

    CopyComTag                        tag;
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

            tag.box = bx;

            if (dst_owner == MyProc)
            {
                tag.fabIndex = i;

                if (src_owner == MyProc)
                {
                    tag.srcIndex = k;

                    TheCPC.m_LocTags->push_back(tag);
                }
                else
                {
                    const int vol = bx.numPts();

                    FabArrayBase::CopyComTag::SetRecvTag(*TheCPC.m_RcvTags,src_owner,tag,*TheCPC.m_RcvVols,vol);
                }
            }
            else if (src_owner == MyProc)
            {
                tag.fabIndex = k;

                const int vol = bx.numPts();

                FabArrayBase::CopyComTag::SetSendTag(*TheCPC.m_SndTags,dst_owner,tag,*TheCPC.m_SndVols,vol);
            }
        }
    }

    if (TheCPC.m_LocTags->empty() && TheCPC.m_SndTags->empty() && TheCPC.m_RcvTags->empty())
    {
        //
        // This MPI proc has no work to do.  Don't store in the cache.
        //
        m_TheCopyCache.erase(cache_it);

        return m_TheCopyCache.end();
    }

    return cache_it;
}

void
FabArrayBase::CPC::FlushCache ()
{
    int stats[3] = {0,0,0}; // size, reused, bytes

    stats[0] = m_TheCopyCache.size();

    for (CPCCacheIter it = m_TheCopyCache.begin(), End = m_TheCopyCache.end();
         it != End;
         ++it)
    {
        stats[2] += it->second.bytes();
        if (it->second.m_reused)
            stats[1]++;
    }

    if (FabArrayBase::verbose)
    {
        ParallelDescriptor::ReduceIntMax(&stats[0], 3, ParallelDescriptor::IOProcessorNumber());

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

int
FabArrayBase::SI::bytes () const
{
    //
    // Get a estimate of number of bytes used by an SI.
    // This doesn't count any "overhead" in the STL containers.
    //
    int cnt = 0;

    cnt += m_LocTags->size()*sizeof(CopyComTag);

    for (CPC::MapOfCopyComTagContainers::const_iterator it = m_SndTags->begin(),
             m_End = m_SndTags->end();
         it != m_End;
         ++it)
    {
        cnt += it->second.size()*sizeof(CopyComTag);
    }

    for (CPC::MapOfCopyComTagContainers::const_iterator it = m_RcvTags->begin(),
             m_End = m_RcvTags->end();
         it != m_End;
         ++it)
    {
        cnt += it->second.size()*sizeof(CopyComTag);
    }

    cnt += 2*m_SndVols->size()*sizeof(int);
    cnt += 2*m_RcvVols->size()*sizeof(int);

    return cnt;
}

FabArrayBase::FBCache FabArrayBase::m_TheFBCache;

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
    int stats[3] = {0,0,0}; // size, reused, bytes

    stats[0] = m_TheFBCache.size();

    for (FBCacheIter it = m_TheFBCache.begin(), End = m_TheFBCache.end();
         it != End;
         ++it)
    {
        stats[2] += it->second.bytes();
        if (it->second.m_reused)
            stats[1]++;
    }

    if (FabArrayBase::verbose)
    {
        ParallelDescriptor::ReduceIntMax(&stats[0], 3, ParallelDescriptor::IOProcessorNumber());

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

FabArrayBase::FBCacheIter
FabArrayBase::TheFB (bool                cross,
                     const FabArrayBase& mf)
{
    const FabArrayBase::SI si(mf.boxArray(), mf.DistributionMap(), mf.nGrow(), cross);

    const int Key = mf.size() + mf.nGrow() + cross;

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
        //
        for (FBCacheIter it = m_TheFBCache.begin(); it != m_TheFBCache.end(); )
        {
            if (!it->second.m_reused)
            {
                m_TheFBCache.erase(it++);

                if (m_TheFBCache.size() < fb_cache_max_size)
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

        if (m_TheFBCache.size() >= fb_cache_max_size && !m_TheFBCache.empty())
            //
            // Get rid of first entry which is the one with the smallest key.
            //
            m_TheFBCache.erase(m_TheFBCache.begin());
    }
    //
    // Got to insert one & then build it.
    //
    FBCacheIter                cache_it = m_TheFBCache.insert(std::make_pair(Key,si));
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

    CopyComTag                        tag;
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
                const Box& bx        = isects[j].second;
                const int  k         = isects[j].first;
                const int  src_owner = dm[k];

                if (dst_owner != MyProc && src_owner != MyProc) continue;

                if (k == i) continue;

                tag.box = bx;

                if (dst_owner == MyProc)
                {
                    tag.fabIndex = i;

                    if (src_owner == MyProc)
                    {
                        tag.srcIndex = k;

                        TheFB.m_LocTags->push_back(tag);
                    }
                    else
                    {
                        const int vol = bx.numPts();

                        FabArrayBase::CopyComTag::SetRecvTag(*TheFB.m_RcvTags,src_owner,tag,*TheFB.m_RcvVols,vol);
                    }
                }
                else if (src_owner == MyProc)
                {
                    tag.fabIndex = k;

                    const int vol = bx.numPts();

                    FabArrayBase::CopyComTag::SetSendTag(*TheFB.m_SndTags,dst_owner,tag,*TheFB.m_SndVols,vol);
                }
            }
        }
    }

    if (TheFB.m_LocTags->empty() && TheFB.m_SndTags->empty() && TheFB.m_RcvTags->empty())
    {
        //
        // This MPI proc has no work to do.  Don't store in the cache.
        //
        m_TheFBCache.erase(cache_it);

        return m_TheFBCache.end();
    }

    return cache_it;
}
