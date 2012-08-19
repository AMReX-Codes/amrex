#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <iosfwd>
#include <map>

#include <CoordSys.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <RealBox.H>

//
// Rectangular problem domain geometry.
//
// This class describes problem domain and coordinate system for
// RECTANGULAR problem domains.  Since the problem domain is RECTANGULAR,
// periodicity is meaningful.
//

class Geometry
    :
    public CoordSys
{
public:
    //
    // The default constructor.
    //
    Geometry ();
    //
    // Constructor taking the rectangular domain.
    //
    Geometry (const Box&     dom,
              const RealBox* rb     =  0,
              int            coord  = -1,
              int*           is_per =  0);
    //
    // The copy constructor.
    //
    Geometry (const Geometry& g);
    //
    // The destructor.
    //
    ~Geometry();
    
    static void Finalize ();
    //
    // Read static values from ParmParse database.
    //
    static void Setup (const RealBox* rb = 0, int coord = -1, int* is_per = 0);
    //
    // Set the rectangular domain after using default constructor.
    //
    void define (const Box& dom, const RealBox* rb = 0, int coord = -1, int* is_per = 0);
    //
    // Returns the problem domain.
    //
    static const RealBox& ProbDomain ();
    //
    // Sets the problem domain.
    //
    static void ProbDomain (const RealBox& rb);
    //
    // Returns the lo end of the problem domain in each dimension.
    //
    static const Real* ProbLo ();
    //
    // Returns the hi end of the problem domain in each dimension.
    //
    static const Real* ProbHi ();
    //
    // Returns the lo end of the problem domain in specified direction.
    //
    static Real ProbLo (int dir);
    //
    // Returns the hi end of the problem domain in specified direction.
    //
    static Real ProbHi (int dir);
    //
    // Returns the overall size of the domain by multiplying the ProbLength's together
    //
    static Real ProbSize ();
    //
    // Returns length of problem domain in specified dimension.
    //
    static Real ProbLength (int dir);
    //
    // Returns our rectangular domain.
    //
    const Box& Domain () const;
    //
    // Sets our rectangular domain.
    //
    void Domain (const Box& bx);
    //
    // Define a multifab of areas and volumes with given grow factor.
    //
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    int             grow) const;

    void GetVolume (FArrayBox&       vol,
                    const BoxArray& grds,
                    int             idx,
                    int             grow) const;
    //
    // Compute d(log(A))/dr at cell centers in given region and
    //           stuff the results into the passed MultiFab.
    //
    void GetDLogA (MultiFab&       dloga,
                   const BoxArray& grds,
                   int             dir,
                   int             grow) const;
    //
    // Compute area of cell faces in given region and stuff
    // stuff the results into the passed MultiFab.
    //
    void GetFaceArea (MultiFab&       area,
                      const BoxArray& grds,
                      int             dir,
                      int             grow) const;

    void GetFaceArea (FArrayBox&      area,
                      const BoxArray& grds,
                      int             idx,
                      int             dir,
                      int             grow) const;
    //
    // Is the domain periodic in the specified direction?
    //
    static bool isPeriodic (int dir);
    //
    // Is domain periodic in any direction?
    //
    static bool isAnyPeriodic ();
    //
    // Is domain periodic in all directions?
    //
    static bool isAllPeriodic ();
    //
    // What's period in specified direction?
    //
    int period (int dir) const;
    //
    // Compute Array of shifts which will translate src so that it will
    // intersect target with non-zero intersection.  the array will be
    // resized internally, so anything previously there will be gone
    // DO NOT return non-periodic shifts, even if the box's do
    // intersect without shifting.  The logic is that you will only do
    // this as a special case if there is some periodicity.
    //
    void periodicShift (const Box&      target,
                        const Box&      src, 
                        Array<IntVect>& out) const;
    //
    // Fill ghost cells of all components with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               bool      do_corners = false,
                               bool      local      = false) const;
    //
    // Fill ghost cells of selected components with periodic data.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               int       src_comp,
                               int       num_comp,
                               bool      do_corners = false,
                               bool      local      = false) const;
    //
    // Sums the values in ghost cells, that can be shifted periodically
    // into valid region, into the corresponding cells in the valid
    // region.  The first routine here does all components while the latter
    // does the specified components.
    //
    void SumPeriodicBoundary (MultiFab& mf) const;

    void SumPeriodicBoundary (MultiFab& mf,
                              int       src_comp,
                              int       num_comp) const;

    void SumPeriodicBoundary (MultiFab&       dstmf,
                              const MultiFab& srcmf) const;

    void SumPeriodicBoundary (MultiFab&       dstmf,
                              const MultiFab& srcmf,
                              int             dcomp,
                              int             scomp,
                              int             ncomp) const;
    //
    // Flush the cache of PIRM information.
    //
    static void FlushPIRMCache ();
    //
    // The size of the PIRM cache.
    //
    static int PIRMCacheSize ();
    //
    // Used by SumPeriodicBoundary().
    //
    struct PIRec
    {
        PIRec ()
            :
            mfid(-1), srcId(-1) {}

        PIRec (int _mfid, int _srcId, const Box& _srcBox, const Box& _dstBox)
            :
            mfid(_mfid), srcId(_srcId), srcBox(_srcBox), dstBox(_dstBox) {}

        int       mfid;
        int       srcId;
        Box       srcBox;
        Box       dstBox;
        FillBoxId fbid;
    };
    //
    // A handy typedef.
    //
    typedef std::vector<PIRec> PIRMVector;
    //
    // Used by FillPeriodicBoundary().
    //
    struct FPBComTag
    {
        Box sbox;
        Box dbox;
        int srcIndex;
        int dstIndex;
    };
    //
    // Used in caching FillPeriodicBoundary().
    //
    struct FPB
    {
        FPB ()
            :
            m_ngrow(-1),
            m_do_corners(false),
            m_reused(false),
            m_LocTags(0),
            m_SndTags(0),
            m_RcvTags(0),
            m_SndVols(0),
            m_RcvVols(0)
            {}

        FPB (const BoxArray&            ba,
             const DistributionMapping& dm,
             const Box&                 domain,
             int                        ngrow,
             bool                       do_corners)
            :
            m_ba(ba),
            m_dm(dm),
            m_domain(domain),
            m_ngrow(ngrow),
            m_do_corners(do_corners),
            m_reused(false),
            m_LocTags(0),
            m_SndTags(0),
            m_RcvTags(0),
            m_SndVols(0),
            m_RcvVols(0)
            {
                BL_ASSERT(ngrow >= 0);
                BL_ASSERT(domain.ok());
            }

        ~FPB ();

        bool operator== (const FPB& rhs) const;

        bool operator!= (const FPB& rhs) const { return !operator==(rhs); }

        int bytes () const;

        BoxArray            m_ba;
        DistributionMapping m_dm;
        Box                 m_domain;
        int                 m_ngrow;
        bool                m_do_corners;
        bool                m_reused;
        //
        // Some useful typedefs.
        //
        typedef std::deque<FPBComTag> FPBComTagsContainer;

        typedef std::map<int,FPBComTagsContainer> MapOfFPBComTagContainers;
        //
        // The cache of local and send/recv info per FabArray::copy().
        //
        FPBComTagsContainer*      m_LocTags;
        MapOfFPBComTagContainers* m_SndTags;
        MapOfFPBComTagContainers* m_RcvTags;
        std::map<int,int>*        m_SndVols;
        std::map<int,int>*        m_RcvVols;
    };
    //
    // Some useful typedefs for the FPB cache.
    //
    typedef std::multimap<int,Geometry::FPB> FPBMMap;

    typedef FPBMMap::iterator FPBMMapIter;

    static FPBMMap m_FPBCache;

    static int fpb_cache_max_size;
    //
    // See if we're got an approprite FPB cached.
    //
    static Geometry::FPBMMapIter GetFPB (const Geometry&      geom,
                                         const Geometry::FPB& fpb,
                                         int                  scomp,
                                         int                  ncomp);
private:
    //
    // Helper functions.
    //
    void read_params ();
    //
    // Static data.
    //
    static int     spherical_origin_fix;
    static bool    is_periodic[BL_SPACEDIM]; // 0 means not periodic
    static RealBox prob_domain;
    //
    // Non-static data.
    //
    Box domain;
};
//
// Nice ASCII output.
//
std::ostream& operator<< (std::ostream&, const Geometry&);
//
// Nice ASCII input.
//
std::istream& operator>> (std::istream&, Geometry&);

inline
const RealBox&
Geometry::ProbDomain ()
{
    return prob_domain;
}

inline
void
Geometry::ProbDomain (const RealBox& rb)
{
    prob_domain = rb;
}

inline
const Box&
Geometry::Domain () const
{
    return domain;
}

inline
void
Geometry::Domain (const Box& bx)
{
    domain = bx;
}

inline
bool
Geometry::isPeriodic (int dir)
{
    return is_periodic[dir] != 0;
}

inline
bool
Geometry::isAnyPeriodic ()
{
    return isPeriodic(0)
#if BL_SPACEDIM>1
        ||   isPeriodic(1)
#endif
#if BL_SPACEDIM>2
        ||   isPeriodic(2)
#endif
        ;
}

inline
bool
Geometry::isAllPeriodic ()
{
#if BL_SPACEDIM==1
    return isPeriodic(0);
#elif BL_SPACEDIM==2
    return (isPeriodic(0) && isPeriodic(1));
#elif BL_SPACEDIM==3
    return (isPeriodic(0) && isPeriodic(1) && isPeriodic(2));
#endif
}

inline
int
Geometry::period (int dir) const
{
    BL_ASSERT(is_periodic[dir]);
    return domain.length(dir);
}

inline
const Real*
Geometry::ProbLo ()
{
    return prob_domain.lo();
}

inline
const Real*
Geometry::ProbHi ()
{
    return prob_domain.hi();
}

inline
Real
Geometry::ProbLo (int dir)
{
    return prob_domain.lo(dir);
}

inline
Real
Geometry::ProbHi (int dir)
{
    return prob_domain.hi(dir);
}

inline
Real
Geometry::ProbLength (int dir)
{
    return prob_domain.length(dir);
}

inline
Real
Geometry::ProbSize ()
{
#if (BL_SPACEDIM == 1) 
    return prob_domain.length(0);
#elif (BL_SPACEDIM == 2) 
    return prob_domain.length(0)*prob_domain.length(1);
#elif (BL_SPACEDIM == 3) 
    return prob_domain.length(0)*prob_domain.length(1)*prob_domain.length(2);
#endif
}

namespace BoxLib
{
    template <class FAB>
    void
    FillPeriodicBoundary (const Geometry& geom,
                          FabArray<FAB>&  mf,
                          int             scomp,
                          int             ncomp,
                          bool            corners=false)
    {
        if (!geom.isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

        Box TheDomain = geom.Domain();
        for (int n = 0; n < BL_SPACEDIM; n++)
            if (mf.boxArray()[0].ixType()[n] == IndexType::NODE)
                TheDomain.surroundingNodes(n);

        if (TheDomain.contains(BoxLib::grow(mf.boxArray().minimalBox(), mf.nGrow()))) return;

#ifndef NDEBUG
        //
        // Don't let folks ask for more grow cells than they have valid region.
        //
        for (int n = 0; n < BL_SPACEDIM; n++)
            if (geom.isPeriodic(n))
                BL_ASSERT(mf.nGrow() <= geom.Domain().length(n));
#endif
        const Geometry::FPB fpb(mf.boxArray(),mf.DistributionMap(),geom.Domain(),mf.nGrow(),corners);

        Geometry::FPBMMapIter cache_it = Geometry::GetFPB(geom,fpb,scomp,ncomp);

        if (cache_it != Geometry::m_FPBCache.end())
        {
            Geometry::FPB& TheFPB = cache_it->second;
            //
            // Do local work.
            //
            for (Geometry::FPB::FPBComTagsContainer::const_iterator it = TheFPB.m_LocTags->begin(),
                     End = TheFPB.m_LocTags->end();
                 it != End;
                 ++it)
            {
                const Geometry::FPBComTag& tag = *it;

                BL_ASSERT(mf.DistributionMap()[tag.dstIndex] == ParallelDescriptor::MyProc());
                BL_ASSERT(mf.DistributionMap()[tag.srcIndex] == ParallelDescriptor::MyProc());

                mf[tag.dstIndex].copy(mf[tag.srcIndex],tag.sbox,scomp,tag.dbox,scomp,ncomp);
            }
        }
    
        if (ParallelDescriptor::NProcs() == 1) return;

#ifdef BL_USE_MPI
        //
        // Do this before prematurely exiting if running in parallel.
        // Otherwise sequence numbers will not match across MPI processes.
        //
        const int SeqNum = ParallelDescriptor::SeqNum();

        if (cache_it == Geometry::m_FPBCache.end())
            //
            // No parallel work to do.
            //
            return;

        Geometry::FPB& TheFPB = cache_it->second;

        typedef typename FAB::value_type value_type;

        Array<MPI_Status>  stats;
        Array<int>         recv_from, index;
        Array<value_type*> recv_data, send_data;
        Array<MPI_Request> recv_reqs, send_reqs;
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        int total_rcvs_volume = 0;

        for (std::map<int,int>::const_iterator it = TheFPB.m_RcvVols->begin(),
                 End = TheFPB.m_RcvVols->end();
             it != End;
             ++it)
        {
            total_rcvs_volume += it->second;
        }

        total_rcvs_volume *= ncomp;

        BL_ASSERT((total_rcvs_volume*sizeof(value_type)) < std::numeric_limits<int>::max());

        value_type* the_recv_data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(total_rcvs_volume*sizeof(value_type)));



        int offset = 0;

        for (Geometry::FPB::MapOfFPBComTagContainers::const_iterator m_it = TheFPB.m_RcvTags->begin(),
                 m_End = TheFPB.m_RcvTags->end();
             m_it != m_End;
             ++m_it)
        {
            const int N = (*TheFPB.m_RcvVols)[m_it->first]*ncomp;

            BL_ASSERT(N < std::numeric_limits<int>::max());

            recv_data.push_back(&the_recv_data[offset]);
            recv_from.push_back(m_it->first);
            recv_reqs.push_back(ParallelDescriptor::Arecv(recv_data.back(),N,m_it->first,SeqNum).req());

            offset += N;
        }
        //
        // Send the FAB data.
        //
        FAB fab;

        for (Geometry::FPB::MapOfFPBComTagContainers::const_iterator m_it = TheFPB.m_SndTags->begin(),
                 m_End = TheFPB.m_SndTags->end();
             m_it != m_End;
             ++m_it)
        {
            const int N = (*TheFPB.m_SndVols)[m_it->first]*ncomp;

            BL_ASSERT(N < std::numeric_limits<int>::max());

            value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
            value_type* dptr = data;

            for (Geometry::FPB::FPBComTagsContainer::const_iterator it = m_it->second.begin(),
                     End = m_it->second.end();
                 it != End;
                 ++it)
            {
                const Box& bx = it->sbox;
                fab.resize(bx, ncomp);
                fab.copy(mf[it->srcIndex],bx,scomp,bx,0,ncomp);
                const int Cnt = bx.numPts()*ncomp;
                memcpy(dptr, fab.dataPtr(), Cnt*sizeof(value_type));
                dptr += Cnt;
            }
            BL_ASSERT(data+N == dptr);

            if (FabArrayBase::do_async_sends)
            {
                send_data.push_back(data);
                send_reqs.push_back(ParallelDescriptor::Asend(data,N,m_it->first,SeqNum).req());
            }
            else
            {
                ParallelDescriptor::Send(data,N,m_it->first,SeqNum);
                BoxLib::The_Arena()->free(data);
            }
        }
        //
        // Now receive and unpack FAB data as it becomes available.
        //
        Geometry::FPB::MapOfFPBComTagContainers::const_iterator m_it;

        const int N_rcvs = TheFPB.m_RcvTags->size();

        index.resize(N_rcvs);
        stats.resize(N_rcvs);

        for (int NWaits = N_rcvs, completed; NWaits > 0; NWaits -= completed)
        {
            ParallelDescriptor::Waitsome(recv_reqs, completed, index, stats);

            for (int k = 0; k < completed; k++)
            {
                const value_type* dptr = recv_data[index[k]];

                BL_ASSERT(dptr != 0);

                m_it = TheFPB.m_RcvTags->find(recv_from[index[k]]);

                BL_ASSERT(m_it != TheFPB.m_RcvTags->end());

                for (Geometry::FPB::FPBComTagsContainer::const_iterator it = m_it->second.begin(),
                         End = m_it->second.end();
                     it != End;
                     ++it)
                {
                    const Box& bx = it->dbox;
                    fab.resize(bx,ncomp);
                    const int Cnt = bx.numPts()*ncomp;
                    memcpy(fab.dataPtr(), dptr, Cnt*sizeof(value_type));
                    mf[it->dstIndex].copy(fab,bx,0,bx,scomp,ncomp);
                    dptr += Cnt;
                }
            }
        }

        BoxLib::The_Arena()->free(the_recv_data);

        if (FabArrayBase::do_async_sends && !TheFPB.m_SndTags->empty())
        {
            //
            // Now grok the asynchronous send buffers & free up send buffer space.
            //
            const int N_snds = TheFPB.m_SndTags->size();

            stats.resize(N_snds);

            BL_MPI_REQUIRE( MPI_Waitall(N_snds, send_reqs.dataPtr(), stats.dataPtr()) );

            for (int i = 0; i < N_snds; i++)
                BoxLib::The_Arena()->free(send_data[i]);
        }

#endif /*BL_USE_MPI*/
    }
}

#endif /*_GEOMETRY_H_*/
