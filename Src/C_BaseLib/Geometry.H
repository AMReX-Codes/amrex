#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <iosfwd>
#include <list>
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
    // Returns length of problem domain in each dimension.
    //
    static const Real* ProbLength ();
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
    // Flush the cache of PIRM information.
    //
    static void FlushPIRMCache ();
    //
    // The size of the PIRM cache.
    //
    static int PIRMCacheSize ();
    //
    // Src/dest box pairs to copy for FillPeriodicBoundary().
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
    typedef std::vector<PIRec> PIRMMap;
    //
    // Used in caching PIRMMaps && FabArrayBase::CommData.
    //
    struct FPB
    {
        FPB ();

        FPB (const BoxArray&            ba,
             const DistributionMapping& dm,
             const Box&                 domain,
             int                        ngrow,
             bool                       do_corners);

        FPB (const FPB& rhs);

        ~FPB ();

        bool operator== (const FPB& rhs) const;
        bool operator!= (const FPB& rhs) const;

        Array<int>                  m_cache;    // Snds cached for CollectData().
        FabArrayBase::CommDataCache m_commdata; // Yet another cache for CollectData().
        PIRMMap                     m_pirm;
        BoxArray                    m_ba;
        DistributionMapping         m_dm;
        Box                         m_domain;
        int                         m_ngrow;
        bool                        m_do_corners;
        bool                        m_reused;
    };

    typedef std::multimap<int,Geometry::FPB> FPBMMap;

    typedef FPBMMap::iterator FPBMMapIter;

    static FPBMMap m_FPBCache;

    static int fpb_cache_max_size;

protected:
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

namespace BoxLib
{
    //
    // Template functions used in implementing FillPeriodicBoundary().
    //
    template <class FAB>
    Geometry::FPB&
    BuildFPB (FabArray<FAB>&       mf,
              const Geometry&      geom,
              const Geometry::FPB& fpb)
    {
        const int key = mf.nGrow() + mf.size();

        Geometry::FPBMMapIter it = Geometry::m_FPBCache.insert(std::make_pair(key,fpb));

        Geometry::PIRMMap& pirm = it->second.m_pirm;

        Array<IntVect> pshifts(27);

        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            const Box& dest = mf[mfi].box();

            BL_ASSERT(dest == BoxLib::grow(mfi.validbox(), mf.nGrow()));

            Box TheDomain = geom.Domain();
            for (int n = 0; n < BL_SPACEDIM; n++)
                if (dest.ixType()[n] == IndexType::NODE)
                    TheDomain.surroundingNodes(n);

            if (!TheDomain.contains(dest))
            {
                const BoxArray& grids = mf.boxArray();

                for (int j = 0; j < grids.size(); j++)
                {
                    Box src = grids[j] & TheDomain;

                    if (fpb.m_do_corners)
                    {
                        for (int i = 0; i < BL_SPACEDIM; i++)
                        {
                            if (!geom.isPeriodic(i))
                            {
                                if (src.smallEnd(i) == geom.Domain().smallEnd(i))
                                    src.growLo(i,mf.nGrow());
                                if (src.bigEnd(i) == geom.Domain().bigEnd(i))
                                    src.growHi(i,mf.nGrow());
                            }
                        }
                    }

                    geom.periodicShift(dest, src, pshifts);

                    for (int i = 0; i < pshifts.size(); i++)
                    {
                        Box shftbox = src + pshifts[i];
                        Box dbx     = dest & shftbox;
                        Box sbx     = dbx - pshifts[i];

                        pirm.push_back(Geometry::PIRec(mfi.index(),j,sbx,dbx));
                    }
                }
            }
        }

        return it->second;
    }

    template <class FAB>
    Geometry::FPB&
    GetFPB (FabArray<FAB>&       mf,
            const Geometry&      geom,
            const Geometry::FPB& fpb,
            int                  scomp,
            int                  ncomp)
    {
        const int key = mf.nGrow() + mf.size();

        std::pair<Geometry::FPBMMapIter,Geometry::FPBMMapIter> er_it = Geometry::m_FPBCache.equal_range(key);
    
        for (Geometry::FPBMMapIter it = er_it.first; it != er_it.second; ++it)
        {
            if (it->second == fpb)
            {
                it->second.m_reused = true;
                //
                // Adjust the ncomp & scomp in CommData.
                //
                Array<ParallelDescriptor::CommData>& cd = it->second.m_commdata.theCommData();

                for (int i = 0; i < cd.size(); i++)
                {
                    cd[i].nComp(ncomp);
                    cd[i].srcComp(scomp);
                }

                return it->second;
            }
        }

        if (Geometry::m_FPBCache.size() >= Geometry::fpb_cache_max_size && Geometry::fpb_cache_max_size != -1)
        {
            //
            // Don't let the size of the cache get too big.
            //
            for (Geometry::FPBMMapIter it = Geometry::m_FPBCache.begin(); it != Geometry::m_FPBCache.end(); )
            {
                if (!it->second.m_reused)
                {
                    Geometry::m_FPBCache.erase(it++);
                    //
                    // Only delete enough entries to stay under limit.
                    //
                    if (Geometry::m_FPBCache.size() < Geometry::fpb_cache_max_size) break;
                }
                else
                {
                    ++it;
                }
            }

            if (Geometry::m_FPBCache.size() >= Geometry::fpb_cache_max_size)
                //
                // Get rid of first entry which is the one with the smallest key.
                //
                Geometry::m_FPBCache.erase(Geometry::m_FPBCache.begin());
        }

        return BoxLib::BuildFPB(mf,geom,fpb);
    }

    template <class FAB>
    void
    FillPeriodicBoundary (const Geometry& geom,
                          FabArray<FAB>&  mf,
                          int             scomp,
                          int             ncomp,
                          bool            corners=false)
    {
        if (!geom.isAnyPeriodic()) return;

        FabArrayCopyDescriptor<FAB> mfcd;

        Geometry::FPB TheFPB(mf.boxArray(),mf.DistributionMap(),geom.Domain(),mf.nGrow(),corners);

        const FabArrayId   mfid = mfcd.RegisterFabArray(&mf);
        Geometry::FPB&     fpb  = BoxLib::GetFPB(mf,geom,TheFPB,scomp,ncomp);
        Geometry::PIRMMap& pirm = fpb.m_pirm;
        //
        // Add boxes we need to collect.
        //
        for (int i = 0; i < pirm.size(); i++)
        {
            pirm[i].fbid = mfcd.AddBox(mfid,
                                       pirm[i].srcBox,
                                       0,
                                       pirm[i].srcId,
                                       scomp,
                                       scomp,
                                       ncomp,
                                       !corners);
        }

        mfcd.CollectData(&fpb.m_cache,&fpb.m_commdata);

        for (int i = 0; i < pirm.size(); i++)
        {
            BL_ASSERT(pirm[i].fbid.box() == pirm[i].srcBox);
            BL_ASSERT(pirm[i].srcBox.sameSize(pirm[i].dstBox));
            BL_ASSERT(mf.DistributionMap()[pirm[i].mfid] == ParallelDescriptor::MyProc());

            mfcd.FillFab(mfid, pirm[i].fbid, mf[pirm[i].mfid], pirm[i].dstBox);
        }
    }
}

#endif /*_GEOMETRY_H_*/
