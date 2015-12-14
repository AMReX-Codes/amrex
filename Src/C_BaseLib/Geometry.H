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
    static const RealBox& ProbDomain () { return prob_domain; }
    //
    // Sets the problem domain.
    //
    static void ProbDomain (const RealBox& rb) { prob_domain = rb; }
    //
    // Returns the lo end of the problem domain in each dimension.
    //
    static const Real* ProbLo ()  { return prob_domain.lo(); }
    //
    // Returns the hi end of the problem domain in each dimension.
    //
    static const Real* ProbHi () { return prob_domain.hi(); }
    //
    // Returns the lo end of the problem domain in specified direction.
    //
    static Real ProbLo (int dir) { return prob_domain.lo(dir); }
    //
    // Returns the hi end of the problem domain in specified direction.
    //
    static Real ProbHi (int dir) { return prob_domain.hi(dir); }
    //
    // Returns the overall size of the domain by multiplying the ProbLength's together
    //
    static Real ProbSize ()
    {
        return D_TERM(prob_domain.length(0),*prob_domain.length(1),*prob_domain.length(2));
    }
    //
    // Returns length of problem domain in specified dimension.
    //
    static Real ProbLength (int dir) { return prob_domain.length(dir); }
    //
    // Returns our rectangular domain.
    //
    const Box& Domain () const { return domain; }
    //
    // Sets our rectangular domain.
    //
    void Domain (const Box& bx) { domain = bx; }
    //
    // Define a multifab of areas and volumes with given grow factor.
    //
    void GetVolume (MultiFab&       vol,
                    const BoxArray& grds,
                    int             grow) const;
    // Fill the pre-built multifab with volume
    void GetVolume (MultiFab&       vol) const;

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
    // Fill the pre-built multifab with area
    void GetFaceArea (MultiFab&       area,
                      int             dir) const;

    void GetFaceArea (FArrayBox&      area,
                      const BoxArray& grds,
                      int             idx,
                      int             dir,
                      int             grow) const;
    //
    // Is the domain periodic in the specified direction?
    //
    static bool isPeriodic (int dir) { return is_periodic[dir] != 0; }
    //
    // Is domain periodic in any direction?
    //
    static bool isAnyPeriodic ()
    {
        return D_TERM(isPeriodic(0),||isPeriodic(1),||isPeriodic(2));
    }
    //
    // Is domain periodic in all directions?
    //
    static bool isAllPeriodic ()
    {
        return D_TERM(isPeriodic(0),&&isPeriodic(1),&&isPeriodic(2));
    }
    //
    // What's period in specified direction?
    //
    int period (int dir) const { BL_ASSERT(is_periodic[dir]); return domain.length(dir); }
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
        FPB ();

        FPB (const BoxArray&            ba,
             const DistributionMapping& dm,
             const Box&                 domain,
             int                        ngrow,
             bool                       do_corners);

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
	bool                m_threadsafe_loc;
	bool                m_threadsafe_rcv;
        //
        // Some useful typedefs.
        //
        typedef std::vector<FPBComTag> FPBComTagsContainer;

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
    // See if we've got an approprite FPB cached.
    //
    static Geometry::FPBMMapIter GetFPB (const Geometry&      geom,
                                         const Geometry::FPB& fpb,
                                         const FabArrayBase&  mf);
    //
    // Send Geometry to sidecar procs.
    //
#ifdef BL_USE_MPI
    static void SendGeometryToSidecars (Geometry *geom);
#endif
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
            if (mf.boxArray().ixType()[n] == IndexType::NODE)
                TheDomain.surroundingNodes(n);
#ifndef NDEBUG
        //
        // Don't let folks ask for more grow cells than they have valid region.
        //
        for (int n = 0; n < BL_SPACEDIM; n++)
            if (geom.isPeriodic(n))
                BL_ASSERT(mf.nGrow() <= geom.Domain().length(n));
#endif
        const Geometry::FPB fpb(mf.boxArray(),mf.DistributionMap(),geom.Domain(),mf.nGrow(),corners);

        Geometry::FPBMMapIter cache_it = Geometry::GetFPB(geom,fpb,mf);

        BL_ASSERT(cache_it != Geometry::m_FPBCache.end());

        const Geometry::FPB& TheFPB = cache_it->second;

        if (ParallelDescriptor::NProcs() == 1)
        {
            //
            // There can only be local work to do.
            //
	    int N_loc = (*TheFPB.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (TheFPB.m_threadsafe_loc)
#endif
	    for (int i=0; i<N_loc; ++i)
	    {
                const Geometry::FPBComTag& tag = (*TheFPB.m_LocTags)[i];
                mf[tag.dstIndex].copy(mf[tag.srcIndex],tag.sbox,scomp,tag.dbox,scomp,ncomp);
            }

            return;
        }

#ifdef BL_USE_MPI
        //
        // Do this before prematurely exiting if running in parallel.
        // Otherwise sequence numbers will not match across MPI processes.
        //
        const int SeqNum = ParallelDescriptor::SeqNum();

        if (TheFPB.m_LocTags->empty() && TheFPB.m_RcvTags->empty() && TheFPB.m_SndTags->empty())
            //
            // No work to do.
            //
            return;

        typedef typename FAB::value_type value_type;

        Array<MPI_Status>  stats;
        Array<int>         recv_from;
        Array<value_type*> recv_data;
        Array<MPI_Request> recv_reqs;
        //
        // Post rcvs. Allocate one chunk of space to hold'm all.
        //
        value_type* the_recv_data = 0;

        FabArrayBase::PostRcvs(*TheFPB.m_RcvTags,*TheFPB.m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);

        //
        // Post send's
        //
	const int N_snds = TheFPB.m_SndTags->size();

	Array<value_type*>                     send_data;
	Array<int>                             send_N;
	Array<int>                             send_rank;
	Array<const Geometry::FPB::FPBComTagsContainer*> send_fctc;

	send_data.reserve(N_snds);
	send_N   .reserve(N_snds);
	send_rank.reserve(N_snds);
	send_fctc.reserve(N_snds);

        for (Geometry::FPB::MapOfFPBComTagContainers::const_iterator m_it = TheFPB.m_SndTags->begin(),
                 m_End = TheFPB.m_SndTags->end();
             m_it != m_End;
             ++m_it)
        {
            std::map<int,int>::const_iterator vol_it = TheFPB.m_SndVols->find(m_it->first);

            BL_ASSERT(vol_it != TheFPB.m_SndVols->end());

            const int N = vol_it->second*ncomp;

            BL_ASSERT(N < std::numeric_limits<int>::max());

            value_type* data = static_cast<value_type*>(BoxLib::The_Arena()->alloc(N*sizeof(value_type)));

	    send_data.push_back(data);
	    send_N   .push_back(N);
	    send_rank.push_back(m_it->first);
	    send_fctc.push_back(&(m_it->second));
	}

#ifdef _OPENMP
#pragma omp parallel for
#endif
	for (int i=0; i<N_snds; ++i)
	{
	    value_type* dptr = send_data[i];
	    BL_ASSERT(dptr != 0);

	    const Geometry::FPB::FPBComTagsContainer& fctc = *send_fctc[i];

            for (Geometry::FPB::FPBComTagsContainer::const_iterator it = fctc.begin();
                 it != fctc.end(); ++it)
            {
                const Box& bx = it->sbox;
                mf[it->srcIndex].copyToMem(bx,scomp,ncomp,dptr);
                const int Cnt = bx.numPts()*ncomp;
                dptr += Cnt;
            }
	}

	Array<MPI_Request> send_reqs;

	if (FabArrayBase::do_async_sends)
	{
	    send_reqs.reserve(N_snds);
	    for (int i=0; i<N_snds; ++i) {
                send_reqs.push_back(ParallelDescriptor::Asend
				    (send_data[i],send_N[i],send_rank[i],SeqNum).req());
            }
	} else {
	    for (int i=0; i<N_snds; ++i) {
                ParallelDescriptor::Send(send_data[i],send_N[i],send_rank[i],SeqNum);
                BoxLib::The_Arena()->free(send_data[i]);
            }
        }

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
	int N_loc = (*TheFPB.m_LocTags).size();
#ifdef _OPENMP
#pragma omp parallel for if (TheFPB.m_threadsafe_loc)
#endif
	for (int i=0; i<N_loc; ++i)
	{
            const Geometry::FPBComTag& tag = (*TheFPB.m_LocTags)[i];

            BL_ASSERT(mf.DistributionMap()[tag.dstIndex] == ParallelDescriptor::MyProc());
            BL_ASSERT(mf.DistributionMap()[tag.srcIndex] == ParallelDescriptor::MyProc());

            mf[tag.dstIndex].copy(mf[tag.srcIndex],tag.sbox,scomp,tag.dbox,scomp,ncomp);
        }

	//
	// wait and unpack

        const int N_rcvs = TheFPB.m_RcvTags->size();

	if (N_rcvs > 0)
	{
	    Array<const Geometry::FPB::FPBComTagsContainer*> recv_fctc;
	    recv_fctc.reserve(N_rcvs);

	    for (int k = 0; k < N_rcvs; k++)
	    {
		Geometry::FPB::MapOfFPBComTagContainers::const_iterator m_it = TheFPB.m_RcvTags->find(recv_from[k]);
                BL_ASSERT(m_it != TheFPB.m_RcvTags->end());
		
		recv_fctc.push_back(&(m_it->second));
	    }

	    stats.resize(N_rcvs);
	    BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
	    
#ifdef _OPENMP
#pragma omp parallel for if (TheFPB.m_threadsafe_rcv)
#endif
	    for (int k = 0; k < N_rcvs; k++) 
	    {
		value_type*  dptr = recv_data[k];
		BL_ASSERT(dptr != 0);
		
		const Geometry::FPB::FPBComTagsContainer& fctc = *recv_fctc[k];

                for (Geometry::FPB::FPBComTagsContainer::const_iterator it = fctc.begin();
                     it != fctc.end(); ++it)
                {
                    const Box& bx  = it->dbox;
                    const int  Cnt = bx.numPts()*ncomp;
                    mf[it->dstIndex].copyFromMem(bx,scomp,ncomp,dptr);
                    dptr += Cnt;
                }
            }
        }

        BoxLib::The_Arena()->free(the_recv_data);

        if (FabArrayBase::do_async_sends && !TheFPB.m_SndTags->empty())
            FabArrayBase::GrokAsyncSends(TheFPB.m_SndTags->size(),send_reqs,send_data,stats);

#endif /*BL_USE_MPI*/
    }
}

#endif /*_GEOMETRY_H_*/
