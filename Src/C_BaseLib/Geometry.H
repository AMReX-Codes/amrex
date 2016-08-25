#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include <iosfwd>
#include <map>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef BL_USE_UPCXX
#include <BLPgas.H>
#endif

#include <CoordSys.H>
#include <MultiFab.H>
#include <ParallelDescriptor.H>
#include <RealBox.H>
#include <Periodicity.H>

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

    Periodicity periodicity () const {
	return Periodicity(IntVect(D_DECL(domain.length(0) * is_periodic[0],
					  domain.length(1) * is_periodic[1],
					  domain.length(2) * is_periodic[2])));
    }

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
    // do_corners: the corner here refers to physical DOMAIN corners
    // that are partially but not all periodic.  If do_corners=true,
    // these corners will be filled with periodically shifted GHOST
    // cell data, because no valid cell data are available for the
    // role.  Note that this does not refer to corner cells of boxes
    // that are not at DOMAIN corners.  Corner cells of a box will
    // always be filled if valid cell data are available periodically.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               bool      do_corners = false,
                               bool      local      = false) const;
    //
    // Fill ghost cells of selected components with periodic data.
    //
    // do_corners: see comments above.
    //
    void FillPeriodicBoundary (MultiFab& mf,
                               int       src_comp,
                               int       num_comp,
                               bool      do_corners = false,
                               bool      local      = false) const;
  
    void FillPeriodicBoundary_local (MultiFab& mf,
				     int       src_comp,
				     int       num_comp,
				     bool      do_corners = false) const;
    void FillPeriodicBoundary_nowait (MultiFab& mf,
				      bool      do_corners = false,
				      bool      local      = false) const;
    void FillPeriodicBoundary_nowait (MultiFab& mf,
				      int       src_comp,
				      int       num_comp,
				      bool      do_corners = false,
				      bool      local      = false) const;
    void FillPeriodicBoundary_finish (MultiFab& mf) const;
    // 
    // Copy periodically shifted data from srcmf to dstmf.
    //
    void PeriodicCopy (MultiFab&       dstmf,
		       const MultiFab& srcmf) const;
    void PeriodicCopy (MultiFab&       dstmf,
		       const MultiFab& srcmf,
		       int             dcomp,
		       int             scomp,
		       int             ncomp,
		       int             dstng = 0,
		       int             srcng = 0) const;
    //
    // Used by FillPeriodicBoundary().
    //
    struct FPBComTag
    {
        Box sbox;
        Box dbox;
        int srcIndex;
        int dstIndex;
	FPBComTag () {}
	FPBComTag (const Box& sb, const Box& db, int sidx, int didx)
	    : sbox(sb), dbox(db), srcIndex(sidx), dstIndex(didx) {}
	// FPBComTag needs to be sortable if it is used in remote communication.
	bool operator< (const FPBComTag& rhs) const {
	    return (srcIndex < rhs.srcIndex) || ((srcIndex == rhs.srcIndex) && (
		   (dstIndex < rhs.dstIndex) || ((dstIndex == rhs.dstIndex) && (
                   (IntVect::Compare()(dbox.smallEnd(),rhs.dbox.smallEnd()))))));
	}
    };
    //
    // Used in caching FillPeriodicBoundary().
    //
    struct FPB
    {
        FPB (const Geometry& geom, const FabArrayBase& fa, bool do_corners=false);
        ~FPB ();

        long bytes () const;

	IndexType m_typ;
	int       m_ngrow;
	bool      m_do_corners;
	//
	bool                m_threadsafe_loc;
	bool                m_threadsafe_rcv;
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
	//
	int m_nuse;

    private:
	long bytesOfMapOfComTagContainers (const Geometry::FPB::MapOfFPBComTagContainers& m) const;
    };
    //
    typedef std::multimap<FabArrayBase::BDKey,FPB*> FPBCache;
    typedef FPBCache::iterator FPBCacheIter;
    //
    static FPBCache                 m_TheFPBCache;
    static FabArrayBase::CacheStats m_FPBC_stats;
    //
    const FPB& getFPB (const FabArrayBase& fa, bool do_corners=false) const;
    //
    static void flushFPB (const FabArrayBase::BDKey& key); // Flushes FPBs related to this key.
    static void flushFPBCache ();                          // Flushes the entire cache.

#ifdef BL_USE_MPI
    //
    // Send Geometry to sidecar procs.
    //
    static void SendGeometryToSidecar (Geometry *geom, int whichSidecar);
    //
    static void BroadcastGeometry (Geometry &geom, int fromProc, MPI_Comm comm);
    static void BroadcastGeometry (Geometry &geom, int fromProc, MPI_Comm comm, bool bcastSource);
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
    FillPeriodicBoundary_nowait (const Geometry& geom,
				 FabArray<FAB>&  mf,
				 int             scomp,
				 int             ncomp,
				 bool            corners=false)
    {
        if (!geom.isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

	const DistributionMapping& dmap = mf.DistributionMap();

#ifndef NDEBUG
        //
        // Don't let folks ask for more grow cells than they have valid region.
        //
        for (int n = 0; n < BL_SPACEDIM; n++)
            if (geom.isPeriodic(n))
                BL_ASSERT(mf.nGrow() <= geom.Domain().length(n));
#endif

        const Geometry::FPB& TheFPB = geom.getFPB(mf,corners);

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

	mf.fpb_corners = corners;
	mf.fpb_scomp   = scomp;
	mf.fpb_ncomp   = ncomp;

#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
	ParallelDescriptor::Mode.set_upcxx_mode();
	ParallelDescriptor::Mode.incr_upcxx();
#endif

#if !defined(BL_USE_MPI3)
	BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif

        //
        // Do this before prematurely exiting if running in parallel.
        // Otherwise sequence numbers will not match across MPI processes.
        //
	int SeqNum;
	{
	    ParallelDescriptor::Color mycolor = mf.color();
	    if (mycolor == ParallelDescriptor::DefaultColor()) {
		SeqNum = ParallelDescriptor::SeqNum();
	    } else if (mycolor == ParallelDescriptor::SubCommColor()) {
		SeqNum = ParallelDescriptor::SubSeqNum();
	    }
	    // else I don't have any data and my SubSeqNum() should not be called.
	}

	const int N_locs = TheFPB.m_LocTags->size();
	const int N_rcvs = TheFPB.m_RcvTags->size();
	const int N_snds = TheFPB.m_SndTags->size();

        if (N_locs == 0 && N_rcvs == 0 && N_snds == 0)
            //
            // No work to do.
            //
            return;

        typedef typename FAB::value_type value_type;

	//
	// Post rcvs. Allocate one chunk of space to hold'm all.
	//

#ifdef BL_USE_MPI3
	MPI_Group tgroup, rgroup, sgroup;
	if (ParallelDescriptor::MPIOneSided())
	    MPI_Comm_group(MPI_COMM_WORLD, &tgroup);
#endif

	if (N_rcvs > 0) {
#ifdef BL_USE_UPCXX
	    FabArrayBase::PostRcvs_PGAS(*TheFPB.m_RcvVols,mf.fpb_the_recv_data,
					mf.fpb_recv_data,mf.fpb_recv_from,
					mf.fpb_ncomp,SeqNum,&BLPgas::fpb_recv_event);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		FabArrayBase::PostRcvs_MPI_Onesided(*TheFPB.m_RcvVols,mf.fpb_the_recv_data,
						    mf.fpb_recv_data,mf.fpb_recv_from,mf.fpb_recv_reqs,
						    mf.fpb_recv_disp,ncomp,SeqNum,
						    ParallelDescriptor::fpb_win);
		MPI_Group_incl(tgroup, mf.fpb_recv_from.size(), mf.fpb_recv_from.dataPtr(), &rgroup);
		MPI_Win_post(rgroup, 0, ParallelDescriptor::fpb_win);
#endif
            } else {
		FabArrayBase::PostRcvs(*TheFPB.m_RcvVols,mf.fpb_the_recv_data,
				       mf.fpb_recv_data,mf.fpb_recv_from,mf.fpb_recv_reqs,
				       mf.fpb_ncomp,SeqNum);
	    }
#endif
	}

        //
        // Post send's
        //
	if (N_snds > 0)
	{
            Array<value_type*>&                    send_data = mf.fpb_send_data;
	    Array<int>                             send_N;
	    Array<int>                             send_rank;
	    Array<const Geometry::FPB::FPBComTagsContainer*> send_fctc;

	    const int N_snds = TheFPB.m_SndTags->size();

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
		
		value_type* data = static_cast<value_type*>
#ifdef BL_USE_UPCXX
		    (BLPgas::alloc(N*sizeof(value_type)));
#else
		    (BoxLib::The_Arena()->alloc(N*sizeof(value_type)));
#endif
		
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

#ifdef BL_USE_UPCXX

	    BLPgas::fpb_send_counter = 0;

	    for (int i=0; i<N_snds; ++i) {
		BLPgas::Send(upcxx::global_ptr<void>((void *)send_data[i], upcxx::myrank()),
			     send_rank[i],
			     send_N[i]*sizeof(value_type),
			     SeqNum,
			     &BLPgas::fpb_send_event,
			     &BLPgas::fpb_send_counter);
	    }
	    
	    while (BLPgas::fpb_send_counter < N_snds)
		upcxx::advance();

#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		Array<MPI_Request> send_reqs;
		Array<MPI_Aint>    send_disp;
		send_reqs.reserve(N_snds);
		send_disp.resize (N_snds);
		
		for (int i=0; i<N_snds; ++i) {
		    send_reqs.push_back(ParallelDescriptor::Arecv
					(&send_disp[i],1,send_rank[i],SeqNum).req());
		}
		
                MPI_Group_incl(tgroup, send_rank.size(), send_rank.dataPtr(), &sgroup);
		MPI_Win_start(sgroup,0,ParallelDescriptor::fpb_win);

	        int send_counter = 0;	    
	        while (send_counter < N_snds) {
		    MPI_Status status;
		    int index;
		
		    MPI_Waitany(N_snds, send_reqs.dataPtr(), &index, &status);
		    BL_ASSERT(status.MPI_TAG == SeqNum);
		    BL_ASSERT(status.MPI_SOURCE == send_rank[index]);
		    MPI_Put(send_data[index], send_N[index]*sizeof(value_type), MPI_CHAR, send_rank[index],
			    send_disp[index], send_N[index]*sizeof(value_type), MPI_CHAR, 
			    ParallelDescriptor::fpb_win);
		    send_counter++;
	        }
#endif
	    } else {
		mf.fpb_send_reqs.reserve(N_snds);

		for (int i=0; i<N_snds; ++i) {
		    mf.fpb_send_reqs.push_back(ParallelDescriptor::Asend
					       (send_data[i], send_N[i], send_rank[i], SeqNum).req());
		}
	    }
#endif
	}

#ifdef BL_USE_MPI3
        if (ParallelDescriptor::MPIOneSided()) {
	    MPI_Group_free(&tgroup);
	    if (N_rcvs > 0) MPI_Group_free(&rgroup);
	    if (N_snds > 0) MPI_Group_free(&sgroup);
        }
#endif

        //
        // Do the local work.  Hope for a bit of communication/computation overlap.
        //
	int N_loc = (*TheFPB.m_LocTags).size();
	if (ParallelDescriptor::TeamSize() > 1 && TheFPB.m_threadsafe_loc)
	{
#ifdef BL_USE_TEAM
#ifdef _OPENMP
#pragma omp parallel
#endif
	    ParallelDescriptor::team_for(0, N_loc, [&] (int i) 
            {
                const auto& tag = (*TheFPB.m_LocTags)[i];

		BL_ASSERT(ParallelDescriptor::sameTeam(dmap[tag.dstIndex]));
		BL_ASSERT(ParallelDescriptor::sameTeam(dmap[tag.srcIndex]));

		mf[tag.dstIndex].copy(mf[tag.srcIndex],tag.sbox,scomp,tag.dbox,scomp,ncomp);
	    });
#endif
	}
	else
	{
#ifdef _OPENMP
#pragma omp parallel for if (TheFPB.m_threadsafe_loc)
#endif
	    for (int i=0; i<N_loc; ++i)
	    {
		const Geometry::FPBComTag& tag = (*TheFPB.m_LocTags)[i];

		BL_ASSERT(ParallelDescriptor::sameTeam(dmap[tag.dstIndex]));
		BL_ASSERT(ParallelDescriptor::sameTeam(dmap[tag.srcIndex]));

		if (dmap[tag.dstIndex] == ParallelDescriptor::MyProc()) {
		    mf[tag.dstIndex].copy(mf[tag.srcIndex],tag.sbox,scomp,tag.dbox,scomp,ncomp);
		}
	    }
	}
#endif /*BL_USE_MPI*/
    } 
	
    template <class FAB>
    void
    FillPeriodicBoundary_finish (const Geometry& geom,
				 FabArray<FAB>&  mf)
    {
        if (!geom.isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

	if (ParallelDescriptor::NProcs() == 1) return;

#ifdef BL_USE_MPI

#if defined(BL_USE_UPCXX)
    ParallelDescriptor::Mode.set_upcxx_mode();
    ParallelDescriptor::Mode.decr_upcxx();
#endif

#if !defined(BL_USE_MPI3)
        BL_ASSERT(!ParallelDescriptor::MPIOneSided());
#endif
	
        const Geometry::FPB& TheFPB = geom.getFPB(mf,mf.fpb_corners);

        typedef typename FAB::value_type value_type;

        const int N_rcvs = TheFPB.m_RcvTags->size();
        const int N_snds = TheFPB.m_SndTags->size();

#ifdef BL_USE_UPCXX
	if (N_rcvs > 0) BLPgas::fpb_recv_event.wait();
#else
        if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
            if (N_snds > 0) MPI_Win_complete(ParallelDescriptor::fpb_win);
	    if (N_rcvs > 0) MPI_Win_wait    (ParallelDescriptor::fpb_win);
#endif
	} else {
	    if (N_rcvs > 0) {
	        Array<MPI_Status> stats(N_rcvs);
	        BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, mf.fpb_recv_reqs.dataPtr(), stats.dataPtr()) );
            }
	}
#endif

	if (N_rcvs > 0)
	{
	    Array<const Geometry::FPB::FPBComTagsContainer*> recv_fctc;
	    recv_fctc.reserve(N_rcvs);

	    for (int k = 0; k < N_rcvs; k++)
	    {
		Geometry::FPB::MapOfFPBComTagContainers::const_iterator m_it 
		    = TheFPB.m_RcvTags->find(mf.fpb_recv_from[k]);
                BL_ASSERT(m_it != TheFPB.m_RcvTags->end());
		recv_fctc.push_back(&(m_it->second));
	    }
	    
#ifdef _OPENMP
#pragma omp parallel for if (TheFPB.m_threadsafe_rcv)
#endif
	    for (int k = 0; k < N_rcvs; k++) 
	    {
		value_type*  dptr = mf.fpb_recv_data[k];
		BL_ASSERT(dptr != 0);
		
		const Geometry::FPB::FPBComTagsContainer& fctc = *recv_fctc[k];

                for (Geometry::FPB::FPBComTagsContainer::const_iterator it = fctc.begin();
                     it != fctc.end(); ++it)
                {
                    const Box& bx  = it->dbox;
                    const int  Cnt = bx.numPts()*mf.fpb_ncomp;
                    mf[it->dstIndex].copyFromMem(bx,mf.fpb_scomp,mf.fpb_ncomp,dptr);
                    dptr += Cnt;
                }
            }

#ifdef BL_USE_UPCXX
	    BLPgas::free(mf.fpb_the_recv_data);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		MPI_Win_detach(ParallelDescriptor::fpb_win, mf.fpb_the_recv_data);
		BoxLib::The_Arena()->free(mf.fpb_the_recv_data);
		mf.fpb_recv_disp.clear();
#endif
	    } else {
		BoxLib::The_Arena()->free(mf.fpb_the_recv_data);
	    }
#endif

	    mf.fpb_recv_from.clear();
	    mf.fpb_recv_data.clear();
	    mf.fpb_recv_reqs.clear();
	}

	if (N_snds > 0) {
#ifdef BL_USE_UPCXX
	    FabArrayBase::WaitForAsyncSends_PGAS(N_snds,mf.fpb_send_data,
					         &BLPgas::fpb_send_event,
					         &BLPgas::fpb_send_counter);
#else
	    if (ParallelDescriptor::MPIOneSided()) {
#if defined(BL_USE_MPI3)
		for (int i = 0; i < N_snds; ++i)
		    BoxLib::The_Arena()->free(mf.fpb_send_data[i]);
#endif
	    } else {
		Array<MPI_Status> stats;
		FabArrayBase::WaitForAsyncSends(N_snds,mf.fpb_send_reqs,mf.fpb_send_data,stats);
	    }
#endif
	    mf.fpb_send_data.clear();
	    mf.fpb_send_reqs.clear();
	}
	
#ifdef BL_USE_TEAM
	ParallelDescriptor::MyTeam().MemoryBarrier();
#endif

#endif /*BL_USE_MPI*/
    }

    template <class FAB>
    void
    PeriodicCopy (const Geometry&      geom,
		  FabArray<FAB>&       dstmf,
		  const FabArray<FAB>& srcmf,
		  int                  dcomp,
		  int                  scomp,
		  int                  ncomp,
		  int                  dstng,
		  int                  srcng,
		  FabArrayBase::CpOp   op)
    {
	if (!geom.isAnyPeriodic() || dstmf.size() == 0 || srcmf.size() == 0) return;

	BL_PROFILE("Geometry::PeriodicCopy");

	const Box&                 domain = geom.Domain();
	const BoxArray&            dstba  = dstmf.boxArray();
	const BoxArray&            srcba  = srcmf.boxArray();
	const DistributionMapping& srcdm  = srcmf.DistributionMap();
	const IndexType&           srctyp = srcba.ixType();

	BL_ASSERT(srctyp == dstba.ixType());
	BL_ASSERT(dstng <= dstmf.nGrow());

	const Box& bbox_dst  = BoxLib::grow(dstba.minimalBox(), dstng);

	IntVect lo, hi;
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo.setVal(i, bbox_dst.bigEnd(i)  -domain.length(i));
	    hi.setVal(i, bbox_dst.smallEnd(i)+domain.length(i));
	}
	Box bbox_src(lo,hi,srctyp);

	BoxList bl(srctyp);
	Array<int> iproc;
	Array<int> boxid;
	Array<IntVect> boxshft;
	Array<IntVect> pshifts(26);

	for (int i = 0; i < srcba.size(); ++i) 
	{
	    const Box& bx0 = BoxLib::grow(srcba[i],srcng);

	    if (bbox_src.isEmpty() || !bbox_src.strictly_contains(bx0)) 
	    {
		geom.periodicShift(bbox_dst,bx0,pshifts);

		for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
		     it != End; ++it)
		{
		    const IntVect& iv   = *it;
		    const Box&     shft = bx0 + iv;
		    const Box&     bx   = shft & bbox_dst;
		    if (bx.ok()) {
			bl.push_back(bx);
			iproc.push_back(srcdm[i]);
			boxid.push_back(i);
			boxshft.push_back(iv);
		    }
		}
	    }
	}
	
	if (bl.size() > 0)
	{
	    BoxArray shftba(bl);

	    iproc.push_back(ParallelDescriptor::MyProc());
	    DistributionMapping shftdm(iproc);

	    FabArray<FAB> shftmf(shftba, ncomp, 0, shftdm);

#ifdef _OPENMP
#pragma omp parallel
#endif
	    for (MFIter mfi(shftmf); mfi.isValid(); ++mfi)
	    {
		const Box& bx = mfi.validbox();
		int i = mfi.index();
		Box sbx = bx - boxshft[i];
		shftmf[mfi].copy(srcmf[boxid[i]], sbx, scomp, bx, 0, ncomp);
	    }

	    dstmf.copy(shftmf, 0, dcomp, ncomp, 0, dstng, Periodicity::NonPeriodic(), op);
	}
    }
}

#endif /*_GEOMETRY_H_*/
