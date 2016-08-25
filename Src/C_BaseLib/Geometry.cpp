
#include <winstd.H>

#include <iostream>

#include <BoxArray.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <MultiFab.H>
#include <FArrayBox.H>
#include <BLProfiler.H>
#include <Utility.H>
#include <SPACE.H>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

#ifdef BL_MEM_PROFILING
#include <MemProfiler.H>
#endif

//
// The definition of some static data members.
//
int     Geometry::spherical_origin_fix = 0;
RealBox Geometry::prob_domain;
bool    Geometry::is_periodic[BL_SPACEDIM] = {D_DECL(0,0,0)};

Geometry::FPBCache       Geometry::m_TheFPBCache;
FabArrayBase::CacheStats Geometry::m_FPBC_stats("FPBCache");

std::ostream&
operator<< (std::ostream&   os,
            const Geometry& g)
{
    os << (CoordSys&) g << g.ProbDomain() << g.Domain();
    return os;
}

std::istream&
operator>> (std::istream& is,
            Geometry&     g)
{
    Box     bx;
    RealBox rb;

    is >> (CoordSys&) g >> rb >> bx;

    g.Domain(bx);
    Geometry::ProbDomain(rb);

    return is;
}

Geometry::FPB::FPB (const Geometry& geom, const FabArrayBase& fa, bool do_corners)
    : m_typ(fa.boxArray().ixType()), m_ngrow(fa.nGrow()), m_do_corners(do_corners),
      m_threadsafe_loc(false), m_threadsafe_rcv(false),
      m_LocTags(0), m_SndTags(0), m_RcvTags(0), m_SndVols(0),m_RcvVols(0),
      m_nuse(0)
{
    BL_PROFILE("Geometry::FPB::FPB()");

    BL_ASSERT(m_ngrow >= 0);
    BL_ASSERT(fa.boxArray().size() > 0);
    BL_ASSERT(isAnyPeriodic());
    BL_ASSERT(geom.domain.ok());

    const BoxArray&            ba     = fa.boxArray();
    const DistributionMapping& dm     = fa.DistributionMap();
    const Array<int>&          imap   = fa.IndexArray();
    const int                  MyProc = ParallelDescriptor::MyProc();
    
    m_LocTags = new FPB::FPBComTagsContainer;
    m_SndTags = new FPB::MapOfFPBComTagContainers;
    m_RcvTags = new FPB::MapOfFPBComTagContainers;
    m_SndVols = new std::map<int,int>;
    m_RcvVols = new std::map<int,int>;

    if (!imap.empty()) 
    {
	// All workers in the same team will have identical copies of tags for local copy
	// so that they can share work.  But for remote communication, they are all different.

	const int nlocal = imap.size();
	const int ng = m_ngrow;
	std::vector<std::pair<int,Box> > isects;
	Array<IntVect> pshifts(26);
	
	Box TheDomain = geom.Domain();
	TheDomain.convert(ba.ixType());
	const Box& GrownDomain = BoxLib::grow(TheDomain,ng);

	FPB::MapOfFPBComTagContainers send_tags; // temp copy
	
	for (int i = 0; i < nlocal; ++i)
	{
	    const   int k_src = imap[i];
	    const Box& bx_src = ba[k_src];
	    
	    if (TheDomain.contains(BoxLib::grow(bx_src,ng))) continue;
	    
	    geom.periodicShift(GrownDomain, bx_src, pshifts);
	    
	    for (Array<IntVect>::const_iterator pit = pshifts.begin(), pEnd = pshifts.end();
		 pit != pEnd; ++pit)
	    {
		const IntVect& iv   = *pit;
		const Box&     shft = bx_src + iv;
		
		ba.intersections(shft, isects, false, ng);
		
		for (int j = 0, M = isects.size(); j < M; ++j)
		{
		    const int  k_dst     = isects[j].first;
		    const Box& bx_dst    = ba[k_dst];
		    Box        bx        = isects[j].second;
		    const int  dst_owner = dm[k_dst];
		    
		    if (m_do_corners) {
			for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
			    if (!isPeriodic(dir)) {
				if (bx.smallEnd(dir) == TheDomain.smallEnd(dir) 
				    && bx_dst.smallEnd(dir) == TheDomain.smallEnd(dir) )
				{
				    bx.growLo(dir,ng);
				}
				if (bx.bigEnd(dir) == TheDomain.bigEnd(dir)
				    && bx_dst.bigEnd(dir) == TheDomain.bigEnd(dir) )
				{
				    bx.growHi(dir,ng);
				}
			    }
			}
		    }
		    
		    if (ParallelDescriptor::sameTeam(dst_owner)) {
			continue; // local copy will be dealt with later
		    } else if (MyProc == dm[k_src]) {
			const BoxList& bl = BoxLib::boxDiff(bx, bx_dst);
			for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
			    send_tags[dst_owner].push_back(FPBComTag(*lit-iv, *lit, k_src, k_dst));
		    }
		}
	    }
	}

	FPB::MapOfFPBComTagContainers recv_tags; // temp copy

	BaseFab<int> localtouch, remotetouch;
	bool check_local = false, check_remote = false;
#ifdef _OPENMP
	if (omp_get_max_threads() > 1) {
	    check_local = true;
	    check_remote = true;
	}
#endif    

	if (ParallelDescriptor::TeamSize() > 1) {
	    check_local = true;
	}
	
	if ( ba.ixType().cellCentered() ) {
	    m_threadsafe_loc = true;
	    m_threadsafe_rcv = true;
	    check_local = false;
	    check_remote = false;
	}

	for (int i = 0; i < nlocal; ++i)
	{
	    const int   k_dst   = imap[i];
	    const Box& bx_dst   = ba[k_dst];
	    const Box& bx_dst_g = BoxLib::grow(bx_dst, ng);
	    
	    if (TheDomain.contains(bx_dst_g)) continue;

	    if (check_local) {
		localtouch.resize(bx_dst_g);
		localtouch.setVal(0);
	    }

	    if (check_remote) {
		remotetouch.resize(bx_dst_g);
		remotetouch.setVal(0);
	    }
	    
	    geom.periodicShift(TheDomain, bx_dst_g, pshifts);
	    
	    for (Array<IntVect>::const_iterator pit = pshifts.begin(), pEnd = pshifts.end();
		 pit != pEnd; ++pit)
	    {
		const IntVect& iv   = *pit;
		const Box&     shft = bx_dst_g + iv;
		
		ba.intersections(shft, isects);
		
		for (int j = 0, M = isects.size(); j < M; ++j)
		{
		    const int k_src     = isects[j].first;
		    Box       bx        = isects[j].second;
		    const int src_owner = dm[k_src];
		    
		    if (m_do_corners) {
			for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
			    if (!geom.isPeriodic(dir)) {
				if (bx.smallEnd(dir) == TheDomain.smallEnd(dir)
				    && bx_dst.smallEnd(dir) == TheDomain.smallEnd(dir) )
				{
				    bx.growLo(dir,ng);
				}
				if (bx.bigEnd(dir) == TheDomain.bigEnd(dir)
				    && bx_dst.bigEnd(dir) == TheDomain.bigEnd(dir) )
				{
				    bx.growHi(dir,ng);
				}
			    }
			}
		    }

		    const BoxList& bl = BoxLib::boxDiff(bx-iv, bx_dst); // destinatin boxes
		    for (BoxList::const_iterator lit = bl.begin(); lit != bl.end(); ++lit)
		    {
			const Box& blbx = *lit;
			
			if (ParallelDescriptor::sameTeam(src_owner)) { // local copy
			    const BoxList tilelist(blbx, FabArrayBase::comm_tile_size);
			    for (BoxList::const_iterator
				     it_tile  = tilelist.begin(),
				     End_tile = tilelist.end();   it_tile != End_tile; ++it_tile)
			    {
				m_LocTags->push_back(FPBComTag((*it_tile)+iv, *it_tile,
								      k_src, k_dst));
			    }
			    if (check_local) {
				localtouch.plus(1, blbx);
			    }
			} else if (MyProc == dm[k_dst]) {
			    recv_tags[src_owner].push_back(FPBComTag(blbx+iv, blbx, k_src, k_dst));
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
	    FPB::MapOfFPBComTagContainers & Tags    = (ipass == 0) ? *m_SndTags : *m_RcvTags;
	    FPB::MapOfFPBComTagContainers & tmpTags = (ipass == 0) ?  send_tags :  recv_tags;
	    std::map<int,int>             & Vols    = (ipass == 0) ? *m_SndVols : *m_RcvVols;
	    
	    for (FPB::MapOfFPBComTagContainers::iterator
		     it  = tmpTags.begin(),
		     End = tmpTags.end();   it != End; ++it)
	    {
		const int key = it->first;
		std::vector<FPBComTag>& fctv = it->second;
		
		// We need to fix the order so that the send and recv processes match.
		std::sort(fctv.begin(), fctv.end());
		
		std::vector<FPBComTag> new_fctv;
		new_fctv.reserve(fctv.size());
		
		for (std::vector<FPBComTag>::const_iterator
			 it2  = fctv.begin(),
			 End2 = fctv.end();   it2 != End2; ++it2)
		{
		    const Box& sbx = it2->sbox;
		    const Box& dbx = it2->dbox;
		    IntVect diff = sbx.smallEnd() - dbx.smallEnd();

		    Vols[key] += sbx.numPts();
		    
		    const BoxList tilelist(sbx, FabArrayBase::comm_tile_size);
		    for (BoxList::const_iterator 
			     it_tile  = tilelist.begin(), 
			     End_tile = tilelist.end();    it_tile != End_tile; ++it_tile)
		    {
			new_fctv.push_back(FPBComTag(*it_tile, (*it_tile)-diff,
						     it2->srcIndex, it2->dstIndex));
		    }
		}

		Tags[key].swap(new_fctv);
	    } 
	}
    }
}

Geometry::FPB::~FPB ()
{
    delete m_LocTags;
    delete m_SndTags;
    delete m_RcvTags;
    delete m_SndVols;
    delete m_RcvVols;
}

long 
Geometry::FPB::bytesOfMapOfComTagContainers (const Geometry::FPB::MapOfFPBComTagContainers& m) const
{
    long r = sizeof(MapOfFPBComTagContainers);
    for (FPB::MapOfFPBComTagContainers::const_iterator it = m.begin(); it != m.end(); ++it) {
	r += sizeof(it->first) + BoxLib::bytesOf(it->second)
	    + BoxLib::gcc_map_node_extra_bytes;
    }
    return r;
}

long
Geometry::FPB::bytes () const
{
    long cnt = sizeof(Geometry::FPB);

    if (m_LocTags)
        cnt += BoxLib::bytesOf(*m_LocTags);

    if (m_SndTags)
	cnt += bytesOfMapOfComTagContainers(*m_SndTags);

    if (m_RcvTags)
        cnt += bytesOfMapOfComTagContainers(*m_RcvTags);

    if (m_SndVols)
	cnt += BoxLib::bytesOf(*m_SndVols);

    if (m_RcvVols)
	cnt += BoxLib::bytesOf(*m_RcvVols);

    return cnt;
}

void
Geometry::flushFPB (const FabArrayBase::BDKey& key)
{
    std::pair<FPBCacheIter,FPBCacheIter> er_it = m_TheFPBCache.equal_range(key);
    for (FPBCacheIter it = er_it.first; it != er_it.second; ++it)
    {
#ifdef BL_MEM_PROFILING
	m_FPBC_stats.bytes -= it->second->bytes();
#endif
	m_FPBC_stats.recordErase(it->second->m_nuse);
	delete it->second;
    }
    m_TheFPBCache.erase(er_it.first, er_it.second);
}

void
Geometry::flushFPBCache ()
{
    for (FPBCacheIter it = m_TheFPBCache.begin(); it != m_TheFPBCache.end(); ++it)
    {
	m_FPBC_stats.recordErase(it->second->m_nuse);
	delete it->second;
    }
    m_TheFPBCache.clear();
#ifdef BL_MEM_PROFILING
    m_FPBC_stats.bytes = 0L;
#endif    
}

const Geometry::FPB&
Geometry::getFPB (const FabArrayBase& fa, bool do_corners) const
{
    BL_PROFILE("FabArrayBase::getFPB()");

    std::pair<FPBCacheIter,FPBCacheIter> er_it = m_TheFPBCache.equal_range(fa.getBDKey());
    for (FPBCacheIter it = er_it.first; it != er_it.second; ++it)
    {
	if (it->second->m_typ        == fa.boxArray().ixType() &&
	    it->second->m_ngrow      == fa.nGrow()             &&
	    it->second->m_do_corners == do_corners)
	{
	    ++(it->second->m_nuse);
	    m_FPBC_stats.recordUse();
	    return *(it->second);
	}
    }

    // Have to build a new one
    FPB* new_fpb = new FPB(*this, fa, do_corners);

#ifdef BL_PROFILE
    m_FPBC_stats.bytes += new_fpb->bytes();
    m_FPBC_stats.bytes_hwm = std::max(m_FPBC_stats.bytes_hwm, m_FPBC_stats.bytes);
#endif

    new_fpb->m_nuse = 1;
    m_FPBC_stats.recordBuild();
    m_FPBC_stats.recordUse();

    m_TheFPBCache.insert(er_it.second, FPBCache::value_type(fa.getBDKey(),new_fpb));

    return *new_fpb;    
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                bool      do_corners,
                                bool      local) const
{
    FillPeriodicBoundary(mf,0,mf.nComp(),do_corners,local);
}

void
Geometry::FillPeriodicBoundary_nowait (MultiFab& mf,
				       bool      do_corners,
				       bool      local) const
{
    FillPeriodicBoundary_nowait(mf,0,mf.nComp(),do_corners,local);
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                int       scomp,
                                int       ncomp,
                                bool      corners,
                                bool      local) const
{
    BL_PROFILE("Geometry::FillPeriodicBoundary()");

    if ( local )
    {
	FillPeriodicBoundary_local(mf, scomp, ncomp, corners);
    }
    else
    {
	if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;
        BoxLib::FillPeriodicBoundary_nowait(*this, mf, scomp, ncomp, corners);
	BoxLib::FillPeriodicBoundary_finish(*this, mf);
    }
}

void
Geometry::FillPeriodicBoundary_local (MultiFab& mf,
				      int       scomp,
				      int       ncomp,
				      bool      corners) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

    BL_PROFILE("Geometry::FillPeriodicBoundary_local()");

    //
    // Do what you can with the FABs you own.  No parallelism allowed.
    //
    
    Box TheDomain = Domain();
    TheDomain.convert(mf.boxArray().ixType());

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<IntVect> pshifts(26);

        for (MFIter mfidst(mf); mfidst.isValid(); ++mfidst)
        {
            const Box& dst = mf[mfidst].box();

            BL_ASSERT(dst == BoxLib::grow(mfidst.validbox(), mf.nGrow()));

            if (TheDomain.contains(dst)) continue;

            // Turn off sharing among threads because this MFIter is inside another MFIter
	    unsigned char flags = MFIter::AllBoxes || MFIter::NoTeamBarrier;
            for (MFIter mfisrc(mf,flags); mfisrc.isValid(); ++mfisrc)
            {
                Box src = mfisrc.validbox() & TheDomain;

                if (corners)
                {
                    for (int i = 0; i < BL_SPACEDIM; i++)
                    {
                        if (!isPeriodic(i))
                        {
                            if (src.smallEnd(i) == Domain().smallEnd(i))
                                src.growLo(i,mf.nGrow());
                            if (src.bigEnd(i) == Domain().bigEnd(i))
                                src.growHi(i,mf.nGrow());
                        }
                    }
                }

                if (TheDomain.contains(src)) continue;

                periodicShift(dst, src, pshifts);

                for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                     it != End;
                     ++it)
                {
                    const IntVect& iv = *it;

                    const Box& shft = src + iv;
                    const Box& dbx  = dst & shft;
                    const Box& sbx  = dbx - iv;

                    mf[mfidst].copy(mf[mfisrc],sbx,scomp,dbx,scomp,ncomp);
                }
            }
        }
    }
}

void
Geometry::FillPeriodicBoundary_nowait (MultiFab& mf,
				       int       scomp,
				       int       ncomp,
				       bool      corners,
				       bool      local) const
{
    BL_PROFILE("Geometry::FillPeriodicBoundary_nowait()");

    if ( local )
    {
	FillPeriodicBoundary_local(mf, scomp, ncomp, corners);
    }
    else
    {
	if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;
        BoxLib::FillPeriodicBoundary_nowait(*this, mf, scomp, ncomp, corners);
    }
}

void
Geometry::FillPeriodicBoundary_finish (MultiFab& mf) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;
    BoxLib::FillPeriodicBoundary_finish(*this, mf);
}

void
Geometry::PeriodicCopy (MultiFab&       dstmf,
			const MultiFab& srcmf) const
{
    PeriodicCopy(dstmf, srcmf, 0, 0, srcmf.nComp());
}

void 
Geometry::PeriodicCopy (MultiFab&       dstmf,
			const MultiFab& srcmf,
			int             dcomp,
			int             scomp,
			int             ncomp,
			int             dstng,
			int             srcng) const
{
    BoxLib::PeriodicCopy(*this, dstmf, srcmf, dcomp, scomp, ncomp, dstng, srcng,
			 FabArrayBase::COPY);
}

Geometry::Geometry () {}

Geometry::Geometry (const Box&     dom,
                    const RealBox* rb,
                    int            coord,
                    int*           is_per)
{
    define(dom,rb,coord,is_per);
}

Geometry::Geometry (const Geometry& g)
{
    ok     = g.ok;
    domain = g.domain;

    D_TERM(dx[0]=g.dx[0];,dx[1]=g.dx[1];,dx[2]=g.dx[2];)
}

Geometry::~Geometry() {}

void
Geometry::define (const Box&     dom,
                  const RealBox* rb,
                  int            coord,
                  int*           is_per)
{
    if (c_sys == undef)
        Setup(rb,coord,is_per);

    domain = dom;
    ok     = true;

    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
    }
    if (Geometry::spherical_origin_fix == 1)
    {
	if (c_sys == SPHERICAL && prob_domain.lo(0) == 0 && BL_SPACEDIM > 1)
        {
            prob_domain.setLo(0,2*dx[0]);

            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
            }
	}
    } 
}

void
Geometry::Finalize ()
{
    c_sys = undef;

    Geometry::flushFPBCache();
    if (ParallelDescriptor::IOProcessor() && BoxLib::verbose) {
	m_FPBC_stats.print();
    }
}

void
Geometry::Setup (const RealBox* rb, int coord, int* isper)
{
    ParmParse pp("geometry");
    //
    // The default behavior is as before.  If rb and coord come
    // in with default values, we require that user set them through pp.
    // If not, use those coming in, and possibly override them w/pp
    //
    Array<Real> prob_lo(BL_SPACEDIM);
    Array<Real> prob_hi(BL_SPACEDIM);
    if (rb == 0  &&  coord==-1)
    {
        pp.get("coord_sys",coord);
        SetCoord( (CoordType) coord );
        pp.getarr("prob_lo",prob_lo,0,BL_SPACEDIM);
        BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
        pp.getarr("prob_hi",prob_hi,0,BL_SPACEDIM);
        BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
        prob_domain.setLo(prob_lo);
        prob_domain.setHi(prob_hi);
    }
    else
    {
        BL_ASSERT(rb != 0  &&  coord != -1);
        pp.query("coord_sys",coord);
        SetCoord( (CoordType) coord );
        prob_domain.setLo(rb->lo());
        prob_domain.setHi(rb->hi());

        if (pp.countval("prob_lo")>0)
        {
            pp.queryarr("prob_lo",prob_lo,0,BL_SPACEDIM);
            BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
            prob_domain.setLo(prob_lo);
        }
        if (pp.countval("prob_hi")>0)
        {
            pp.queryarr("prob_hi",prob_hi,0,BL_SPACEDIM);
            BL_ASSERT(prob_hi.size() == BL_SPACEDIM);
            prob_domain.setHi(prob_hi);
        }
    }
    pp.query("spherical_origin_fix", Geometry::spherical_origin_fix);
    //
    // Now get periodicity info.
    //
    if (isper == 0)
    {
        Array<int> is_per(BL_SPACEDIM);
        pp.queryarr("is_periodic",is_per,0,BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; n++)  
            is_periodic[n] = is_per[n];
    }
    else
    {
        for (int n = 0; n < BL_SPACEDIM; n++)  
            is_periodic[n] = isper[n];
    }

#ifdef BL_MEM_PROFILING
    MemProfiler::add(m_FPBC_stats.name, std::function<MemProfiler::MemInfo()>
		     ([] () -> MemProfiler::MemInfo {
			 return {m_FPBC_stats.bytes, m_FPBC_stats.bytes_hwm};
		     }));
#endif

    BoxLib::ExecOnFinalize(Geometry::Finalize);
}

void
Geometry::GetVolume (MultiFab&       vol,
                     const BoxArray& grds,
                     int             ngrow) const
{
    vol.define(grds,1,ngrow,Fab_allocate);
    GetVolume(vol);
}

void
Geometry::GetVolume (MultiFab&       vol) const
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(vol,true); mfi.isValid(); ++mfi)
    {
	CoordSys::SetVolume(vol[mfi], mfi.growntilebox());
    }
}

void
Geometry::GetVolume (FArrayBox&      vol,
                     const BoxArray& grds,
                     int             idx,
                     int             ngrow) const
{
    CoordSys::GetVolume(vol, BoxLib::grow(grds[idx],ngrow));
}

#if (BL_SPACEDIM <= 2)
void
Geometry::GetDLogA (MultiFab&       dloga,
                    const BoxArray& grds, 
                    int             dir,
                    int             ngrow) const
{
    dloga.define(grds,1,ngrow,Fab_allocate);
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(dloga,true); mfi.isValid(); ++mfi)
    {
	CoordSys::SetDLogA(dloga[mfi], mfi.growntilebox(), dir);
    }
}
#endif

void
Geometry::GetFaceArea (MultiFab&       area,
                       const BoxArray& grds,
                       int             dir,
                       int             ngrow) const
{
    BoxArray edge_boxes(grds);
    edge_boxes.surroundingNodes(dir);
    area.define(edge_boxes,1,ngrow,Fab_allocate);

    GetFaceArea(area, dir);
}

void
Geometry::GetFaceArea (MultiFab&       area,
                       int             dir) const
{
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter mfi(area,true); mfi.isValid(); ++mfi)
    {
	CoordSys::SetFaceArea(area[mfi],mfi.growntilebox(),dir);
    }
}

void
Geometry::GetFaceArea (FArrayBox&      area,
                       const BoxArray& grds,
                       int             idx,
                       int             dir,
                       int             ngrow) const
{
    CoordSys::GetFaceArea(area, BoxLib::grow(grds[idx],ngrow), dir);
}

void
Geometry::periodicShift (const Box&      target,
                         const Box&      src, 
                         Array<IntVect>& out) const
{
    out.resize(0);

    Box locsrc(src);

    int nist,njst,nkst;
    int niend,njend,nkend;
    nist = njst = nkst = 0;
    niend = njend = nkend = 0;
    D_TERM( nist , =njst , =nkst ) = -1;
    D_TERM( niend , =njend , =nkend ) = +1;

    int ri,rj,rk;
    for (ri = nist; ri <= niend; ri++)
    {
        if (ri != 0 && !is_periodic[0])
            continue;
        if (ri != 0 && is_periodic[0])
            locsrc.shift(0,ri*domain.length(0));

        for (rj = njst; rj <= njend; rj++)
        {
            if (rj != 0 && !is_periodic[1])
                continue;
            if (rj != 0 && is_periodic[1])
                locsrc.shift(1,rj*domain.length(1));

            for (rk = nkst; rk <= nkend; rk++)
            {
                if (rk!=0
#if (BL_SPACEDIM == 3)
                    && !is_periodic[2]
#endif
                    )
                {
                    continue;
                }
                if (rk!=0
#if (BL_SPACEDIM == 3)
                    && is_periodic[2]
#endif
                    )
                {
                    locsrc.shift(2,rk*domain.length(2));
                }

                if (ri == 0 && rj == 0 && rk == 0)
                    continue;
                //
                // If losrc intersects target, then add to "out".
                //
                if (target.intersects(locsrc))
                {
                    out.push_back(IntVect(D_DECL(ri*domain.length(0),
                                                 rj*domain.length(1),
                                                 rk*domain.length(2))));
                }
                if (rk != 0
#if (BL_SPACEDIM == 3)
                    && is_periodic[2]
#endif
                    )
                {
                    locsrc.shift(2,-rk*domain.length(2));
                }
            }
            if (rj != 0 && is_periodic[1])
                locsrc.shift(1,-rj*domain.length(1));
        }
        if (ri != 0 && is_periodic[0])
            locsrc.shift(0,-ri*domain.length(0));
    }
}

#ifdef BL_USE_MPI
void
Geometry::SendGeometryToSidecar (Geometry *geom, int whichSidecar)
{
  int fromProc;

  MPI_Comm commSource = ParallelDescriptor::CommunicatorComp();
  MPI_Comm commInter  = ParallelDescriptor::CommunicatorInter(whichSidecar);
  MPI_Comm comm = commInter;

  bool bcastSource(ParallelDescriptor::Communicator() == commSource);

  if(bcastSource) {
    fromProc = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
    BL_ASSERT(ParallelDescriptor::IOProcessorNumber() == 0);  // ---- because we are assuming this in commDest
  }
  if( ! bcastSource) {
    fromProc = 0;  // ---- really the rank of MPI_ROOT in commSource
  }

  Geometry::BroadcastGeometry(*geom, fromProc, comm, bcastSource);
}



void
Geometry::BroadcastGeometry (Geometry &geom, int fromProc, MPI_Comm comm)
{
  bool bcastSource(ParallelDescriptor::MyProc() == fromProc);
  Geometry::BroadcastGeometry(geom, fromProc, comm, bcastSource);
}



void
Geometry::BroadcastGeometry (Geometry &geom, int fromProc, MPI_Comm comm, bool bcastSource)
{
  int coord;
  int is_periodic[BL_SPACEDIM];
  Real realBox_lo[BL_SPACEDIM];
  Real realBox_hi[BL_SPACEDIM];
  Array<int> baseBoxAI;

  CoordSys::BroadcastCoordSys(geom, fromProc, comm, bcastSource);

  if(bcastSource) {  // ---- initialize the source data
    const RealBox &realBox = geom.ProbDomain();
    for(int n(0); n < BL_SPACEDIM; ++n) {
      realBox_lo[n] = realBox.lo(n);
      realBox_hi[n] = realBox.hi(n);
      is_periodic[n] = geom.isPeriodic(n);
    }
    coord = geom.CoordInt();
    baseBoxAI = BoxLib::SerializeBox(geom.Domain());
  }


  // ---- do the broadcasts
  if( ! bcastSource) {
    baseBoxAI.resize(BoxLib::SerializeBoxSize());
  }
  ParallelDescriptor::Bcast(baseBoxAI.dataPtr(), baseBoxAI.size(), fromProc, comm);

  ParallelDescriptor::Bcast(realBox_lo, BL_SPACEDIM, fromProc, comm);
  ParallelDescriptor::Bcast(realBox_hi, BL_SPACEDIM, fromProc, comm);

  ParallelDescriptor::Bcast(&coord, 1, fromProc, comm);
  ParallelDescriptor::Bcast(is_periodic, BL_SPACEDIM, fromProc, comm);
  ParallelDescriptor::Bcast(&Geometry::spherical_origin_fix, 1, fromProc, comm);


  if( ! bcastSource) {    // ---- define the destination geometry
    Box baseBox(BoxLib::UnSerializeBox(baseBoxAI));
    RealBox realBox(realBox_lo, realBox_hi);

    geom.define(baseBox, &realBox, coord, is_periodic);
  }
}
#endif
