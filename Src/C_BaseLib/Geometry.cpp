
#include <winstd.H>

#include <iostream>

#include <BoxArray.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <MultiFab.H>
#include <FArrayBox.H>
#include <BLProfiler.H>

#ifdef BL_LAZY
#include <Lazy.H>
#endif

//
// The definition of some static data members.
//
int     Geometry::spherical_origin_fix;
RealBox Geometry::prob_domain;
bool    Geometry::is_periodic[BL_SPACEDIM];

namespace
{
    bool verbose;
}

const int fpb_cache_max_size_def = 25;

int Geometry::fpb_cache_max_size = fpb_cache_max_size_def;

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

Geometry::FPB::FPB ()
    :
    m_ngrow(-1),
    m_do_corners(false),
    m_reused(false),
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0) {}

Geometry::FPB::FPB (const BoxArray&            ba,
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
    m_threadsafe_loc(false),
    m_threadsafe_rcv(false),
    m_LocTags(0),
    m_SndTags(0),
    m_RcvTags(0),
    m_SndVols(0),
    m_RcvVols(0)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(domain.ok());
}

Geometry::FPB::~FPB ()
{
    delete m_LocTags;
    delete m_SndTags;
    delete m_RcvTags;
    delete m_SndVols;
    delete m_RcvVols;
}

bool
Geometry::FPB::operator== (const FPB& rhs) const
{
    return
        m_ngrow == rhs.m_ngrow && m_do_corners == rhs.m_do_corners && m_domain == rhs.m_domain && m_ba == rhs.m_ba && m_dm == rhs.m_dm;
}

int
Geometry::FPB::bytes () const
{
    int cnt = sizeof(Geometry::FPB);

    if (m_LocTags)
    {
        cnt += sizeof(FPBComTagsContainer) + m_LocTags->size()*sizeof(FPBComTag);
    }

    if (m_SndTags)
    {
        cnt += sizeof(MapOfFPBComTagContainers);

        cnt += m_SndTags->size()*sizeof(MapOfFPBComTagContainers::value_type);

        for (FPB::MapOfFPBComTagContainers::const_iterator it = m_SndTags->begin(),
                 m_End = m_SndTags->end();
             it != m_End;
             ++it)
        {
            cnt += it->second.size()*sizeof(FPBComTag);
        }
    }

    if (m_RcvTags)
    {
        cnt += sizeof(MapOfFPBComTagContainers);

        cnt += m_SndTags->size()*sizeof(MapOfFPBComTagContainers::value_type);

        for (FPB::MapOfFPBComTagContainers::const_iterator it = m_RcvTags->begin(),
                 m_End = m_RcvTags->end();
             it != m_End;
             ++it)
        {
            cnt += it->second.size()*sizeof(FPBComTag);
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

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                bool      do_corners,
                                bool      local) const
{
    FillPeriodicBoundary(mf,0,mf.nComp(),do_corners,local);
}

void
Geometry::SumPeriodicBoundary (MultiFab& mf) const
{
    SumPeriodicBoundary(mf,0,mf.nComp());
}

void
Geometry::SumPeriodicBoundary (MultiFab&       dstmf,
                               const MultiFab& srcmf) const
{
    BL_ASSERT(dstmf.nComp() >= srcmf.nComp());

    SumPeriodicBoundary(dstmf, srcmf, 0, 0, srcmf.nComp());
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                int       scomp,
                                int       ncomp,
                                bool      corners,
                                bool      local) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

    BL_PROFILE("Geometry::FillPeriodicBoundary()");

    Box TheDomain = Domain();
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (mf.boxArray().ixType()[n] == IndexType::NODE)
            TheDomain.surroundingNodes(n);

    if ( local )
    {
        //
        // Do what you can with the FABs you own.  No parallelism allowed.
        //

#ifdef _OPENMP
#pragma omp parallel
#endif
    {
        Array<IntVect> pshifts(27);

        for (MFIter mfidst(mf); mfidst.isValid(); ++mfidst)
        {
            const Box& dst = mf[mfidst].box();

            BL_ASSERT(dst == BoxLib::grow(mfidst.validbox(), mf.nGrow()));

            if (TheDomain.contains(dst)) continue;

            // Turn off sharing among threas because this MFIter is inside another MFIter
	    unsigned char flags = MFIter::NoSharing;
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
    else
    {
        BoxLib::FillPeriodicBoundary(*this, mf, scomp, ncomp, corners);
    }
}

//
// Some useful typedefs.
//
typedef FabArrayBase::CopyComTag::CopyComTagsContainer CopyComTagsContainer;

typedef FabArrayBase::CopyComTag::MapOfCopyComTagContainers MapOfCopyComTagContainers;

static
void
SumPeriodicBoundaryInnards (MultiFab&       dstmf,
                            const MultiFab& srcmf,
                            const Geometry& geom,
                            int             scomp,
                            int             dcomp,
                            int             ncomp)
{
    BL_PROFILE("SumPeriodicBoundaryInnards()");

#ifndef NDEBUG
    //
    // Don't let folks ask for more grow cells than they have valid region.
    //
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (geom.isPeriodic(n))
            BL_ASSERT(srcmf.nGrow() <= geom.Domain().length(n));
#endif
    //
    // Note that in the usual case dstmf == srcmf.
    //
    Box TheDomain = geom.Domain();
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (srcmf.boxArray().ixType()[n] == IndexType::NODE)
            TheDomain.surroundingNodes(n);

    MapOfCopyComTagContainers  m_SndTags, m_RcvTags;
    std::map<int,int>          m_SndVols, m_RcvVols;
    Array<IntVect>             pshifts(27);
    const int                  ngrow  = srcmf.nGrow();
    const int                  MyProc = ParallelDescriptor::MyProc();
    const BoxArray&            srcba  = srcmf.boxArray();
    const BoxArray&            dstba  = dstmf.boxArray();
    const DistributionMapping& srcdm  = srcmf.DistributionMap();
    const DistributionMapping& dstdm  = dstmf.DistributionMap();

    for (int i = 0, K = dstba.size(); i < K; i++)
    {
        const Box& dst       = dstba[i];
        const int  dst_owner = dstdm[i];

        if (TheDomain.contains(BoxLib::grow(dst,ngrow))) continue;
        //
        // We may overlap with the periodic boundary.  Some other ghost
        // region(s) could be periodically-shiftable into our valid region.
        //
        for (int j = 0, N = srcba.size(); j < N; j++)
        {
            const int src_owner = srcdm[j];

            if (dst_owner != MyProc && src_owner != MyProc) continue;

            const Box& src = BoxLib::grow(srcba[j],ngrow);

            if (TheDomain.contains(src)) continue;

            geom.periodicShift(dst, src, pshifts);

            for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                 it != End;
                 ++it)
            {
                const IntVect& iv  = *it;
                const Box&    shft = src + iv;
                const Box&    dbx  = dst & shft;
                const Box&    sbx  = dbx - iv;

                FabArrayBase::CopyComTag tag;

                if (dst_owner == MyProc)
                {
                    if (src_owner == MyProc)
                    {
                        //
                        // Do the local work right here.
                        //
                        dstmf[i].plus(srcmf[j],sbx,dbx,scomp,dcomp,ncomp);
                    }
                    else
                    {
                        tag.box      = dbx;
                        tag.fabIndex = i;

                        FabArrayBase::SetRecvTag(m_RcvTags,src_owner,tag,m_RcvVols,tag.box);
                    }
                }
                else if (src_owner == MyProc)
                {
                    tag.box      = sbx;
                    tag.fabIndex = j;

                    FabArrayBase::SetSendTag(m_SndTags,dst_owner,tag,m_SndVols,tag.box);
                }
            }
        }
    }

#ifdef BL_USE_MPI
    if (ParallelDescriptor::NProcs() == 1) return;
    //
    // Do this before prematurely exiting if running in parallel.
    // Otherwise sequence numbers will not match across MPI processes.
    //
    const int SeqNum = ParallelDescriptor::SeqNum();

    if (m_SndTags.empty() && m_RcvTags.empty())
        //
        // No parallel work for this MPI process to do.
        //
        return;

    Array<MPI_Status>  stats;
    Array<int>         recv_from;
    Array<Real*>       recv_data;
    Array<MPI_Request> recv_reqs;
    //
    // Post rcvs. Allocate one chunk of space to hold'm all.
    //
    Real* the_recv_data = 0;

    FabArrayBase::PostRcvs(m_RcvTags,m_RcvVols,the_recv_data,recv_data,recv_from,recv_reqs,ncomp,SeqNum);

    //
    // Post send's
    //
    const int N_snds = m_SndTags.size();

    Array<Real*>                       send_data;
    Array<int>                         send_N;
    Array<int>                         send_rank;
    Array<const CopyComTagsContainer*> send_cctc;

    send_data.reserve(N_snds);
    send_N   .reserve(N_snds);
    send_rank.reserve(N_snds);
    send_cctc.reserve(N_snds);

    for (MapOfCopyComTagContainers::const_iterator m_it = m_SndTags.begin(),
             m_End = m_SndTags.end();
         m_it != m_End;
         ++m_it)
    {
        std::map<int,int>::const_iterator vol_it = m_SndVols.find(m_it->first);

        BL_ASSERT(vol_it != m_SndVols.end());

        const int N = vol_it->second*ncomp;

        BL_ASSERT(N < std::numeric_limits<int>::max());

        Real* data = static_cast<Real*>(BoxLib::The_Arena()->alloc(N*sizeof(Real)));

	send_data.push_back(data);
	send_N   .push_back(N);
	send_rank.push_back(m_it->first);
	send_cctc.push_back(&(m_it->second));
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i=0; i<N_snds; ++i)
    {
	Real* dptr = send_data[i];
	BL_ASSERT(dptr != 0);
	
	const CopyComTagsContainer& cctc = *send_cctc[i];
	
	for (CopyComTagsContainer::const_iterator it = cctc.begin();
	     it != cctc.end(); ++it)
	{
	    const Box& bx = it->box;
	    srcmf[it->fabIndex].copyToMem(bx,scomp,ncomp,dptr);
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
    // wait and unpack
    //
    const int N_rcvs = m_RcvTags.size();

    if (N_rcvs > 0)
    {
	Array<const CopyComTagsContainer*> recv_cctc;
	recv_cctc.reserve(N_rcvs);

	for (int k=0; k<N_rcvs; k++)
	{
	    MapOfCopyComTagContainers::const_iterator m_it = m_RcvTags.find(recv_from[k]);
            BL_ASSERT(m_it != m_RcvTags.end());

	    recv_cctc.push_back(&(m_it->second));
	}

	stats.resize(N_rcvs);
	BL_MPI_REQUIRE( MPI_Waitall(N_rcvs, recv_reqs.dataPtr(), stats.dataPtr()) );
	    
	// unsafe to use OMP here
	{
	    FArrayBox fab;

	    for (int k = 0; k < N_rcvs; k++) 
	    {
		Real* dptr = recv_data[k];
		BL_ASSERT(dptr != 0);
		
		const CopyComTagsContainer& cctc = *recv_cctc[k];
		
		for (CopyComTagsContainer::const_iterator it = cctc.begin();
		     it != cctc.end(); ++it)
		{
		    const Box& bx = it->box;
		    fab.resize(bx,ncomp);
		    const int Cnt = bx.numPts()*ncomp;
		    memcpy(fab.dataPtr(),dptr,Cnt*sizeof(Real));
		    dstmf[it->fabIndex].plus(fab,bx,bx,0,dcomp,ncomp);
		    dptr += Cnt;
		}
	    }
        }
    }

    BoxLib::The_Arena()->free(the_recv_data);

    if (FabArrayBase::do_async_sends && !m_SndTags.empty())
        FabArrayBase::GrokAsyncSends(m_SndTags.size(),send_reqs,send_data,stats);

#endif /*BL_USE_MPI*/
}

void
Geometry::SumPeriodicBoundary (MultiFab& mf,
                               int       scomp,
                               int       ncomp) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

    SumPeriodicBoundaryInnards(mf,mf,*this,scomp,scomp,ncomp);
}

void
Geometry::SumPeriodicBoundary (MultiFab&       dstmf,
                               const MultiFab& srcmf,
                               int             dcomp,
                               int             scomp,
                               int             ncomp) const
{
    if (!isAnyPeriodic() || srcmf.nGrow() == 0 || srcmf.size() == 0 || dstmf.size() == 0) return;

    BL_ASSERT(scomp+ncomp <= srcmf.nComp());
    BL_ASSERT(dcomp+ncomp <= dstmf.nComp());
    BL_ASSERT(srcmf.boxArray().ixType() == dstmf.boxArray().ixType());

    SumPeriodicBoundaryInnards(dstmf,srcmf,*this,scomp,dcomp,ncomp);
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

    Geometry::FlushPIRMCache();
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
    //
    // Set default values here!!!
    //
    verbose                        = true;
    Geometry::spherical_origin_fix = 0;
    Geometry::fpb_cache_max_size   = fpb_cache_max_size_def;

    D_EXPR(is_periodic[0]=0, is_periodic[1]=0, is_periodic[2]=0);

    pp.query("verbose",              verbose);
    pp.query("spherical_origin_fix", Geometry::spherical_origin_fix);
    pp.query("fpb_cache_max_size",   Geometry::fpb_cache_max_size);
    //
    // Don't let the cache size get too small.  This simplifies some logic later.
    //
    if (Geometry::fpb_cache_max_size < 1)
        Geometry::fpb_cache_max_size = 1;
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
    dloga.define(grds,1,ngrow,Fab_noallocate);
    for (MFIter mfi(dloga); mfi.isValid(); ++mfi)
    {
        const Box& gbx = BoxLib::grow(grds[mfi.index()],ngrow);
        dloga.setFab(mfi,CoordSys::GetDLogA(gbx,dir));
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

//
// The cache.
//
Geometry::FPBMMap Geometry::m_FPBCache;

Geometry::FPBMMapIter
Geometry::GetFPB (const Geometry&      geom,
                  const Geometry::FPB& fpb,
                  const FabArrayBase&  mf)
{
    BL_PROFILE("Geometry::GetFPB");

    BL_ASSERT(fpb.m_ngrow > 0);
    BL_ASSERT(fpb.m_ba.size() > 0);
    BL_ASSERT(geom.isAnyPeriodic());

    const BoxArray&            ba     = fpb.m_ba;
    const DistributionMapping& dm     = fpb.m_dm;
    const int                  MyProc = ParallelDescriptor::MyProc();
    const IntVect&             Typ    = ba[0].type();
    const int                  Scale  = D_TERM(Typ[0],+3*Typ[1],+5*Typ[2]) + 11;
    const int                  Key    = ba.size() + ba[0].numPts() + Scale + fpb.m_ngrow;

    std::pair<Geometry::FPBMMapIter,Geometry::FPBMMapIter> er_it = m_FPBCache.equal_range(Key);
    
    for (Geometry::FPBMMapIter it = er_it.first; it != er_it.second; ++it)
    {
        if (it->second == fpb)
        {
            it->second.m_reused = true;

            return it;
        }
    }

    if (m_FPBCache.size() >= Geometry::fpb_cache_max_size)
    {
        //
        // Don't let the size of the cache get too big.
        // Get rid of entries with the biggest largest key that haven't been reused.
        // Otherwise just remove the entry with the largest key.
        //
        Geometry::FPBMMapIter End      = m_FPBCache.end();
        Geometry::FPBMMapIter last_it  = End;
        Geometry::FPBMMapIter erase_it = End;

        for (Geometry::FPBMMapIter it = m_FPBCache.begin(); it != End; ++it)
        {
            last_it = it;

            if (!it->second.m_reused)
                erase_it = it;
        }

        if (erase_it != End)
        {
            m_FPBCache.erase(erase_it);
        }
        else if (last_it != End)
        {
            m_FPBCache.erase(last_it);
        }
    }
    //
    // Got to insert one & then build it.
    //
    Geometry::FPBMMapIter cache_it = m_FPBCache.insert(FPBMMap::value_type(Key,fpb));
    FPB&                  TheFPB   = cache_it->second;
    //
    // Here's where we allocate memory for the cache innards.
    // We do this so we don't have to build objects of these types
    // each time we search the cache.  Otherwise we'd be constructing
    // and destroying said objects quite frequently.
    //
    TheFPB.m_LocTags = new FPB::FPBComTagsContainer;
    TheFPB.m_SndTags = new FPB::MapOfFPBComTagContainers;
    TheFPB.m_RcvTags = new FPB::MapOfFPBComTagContainers;
    TheFPB.m_SndVols = new std::map<int,int>;
    TheFPB.m_RcvVols = new std::map<int,int>;

    if (mf.IndexMap().empty())
        //
        // We don't own any of the relevant FABs so can't possibly have any work to do.
        //
        return cache_it;

    Box TheDomain = geom.Domain();
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (ba.ixType()[n] == IndexType::NODE)
            TheDomain.surroundingNodes(n);

    Array<IntVect> pshifts(27);

    for (int i = 0, N = ba.size(); i < N; i++)
    {
        const Box& dst      = BoxLib::grow(ba[i],fpb.m_ngrow);
        const int dst_owner = dm[i];

        if (TheDomain.contains(dst)) continue;

        for (int j = 0, N = ba.size(); j < N; j++)
        {
            const int src_owner = dm[j];

            if (dst_owner != MyProc && src_owner != MyProc) continue;

            Box src = ba[j] & TheDomain;

            if (TheDomain.contains(BoxLib::grow(src,fpb.m_ngrow))) continue;

            if (fpb.m_do_corners)
            {
                for (int i = 0; i < BL_SPACEDIM; i++)
                {
                    if (!geom.isPeriodic(i))
                    {
                        src.growLo(i,fpb.m_ngrow);
                        src.growHi(i,fpb.m_ngrow);
                    }
                }
            }

            geom.periodicShift(dst, src, pshifts);

            for (Array<IntVect>::const_iterator it = pshifts.begin(), End = pshifts.end();
                 it != End;
                 ++it)
            {
                const IntVect& iv   = *it;
                const Box&     shft = src + iv;

                Box dbox_tmp = dst & shft;

		const BoxList tilelist(dbox_tmp, FabArrayBase::comm_tile_size);

		for (BoxList::const_iterator it = tilelist.begin(), End = tilelist.end(); it != End; ++it)
                {
                    FPBComTag tag;
                    tag.dbox     = *it;
                    tag.sbox     = tag.dbox - iv;
                    tag.dstIndex = i;
                    tag.srcIndex = j;

                    if (dst_owner == MyProc)
                    {
                        if (src_owner == MyProc)
                        {
                            TheFPB.m_LocTags->push_back(tag);
                        }
                        else
                        {
                            FabArrayBase::SetRecvTag(*TheFPB.m_RcvTags,src_owner,tag,*TheFPB.m_RcvVols,tag.dbox);
                        }
                    }
                    else if (src_owner == MyProc)
                    {
                        FabArrayBase::SetSendTag(*TheFPB.m_SndTags,dst_owner,tag,*TheFPB.m_SndVols,tag.dbox);
                    }
                }
            }
        }
    }
    //
    // Squeeze out any unused memory ...
    //
    FPB::FPBComTagsContainer tmp(*TheFPB.m_LocTags); 

    TheFPB.m_LocTags->swap(tmp);

    for (FPB::MapOfFPBComTagContainers::iterator it = TheFPB.m_SndTags->begin(), End = TheFPB.m_SndTags->end(); it != End; ++it)
    {
        FPB::FPBComTagsContainer tmp(it->second);

        it->second.swap(tmp);
    }

    for (FPB::MapOfFPBComTagContainers::iterator it = TheFPB.m_RcvTags->begin(), End = TheFPB.m_RcvTags->end(); it != End; ++it)
    {
        FPB::FPBComTagsContainer tmp(it->second);

        it->second.swap(tmp);
    }

    //
    // set thread safety
    //
    if ( ba[0].cellCentered() ) {
	TheFPB.m_threadsafe_loc = true;
	TheFPB.m_threadsafe_rcv = true;
    } else {
	TheFPB.m_threadsafe_loc = false;
	TheFPB.m_threadsafe_rcv = false;
    }

    return cache_it;
}

void
Geometry::FlushPIRMCache ()
{
    long stats[3] = {0,0,0}; // size, reused, bytes

    stats[0] = m_FPBCache.size();

    for (FPBMMapIter it = m_FPBCache.begin(), End = m_FPBCache.end(); it != End; ++it)
    {
        stats[2] += it->second.bytes();
        if (it->second.m_reused)
            stats[1]++;
    }

    if (verbose)
    {
#ifdef BL_LAZY
	Lazy::QueueReduction( [=] () mutable {
#endif
        ParallelDescriptor::ReduceLongMax(&stats[0], 3, ParallelDescriptor::IOProcessorNumber());
        if (stats[0] > 0 && ParallelDescriptor::IOProcessor())
        {
            std::cout << "Geometry::TheFPBCache: max size: "
                      << stats[0]
                      << ", max # reused: "
                      << stats[1]
                      << ", max bytes used: "
                      << stats[2]
                      << std::endl;
        }
#ifdef BL_LAZY
	});
#endif
    }

    m_FPBCache.clear();
}

int
Geometry::PIRMCacheSize ()
{
    return m_FPBCache.size();
}

#ifdef BL_USE_MPI
void
Geometry::SendGeometryToSidecars (Geometry *geom)
{
  const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;

  int coord;
  int is_periodic[BL_SPACEDIM];

  if (ParallelDescriptor::Communicator() == ParallelDescriptor::CommunicatorComp())
  {
    // Data to construct base Box
    const Box& baseBox = geom->Domain();
    const IntVect box_index_type = baseBox.type();
    const int* box_index_type_IV = box_index_type.getVect();
    const int *smallEnd = baseBox.smallEnd().getVect();
    const int *bigEnd = baseBox.bigEnd().getVect();

    // Data to construct RealBox
    const RealBox& realBox = geom->ProbDomain();
    const Real *realBox_lo = realBox.lo();
    const Real *realBox_hi = realBox.hi();

    CoordType coordtype = geom->Coord();
    // UGHHH STOP USING ENUMS IN MPI CODE. IT DESTROYS LIVES.
    int coord;
    switch (coordtype)
    {
        case undef:
            coord = -1;
            break;
        case cartesian:
            coord = 0;
            break;
        case RZ:
            coord = 1;
            break;
        case SPHERICAL:
            coord = 2;
            break;
    }

    for (unsigned int i = 0; i < BL_SPACEDIM; ++i)
    {
      is_periodic[i] = geom->isPeriodic(i);
    }

      // Step 1: send the base Box
      ParallelDescriptor::Bcast(const_cast<int*>(box_index_type_IV), BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(const_cast<int*>(smallEnd)      , BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(const_cast<int*>(bigEnd)        , BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Step 2: send the RealBox
      ParallelDescriptor::Bcast(const_cast<Real*>(realBox_lo), BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
      ParallelDescriptor::Bcast(const_cast<Real*>(realBox_hi), BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Step 3: send the coordinates
      ParallelDescriptor::Bcast(&coord, 1, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());

      // Step 4: send the periodicity flags
      ParallelDescriptor::Bcast(is_periodic, BL_SPACEDIM, MPI_IntraGroup_Broadcast_Rank, ParallelDescriptor::CommunicatorInter());
  }
  else
  {
    int box_index_type[BL_SPACEDIM];
    int smallEnd[BL_SPACEDIM];
    int bigEnd[BL_SPACEDIM];

    Real realBox_lo[BL_SPACEDIM];
    Real realBox_hi[BL_SPACEDIM];

    ParallelDescriptor::Bcast(box_index_type, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
    ParallelDescriptor::Bcast(smallEnd      , BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
    ParallelDescriptor::Bcast(bigEnd        , BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());

    ParallelDescriptor::Bcast(realBox_lo, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());
    ParallelDescriptor::Bcast(realBox_hi, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());

    ParallelDescriptor::Bcast(&coord, 1, 0, ParallelDescriptor::CommunicatorInter());

    ParallelDescriptor::Bcast(is_periodic, BL_SPACEDIM, 0, ParallelDescriptor::CommunicatorInter());

    // Now reconstruct all the parts
    IntVect smallEnd_IV(smallEnd);
    IntVect bigEnd_IV(bigEnd);
    IntVect box_index_type_IV(box_index_type);
    Box baseBox(smallEnd_IV, bigEnd_IV, box_index_type_IV);

    RealBox realBox;
    for (int n = 0; n < BL_SPACEDIM; n++) {
      realBox.setLo(n, realBox_lo[n]);
      realBox.setHi(n, realBox_hi[n]);
    }

    geom->define(baseBox, &realBox, coord, is_periodic);
  }
}
#endif
