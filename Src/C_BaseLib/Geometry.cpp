
#include <winstd.H>

#include <iostream>

#include <BoxArray.H>
#include <Geometry.H>
#include <ParmParse.H>
#include <MultiFab.H>
#include <FArrayBox.H>
//
// The definition of static data members.
//
int     Geometry::spherical_origin_fix;
RealBox Geometry::prob_domain;
bool    Geometry::is_periodic[BL_SPACEDIM];

Geometry::FPBMMap Geometry::m_FPBCache;

namespace
{
    bool verbose;
}

const int fpb_cache_max_size_def = 30; // -1 ==> no maximum size

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

std::ostream&
operator<< (std::ostream&          os,
            const Geometry::PIRec& pir)
{
    os << "mfi: "
       << pir.mfid
       << " from (Box "
       << pir.srcId
       << ") "
       << pir.srcBox
       << " to "
       << pir.dstBox;

    return os;
}

std::ostream&
operator<< (std::ostream&               os,
	    const Geometry::PIRMVector& pirm)
{
    for (Geometry::PIRMVector::const_iterator it = pirm.begin(), End = pirm.end();
         it != End;
         ++it)
    {
        os << *it << '\n';
    }
    return os;
}

Geometry::FPB::FPB ()
    :
    m_ngrow(-1),
    m_do_corners(false),
    m_reused(false)
{}

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
    m_reused(false)
{
    BL_ASSERT(ngrow >= 0);
    BL_ASSERT(domain.ok());
}

Geometry::FPB::FPB (const FPB& rhs)
    :
    m_cache(rhs.m_cache),
    m_commdata(rhs.m_commdata),
    m_pirm(rhs.m_pirm),
    m_ba(rhs.m_ba),
    m_dm(rhs.m_dm),
    m_domain(rhs.m_domain),
    m_ngrow(rhs.m_ngrow),
    m_do_corners(rhs.m_do_corners),
    m_reused(rhs.m_reused)
{}

Geometry::FPB::~FPB () {}

bool
Geometry::FPB::operator== (const FPB& rhs) const
{
    return
        m_ngrow      == rhs.m_ngrow      &&
        m_do_corners == rhs.m_do_corners &&
        m_domain     == rhs.m_domain     &&
        m_ba         == rhs.m_ba         &&
        m_dm         == rhs.m_dm;
}

bool
Geometry::FPB::operator!= (const FPB& rhs) const
{
    return !operator==(rhs);
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

int
Geometry::PIRMCacheSize ()
{
    return m_FPBCache.size();
}

void
Geometry::FlushPIRMCache ()
{
    if (ParallelDescriptor::IOProcessor() && !m_FPBCache.empty() && verbose)
    {
        int reused = 0;

        for (FPBMMapIter it = m_FPBCache.begin(), End = m_FPBCache.end(); it != End; ++it)
            if (it->second.m_reused)
                reused++;

        std::cout << "Geometry::PIRMCache.size() = " << m_FPBCache.size() << ", # reused = " << reused << '\n';
    }
    m_FPBCache.clear();
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                int       scomp,
                                int       ncomp,
                                bool      corners,
                                bool      local) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

    Box TheDomain = Domain();
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (mf.boxArray()[0].ixType()[n] == IndexType::NODE)
            TheDomain.surroundingNodes(n);

    if (TheDomain.contains(BoxLib::grow(mf.boxArray().minimalBox(), mf.nGrow()))) return;

    if ( local )
    {
        //
        // Do what you can with the FABs you own.  No parallelism allowed.
        //
        Array<IntVect> pshifts(27);

        for (MFIter mfidst(mf); mfidst.isValid(); ++mfidst)
        {
            const Box& dst = mf[mfidst].box();

            BL_ASSERT(dst == BoxLib::grow(mfidst.validbox(), mf.nGrow()));

            if (!TheDomain.contains(dst))
            {
                for (MFIter mfisrc(mf); mfisrc.isValid(); ++mfisrc)
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

                    periodicShift(dst, src, pshifts);

                    for (int i = 0, N = pshifts.size(); i < N; i++)
                    {
                        Box shftbox = src + pshifts[i];
                        Box dbx     = dst & shftbox;
                        Box sbx     = dbx - pshifts[i];

                        mf[mfidst].copy(mf[mfisrc], sbx, scomp, dbx, scomp, ncomp);
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

static
void
SumPeriodicBoundaryInnards (const MultiFab&         dstmf,
                            const MultiFab&         srcmf,
                            const Geometry&         geom,
                            MultiFabCopyDescriptor& mfcd,
                            const FabArrayId&       mfid,
                            Geometry::PIRMVector&   pirm,
                            int                     scomp,
                            int                     ncomp)
{
    //
    // Note that in the usual case dstmf == srcmf.
    //
    Array<IntVect> pshifts(27);

    Box TheDomain = geom.Domain();
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (srcmf.boxArray()[0].ixType()[n] == IndexType::NODE)
            TheDomain.surroundingNodes(n);

    if (TheDomain.contains(BoxLib::grow(srcmf.boxArray().minimalBox(), srcmf.nGrow()))) return;

    for (MFIter mfi(dstmf); mfi.isValid(); ++mfi)
    {
        const Box& dest = mfi.validbox();

        if (TheDomain.contains(BoxLib::grow(dest,srcmf.nGrow()))) continue;
        //
        // We may overlap with the periodic boundary.  Some other ghost
        // region(s) could be periodically-shiftable into our valid region.
        //
        const BoxArray& grids = srcmf.boxArray();

        for (int j = 0, N = grids.size(); j < N; j++)
        {
            const Box src = BoxLib::grow(grids[j],srcmf.nGrow());

            if (TheDomain.contains(src)) continue;

            geom.periodicShift(dest, src, pshifts);

            for (int i = 0, M = pshifts.size(); i < M; i++)
            {
                const Box shftbox = src + pshifts[i];
                const Box dbx     = dest & shftbox;
                const Box sbx     = dbx - pshifts[i];

                pirm.push_back(Geometry::PIRec(mfi.index(),j,sbx,dbx));

                pirm.back().fbid = mfcd.AddBox(mfid,
                                               sbx,
                                               0,
                                               j,
                                               scomp,
                                               0,
                                               ncomp,
                                               false);
            }
        }
    }
}

void
Geometry::SumPeriodicBoundary (MultiFab& mf,
                               int       scomp,
                               int       ncomp) const
{
    if (!isAnyPeriodic() || mf.nGrow() == 0 || mf.size() == 0) return;

#ifndef NDEBUG
    //
    // Don't let folks ask for more grow cells than they have valid region.
    //
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (isPeriodic(n))
            BL_ASSERT(mf.nGrow() <= Domain().length(n));
#endif

    PIRMVector             pirm;
    MultiFabCopyDescriptor mfcd;
    const FabArrayId       mfid = mfcd.RegisterFabArray(&mf);

    SumPeriodicBoundaryInnards(mf,mf,*this,mfcd,mfid,pirm,scomp,ncomp);

    int nrecv = pirm.size();
    ParallelDescriptor::ReduceIntMax(nrecv);
    if (nrecv == 0)
        //
        // There's no parallel work to do.
        //
        return;

    mfcd.CollectData();

    FArrayBox fab;

    for (Geometry::PIRMVector::const_iterator it = pirm.begin(), End = pirm.end();
         it != End;
         ++it)
    {
        BL_ASSERT(it->fbid.box() == it->srcBox);
        BL_ASSERT(it->srcBox.sameSize(it->dstBox));
        BL_ASSERT(mf.DistributionMap()[it->mfid] == ParallelDescriptor::MyProc());

        fab.resize(it->srcBox,ncomp);

        mfcd.FillFab(mfid, it->fbid, fab);

        mf[it->mfid].plus(fab, fab.box(), it->dstBox, 0, scomp, ncomp);
    }
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
    BL_ASSERT(srcmf.boxArray()[0].ixType() == dstmf.boxArray()[0].ixType());

#ifndef NDEBUG
    //
    // Don't let folks ask for more grow cells than they have valid region.
    //
    for (int n = 0; n < BL_SPACEDIM; n++)
        if (isPeriodic(n))
            BL_ASSERT(srcmf.nGrow() <= Domain().length(n));
#endif

    PIRMVector             pirm;
    MultiFabCopyDescriptor mfcd;
    const FabArrayId       mfid = mfcd.RegisterFabArray(const_cast<MultiFab*>(&srcmf));

    SumPeriodicBoundaryInnards(dstmf,srcmf,*this,mfcd,mfid,pirm,scomp,ncomp);

    int nrecv = pirm.size();
    ParallelDescriptor::ReduceIntMax(nrecv);
    if (nrecv == 0)
        //
        // There's no parallel work to do.
        //
        return;

    mfcd.CollectData();

    FArrayBox fab;

    for (Geometry::PIRMVector::const_iterator it = pirm.begin(), End = pirm.end();
         it != End;
         ++it)
    {
        BL_ASSERT(it->fbid.box() == it->srcBox);
        BL_ASSERT(it->srcBox.sameSize(it->dstBox));
        BL_ASSERT(dstmf.DistributionMap()[it->mfid] == ParallelDescriptor::MyProc());

        fab.resize(it->srcBox,ncomp);

        mfcd.FillFab(mfid, it->fbid, fab);

        dstmf[it->mfid].plus(fab, fab.box(), it->dstBox, 0, dcomp, ncomp);
    }
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

    Geometry::m_FPBCache.clear();
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
    verbose                        = false;
    Geometry::spherical_origin_fix = 0;
    Geometry::fpb_cache_max_size   = fpb_cache_max_size_def;

    D_EXPR(is_periodic[0]=0, is_periodic[1]=0, is_periodic[2]=0);

    pp.query("verbose",              verbose);
    pp.query("spherical_origin_fix", Geometry::spherical_origin_fix);
    pp.query("fpb_cache_max_size",   Geometry::fpb_cache_max_size);
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
    vol.define(grds,1,ngrow,Fab_noallocate);
    for (MFIter mfi(vol); mfi.isValid(); ++mfi)
    {
        Box gbx = BoxLib::grow(grds[mfi.index()],ngrow);
        vol.setFab(mfi.index(),CoordSys::GetVolume(gbx));
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
        Box gbx = BoxLib::grow(grds[mfi.index()],ngrow);
        dloga.setFab(mfi.index(),CoordSys::GetDLogA(gbx,dir));
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
    area.define(edge_boxes,1,ngrow,Fab_noallocate);
    for (MFIter mfi(area); mfi.isValid(); ++mfi)
    {
        Box gbx = BoxLib::grow(grds[mfi.index()],ngrow);
        area.setFab(mfi.index(),CoordSys::GetFaceArea(gbx,dir));
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
    Box locsrc(src);
    out.resize(0);

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
                    IntVect sh;
                    D_TERM(sh.setVal(0,ri*domain.length(0));,
                           sh.setVal(1,rj*domain.length(1));,
                           sh.setVal(2,rk*domain.length(2));)
                    out.resize(out.size()+1); 
                    out[out.size()-1] = sh;
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
