//BL_COPYRIGHT_NOTICE

//
// $Id: Geometry.cpp,v 1.19 1998-06-13 17:28:30 lijewski Exp $
//

#include <Geometry.H>
#include <ParmParse.H>
#include <BoxArray.H>

//
// The definition of static data members.
//
RealBox Geometry::prob_domain;

bool Geometry::is_periodic[BL_SPACEDIM];

list<Geometry::FPB> Geometry::m_Cache;

ostream&
operator<< (ostream&        os,
            const Geometry& g)
{
    os << (CoordSys&) g;
    os << g.prob_domain;
    os << g.domain;
    return os;
}

istream&
operator>> (istream&  is,
            Geometry& g)
{
    is >> (CoordSys&) g;
    is >> g.prob_domain;
    is >> g.domain;
    return is;
}

ostream&
operator << (ostream&               os,
	     const Geometry::PIRec& pir)
{
    os << "  From (Box "
       << pir.srcId
       << ") "
       << pir.srcBox
       << " to "
       << pir.dstBox;
    return os;
}

ostream&
operator<< (ostream&                 os,
	    const Geometry::PIRMMap& pirm)
{
    if (pirm.size() > 0)
    {
	Geometry::PIRMMap::const_iterator it = pirm.begin();
	Geometry::PIRMMap::key_type key = (*it).first;
	os << "Key: " << key << '\n';
	for ( ; it != pirm.end(); ++it )
	{
	    Geometry::PIRMMap::key_type key1 = (*it).first;
	    if (key != key1)
	    {
		key = key1;
		os << "Key: " << key << '\n';
	    }
	    os << (*it).second << '\n';
	}   
    }
    return os;
}

bool
Geometry::FPB::operator== (const FPB& rhs) const
{
    return m_ba == rhs.m_ba && m_ngrow == rhs.m_ngrow && m_noovlp == rhs.m_noovlp;
}

void
Geometry::FlushPIRMCache ()
{
    for (list<FPB>::iterator it = m_Cache.begin(); it != m_Cache.end(); ++it)
    {
        delete (*it).m_pirm;
    }
    m_Cache.clear();
}

Geometry::PIRMMap&
Geometry::buildPIRMMap (const BoxArray& grids,
                        int             nGrow,
                        bool            no_ovlp) const
{
    assert(isAnyPeriodic());
    //
    // Add new FPBs to the front of the cache.
    //
    m_Cache.push_front(FPB(new PIRMMap,grids,nGrow,no_ovlp));

    PIRMMap& pirm = *m_Cache.front().m_pirm;

    const int len = grids.length();
    //
    // Make a junk multifab to access its iterator.
    // Don't allocate any mem for it.
    //
    MultiFab mf(grids, 1, nGrow, Fab_noallocate);

    Array<IntVect> pshifts(27);
    //
    // Do only those I own.
    //
    for (ConstMultiFabIterator mfmfi(mf); mfmfi.isValid(); ++mfmfi)
    {
        Box dest = ::grow(mfmfi.validbox(), nGrow);

        assert(dest.ixType().cellCentered() || dest.ixType().nodeCentered());

        if (no_ovlp)
        {
            const Box& validbox = mfmfi.validbox();

            for (int idir = 0; idir < BL_SPACEDIM; idir++)
            {
                //
                // Shrink box if the +/- direction is not physical boundary.
                //
                if (validbox.smallEnd(idir) != Domain().smallEnd(idir))
                    dest.growLo(idir,-nGrow);
                if (validbox.bigEnd(idir) != Domain().bigEnd(idir))
                    dest.growHi(idir,-nGrow);
            }
        }

        bool doit = false;

        if (dest.ixType().cellCentered())
        {
            doit = !Domain().contains(dest);
        }
        else
        {
            doit = !::grow(::surroundingNodes(Domain()),-1).contains(dest);
        }

        if (doit)
        {
            for (int j = 0; j < len; j++)
            {
                periodicShift(dest, grids[j], pshifts);

                for (int iiv = 0; iiv < pshifts.length(); iiv++)
                {
                    Box shbox(grids[j]);
                    shbox.shift(pshifts[iiv]);
                    Box srcBox = dest & shbox;
                    assert(srcBox.ok());
                    //
                    // OK, we got an intersection.
                    //
                    Box dstBox = srcBox;
                    dstBox.shift(-pshifts[iiv]);

                    pirm.insert(PIRMMap::value_type(mfmfi.index(),
                                                    PIRec(j,dstBox,srcBox)));
                }
            }
        }
    }

    return pirm;
}

Geometry::PIRMMap&
Geometry::getPIRMMap (const BoxArray& grids,
                      int             nGrow,
                      bool            no_ovlp) const
{
    assert(isAnyPeriodic());
    //
    // Have we already made one with appropriate characteristics?
    //
    FPB fpb(grids,nGrow,no_ovlp);

    list<FPB>::iterator find_it = m_Cache.begin();

    for ( ; find_it != m_Cache.end(); ++find_it)
    {
        if (*find_it == fpb)
            break;
    }

    if (!(find_it == m_Cache.end()))
    {
        return *(*find_it).m_pirm;
    }
    else
    {
        return buildPIRMMap(grids,nGrow,no_ovlp);
    }
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                int       src_comp,
                                int       num_comp,
                                bool      no_ovlp) const
{
    if (!isAnyPeriodic())
        return;

    RunStats fpb_stats("fill_periodic_bndry");

    fpb_stats.start();

    PIRMMap& pirm = getPIRMMap(mf.boxArray(),mf.nGrow(),no_ovlp);

    MultiFabCopyDescriptor& mfcd = mf.theFPBmfcd(*this,
                                                 src_comp,
                                                 num_comp,
                                                 no_ovlp);
    //
    // Gather/scatter distributed data to (local) internal buffers.
    //
    mfcd.CollectData();

    const MultiFabId TheFPBMultiFabId = 0;

    typedef PIRMMap::iterator PIRMMapIt;
    //
    // Loop over my receiving fabs, copy periodic regions from buffered data.
    //
    for (MultiFabIterator mfi(mf); mfi.isValid(false); ++mfi)
    {
        pair<PIRMMapIt,PIRMMapIt> range = pirm.equal_range(mfi.index());
        //
	// For each PIRec on this fab box ...
        //
	for (PIRMMapIt p_it = range.first; p_it != range.second; ++p_it)
	{
	    const FillBoxId& fbid = (*p_it).second.fbid;

	    assert(fbid.box() == (*p_it).second.srcBox);

	    mfcd.FillFab(TheFPBMultiFabId, fbid, mfi(), (*p_it).second.dstBox);
        }
    }

    fpb_stats.end();
}

Geometry::Geometry () {}

Geometry::Geometry (const Box& dom)
{
    define(dom);
}

Geometry::Geometry (const Geometry& g)
{
    define(g.domain);
}

Geometry::~Geometry() {}

void
Geometry::define (const Box& dom)
{
    if (c_sys == undef)
        Setup();
    domain = dom;
    ok = true;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
    }
}

void
Geometry::Setup ()
{
    ParmParse pp("geometry");

    int coord;
    pp.get("coord_sys",coord);
    SetCoord( (CoordType) coord );

    Array<Real> prob_lo(BL_SPACEDIM);
    pp.getarr("prob_lo",prob_lo,0,BL_SPACEDIM);
    assert(prob_lo.length() == BL_SPACEDIM);
    Array<Real> prob_hi(BL_SPACEDIM);
    pp.getarr("prob_hi",prob_hi,0,BL_SPACEDIM);
    assert(prob_lo.length() == BL_SPACEDIM);
    prob_domain.setLo(prob_lo);
    prob_domain.setHi(prob_hi);
    //
    // Now get periodicity info.
    //
    D_EXPR( is_periodic[0]=0, is_periodic[1]=0, is_periodic[2]=0);
    if (pp.contains("period_0"))
    {
        is_periodic[0] = 1;
    }
#if BL_SPACEDIM>1
    if (pp.contains("period_1"))
    {
        is_periodic[1] = 1;
    }
#endif
#if BL_SPACEDIM>2
    if (pp.contains("period_2"))
    {
        is_periodic[2] = 1;
    }
#endif
}

void
Geometry::GetVolume (MultiFab&       vol,
                     const BoxArray& grds,
                     int             ngrow) const
{
    vol.define(grds,1,ngrow,Fab_noallocate);
    for (MultiFabIterator mfi(vol); mfi.isValid(); ++mfi)
    {
        Box gbx(grow(grds[mfi.index()],ngrow));
        vol.setFab(mfi.index(),CoordSys::GetVolume(gbx));
    }
}

#if (BL_SPACEDIM == 2)
void
Geometry::GetDLogA (MultiFab&       dloga,
                    const BoxArray& grds, 
                    int             dir,
                    int             ngrow) const
{
    dloga.define(grds,1,ngrow,Fab_noallocate);
    for (MultiFabIterator mfi(dloga); mfi.isValid(); ++mfi)
    {
        Box gbx(grow(grds[mfi.index()],ngrow));
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
    for (MultiFabIterator mfi(area); mfi.isValid(); ++mfi)
    {
        Box gbx(grow(grds[mfi.index()],ngrow));
        area.setFab(mfi.index(),CoordSys::GetFaceArea(gbx,dir));
    }
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
                        out.resize(out.length()+1); 
                        out[out.length()-1] = sh;
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
