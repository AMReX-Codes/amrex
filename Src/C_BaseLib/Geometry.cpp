//BL_COPYRIGHT_NOTICE

//
// $Id: Geometry.cpp,v 1.4 1998-04-22 23:39:32 marc Exp $
//

#include <Geometry.H>
#include <ParmParse.H>
#include <BoxArray.H>

//
// Static data members.
//
RealBox Geometry::prob_domain;

bool Geometry::is_periodic[BL_SPACEDIM];

std::ostream&
operator << (std::ostream&          os,
	     const Geometry::PIRec& pir)
{
    os << "  From (Box " << pir.srcId << ") " << pir.srcBox << " to " << pir.dstBox;
    return os;
}

std::ostream&
operator << (std::ostream&            os,
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

#ifdef __GNUG__
bool
Geometry::PIR::operator == (const Geometry::PIR& rhs) const
{
    return ( (srcno   == rhs.srcno ) &&
	     (destno  == rhs.destno )&&
	     (srcbox  == rhs.srcbox ) &&
	     (destbox == rhs.destbox ) ); 
}
#endif

Geometry::PIRMMap
Geometry:: computePIRMMapForMultiFab(const BoxArray& grids,
				     int             nGrow) const
{
    // Build a MultiMap of <i,PIRec> pairs, where i is the index of
    // a "local" box in grids, and the PIRec contains srcBox/dstBox pairs
    // mapping valid data to grow regions over periodic boundaries
    PIRMMap pirmmap;
    if( ! isAnyPeriodic() ) return pirmmap;

    const Box& domain = Domain();
    int len = grids.length();

    // Make a junk multifab to access its iterator, dont allocate any mem for it
    int nComp = 1;
    MultiFab mf(grids, nComp, nGrow, Fab_noallocate);

    // Do only those I own
    for (ConstMultiFabIterator mfmfi(mf); mfmfi.isValid(); ++mfmfi)
    {
	Box dest = Box(mfmfi.validbox()).grow(nGrow);
	if( ! domain.contains(dest) )
	{
	    for( int j=0; j<len; j++ )
	    {
		Box src = grids[j];
		Array<IntVect> pshifts;
		periodicShift( dest, src, pshifts );
		for( int iiv=0; iiv<pshifts.length(); iiv++ )
		{
		    IntVect iv = pshifts[iiv];
		    Box shbox( src );
		    D_TERM( shbox.shift(0,iv[0]);,
			    shbox.shift(1,iv[1]);,
			    shbox.shift(2,iv[2]); );
		    
		    Box intbox = dest & shbox;
		    assert( intbox.ok() );
		    // ok, we got an intersection
		    Box srcBox = intbox;
		    D_TERM( intbox.shift(0,-iv[0]);,
			    intbox.shift(1,-iv[1]);,
			    intbox.shift(2,-iv[2]); );
		    Box dstBox = intbox;
		    pirmmap.insert(PIRMMap::value_type(mfmfi.index(),
						       PIRec(j,dstBox,srcBox)));
		}
	    }
	}
    }

    return pirmmap;
}

void
Geometry::FillPeriodicFabArray (FabArray<Real,FArrayBox>& fa,
				PIRMMap&                  pirm,
				int                       sComp,
				int                       nComp) const
{
    // Assumes PIRec MultiMap built correctly (i.e. each entry indexed on
    // local fab id, contains srcBox/dstBox pairs to copy valid data from
    // other fabs in the array).  Assumes box-pairs constructed so that
    // all data is "fillable" from the valid region of "fa".
    FabArrayCopyDescriptor<Real,FArrayBox> facd;
    FabArrayId faid = facd.RegisterFabArray(&fa);
    
    typedef PIRMMap::iterator PIRMMapIt;
    BoxList unfilledBoxes((*pirm.begin()).second.srcBox.ixType());

    // Register boxes in copy decriptor (should be no unfilled boxes when finished)
    for (PIRMMapIt p_it = pirm.begin(); p_it != pirm.end(); ++p_it)
    {
	(*p_it).second.fbid = 
	    facd.AddBox(faid, (*p_it).second.srcBox, unfilledBoxes,
			(*p_it).second.srcId, sComp, sComp, nComp);
    }
    assert(unfilledBoxes.length() == 0);

    // Gather/scatter distributed data to (local) internal buffers
    facd.CollectData();

    // Loop over my receiving fabs, copy periodic regions from buffered data
    for (FabArrayIterator<Real,FArrayBox> fai(fa); fai.isValid(false); ++fai)
    {
        std::pair<PIRMMapIt,PIRMMapIt> range = pirm.equal_range(fai.index());

	// For each PIRec on this fab box...
	for (PIRMMapIt p_it = range.first; p_it != range.second; ++p_it)
	{
	    const FillBoxId& fbid = (*p_it).second.fbid;
	    assert(fbid.box() == (*p_it).second.srcBox);
	    FArrayBox overlapFab(fbid.box(), nComp);

	    // Fill fab with buffered data, copy into destination
	    facd.FillFab(faid, fbid, overlapFab);
	    fai().copy(overlapFab, overlapFab.box(), sComp,
		       (*p_it).second.dstBox, sComp, nComp);
        }
    }
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf) const
{
    // Get list of intersection boxes
    PIRMMap pirm = computePIRMMapForMultiFab(mf.boxArray(), mf.nGrow());
    const int sComp = 0;
    const int nComp = mf.nComp();

    // Fill intersection list
    FillPeriodicFabArray(mf, pirm, sComp, nComp);
}

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
