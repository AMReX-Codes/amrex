//BL_COPYRIGHT_NOTICE

//
// $Id: Geometry.cpp,v 1.3 1997-12-11 23:27:00 lijewski Exp $
//

#include <Geometry.H>
#include <ParmParse.H>
#include <BoxArray.H>
#include <MultiFab.H>

//
// Static data members.
//
RealBox Geometry::prob_domain;

bool Geometry::is_periodic[BL_SPACEDIM];

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
