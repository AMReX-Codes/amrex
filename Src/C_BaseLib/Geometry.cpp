
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
int Geometry::spherical_origin_fix = 0;

RealBox Geometry::prob_domain;

bool Geometry::is_periodic[BL_SPACEDIM];

Geometry::FPBMMap Geometry::m_FPBCache;

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
operator<< (std::ostream&            os,
	    const Geometry::PIRMMap& pirm)
{
    for (int i = 0; i < pirm.size(); i++)
        os << pirm[i] << '\n';
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

const RealBox&
Geometry::ProbDomain ()
{
    return prob_domain;
}

void
Geometry::ProbDomain (const RealBox& rb)
{
    prob_domain = rb;
}

const Box&
Geometry::Domain () const
{
    return domain;
}

void
Geometry::Domain (const Box& bx)
{
    domain = bx;
}

bool
Geometry::isPeriodic (int dir)
{
    return is_periodic[dir] != 0;
}

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

int
Geometry::period (int dir) const
{
    BL_ASSERT(is_periodic[dir]);
    return domain.length(dir);
}

const Real*
Geometry::ProbLo ()
{
    return prob_domain.lo();
}

const Real*
Geometry::ProbHi ()
{
    return prob_domain.hi();
}

Real
Geometry::ProbLo (int dir)
{
    return prob_domain.lo(dir);
}

Real
Geometry::ProbHi (int dir)
{
    return prob_domain.hi(dir);
}

const Real*
Geometry::ProbLength ()
{
    return prob_domain.length();
}

Real
Geometry::ProbLength (int dir)
{
    return prob_domain.length(dir);
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                bool do_corners,
                                bool local) const
{
    FillPeriodicBoundary(mf,0,mf.nComp(),do_corners,local);
}

int
Geometry::PIRMCacheSize ()
{
    return m_FPBCache.size();
}

void
Geometry::FlushPIRMCache ()
{
    if (ParallelDescriptor::IOProcessor() && m_FPBCache.size())
        std::cout << "Geometry::PIRMCacheSize() = " << m_FPBCache.size() << '\n';
    m_FPBCache.clear();
}

void
Geometry::FillPeriodicBoundary (MultiFab& mf,
                                int       scomp,
                                int       ncomp,
                                bool      corners,
                                bool      local) const
{
    BL_PROFILE(BL_PROFILE_THIS_NAME() + "::FillPeriodicBoundary(mf)");

    if (!isAnyPeriodic()) return;

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

            Box TheDomain = Domain();
            for (int n = 0; n < BL_SPACEDIM; n++)
                if (dst.ixType()[n] == IndexType::NODE)
                    TheDomain.surroundingNodes(n);

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

                    for (int i = 0; i < pshifts.size(); i++)
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

Geometry::Geometry () {}

Geometry::Geometry (const Box&     dom,
                    const RealBox* rb,
                    int            coord)
{
    define(dom,rb,coord);
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
                  int            coord)
{
    if (c_sys == undef)
        Setup(rb,coord);
    domain = dom;
    ok     = true;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
    }
    if (spherical_origin_fix == 1)
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
Geometry::Setup (const RealBox* rb, int coord)
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

    spherical_origin_fix = 0;
    pp.query("spherical_origin_fix",spherical_origin_fix);
    //
    // Now get periodicity info.
    //
    D_EXPR(is_periodic[0]=0, is_periodic[1]=0, is_periodic[2]=0);

    Array<int> is_per(BL_SPACEDIM);
    pp.queryarr("is_periodic",is_per,0,BL_SPACEDIM);
    for (int n = 0; n < BL_SPACEDIM; n++)  
      is_periodic[n] = is_per[n];
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
