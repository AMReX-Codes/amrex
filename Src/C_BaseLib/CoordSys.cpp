//BL_COPYRIGHT_NOTICE

//
// $Id: CoordSys.cpp,v 1.6 1998-02-17 23:02:16 car Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <cmath>
#else
#include <math.h>
#endif

#include <Misc.H>
#include <CoordSys.H>
#include <COORDSYS_F.H>
#include <FArrayBox.H>

#if (BL_SPACEDIM==2)
const double RZFACTOR = 2*M_PI;
#endif

CoordSys::CoordType CoordSys::c_sys = CoordSys::undef;

Real CoordSys::offset[BL_SPACEDIM];

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

IntVect
CoordSys::CellIndex (const Real* point) const
{
    assert(ok);
    assert(point != 0);
    IntVect ix;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

IntVect
CoordSys::LowerIndex (const Real* point) const
{
    assert(ok);
    assert(point != 0);
    IntVect ix;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

IntVect
CoordSys::UpperIndex(const Real* point) const
{
    assert(ok);
    assert(point != 0);
    IntVect ix;
    for (int k = 0; k < BL_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

FArrayBox*
CoordSys::GetVolume (const Box& region) const 
{
    FArrayBox* vol = new FArrayBox();
    GetVolume(*vol,region);
    return vol;
}

void
CoordSys::GetVolume (FArrayBox& vol,
                     const Box& region) const 
{
    assert(ok);
    assert(region.cellCentered());
    vol.resize(region,1);
    DEF_LIMITS(vol,vol_dat,vlo,vhi);
    int coord = (int) c_sys;
    FORT_SETVOL(vol_dat,ARLIM(vlo),ARLIM(vhi),offset,dx,&coord);
}

#if (BL_SPACEDIM == 2)
FArrayBox*
CoordSys::GetDLogA (const Box& region,
                    int        dir) const
{
    FArrayBox* dloga = new FArrayBox();
    GetDLogA(*dloga,region,dir);
    return dloga;
}

void
CoordSys::GetDLogA (FArrayBox& dloga,
                    const Box& region,
                    int        dir) const
{
    assert(ok);
    assert(region.cellCentered());
    dloga.resize(region,1);
    DEF_LIMITS(dloga,dloga_dat,dlo,dhi);
    int coord = (int) c_sys;
    FORT_SETDLOGA(dloga_dat,ARLIM(dlo),ARLIM(dhi),offset,dx,&dir,&coord);
}
#endif /*BL_SPACEDIM == 2*/

FArrayBox*
CoordSys::GetFaceArea (const Box& region,
                       int        dir) const 
{
    FArrayBox* area = new FArrayBox();
    GetFaceArea(*area,region,dir);
    return area;
}

void
CoordSys::GetFaceArea (FArrayBox& area, 
                       const Box& region,
                       int        dir) const
{
    assert(ok);
    assert(region.cellCentered());
    Box reg(region);
    reg.surroundingNodes(dir);
    area.resize(reg,1);
    DEF_LIMITS(area,area_dat,lo,hi)
    int coord = (int) c_sys;
    FORT_SETAREA(area_dat,ARLIM(lo),ARLIM(hi),offset,dx,&dir,&coord);
}

void
CoordSys::GetEdgeLoc (Array<Real>& loc, 
                      const Box&   region,
                      int          dir) const 
{
    assert(ok);
    assert(region.cellCentered());
    const int* lo = region.loVect();
    const int* hi = region.hiVect();
    int len = hi[dir] - lo[dir] + 2;
    loc.resize(len);
    Real off = offset[dir] + dx[dir]*lo[dir];
    for (int i = 0; i < len; i++)
    {
        loc[i] = off + dx[dir]*i;
    }
}

void
CoordSys::GetCellLoc (Array<Real>& loc, 
                      const Box&   region,
                      int          dir) const
{
    assert(ok);
    assert(region.cellCentered());
    const int* lo = region.loVect();
    const int* hi = region.hiVect();
    int len = hi[dir] - lo[dir] + 1;
    loc.resize(len);
    Real off = offset[dir] + dx[dir]*(0.5 + (Real)lo[dir]);
    for (int i = 0; i < len; i++)
    {
        loc[i] = off + dx[dir]*i;
    }
}

void
CoordSys::GetEdgeVolCoord (Array<Real>& vc,
                           const Box&   region,
                           int          dir) const
{
    //
    // In cartesian and Z direction of RZ volume coordinates
    // are idential to physical distance from axis.
    //
    GetEdgeLoc(vc,region,dir);
    //
    // In R direction of RZ, vol coord = (r^2)/2
    //
#if (BL_SPACEDIM == 2)
    if (dir == 0 && c_sys == RZ)
    {
        int len = vc.length();
        for (int i = 0; i < len; i++)
        {
            Real r = vc[i];
            vc[i] = 0.5*r*r;
        }
    }
#endif    
}

void
CoordSys::GetCellVolCoord (Array<Real>& vc,
                           const Box&   region,
                           int          dir) const
{
    //
    // In cartesian and Z direction of RZ volume coordinates
    // are idential to physical distance from axis.
    //
    GetCellLoc(vc,region,dir);
    //
    // In R direction of RZ, vol coord = (r^2)/2.
    //
#if (BL_SPACEDIM == 2)
    if (dir == 0 && c_sys == RZ)
    {
        int len = vc.length();
        for (int i = 0; i < len; i++)
        {
            Real r = vc[i];
            vc[i] = 0.5*r*r;
        }
    }
#endif    
}

ostream&
operator<< (ostream&        os,
            const CoordSys& c)
{
    os << '(' << (int) c.c_sys << ' ';
    os << D_TERM( '(' << c.offset[0] , <<
                  ',' << c.offset[1] , <<
                  ',' << c.offset[2])  << ')';
    os << D_TERM( '(' << c.dx[0] , <<
                  ',' << c.dx[1] , <<
                  ',' << c.dx[2])  << ')';
    os << ' ' << int(c.ok) << ")\n";
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

istream&
operator>> (istream&  is,
            CoordSys& c)
{
    int coord;
    is.ignore(BL_IGNORE_MAX, '(') >> coord;
    c.c_sys = (CoordSys::CoordType) coord;
    D_EXPR(is.ignore(BL_IGNORE_MAX, '(') >> c.offset[0],
           is.ignore(BL_IGNORE_MAX, ',') >> c.offset[1],
           is.ignore(BL_IGNORE_MAX, ',') >> c.offset[2]);
    is.ignore(BL_IGNORE_MAX, ')');
    D_EXPR(is.ignore(BL_IGNORE_MAX, '(') >> c.dx[0],
           is.ignore(BL_IGNORE_MAX, ',') >> c.dx[1],
           is.ignore(BL_IGNORE_MAX, ',') >> c.dx[2]);
    is.ignore(BL_IGNORE_MAX, ')');
    int tmp;
    is >> tmp;
    c.ok = tmp?true:false;
    is.ignore(BL_IGNORE_MAX, '\n');
    return is;
}

Real
CoordSys::Volume (const IntVect& point) const
{
    Real xhi[BL_SPACEDIM];
    Real xlo[BL_SPACEDIM];
    HiNode(point,xhi);
    LoNode(point,xlo);
    return Volume(xlo,xhi);
}

Real 
CoordSys::Volume (const Real xlo[BL_SPACEDIM], 
                  const Real xhi[BL_SPACEDIM]) const
{
    switch (c_sys)
    {
    case cartesian:
        return (xhi[0]-xlo[0])
#if (BL_SPACEDIM>=2)                       
            *(xhi[1]-xlo[1])
#endif
#if (BL_SPACEDIM>=3)                       
            *(xhi[2]-xlo[2])
#endif
            ;
#if (BL_SPACEDIM==2)
    case RZ:
        return (0.5*RZFACTOR)*(xhi[1]-xlo[1])*(xhi[0]*xhi[0]-xlo[0]*xlo[0]);
#endif
    default:
        assert(0);
    }
    return 0;
}                      

Real
CoordSys::AreaLo (const IntVect& point,
                  int            dir) const
{
#if (BL_SPACEDIM==2)
    Real xlo[BL_SPACEDIM];
    switch (c_sys)
    {
    case cartesian:
        switch (dir)
        {
        case 0: return dx[1];
        case 1: return dx[0];
        }
    case RZ:
        LoNode(point,xlo);
        switch (dir)
        {
        case 0: return RZFACTOR*dx[1]*xlo[0];
        case 1: return ((xlo[0]+dx[0])*(xlo[0]+dx[0])-xlo[0]*xlo[0])*(0.5*RZFACTOR);
        }
    default:
        assert(0);
    }
#endif
#if (BL_SPACEDIM==3)
    switch (dir)
    {
    case 0: return dx[1]*dx[2];
    case 1: return dx[0]*dx[2];
    case 2: return dx[1]*dx[0];
    }
#endif
    return 0;
}

Real
CoordSys::AreaHi (const IntVect& point,
                  int            dir) const
{
#if (BL_SPACEDIM==2)
    Real xhi[BL_SPACEDIM];
    switch (c_sys)
    {
    case cartesian:
        switch (dir)
        {
        case 0: return dx[1];
        case 1: return dx[0];
        }
    case RZ:
        HiNode(point,xhi);
        switch (dir)
        {
        case 0: return RZFACTOR*dx[1]*xhi[0];
        case 1: return (xhi[0]*xhi[0]-(xhi[0]-dx[0])*(xhi[0]-dx[0]))*(RZFACTOR*0.5);
        }
    default:
        assert(0);
    }
#endif
#if (BL_SPACEDIM==3)
    switch (dir)
    {
    case 0: return dx[1]*dx[2];
    case 1: return dx[0]*dx[2];
    case 2: return dx[1]*dx[0];
    }
#endif
    return 0;
}
