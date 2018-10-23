
#include <iostream>
#include <limits>

#include <AMReX_CoordSys.H>
#include <AMReX_COORDSYS_F.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParallelDescriptor.H>

#if (AMREX_SPACEDIM==2)
namespace {
    constexpr double RZFACTOR = 2*3.14159265358979323846264338327950288;
}
#endif

//
// A handy macro.
//
#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

namespace amrex {

//
// The definition of static data members.
//
CoordSys::CoordType CoordSys::c_sys = CoordSys::undef;

Real CoordSys::offset[AMREX_SPACEDIM];

CoordSys::CoordType
CoordSys::Coord ()
{
    return c_sys;
}

int
CoordSys::CoordInt ()
{
    switch (c_sys)
    {
        case undef:
            return -1;
        case cartesian:
            return 0;
        case RZ:
            return 1;
        case SPHERICAL:
            return 2;
        default:
            return -1;
    }
}

void
CoordSys::SetOffset (const Real* x_lo)
{
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        offset[k] = x_lo[k];
    }
}

CoordSys::CoordSys ()
{
    AMREX_D_TERM(dx[0]=0;,dx[1]=0;,dx[2]=0;)
    AMREX_D_TERM(inv_dx[0]=std::numeric_limits<Real>::infinity();,
           inv_dx[1]=std::numeric_limits<Real>::infinity();,
           inv_dx[2]=std::numeric_limits<Real>::infinity();)
    ok = false;
}

void
CoordSys::define (const Real* cell_dx)
{
    BL_ASSERT(c_sys != undef);
    ok = true;
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        dx[k] = cell_dx[k];
	inv_dx[k] = 1.0/dx[k];
    }
}

CoordSys::CoordSys (const Real* cell_dx)
{
    define(cell_dx);
}

void
CoordSys::CellCenter (const IntVect& point,
                      Real*          loc) const
{
    BL_ASSERT(ok);
    BL_ASSERT(loc != 0);
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        loc[k] = offset[k] + dx[k]*(0.5+ (Real)point[k]);
    }
}

void
CoordSys::CellCenter (const IntVect& point,
                      Vector<Real>&   loc) const
{
    BL_ASSERT(ok);
    loc.resize(AMREX_SPACEDIM);
    CellCenter(point, loc.dataPtr());
}

void
CoordSys::LoFace (const IntVect& point,
                  int            dir,
                  Real*          loc) const
{
    BL_ASSERT(ok);
    BL_ASSERT(loc != 0);
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        Real off = (k == dir) ? 0.0 : 0.5;
        loc[k] = offset[k] + dx[k]*(off + (Real)point[k]);
    }
}

void
CoordSys::LoFace (const IntVect& point,
                  int            dir,
                  Vector<Real>&   loc) const
{
    loc.resize(AMREX_SPACEDIM);
    LoFace(point,dir, loc.dataPtr());
}

void
CoordSys::HiFace (const IntVect& point,
                  int            dir,
                  Real*          loc) const
{
    BL_ASSERT(ok);
    BL_ASSERT(loc != 0);
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        Real off = (k == dir) ? 1.0 : 0.5;
        loc[k] = offset[k] + dx[k]*(off + (Real)point[k]);
    }
}

void
CoordSys::HiFace (const IntVect& point,
                  int            dir,
                  Vector<Real>&   loc) const
{
    loc.resize(AMREX_SPACEDIM);
    HiFace(point,dir, loc.dataPtr());
}

void
CoordSys::LoNode (const IntVect& point,
                  Real*          loc) const
{
    BL_ASSERT(ok);
    BL_ASSERT(loc != 0);
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        loc[k] = offset[k] + dx[k]*point[k];
    }
}

void
CoordSys::LoNode (const IntVect& point,
                  Vector<Real>&   loc) const
{
    loc.resize(AMREX_SPACEDIM);
    LoNode(point, loc.dataPtr());
}

void
CoordSys::HiNode (const IntVect& point,
                  Real*          loc) const
{
    BL_ASSERT(ok);
    BL_ASSERT(loc != 0);
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        loc[k] = offset[k] + dx[k]*(point[k] + 1);
    }
}

void
CoordSys::HiNode (const IntVect& point,
                  Vector<Real>&   loc) const
{
    loc.resize(AMREX_SPACEDIM);
    HiNode(point, loc.dataPtr());
}

IntVect
CoordSys::CellIndex (const Real* point) const
{
    BL_ASSERT(ok);
    BL_ASSERT(point != 0);
    IntVect ix;
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

IntVect
CoordSys::LowerIndex (const Real* point) const
{
    BL_ASSERT(ok);
    BL_ASSERT(point != 0);
    IntVect ix;
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        ix[k] = (int) ((point[k]-offset[k])/dx[k]);
    }
    return ix;
}

IntVect
CoordSys::UpperIndex(const Real* point) const
{
    BL_ASSERT(ok);
    BL_ASSERT(point != 0);
    IntVect ix;
    for (int k = 0; k < AMREX_SPACEDIM; k++)
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
    BL_ASSERT(ok);
    BL_ASSERT(region.cellCentered());
    vol.resize(region,1);
    DEF_LIMITS(vol,vol_dat,vlo,vhi);
    int coord = (int) c_sys;
    const int* rlo = region.loVect();
    const int* rhi = region.hiVect();
    amrex_setvol(AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),vol_dat,AMREX_ARLIM(vlo),AMREX_ARLIM(vhi),offset,dx,&coord);
}

void
CoordSys::SetVolume (FArrayBox& vol,
                     const Box& region) const 
{
    BL_ASSERT(ok);
    BL_ASSERT(region.cellCentered());
    DEF_LIMITS(vol,vol_dat,vlo,vhi);
    int coord = (int) c_sys;
    const int* rlo = region.loVect();
    const int* rhi = region.hiVect();
    amrex_setvol(AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),vol_dat,AMREX_ARLIM(vlo),AMREX_ARLIM(vhi),offset,dx,&coord);
}

void
CoordSys::GetDLogA (FArrayBox& dloga,
                    const Box& region,
                    int        dir) const
{
    BL_ASSERT(ok);
    BL_ASSERT(region.cellCentered());
    dloga.resize(region,1);
    DEF_LIMITS(dloga,dloga_dat,dlo,dhi);
    int coord = (int) c_sys;
    amrex_setdloga(dloga_dat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),offset,dx,&dir,&coord);
}

FArrayBox*
CoordSys::GetDLogA (const Box& region,
                    int        dir) const
{
    FArrayBox* dloga = new FArrayBox();
    GetDLogA(*dloga,region,dir);
    return dloga;
}

void
CoordSys::SetDLogA (FArrayBox& dloga,
                    const Box& region,
                    int        dir) const
{
    BL_ASSERT(ok);
    BL_ASSERT(region.cellCentered());
    DEF_LIMITS(dloga,dloga_dat,dlo,dhi);
    int coord = (int) c_sys;
    amrex_setdloga(dloga_dat,AMREX_ARLIM(dlo),AMREX_ARLIM(dhi),offset,dx,&dir,&coord);
}

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
    BL_ASSERT(ok);
    BL_ASSERT(region.cellCentered());
    Box reg(region);
    reg.surroundingNodes(dir);
    area.resize(reg,1);
    DEF_LIMITS(area,area_dat,lo,hi)
    int coord = (int) c_sys;
    const int* rlo = reg.loVect();
    const int* rhi = reg.hiVect();
    amrex_setarea(AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),area_dat,AMREX_ARLIM(lo),AMREX_ARLIM(hi),offset,dx,&dir,&coord);
}

void
CoordSys::SetFaceArea (FArrayBox& area, 
                       const Box& region,
                       int        dir) const
{
    BL_ASSERT(ok);
    DEF_LIMITS(area,area_dat,lo,hi)
    int coord = (int) c_sys;
    const int* rlo = region.loVect();
    const int* rhi = region.hiVect();
    amrex_setarea(AMREX_ARLIM(rlo),AMREX_ARLIM(rhi),area_dat,AMREX_ARLIM(lo),AMREX_ARLIM(hi),offset,dx,&dir,&coord);
}

void
CoordSys::GetEdgeLoc (Vector<Real>& loc, 
                      const Box&   region,
                      int          dir) const 
{
    BL_ASSERT(ok);
    BL_ASSERT(region.cellCentered());
    const int* lo = region.loVect();
    const int* hi = region.hiVect();
    int len       = hi[dir] - lo[dir] + 2;
    Real off      = offset[dir] + dx[dir]*lo[dir];
    loc.resize(len);
    for (int i = 0; i < len; i++)
    {
        loc[i] = off + dx[dir]*i;
    }
}

void
CoordSys::GetCellLoc (Vector<Real>& loc, 
                      const Box&   region,
                      int          dir) const
{
    BL_ASSERT(ok);
    BL_ASSERT(region.cellCentered());
    const int* lo = region.loVect();
    const int* hi = region.hiVect();
    int len       = hi[dir] - lo[dir] + 1;
    Real off = offset[dir] + dx[dir]*(0.5 + (Real)lo[dir]);
    loc.resize(len);
    for (int i = 0; i < len; i++)
    {
        loc[i] = off + dx[dir]*i;
    }
}

void
CoordSys::GetEdgeVolCoord (Vector<Real>& vc,
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
#if (AMREX_SPACEDIM == 2)
    if (dir == 0 && c_sys == RZ)
    {
        int len = vc.size();
        for (int i = 0; i < len; i++)
        {
            Real r = vc[i];
            vc[i] = 0.5*r*r;
        }
    }
    else
    {
        if(dir == 0 && c_sys == SPHERICAL) {
            int len = vc.size();
            int i;
            for (i = 0; i < len; i++) {
                Real r = vc[i];
                vc[i] = 0.3*r*r*r;
            }
        }
    }
#endif    
}

void
CoordSys::GetCellVolCoord (Vector<Real>& vc,
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
#if (AMREX_SPACEDIM == 2)
    if (dir == 0 && c_sys == RZ)
    {
        int len = vc.size();
        for (int i = 0; i < len; i++)
        {
            Real r = vc[i];
            vc[i] = 0.5*r*r;
        }
    }
    else
    {
        if(dir == 0 && c_sys == SPHERICAL) {
            int len = vc.size();
            int i;
            for (i = 0; i < len; i++) {
                Real r = vc[i];
                vc[i] = 0.3*r*r*r;
            }
        }
    }
#endif    
}

std::ostream&
operator<< (std::ostream&   os,
            const CoordSys& c)
{
    os << '(' << (int) c.Coord() << ' ';
    os << AMREX_D_TERM( '(' << c.Offset(0) , <<
                  ',' << c.Offset(1) , <<
                  ',' << c.Offset(2))  << ')';
    os << AMREX_D_TERM( '(' << c.CellSize(0) , <<
                  ',' << c.CellSize(1) , <<
                  ',' << c.CellSize(2))  << ')';
    os << ' ' << int(c.ok) << ")\n";
    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            CoordSys&     c)
{
    int coord;
    is.ignore(BL_IGNORE_MAX, '(') >> coord;
    c.c_sys = (CoordSys::CoordType) coord;
    AMREX_D_EXPR(is.ignore(BL_IGNORE_MAX, '(') >> c.offset[0],
           is.ignore(BL_IGNORE_MAX, ',') >> c.offset[1],
           is.ignore(BL_IGNORE_MAX, ',') >> c.offset[2]);
    is.ignore(BL_IGNORE_MAX, ')');
    AMREX_D_EXPR(is.ignore(BL_IGNORE_MAX, '(') >> c.dx[0],
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
    Real xhi[AMREX_SPACEDIM];
    Real xlo[AMREX_SPACEDIM];
    HiNode(point,xhi);
    LoNode(point,xlo);
    return Volume(xlo,xhi);
}

Real 
CoordSys::Volume (const Real xlo[AMREX_SPACEDIM], 
                  const Real xhi[AMREX_SPACEDIM]) const
{
    switch (c_sys)
    {
    case cartesian:
        return AMREX_D_TERM((xhi[0]-xlo[0]),
                      *(xhi[1]-xlo[1]),
                      *(xhi[2]-xlo[2]));
#if (AMREX_SPACEDIM==2)
    case RZ:
        return (0.5*RZFACTOR)*(xhi[1]-xlo[1])*(xhi[0]*xhi[0]-xlo[0]*xlo[0]);
#endif
    default:
        BL_ASSERT(0);
    }
    return 0;
}                      

Real
CoordSys::AreaLo (const IntVect& point,
                  int            dir) const
{
#if (AMREX_SPACEDIM==2)
    Real xlo[AMREX_SPACEDIM];
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
        BL_ASSERT(0);
    }
#endif
#if (AMREX_SPACEDIM==3)
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
#if (AMREX_SPACEDIM==2)
    Real xhi[AMREX_SPACEDIM];
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
        BL_ASSERT(0);
    }
#endif
#if (AMREX_SPACEDIM==3)
    switch (dir)
    {
    case 0: return dx[1]*dx[2];
    case 1: return dx[0]*dx[2];
    case 2: return dx[1]*dx[0];
    }
#endif
    return 0;
}


#ifdef BL_USE_MPI

void
CoordSys::BroadcastCoordSys (CoordSys &cSys, int fromProc, MPI_Comm comm)
{
  bool bcastSource(ParallelDescriptor::MyProc() == fromProc);
  CoordSys::BroadcastCoordSys(cSys, fromProc, comm, bcastSource);
}



void
CoordSys::BroadcastCoordSys (CoordSys &cSys, int fromProc, MPI_Comm comm, bool bcastSource)
{
  int coord, iOK;

  if(bcastSource) {  // ---- initialize the source data
    coord = cSys.CoordInt();
    iOK = cSys.ok;
  }

  ParallelDescriptor::Bcast(&coord, 1, fromProc, comm);
  ParallelDescriptor::Bcast(cSys.offset, AMREX_SPACEDIM, fromProc, comm);
  ParallelDescriptor::Bcast(cSys.dx, AMREX_SPACEDIM, fromProc, comm);
  ParallelDescriptor::Bcast(cSys.inv_dx, AMREX_SPACEDIM, fromProc, comm);
  ParallelDescriptor::Bcast(&iOK, 1, fromProc, comm);

  if( ! bcastSource) {
    cSys.c_sys = (CoordSys::CoordType) coord;
    cSys.ok = iOK;
  }
}
#endif

}
