

#include <iostream>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_SPACE.H>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

//
// The definition of some static data members.
//
int     Geometry::spherical_origin_fix = 0;
RealBox Geometry::prob_domain;
bool    Geometry::is_periodic[BL_SPACEDIM] = {AMREX_D_DECL(0,0,0)};

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

Geometry::Geometry () {}

Geometry::Geometry (const Box&     dom,
                    const RealBox* rb,
                    int            coord,
                    int*           is_per)
{
    define(dom,rb,coord,is_per);
}

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
	inv_dx[k] = 1.0/dx[k];
    }
    if (Geometry::spherical_origin_fix == 1)
    {
	if (c_sys == SPHERICAL && prob_domain.lo(0) == 0 && BL_SPACEDIM > 1)
        {
            prob_domain.setLo(0,2*dx[0]);

            for (int k = 0; k < BL_SPACEDIM; k++)
            {
                dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
		inv_dx[k] = 1.0/dx[k];
            }
	}
    }

#ifdef AMREX_USE_DEVICE
    init_device();
    set_device();
#endif
}

void
Geometry::Finalize ()
{
    c_sys = undef;
}

void
Geometry::Setup (const RealBox* rb, int coord, int* isper)
{
#ifdef _OPENMP
    BL_ASSERT(!omp_in_parallel());
#endif

    static bool first = true;
    if (first) {
	first = false;
    } else {
	return;
    }

    ParmParse pp("geometry");

    if (coord >=0 && coord <= 2) {
        SetCoord( (CoordType) coord );        
    } else {
        coord = 0;  // default is Cartesian coordinates
        pp.query("coord_sys",coord);
        SetCoord( (CoordType) coord );        
    }

    if (rb == nullptr) {
        Vector<Real> prob_lo(BL_SPACEDIM);
        Vector<Real> prob_hi(BL_SPACEDIM);
        pp.getarr("prob_lo",prob_lo,0,BL_SPACEDIM);
        BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
        pp.getarr("prob_hi",prob_hi,0,BL_SPACEDIM);
        BL_ASSERT(prob_lo.size() == BL_SPACEDIM);
        prob_domain.setLo(prob_lo);
        prob_domain.setHi(prob_hi);
        SetOffset(prob_lo.data());
    } else {
        prob_domain.setLo(rb->lo());
        prob_domain.setHi(rb->hi());
        SetOffset(rb->lo());
    }

    pp.query("spherical_origin_fix", Geometry::spherical_origin_fix);

    //
    // Now get periodicity info.
    //
    if (isper == nullptr)
    {
        Vector<int> is_per(BL_SPACEDIM,0);
        pp.queryarr("is_periodic",is_per,0,BL_SPACEDIM);
        for (int n = 0; n < BL_SPACEDIM; n++)  
            is_periodic[n] = is_per[n];
    }
    else
    {
        for (int n = 0; n < BL_SPACEDIM; n++)  
            is_periodic[n] = isper[n];
    }

    amrex::ExecOnFinalize(Geometry::Finalize);
}

void
Geometry::GetVolume (MultiFab&       vol,
                     const BoxArray& grds,
		     const DistributionMapping& dm,
                     int             ngrow) const
{
    vol.define(grds,dm,1,ngrow,MFInfo(),FArrayBoxFactory());
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
#ifdef AMREX_USE_DEVICE
        Device::synchronize();
#endif
	CoordSys::SetVolume(vol[mfi], mfi.growntilebox());
    }
}

void
Geometry::GetVolume (FArrayBox&      vol,
                     const BoxArray& grds,
                     int             idx,
                     int             ngrow) const
{
    CoordSys::GetVolume(vol, amrex::grow(grds[idx],ngrow));
}

#if (BL_SPACEDIM <= 2)
void
Geometry::GetDLogA (MultiFab&       dloga,
                    const BoxArray& grds, 
                    const DistributionMapping& dm,
                    int             dir,
                    int             ngrow) const
{
    dloga.define(grds,dm,1,ngrow,MFInfo(),FArrayBoxFactory());
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
		       const DistributionMapping& dm,
                       int             dir,
                       int             ngrow) const
{
    BoxArray edge_boxes(grds);
    edge_boxes.surroundingNodes(dir);
    area.define(edge_boxes,dm,1,ngrow,MFInfo(),FArrayBoxFactory());

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
    CoordSys::GetFaceArea(area, amrex::grow(grds[idx],ngrow), dir);
}

void
Geometry::periodicShift (const Box&      target,
                         const Box&      src, 
                         Vector<IntVect>& out) const
{
    out.resize(0);

    Box locsrc(src);

    int nist,njst,nkst;
    int niend,njend,nkend;
    nist = njst = nkst = 0;
    niend = njend = nkend = 0;
    AMREX_D_TERM( nist , =njst , =nkst ) = -1;
    AMREX_D_TERM( niend , =njend , =nkend ) = +1;

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
                    out.push_back(IntVect(AMREX_D_DECL(ri*domain.length(0),
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
  Vector<int> baseBoxAI;

  CoordSys::BroadcastCoordSys(geom, fromProc, comm, bcastSource);

  if(bcastSource) {  // ---- initialize the source data
    const RealBox &realBox = geom.ProbDomain();
    for(int n(0); n < BL_SPACEDIM; ++n) {
      realBox_lo[n] = realBox.lo(n);
      realBox_hi[n] = realBox.hi(n);
      is_periodic[n] = geom.isPeriodic(n);
    }
    coord = geom.CoordInt();
    baseBoxAI = amrex::SerializeBox(geom.Domain());
  }


  // ---- do the broadcasts
  if( ! bcastSource) {
    baseBoxAI.resize(amrex::SerializeBoxSize());
  }
  ParallelDescriptor::Bcast(baseBoxAI.dataPtr(), baseBoxAI.size(), fromProc, comm);

  ParallelDescriptor::Bcast(realBox_lo, BL_SPACEDIM, fromProc, comm);
  ParallelDescriptor::Bcast(realBox_hi, BL_SPACEDIM, fromProc, comm);

  ParallelDescriptor::Bcast(&coord, 1, fromProc, comm);
  ParallelDescriptor::Bcast(is_periodic, BL_SPACEDIM, fromProc, comm);
  ParallelDescriptor::Bcast(&Geometry::spherical_origin_fix, 1, fromProc, comm);


  if( ! bcastSource) {    // ---- define the destination geometry
    Box baseBox(amrex::UnSerializeBox(baseBoxAI));
    RealBox realBox(realBox_lo, realBox_hi);

    geom.define(baseBox, &realBox, coord, is_periodic);
  }
}
#endif

}
