

#include <iostream>

#include <AMReX_Algorithm.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_SPACE.H>

#include <AMReX_OpenMP.H>

namespace amrex {

std::ostream&
operator<< (std::ostream&   os,
            const Geometry& g)
{
    os << (CoordSys&) g << g.ProbDomain() << g.Domain() << 'P' << IntVect(g.isPeriodic());
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
    g.ProbDomain(rb);

    int ic = is.peek();
    if (ic == static_cast<int>('P')) {
        char c;
        is >> c;
        IntVect is_per;
        is >> is_per;
        g.setPeriodicity({{AMREX_D_DECL(is_per[0],is_per[1],is_per[2])}});
    } else {
        g.setPeriodicity(DefaultGeometry().isPeriodic());
    }

    return is;
}

Geometry::Geometry () noexcept
{
    if (!AMReX::empty()) *this = DefaultGeometry();
}

Geometry::Geometry (const Box& dom, const RealBox* rb, int coord,
                    int const* is_per) noexcept
{
    define(dom,rb,coord,is_per);
}

Geometry::Geometry (const Box& dom, const RealBox& rb, int coord,
                    Array<int,AMREX_SPACEDIM> const& is_per) noexcept
{
    define(dom,rb,coord,is_per);
}

void
Geometry::define (const Box& dom, const RealBox& rb, int coord,
                  Array<int,AMREX_SPACEDIM> const& is_per) noexcept
{
    define(dom, &rb, coord, is_per.data());
}

void
Geometry::define (const Box& dom, const RealBox* rb, int coord,
                  int const* is_per) noexcept
{
    Setup(rb,coord,is_per);

    Geometry* gg = AMReX::top()->getDefaultGeometry();

    if (coord == -1) {
        c_sys = gg->Coord();
    } else {
        c_sys = static_cast<CoordType>(coord);
    }

    if (is_per == nullptr) {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            is_periodic[idim] = gg->isPeriodic(idim);
        }
    } else {
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            is_periodic[idim] = is_per[idim];
        }
    }

    if (rb == nullptr) {
        prob_domain = gg->ProbDomain();
    } else {
        prob_domain = *rb;
    }

    domain = dom;
    ok     = true;

    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        offset[k] = prob_domain.lo(k);
        dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
	inv_dx[k] = 1.0/dx[k];
    }

    computeRoundoffDomain();
}

void
Geometry::Setup (const RealBox* rb, int coord, int const* isper) noexcept
{
    Geometry* gg = AMReX::top()->getDefaultGeometry();

    if (gg->ok) return;

    BL_ASSERT(!OpenMP::in_parallel());

    ParmParse pp("geometry");

    if (coord >=0 && coord <= 2) {
        gg->SetCoord( (CoordType) coord );        
    } else {
        coord = 0;  // default is Cartesian coordinates
        pp.query("coord_sys",coord);
        gg->SetCoord( (CoordType) coord );        
    }

    if (rb == nullptr) {
        Vector<Real> prob_lo(AMREX_SPACEDIM);
        Vector<Real> prob_hi(AMREX_SPACEDIM);
        pp.getarr("prob_lo",prob_lo,0,AMREX_SPACEDIM);
        BL_ASSERT(prob_lo.size() == AMREX_SPACEDIM);
        pp.getarr("prob_hi",prob_hi,0,AMREX_SPACEDIM);
        BL_ASSERT(prob_hi.size() == AMREX_SPACEDIM);
        gg->prob_domain.setLo(prob_lo);
        gg->prob_domain.setHi(prob_hi);
        gg->SetOffset(prob_lo.data());
    } else {
        gg->prob_domain.setLo(rb->lo());
        gg->prob_domain.setHi(rb->hi());
        gg->SetOffset(rb->lo());
    }

    //
    // Now get periodicity info.
    //
    if (isper == nullptr)
    {
        Vector<int> is_per(AMREX_SPACEDIM,0);
        pp.queryarr("is_periodic",is_per,0,AMREX_SPACEDIM);
        for (int n = 0; n < AMREX_SPACEDIM; n++) {
            gg->is_periodic[n] = is_per[n];
        }
    }
    else
    {
        for (int n = 0; n < AMREX_SPACEDIM; n++) {
            gg->is_periodic[n] = isper[n];
        }
    }

    gg->ok = true;
}

void
Geometry::ResetDefaultProbDomain (const RealBox& rb) noexcept
{
    Geometry* gg = AMReX::top()->getDefaultGeometry();
    gg->prob_domain.setLo(rb.lo());
    gg->prob_domain.setHi(rb.hi());
    gg->SetOffset(rb.lo());
}

void
Geometry::ResetDefaultPeriodicity (const Array<int,AMREX_SPACEDIM>& is_per) noexcept
{
    Geometry* gg = AMReX::top()->getDefaultGeometry();
    gg->setPeriodicity(is_per);
}

void
Geometry::ResetDefaultCoord (int coord) noexcept
{
    AMREX_ASSERT(coord >= -1 && coord <= 2);
    Geometry* gg = AMReX::top()->getDefaultGeometry();
    gg->SetCoord(static_cast<CoordType>(coord));
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(vol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
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
    CoordSys::GetVolume(vol, amrex::grow(grds[idx],ngrow));
}

#if (AMREX_SPACEDIM <= 2)
void
Geometry::GetDLogA (MultiFab&       dloga,
                    const BoxArray& grds, 
                    const DistributionMapping& dm,
                    int             dir,
                    int             ngrow) const
{
    dloga.define(grds,dm,1,ngrow,MFInfo(),FArrayBoxFactory());
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(dloga,TilingIfNotGPU()); mfi.isValid(); ++mfi)
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
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(area,TilingIfNotGPU()); mfi.isValid(); ++mfi)
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
                         Vector<IntVect>& out) const noexcept
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
            if (rj != 0
#if (AMREX_SPACEDIM > 1)
                && !is_periodic[1]
#endif
                )
            {
                continue;
            }
            if (rj != 0
#if (AMREX_SPACEDIM > 1)
                && is_periodic[1]
#endif
                )
            {
                locsrc.shift(1,rj*domain.length(1));
            }

            for (rk = nkst; rk <= nkend; rk++)
            {
                if (rk!=0
#if (AMREX_SPACEDIM == 3)
                    && !is_periodic[2]
#endif
                    )
                {
                    continue;
                }
                if (rk!=0
#if (AMREX_SPACEDIM == 3)
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
#if (AMREX_SPACEDIM == 3)
                    && is_periodic[2]
#endif
                    )
                {
                    locsrc.shift(2,-rk*domain.length(2));
                }
            }
            if (rj != 0
#if (AMREX_SPACEDIM > 1)
                && is_periodic[1]
#endif
                )
            {
                locsrc.shift(1,-rj*domain.length(1));
            }
        }
        if (ri != 0 && is_periodic[0])
            locsrc.shift(0,-ri*domain.length(0));
    }
}

Box
Geometry::growNonPeriodicDomain (int ngrow) const noexcept
{
    Box b = Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (!isPeriodic(idim)) {
            b.grow(idim,ngrow);
        }
    }
    return b;
}

Box
Geometry::growPeriodicDomain (int ngrow) const noexcept
{
    Box b = Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (isPeriodic(idim)) {
            b.grow(idim,ngrow);
        }
    }
    return b;
}

void
Geometry::computeRoundoffDomain ()
{
    roundoff_domain = prob_domain;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        int ilo = Domain().smallEnd(idim);
        int ihi = Domain().bigEnd(idim);
        Real plo = ProbLo(idim);
        Real phi = ProbHi(idim);
        Real idx = InvCellSize(idim);
        Real deltax = CellSize(idim);

#ifdef AMREX_SINGLE_PRECISION_PARTICLES
        Real tolerance = std::max(1.e-4*deltax, 1.e-10*phi);
#else
        Real tolerance = std::max(1.e-8*deltax, 1.e-14*phi);
#endif
        // bisect the point at which the cell no longer maps to inside the domain
        Real lo = static_cast<Real>(phi) - Real(0.5)*static_cast<Real>(deltax);
        Real hi = static_cast<Real>(phi) + Real(0.5)*static_cast<Real>(deltax);

        Real mid = bisect(lo, hi,
                          [=] AMREX_GPU_HOST_DEVICE (Real x) -> Real
                          {
                              int i = int(Math::floor((x - plo)*idx)) + ilo;
                              bool inside = i >= ilo and i <= ihi;
                              return static_cast<Real>(inside) - Real(0.5);
                          }, tolerance);
        roundoff_domain.setHi(idim, mid - tolerance);
    }
}

bool
Geometry::outsideRoundoffDomain (AMREX_D_DECL(Real x, Real y, Real z)) const
{
    bool outside = AMREX_D_TERM(x <  roundoff_domain.lo(0)
                             || x >= roundoff_domain.hi(0),
                             || y <  roundoff_domain.lo(1)
                             || y >= roundoff_domain.hi(1),
                             || z <  roundoff_domain.lo(2)
                             || z >= roundoff_domain.hi(2));
    return outside;
}

bool
Geometry::insideRoundoffDomain (AMREX_D_DECL(Real x, Real y, Real z)) const
{
    return !outsideRoundoffDomain(AMREX_D_DECL(x, y, z));
}

}
