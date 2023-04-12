
#include <AMReX_Algorithm.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_SPACE.H>
#include <AMReX_COORDSYS_C.H>

#include <AMReX_OpenMP.H>

#include <iostream>

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

    computeRoundoffDomain();
}

void
Geometry::Setup (const RealBox* rb, int coord, int const* isper) noexcept
{
    Geometry* gg = AMReX::top()->getDefaultGeometry();

    if (gg->ok) return;

    AMREX_ASSERT(!OpenMP::in_parallel());

    ParmParse pp("geometry");

    if (coord >=0 && coord <= 2) {
        gg->SetCoord( (CoordType) coord );
    } else {
        coord = 0;  // default is Cartesian coordinates
        pp.queryAdd("coord_sys",coord);
        gg->SetCoord( (CoordType) coord );
    }

    if (rb == nullptr) {
        Vector<Real> prob_lo(AMREX_SPACEDIM);
        Vector<Real> prob_hi(AMREX_SPACEDIM);
        Vector<Real> prob_extent(AMREX_SPACEDIM);

        for (int i = 0; i < AMREX_SPACEDIM; i++) { prob_lo[i] = 0.; }
        pp.queryAdd("prob_lo", prob_lo, AMREX_SPACEDIM);
        AMREX_ASSERT(prob_lo.size() == AMREX_SPACEDIM);

        bool read_prob_hi = pp.queryarr("prob_hi",prob_hi,0,AMREX_SPACEDIM);
        AMREX_ASSERT(prob_hi.size() == AMREX_SPACEDIM);

        bool read_prob_extent = pp.queryarr("prob_extent",prob_extent,0,AMREX_SPACEDIM);
        AMREX_ASSERT(prob_extent.size() == AMREX_SPACEDIM);

        // We enforce that one and only one of prob_hi vs prob_extent is input
        AMREX_ALWAYS_ASSERT(  read_prob_hi || read_prob_extent);
        AMREX_ALWAYS_ASSERT(!(read_prob_hi && read_prob_extent));

        if (read_prob_extent)
        {
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                prob_hi[i] = prob_lo[i] + prob_extent[i];
        }

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
        pp.queryAdd("is_periodic", is_per);
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
Geometry::GetVolume (MultiFab& vol) const
{
    const auto a_dx = CellSizeArray();
    if (IsCartesian()) {
        vol.setVal(AMREX_D_TERM(a_dx[0],*a_dx[1],*a_dx[2]), 0, 1, vol.nGrowVect());
    } else {
#if (AMREX_SPACEDIM == 3)
        amrex::Abort("Geometry::GetVolume: for 3d, only Cartesian is supported");
#else
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            GpuArray<Real,AMREX_SPACEDIM> a_offset{{AMREX_D_DECL(offset[0],offset[1],offset[2])}};
            int coord = (int) c_sys;
            auto const& ma = vol.arrays();
            ParallelFor(vol, vol.nGrowVect(),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                amrex_setvol(makeSingleCellBox(i,j,k), ma[box_no], a_offset, a_dx, coord);
            });
            Gpu::streamSynchronize();
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(vol,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                CoordSys::SetVolume(vol[mfi], mfi.growntilebox());
            }
        }
#endif
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

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& ma = dloga.arrays();
        GpuArray<Real,AMREX_SPACEDIM> a_offset{AMREX_D_DECL(offset[0],offset[1],offset[2])};
        GpuArray<Real,AMREX_SPACEDIM> a_dx    {AMREX_D_DECL(    dx[0],    dx[1],    dx[2])};
        int coord = (int) c_sys;
        ParallelFor(dloga, IntVect(ngrow),
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            amrex_setdloga(makeSingleCellBox(i,j,k), ma[box_no], a_offset, a_dx, dir, coord);
        });
        Gpu::streamSynchronize();
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dloga,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            CoordSys::SetDLogA(dloga[mfi], mfi.growntilebox(), dir);
        }
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
    const auto a_dx = CellSizeArray();

    if (IsCartesian()) {
#if (AMREX_SPACEDIM == 1)
        amrex::ignore_unused(a_dx);
        const Real a0 = 1._rt;
#elif (AMREX_SPACEDIM == 2)
        const Real a0 = (dir == 0) ? a_dx[1] : a_dx[0];
#else
        const Real a0 = (dir == 0) ? a_dx[1]*a_dx[2]
            : ((dir == 1) ? a_dx[0]*a_dx[2] : a_dx[0]*a_dx[1]);
#endif
        area.setVal(a0, 0, 1, area.nGrowVect());
    } else {
#if (AMREX_SPACEDIM == 3)
        amrex::Abort("Geometry::GetFaceArea:: for 3d, only Cartesian is supported");
#else
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            GpuArray<Real,AMREX_SPACEDIM> a_offset{AMREX_D_DECL(offset[0],offset[1],offset[2])};
            int coord = (int) c_sys;
            auto const& ma = area.arrays();
            ParallelFor(area, area.nGrowVect(),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                amrex_setarea(makeSingleCellBox(i,j,k), ma[box_no], a_offset, a_dx, dir, coord);
            });
            Gpu::streamSynchronize();
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(area,TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                CoordSys::SetFaceArea(area[mfi],mfi.growntilebox(),dir);
            }
        }
#endif
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

    int nist  = -1;
    int niend =  1;
#if (AMREX_SPACEDIM > 1)
    int njst  = -1;
    int njend =  1;
#else
    int njst  = 0;
    int njend = 0;
#endif
#if (AMREX_SPACEDIM > 2)
    int nkst  = -1;
    int nkend =  1;
#else
    int nkst  = 0;
    int nkend = 0;
#endif

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
Geometry::growNonPeriodicDomain (IntVect const& ngrow) const noexcept
{
    Box b = Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (!isPeriodic(idim)) {
            b.grow(idim,ngrow[idim]);
        }
    }
    return b;
}

Box
Geometry::growPeriodicDomain (IntVect const& ngrow) const noexcept
{
    Box b = Domain();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (isPeriodic(idim)) {
            b.grow(idim,ngrow[idim]);
        }
    }
    return b;
}

Box
Geometry::growNonPeriodicDomain (int ngrow) const noexcept
{
    return growNonPeriodicDomain(IntVect(ngrow));
}

Box
Geometry::growPeriodicDomain (int ngrow) const noexcept
{
    return growPeriodicDomain(IntVect(ngrow));
}

void
Geometry::computeRoundoffDomain ()
{
    for (int k = 0; k < AMREX_SPACEDIM; k++)
    {
        offset[k] = prob_domain.lo(k);
        dx[k] = prob_domain.length(k)/(Real(domain.length(k)));
        inv_dx[k] = 1.0_rt/dx[k];
    }

    constexpr int maxiters = 200;

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        int ilo = Domain().smallEnd(idim);
        int ihi = Domain().bigEnd(idim);
        int ncells = ihi-ilo+1;
        Real plo = ProbLo(idim);
        Real phi = ProbHi(idim);
        Real dxinv = InvCellSize(idim);

        // Check that the grid is well formed and that deltax > roundoff
        AMREX_ASSERT((plo + ihi*CellSize(idim)) < (plo + (ihi + 1)*CellSize(idim)));

        // roundoff_lo will be the lowest value that will be inside the domain
        // roundoff_hi will be the highest value that will be inside the domain
        roundoff_lo[idim] = static_cast<ParticleReal>(plo);
        roundoff_hi[idim] = static_cast<ParticleReal>(phi);
        auto& rlo = roundoff_lo[idim];
        auto& rhi = roundoff_hi[idim];

        auto is_outside = [=] (auto x) -> bool
        {
            auto idx = int(std::floor((x - plo)*dxinv));
            return (idx < 0) || (idx >= ncells);
        };

        auto is_inside = [=] (auto x) -> bool
        {
            return !is_outside(x);
        };

        ParticleReal rlo_out;
        if (is_inside(rlo))
        {
            auto eps = std::numeric_limits<ParticleReal>::epsilon() * (rhi-rlo);
            int iters = 0;
            rlo_out = rlo - eps;
            while (is_inside(rlo_out) && iters < maxiters) {
                eps *= ParticleReal(2.);
                rlo_out = rlo - eps;
                ++iters;
            }
            // The assertion on rlo_out makes sure the compiler cannot optimize it away.
            AMREX_ALWAYS_ASSERT(rlo_out > std::numeric_limits<ParticleReal>::lowest()
                                && iters < maxiters);
        }
        else
        {
            auto eps = std::numeric_limits<ParticleReal>::epsilon() * (rhi-rlo);
            int iters = 0;
            auto rtmp = rlo + eps;
            while (is_outside(rtmp) && iters < maxiters) {
                eps *= ParticleReal(2.);
                rtmp = rlo + eps;
                ++iters;
            }
            rlo_out = rlo;
            rlo = rtmp;
            // The assertion on rtmp makes sure the compiler cannot optimize it away.
            AMREX_ALWAYS_ASSERT(rtmp > std::numeric_limits<ParticleReal>::lowest()
                                && iters < maxiters);
        }

        {
            // rlo_out is outside
            // rlo     is inside
            int iters = 0;
            auto epsilon = std::numeric_limits<ParticleReal>::epsilon()
                * std::max(ParticleReal(CellSize(idim)),std::abs(rlo))
                * ParticleReal(2.0);
            auto rlo_minus = rlo-epsilon;
            bool rlo_minus_is_inside = is_inside(rlo_minus);
            while (rlo_minus_is_inside && iters < maxiters) {
                auto rmid = ParticleReal(0.5)*(rlo_out+rlo);
                if (rmid == rlo_out || rmid == rlo) {
                    break;
                }
                if (is_inside(rmid)) {
                    rlo = rmid;
                    epsilon = std::numeric_limits<ParticleReal>::epsilon()
                        * std::max(ParticleReal(CellSize(idim)),std::abs(rlo))
                        * ParticleReal(2.0);
                    rlo_minus = rlo - epsilon;
                    rlo_minus_is_inside = is_inside(rlo_minus);
                } else {
                    rlo_out = rmid;
                }
                ++iters;
            }
            // The assertion on rlo_minus makes sure the compiler cannot optimize it away.
            AMREX_ALWAYS_ASSERT(rlo_minus > std::numeric_limits<ParticleReal>::lowest()
                                && iters < maxiters);
        }

        ParticleReal rhi_out;
        if (is_inside(rhi))
        {
            auto eps = std::numeric_limits<ParticleReal>::epsilon() * (rhi-rlo);
            int iters = 0;
            rhi_out = rhi + eps;
            while (is_inside(rhi_out) && iters < maxiters) {
                eps *= ParticleReal(2.);
                rhi_out = rhi + eps;
                ++iters;
            }
            // The assertion on rhi_out makes sure the compiler cannot optimize it away.
            AMREX_ALWAYS_ASSERT(rhi_out > std::numeric_limits<ParticleReal>::lowest()
                                && iters < maxiters);
        }
        else
        {
            auto eps = std::numeric_limits<ParticleReal>::epsilon() * (rhi-rlo);
            int iters = 0;
            // Yes, we have to write it this way for Intel compiler.
            // is_outside(rhi-eps) could be different from is_outside(rtmp),
            // where rtmp = rhs-eps.
            auto rtmp = rhi - eps;
            while (is_outside(rtmp) && iters < maxiters) {
                eps *= ParticleReal(2.);
                rtmp = rhi - eps;
                ++iters;
            }
            rhi_out = rhi;
            rhi = rtmp;
            // The assertion on rtmp makes sure the compiler cannot optimize it away.
            AMREX_ALWAYS_ASSERT(rtmp > std::numeric_limits<ParticleReal>::lowest()
                                && iters < maxiters);
        }

        {
            // rhi     is inside
            // rhi_out is outside
            int iters = 0;
            auto epsilon = std::numeric_limits<ParticleReal>::epsilon()
                * std::max(ParticleReal(CellSize(idim)),std::abs(rhi))
                * ParticleReal(2.0);
            auto rhi_plus = rhi+epsilon;
            bool rhi_plus_is_inside = is_inside(rhi_plus);
            while (rhi_plus_is_inside && iters < maxiters) {
                auto rmid = ParticleReal(0.5)*(rhi+rhi_out);
                if (rmid == rhi || rmid == rhi_out) {
                    break;
                }
                if (is_inside(rmid)) {
                    rhi = rmid;
                    epsilon = std::numeric_limits<ParticleReal>::epsilon()
                        * std::max(ParticleReal(CellSize(idim)),std::abs(rhi))
                        * ParticleReal(2.0);
                    rhi_plus = rhi + epsilon;
                    rhi_plus_is_inside = is_inside(rhi_plus);
                } else {
                    rhi_out = rmid;
                }
                ++iters;
            }
            // The assertion on rhi_plus makes sure the compiler cannot optimize it away.
            AMREX_ALWAYS_ASSERT(rhi_plus > std::numeric_limits<ParticleReal>::lowest()
                                && iters < maxiters);
        }
    }
}

bool
Geometry::outsideRoundoffDomain (AMREX_D_DECL(ParticleReal x, ParticleReal y, ParticleReal z)) const
{
    bool outside = AMREX_D_TERM(x < roundoff_lo[0]
                             || x > roundoff_hi[0],
                             || y < roundoff_lo[1]
                             || y > roundoff_hi[1],
                             || z < roundoff_lo[2]
                             || z > roundoff_hi[2]);
    return outside;
}

bool
Geometry::insideRoundoffDomain (AMREX_D_DECL(ParticleReal x, ParticleReal y, ParticleReal z)) const
{
    return !outsideRoundoffDomain(AMREX_D_DECL(x, y, z));
}

}
