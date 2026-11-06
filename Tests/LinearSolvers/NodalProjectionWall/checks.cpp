// The checks of MLNodeLaplacian::compRHS without a solve (checks = 1): the nodal divergence
// at walls with and without MLNodeLaplacian::setWallDivuCorrection, its consistency with a
// cell-centred source, its volume conservation, and the wall nodes of a two-level hierarchy.
// The checks with a known answer stop the test (amrex::Abort in main) if they fail.

#include "NodalProjectionWall.H"

#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

using namespace amrex;

namespace {

// A pseudo-random number in [-1, 1] for cell (i,j,k) and component c, the same on every run
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real hash_unit (int i, int j, int k, int c)
{
    unsigned s = (unsigned(i)*73856093u) ^ (unsigned(j)*19349663u) ^ (unsigned(k)*83492791u)
                 ^ (unsigned(c)*2654435761u);
    s ^= s >> 13;
    s *= 0x5bd1e995u;
    s ^= s >> 15;
    return Real(s & 0xffffffu) / Real(0xffffff) * Real(2) - Real(1);
}

// Direct check of the nodal divergence stencil at a wall, with no solve: u = (x y, -y^2/2) is
// divergence-free with u_n'' = -1 on the bottom wall.  compRHS should give 0 at interior
// nodes (the stencil is exact for quadratics there).  On the bottom-wall nodes, if the
// tangential derivative is taken at the cell centres (y = h/2) instead of averaged over the
// half-cell control volume (mean h/4), the divergence comes out as h/4 instead of 0, i.e.
// (h/4) |u_n''|.  Returns max |rhs| at interior nodes, and the mean, min and max of rhs/(h/4) on
// the wall nodes away from the corners.
std::array<Real,4> divu_wall_check (int n, int max_grid_size)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), Real(nz)*h)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, 1)};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
    vel.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const& v = vel.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = (Real(i)+Real(0.5))*h;
            Real const y = (Real(j)+Real(0.5))*h;
            v(i,j,k,0) = x*y;
            v(i,j,k,1) = -Real(0.5)*y*y;
        });
    }
    vel.FillBoundary(geom.periodicity());

    Array<LinOpBCType,AMREX_SPACEDIM> const bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                            LinOpBCType::Neumann,
                                                            LinOpBCType::Periodic)};
    MLNodeLaplacian linop({geom}, {ba}, {dm});
    linop.setDomainBC(bc, bc);
    MultiFab sigma(ba, dm, 1, 0);
    sigma.setVal(1.0);
    linop.setSigma(0, sigma);
    linop.setWallDivuCorrection(wall_corr);
    BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
    MultiFab rhs(nba, dm, 1, 0);
    linop.compRHS({&rhs}, {&vel}, {}, {});
    // Host copy for the loops below (rhs is in device memory on GPU builds)
    MultiFab rhs_h(nba, dm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    rhs_h.ParallelCopy(rhs, 0, 0, 1);
    Gpu::streamSynchronize();

    // Interior nodes: max |rhs|.  Bottom-wall nodes (j = 0), corners left out: mean and
    // spread of rhs / (h/4).
    Real max_int = 0, sum_wall = 0, min_wall = 1.e30, max_wall = -1.e30;
    int nwall = 0;
    int const e = n/8;
    for (MFIter mfi(rhs_h); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const& r = rhs_h.const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
            if (i > 0 && i < n && j > 0 && j < n) {
                max_int = std::max(max_int, std::abs(r(i,j,k)));
            } else if (j == 0 && i >= e && i <= n-e
#if (AMREX_SPACEDIM == 3)
                       && k == 0
#endif
                       ) {
                Real const q = r(i,j,k) / (Real(0.25)*h);
                sum_wall += q; min_wall = std::min(min_wall, q); max_wall = std::max(max_wall, q);
                ++nwall;
            }
        });
    }
    ParallelDescriptor::ReduceRealMax(max_int);
    ParallelDescriptor::ReduceRealSum(sum_wall);
    ParallelDescriptor::ReduceRealMin(min_wall);
    ParallelDescriptor::ReduceRealMax(max_wall);
    ParallelDescriptor::ReduceIntSum(nwall);
    return {max_int, sum_wall/Real(nwall), min_wall, max_wall};
}

// Consistency of the nodal divergence with a cell-centred source, no solve: u = (x (y+y^2),
// y (x+x^2)) has div u = S = x + x^2 + y + y^2, so compRHS({u}, rhcc = -S) is zero exactly;
// what it returns is the discretisation's inconsistency between how it averages div u and
// S over each node's control volume.  The tangential velocity grows linearly (then
// quadratically) away from each wall, as a no-slip profile does.  Prints max |rhs| / h^2
// at interior nodes, on the x = 0 and y = 0 wall nodes away from the other walls, and at
// the (0,0) corner node: O(1) in these units is second order, O(1/h) is first.  Returns them.
std::array<Real,3> source_wall_check (int n, int max_grid_size)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), Real(nz)*h)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, 1)};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
    MultiFab src(ba, dm, 1, 1);
    vel.setVal(0.0);
    src.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const& v = vel.array(mfi);
        auto const& sa = src.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = (Real(i)+Real(0.5))*h;
            Real const y = (Real(j)+Real(0.5))*h;
            v(i,j,k,0) = x*(y + y*y);
            v(i,j,k,1) = y*(x + x*x);
            sa(i,j,k) = -(x + x*x + y + y*y);
        });
    }
    vel.FillBoundary(geom.periodicity());
    src.FillBoundary(geom.periodicity());

    Array<LinOpBCType,AMREX_SPACEDIM> const bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                            LinOpBCType::Neumann,
                                                            LinOpBCType::Periodic)};
    MLNodeLaplacian linop({geom}, {ba}, {dm});
    linop.setDomainBC(bc, bc);
    MultiFab sigma(ba, dm, 1, 0);
    sigma.setVal(1.0);
    linop.setSigma(0, sigma);
    linop.setWallDivuCorrection(wall_corr);
    BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
    MultiFab rhs(nba, dm, 1, 0);
    linop.compRHS({&rhs}, {&vel}, {}, {&src});
    MultiFab rhs_h(nba, dm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    rhs_h.ParallelCopy(rhs, 0, 0, 1);
    Gpu::streamSynchronize();

    Real max_int = 0, max_wall = 0, corner = 0;
    int const e = n/8;
    for (MFIter mfi(rhs_h); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const& r = rhs_h.const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
#if (AMREX_SPACEDIM == 3)
            if (k != 0) { return; }
#endif
            Real const a = std::abs(r(i,j,k)) / (h*h);
            if (i > 0 && i < n && j > 0 && j < n) {
                max_int = std::max(max_int, a);
            } else if ((j == 0 && i >= e && i <= n-e) || (i == 0 && j >= e && j <= n-e)) {
                max_wall = std::max(max_wall, a);
            } else if (i == 0 && j == 0) {
                corner = a;
            }
        });
    }
    ParallelDescriptor::ReduceRealMax(max_int);
    ParallelDescriptor::ReduceRealMax(max_wall);
    ParallelDescriptor::ReduceRealMax(corner);
    return {max_int, max_wall, corner};
}

// Where the divergence error sits near a corner: the Burggraf velocity is divergence-free,
// so compRHS should return zero; prints rhs / h^2 on the 3x3 nodes at the (0,0) corner
// (node (0,0) is the corner, (i,0) are on the y = 0 wall, (0,j) on the x = 0 wall), and at the
// lid corner (0,n).
void burggraf_corner_rhs (int n, int max_grid_size)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), Real(nz)*h)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, 1)};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
    vel.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        auto const& v = vel.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = (Real(i)+Real(0.5))*h;
            Real const y = (Real(j)+Real(0.5))*h;
            v(i,j,k,0) =  Real(8) * bf(x) * bg1(y);
            v(i,j,k,1) = -Real(8) * bf1(x) * bg(y);
        });
    }
    vel.FillBoundary(geom.periodicity());

    Array<LinOpBCType,AMREX_SPACEDIM> const bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                            LinOpBCType::Neumann,
                                                            LinOpBCType::Periodic)};
    MLNodeLaplacian linop({geom}, {ba}, {dm});
    linop.setDomainBC(bc, bc);
    MultiFab sigma(ba, dm, 1, 0);
    sigma.setVal(1.0);
    linop.setSigma(0, sigma);
    linop.setWallDivuCorrection(wall_corr);
    BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
    MultiFab rhs(nba, dm, 1, 0);
    linop.compRHS({&rhs}, {&vel}, {}, {});

    // Gather the corner blocks on the I/O rank: (0,0) and the lid corner (0,n)
    for (int cj0 : {0, n-2}) {
    Box const cbox(IntVect(AMREX_D_DECL(0, cj0, 0)), IntVect(AMREX_D_DECL(2, cj0+2, 0)));
    BoxArray const cba(cbox);
    DistributionMapping const cdm(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab cmf(cba, cdm, 1, 0, MFInfo().SetArena(The_Pinned_Arena()));
    cmf.ParallelCopy(rhs, 0, 0, 1);
    Gpu::streamSynchronize();
    for (MFIter mfi(cmf); mfi.isValid(); ++mfi) {
        auto const& r = cmf.const_array(mfi);
        amrex::Print() << "  n = " << n << ((cj0 == 0) ? "  corner (0,0)" : "  lid corner (0,n)") << "\n";
        for (int j = cj0+2; j >= cj0; --j) {
            char line[256];
            std::snprintf(line, sizeof(line), "    j=%-4d  %12.4f  %12.4f  %12.4f\n", j,
                          double(r(0,j,0)/(h*h)), double(r(1,j,0)/(h*h)), double(r(2,j,0)/(h*h)));
            amrex::Print() << line;
        }
        amrex::Print() << "              i=0           i=1           i=2\n";
    }
    }
}

// Volume conservation of compRHS, with no solve.  The cell velocity is random, so it changes by
// O(1) from cell to cell, as it does next to a corner where the boundary data do not match.
// Walls at y = 0 and 1, and at x = 0 and 1 either walls too (closed = 1) or inflow faces with
// u = 1 and 0.7 in the ghost cells.  The sum of rhs over the nodes' control volumes equals the
// flux through the domain boundary, 0 or 0.7 - 1 per unit height, for any velocity, if the
// fluxes through the faces of the control volumes telescope.  The original stencil and the
// correction both do; returns the difference.  In 3D z is periodic, or with zwall = 1 walls too, so that the box
// has nodes on three walls (or two walls and an inflow face).
Real volume_conservation_check (int n, int max_grid_size, int closed, int zwall)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
    Real const lz = Real(nz)*h;
#else
    Real const lz = Real(1);
    amrex::ignore_unused(zwall);
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), lz)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, (zwall != 0) ? 0 : 1)};
    Geometry const geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    DistributionMapping const dm(ba);

    Real const ulo = (closed != 0) ? Real(0) : Real(1);
    Real const uhi = (closed != 0) ? Real(0) : Real(0.7);
    MultiFab vel(ba, dm, AMREX_SPACEDIM, 1);
    vel.setVal(0.0);
    for (MFIter mfi(vel); mfi.isValid(); ++mfi) {
        Box const& vbx = mfi.validbox();
        auto const& v = vel.array(mfi);
        ParallelFor(vbx, AMREX_SPACEDIM, [=] AMREX_GPU_DEVICE (int i, int j, int k, int c) noexcept
        {
            v(i,j,k,c) = hash_unit(i,j,k,c);
        });
        if (closed != 0) { continue; }
        // The inflow ghost cells, including those at the box corners: u on the face, and a
        // random tangential velocity, which compRHS must not see
        for (int side = 0; side < 2; ++side) {
            if ((side == 0 && vbx.smallEnd(0) != 0) || (side == 1 && vbx.bigEnd(0) != n-1)) {
                continue;
            }
            Box const gbx = amrex::grow((side == 0) ? amrex::adjCellLo(vbx, 0)
                                                    : amrex::adjCellHi(vbx, 0),
                                        IntVect(AMREX_D_DECL(0,1,1)));
            Real const un = (side == 0) ? ulo : uhi;
            ParallelFor(gbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (j < 0 || j >= n) { return; }
#if (AMREX_SPACEDIM == 3)
                if (zwall != 0 && (k < 0 || k >= nz)) { return; }
#endif
                v(i,j,k,0) = un;
                for (int c = 1; c < AMREX_SPACEDIM; ++c) { v(i,j,k,c) = hash_unit(i,j,k,c); }
            });
        }
    }
    vel.FillBoundary(geom.periodicity());

    LinOpBCType const xbc = (closed != 0) ? LinOpBCType::Neumann : LinOpBCType::inflow;
    Array<LinOpBCType,AMREX_SPACEDIM> const bc{AMREX_D_DECL(xbc, LinOpBCType::Neumann,
                                                            (zwall != 0) ? LinOpBCType::Neumann
                                                                         : LinOpBCType::Periodic)};
    MLNodeLaplacian linop({geom}, {ba}, {dm});
    linop.setDomainBC(bc, bc);
    MultiFab sigma(ba, dm, 1, 0);
    sigma.setVal(1.0);
    linop.setSigma(0, sigma);
    linop.setWallDivuCorrection(wall_corr);
    BoxArray const nba = amrex::convert(ba, IntVect::TheNodeVector());
    MultiFab rhs(nba, dm, 1, 0);
    linop.compRHS({&rhs}, {&vel}, {}, {});

    DistributionMapping const dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    MultiFab rhs1(BoxArray(amrex::surroundingNodes(domain)), dm1, 1, 0,
                  MFInfo().SetArena(The_Pinned_Arena()));
    rhs1.ParallelCopy(rhs, 0, 0, 1);
    Gpu::streamSynchronize();
    Real sum = 0;
    for (MFIter mfi(rhs1); mfi.isValid(); ++mfi) {
        auto const& a = rhs1.const_array(mfi);
        Box bx = mfi.validbox();
#if (AMREX_SPACEDIM == 3)
        if (zwall == 0) { bx.setBig(2, nz-1); }   // z is periodic: node nz is node 0
#endif
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
        {
            Real vol = AMREX_D_TERM(h, *h, *h);
            if (i == 0 || i == n) { vol *= Real(0.5); }
            if (j == 0 || j == n) { vol *= Real(0.5); }
#if (AMREX_SPACEDIM == 3)
            if (zwall != 0 && (k == 0 || k == nz)) { vol *= Real(0.5); }
#endif
            sum += a(i,j,k) * vol;
        });
    }
    ParallelDescriptor::ReduceRealSum(sum);
    return sum - (uhi - ulo) * lz;
}

// The wall correction with two AMR levels, no solve.  Walls on all four sides (z periodic in
// 3D), and a fine level (ratio 2) over the coarse cells [3n/8, 5n/8) x [0, n/4), touching the
// bottom wall.  The rhs at a node on the coarse/fine boundary is rebuilt after the single-level
// pass, from the coarse cells and the fine level (mlndlap_divu_cf_contrib,
// mlndlap_divu_fine_contrib), with no wall correction.  Two divergence-free velocities:
//   field 0: u = (x y, -y^2/2), as in divu_wall_check: u_n'' = -1 on the bottom wall, so an
//            uncorrected wall node gives rhs = h/4 and a corrected one 0;
//   field 1: u = (y, 0), a linear shear, with u_n'' = 0: every wall node gives 0.
// Returns, at k = 0:
//   [0] rhs / (h/4), mean over the coarse bottom-wall nodes i = n/8 .. 3n/8 - 1;
//   [1] rhs / (h/4), h the fine cell size, mean over the fine bottom-wall nodes two fine cells
//       or more from the ends of the fine level;
//   [2], [3] rhs at the coarse/fine nodes on the bottom wall, i = 3n/8 and 5n/8;
//   [4] rhs at the corner of the fine level away from the wall, (3n/8, n/4);
//   [5] max |rhs| / (h/4) over the other coarse/fine nodes, on the sides and the top of the
//       fine level.
// [2] to [4] are corners of the coarse/fine boundary, where the stencil is not exact even for a
// linear shear: it takes the velocity in each cell as constant, and the error that leaves, which
// cancels between the four cells of a uniform node, does not where they differ in size.
std::array<Real,6> amr_wall_check (int n, int max_grid_size, int field)
{
    Real const h = Real(1)/Real(n);
#if (AMREX_SPACEDIM == 3)
    int const nz = 8;
#endif
    Box const domain(IntVect(0), IntVect(AMREX_D_DECL(n-1, n-1, nz-1)));
    RealBox const rb({AMREX_D_DECL(Real(0), Real(0), Real(0))},
                     {AMREX_D_DECL(Real(1), Real(1), Real(nz)*h)});
    Array<int,AMREX_SPACEDIM> const is_periodic{AMREX_D_DECL(0, 0, 1)};
    Vector<Geometry> geom(2);
    geom[0].define(domain, rb, CoordSys::cartesian, is_periodic);
    geom[1].define(amrex::refine(domain, 2), rb, CoordSys::cartesian, is_periodic);

    int const ilo = 3*n/8, ihi = 5*n/8, jhi = n/4;
    Vector<BoxArray> ba(2);
    ba[0] = BoxArray(domain);
    ba[0].maxSize(max_grid_size);
    ba[1] = BoxArray(amrex::refine(Box(IntVect(AMREX_D_DECL(ilo, 0, 0)),
                                       IntVect(AMREX_D_DECL(ihi-1, jhi-1, nz-1))), 2));
    ba[1].maxSize(max_grid_size);
    Vector<DistributionMapping> dm{DistributionMapping(ba[0]), DistributionMapping(ba[1])};

    // Cell-centre values; ghost cells inside the domain too (as a fill from the coarse level
    // would give them), and 0 outside it (the wall velocity)
    Vector<MultiFab> vel(2);
    for (int lev = 0; lev < 2; ++lev) {
        Real const hl = h / Real(lev+1);
        Box const ldom = geom[lev].growPeriodicDomain(1);
        vel[lev].define(ba[lev], dm[lev], AMREX_SPACEDIM, 1);
        for (MFIter mfi(vel[lev]); mfi.isValid(); ++mfi) {
            auto const& v = vel[lev].array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real const x = (Real(i)+Real(0.5))*hl;
                Real const y = (Real(j)+Real(0.5))*hl;
                bool const in = ldom.contains(IntVect(AMREX_D_DECL(i,j,k)));
                v(i,j,k,0) = in ? ((field == 0) ? x*y : y) : Real(0);
                v(i,j,k,1) = in ? ((field == 0) ? -Real(0.5)*y*y : Real(0)) : Real(0);
#if (AMREX_SPACEDIM == 3)
                v(i,j,k,2) = Real(0);
#endif
            });
        }
    }

    Array<LinOpBCType,AMREX_SPACEDIM> const bc{AMREX_D_DECL(LinOpBCType::Neumann,
                                                            LinOpBCType::Neumann,
                                                            LinOpBCType::Periodic)};
    MLNodeLaplacian linop(geom, ba, dm);
    linop.setDomainBC(bc, bc);
    Vector<MultiFab> sigma(2);
    Vector<MultiFab> rhs(2);
    for (int lev = 0; lev < 2; ++lev) {
        sigma[lev].define(ba[lev], dm[lev], 1, 0);
        sigma[lev].setVal(1.0);
        linop.setSigma(lev, sigma[lev]);
        rhs[lev].define(amrex::convert(ba[lev], IntVect::TheNodeVector()), dm[lev], 1, 0);
    }
    linop.setWallDivuCorrection(wall_corr);
    linop.compRHS({&rhs[0], &rhs[1]}, {&vel[0], &vel[1]}, {}, {});

    DistributionMapping const dm1(Vector<int>{ParallelDescriptor::IOProcessorNumber()});
    Vector<MultiFab> rhs1(2);
    for (int lev = 0; lev < 2; ++lev) {
        rhs1[lev].define(BoxArray(amrex::surroundingNodes(ba[lev].minimalBox())), dm1, 1, 0,
                         MFInfo().SetArena(The_Pinned_Arena()));
        rhs1[lev].ParallelCopy(rhs[lev], 0, 0, 1);
    }
    Gpu::streamSynchronize();

    std::array<Real,6> r{};
    if (ParallelDescriptor::IOProcessor()) {
        auto const& c = rhs1[0][0].const_array();
        auto const& f = rhs1[1][0].const_array();
        Real const qc = Real(0.25)*h;
        Real const qf = Real(0.25)*h/Real(2);
        for (int i = n/8; i < ilo; ++i) { r[0] += c(i,0,0) / qc; }
        r[0] /= Real(ilo - n/8);
        int nf = 0;
        for (int i = 2*ilo+2; i <= 2*ihi-2; ++i) { r[1] += f(i,0,0) / qf; ++nf; }
        r[1] /= Real(nf);
        r[2] = c(ilo,0,0);
        r[3] = c(ihi,0,0);
        r[4] = c(ilo,jhi,0);
        for (int j = 1; j < jhi; ++j) {
            r[5] = std::max({r[5], std::abs(c(ilo,j,0)/qc), std::abs(c(ihi,j,0)/qc)});
        }
        for (int i = ilo+1; i < ihi; ++i) { r[5] = std::max(r[5], std::abs(c(i,jhi,0)/qc)); }
    }
    ParallelDescriptor::Bcast(r.data(), r.size(), ParallelDescriptor::IOProcessorNumber());
    return r;
}

}

int run_checks (Vector<int> const& n_cell, int max_grid_size)
{
    int nfail = 0;
    auto check = [&] (bool ok, char const* what, int n)
    {
        if (!ok) {
            ++nfail;
            amrex::Print() << "  FAILED: " << what << ", n = " << n << "\n";
        }
    };
    // Rounding errors of rhs relative to its terms, which are O(1/h)
    Real const tol = Real(1.e-10);
    // What the wall nodes of u = (x y, -y^2/2) give, rhs / (h/4)
    Real const wall_expected = wall_corr ? Real(0) : Real(1);
    char line[256];

    amrex::Print() << "\nNodal divergence of u = (x y, -y^2/2) (div u = 0, u_n'' = -1 on the bottom\n"
                   << "wall): max |rhs| at interior nodes, and rhs / (h/4) on the bottom-wall\n"
                   << "nodes away from the corners (mean, min, max): 1 without the correction, 0 with it\n"
                   << "     n   max interior         mean          min          max\n";
    for (int n : n_cell) {
        auto const q = divu_wall_check(n, max_grid_size);
        std::snprintf(line, sizeof(line), "  %4d   %12.3e   %10.4f   %10.4f   %10.4f\n",
                      n, double(q[0]), double(q[1]), double(q[2]), double(q[3]));
        amrex::Print() << line;
        check(q[0] < tol, "interior divergence of u = (x y, -y^2/2)", n);
        check(std::abs(q[2]-wall_expected) < tol && std::abs(q[3]-wall_expected) < tol,
              "wall divergence of u = (x y, -y^2/2)", n);
    }

    amrex::Print() << "\nNodal divergence with a cell-centred source, u = (x(y+y^2), y(x+x^2)),\n"
                   << "S = x+x^2+y+y^2 (div u - S = 0): max |rhs| / h^2 at interior nodes, on the\n"
                   << "x = 0 and y = 0 wall nodes away from the other walls, and at the (0,0) corner\n"
                   << "     n       interior           wall         corner\n";
    for (int n : n_cell) {
        auto const q = source_wall_check(n, max_grid_size);
        std::snprintf(line, sizeof(line), "  %4d   %12.4f   %12.4f   %12.4f\n",
                      n, double(q[0]), double(q[1]), double(q[2]));
        amrex::Print() << line;
        Real const h2 = Real(1)/(Real(n)*Real(n));
        check(std::max({q[0], q[1], q[2]})*h2 < tol, "velocity and source treated alike", n);
    }

    amrex::Print() << "\nVolume conservation, random cell velocity: sum of rhs over the control volumes\n"
                   << "minus the flux through the boundary (0 to round-off if the face fluxes\n"
                   << "telescope); closed box, and walls with inflow faces at both x ends\n"
#if (AMREX_SPACEDIM == 3)
                   << "(z periodic, then z walls too)\n"
                   << "     n     closed box    with inflow faces     closed, z walls   inflow, z walls\n";
#else
                   << "     n     closed box    with inflow faces\n";
#endif
    for (int n : n_cell) {
#if (AMREX_SPACEDIM == 3)
        std::array<Real,4> const q{volume_conservation_check(n, max_grid_size, 1, 0),
                                   volume_conservation_check(n, max_grid_size, 0, 0),
                                   volume_conservation_check(n, max_grid_size, 1, 1),
                                   volume_conservation_check(n, max_grid_size, 0, 1)};
        std::snprintf(line, sizeof(line), "  %4d   %12.3e   %12.3e         %12.3e      %12.3e\n", n,
                      double(q[0]), double(q[1]), double(q[2]), double(q[3]));
#else
        std::array<Real,2> const q{volume_conservation_check(n, max_grid_size, 1, 0),
                                   volume_conservation_check(n, max_grid_size, 0, 0)};
        std::snprintf(line, sizeof(line), "  %4d   %12.3e   %12.3e\n", n, double(q[0]), double(q[1]));
#endif
        amrex::Print() << line;
        for (auto const v : q) {
            check(std::abs(v) < Real(1.e-12), "volume conservation", n);
        }
    }

    amrex::Print() << "\nTwo levels, fine level (ratio 2) on the bottom wall over the coarse cells\n"
                   << "[3n/8, 5n/8) x [0, n/4).  rhs / (h/4), h of the node's level, on the coarse\n"
                   << "and fine wall nodes away from the coarse/fine boundary (1 uncorrected, 0\n"
                   << "corrected, for u = (x y, -y^2/2)); rhs at the corners of the coarse/fine\n"
                   << "boundary: on the wall, and at (3n/8, n/4); max |rhs| / (h/4) at the other\n"
                   << "coarse/fine nodes\n"
                   << "               wall / (h/4)              rhs at the c/f corners             other c/f\n"
                   << "     n      coarse       fine     wall i=3n/8  wall i=5n/8  (3n/8,n/4)     / (h/4)\n";
    for (int field = 0; field < 2; ++field) {
        amrex::Print() << ((field == 0) ? "  u = (x y, -y^2/2)\n" : "  u = (y, 0)\n");
        Real const expected = (field == 0) ? wall_expected : Real(0);
        for (int n : n_cell) {
            auto const q = amr_wall_check(n, max_grid_size, field);
            std::snprintf(line, sizeof(line), "  %4d   %9.4f  %9.4f     %10.3e   %10.3e   %10.3e   %9.4f\n",
                          n, double(q[0]), double(q[1]), double(q[2]), double(q[3]),
                          double(q[4]), double(q[5]));
            amrex::Print() << line;
            check(std::abs(q[0]-expected) < tol && std::abs(q[1]-expected) < tol,
                  "wall divergence on the coarse and fine levels", n);
        }
    }

    amrex::Print() << "\nBurggraf (div u = 0): rhs / h^2 on the 3x3 nodes at the (0,0) corner and at the\n"
                   << "lid corner (0,n)\n";
    for (int n : n_cell) {
        burggraf_corner_rhs(n, max_grid_size);
    }
    amrex::Print() << '\n';
    return nfail;
}
