#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MLNodeLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

using namespace amrex;

namespace {

struct Parameters
{
    std::string mode = "rz_mms";
    std::vector<int> n_cell = {32, 64, 128};
    int max_grid_size = 16;
    int max_coarsening_level = 100;
    int verbose = 0;
    int bottom_verbose = 0;
    int max_iter = 100;
    int max_fmg_iter = 0;
    int do_assert = 1;
    int min_mg_levels = 1;
    Real reltol = 1.e-11;
    Real abstol = 0.0;
    Real min_l2_order = 1.8;
    Real min_linf_order = 1.8;
    Real max_finest_l2 = 5.e-5;
    Real max_finest_linf = 1.e-4;
};

struct SolveResult
{
    int n_cell = 0;
    int nprocs = 1;
    int mg_levels = 0;
    int iterations = 0;
    Real residual = 0.0;
    Real l2 = 0.0;
    Real linf = 0.0;
    Real setup_time = 0.0;
    Real solve_time = 0.0;
    Real total_time = 0.0;
};

Parameters read_parameters ()
{
    Parameters p;
    ParmParse pp;
    pp.query("mode", p.mode);
    pp.queryarr("n_cell", p.n_cell);
    pp.query("max_grid_size", p.max_grid_size);
    pp.query("max_coarsening_level", p.max_coarsening_level);
    pp.query("verbose", p.verbose);
    pp.query("bottom_verbose", p.bottom_verbose);
    pp.query("max_iter", p.max_iter);
    pp.query("max_fmg_iter", p.max_fmg_iter);
    pp.query("min_mg_levels", p.min_mg_levels);
    pp.query("reltol", p.reltol);
    pp.query("abstol", p.abstol);
    pp.query("do_assert", p.do_assert);
    pp.query("min_l2_order", p.min_l2_order);
    pp.query("min_linf_order", p.min_linf_order);
    pp.query("max_finest_l2", p.max_finest_l2);
    pp.query("max_finest_linf", p.max_finest_linf);
    return p;
}

Geometry make_rz_geometry (int n_cell)
{
    Box domain(IntVect(AMREX_D_DECL(0, 0, 0)),
               IntVect(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)));
    RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)},
               {AMREX_D_DECL(1.0, 1.0, 1.0)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, rb, CoordSys::RZ, is_periodic);
}

Geometry make_cartesian_geometry (int n_cell)
{
    Box domain(IntVect(AMREX_D_DECL(0, 0, 0)),
               IntVect(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)));
    RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)},
               {AMREX_D_DECL(1.0, 1.0, 1.0)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    return Geometry(domain, rb, CoordSys::cartesian, is_periodic);
}

void fill_rz_sigma (MultiFab& sigma, Geometry const& geom)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sigma, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const sig = sigma.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            sig(i,j,k) = plo[0] + (static_cast<Real>(i) + Real(0.5)) * dx[0];
        });
    }
}

void fill_rz_radial_mms (MultiFab& exact, MultiFab& rhs, MultiFab& sigma,
                         Geometry const& geom)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real a = pi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const ex = exact.array(mfi);
        Array4<Real> const rh = rhs.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + static_cast<Real>(i) * dx[0];
            Real const s = std::sin(a*r);
            Real const c = std::cos(a*r);
            ex(i,j,k) = s;
            rh(i,j,k) = a*c - r*a*a*s;
        });
    }

    fill_rz_sigma(sigma, geom);
}

void fill_rz_full_mms (MultiFab& exact, MultiFab& rhs, MultiFab& sigma,
                       Geometry const& geom)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real a = pi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const ex = exact.array(mfi);
        Array4<Real> const rh = rhs.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + static_cast<Real>(i) * dx[0];
            Real const z = plo[1] + static_cast<Real>(j) * dx[1];
            Real const sr = std::sin(a*r);
            Real const cr = std::cos(a*r);
            Real const sz = std::sin(a*z);
            ex(i,j,k) = sr * sz;
            rh(i,j,k) = (a*cr - Real(2.0)*r*a*a*sr) * sz;
        });
    }

    fill_rz_sigma(sigma, geom);
}

void fill_cartesian_mms (MultiFab& exact, MultiFab& rhs, MultiFab& sigma,
                         Geometry const& geom)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real a = pi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const ex = exact.array(mfi);
        Array4<Real> const rh = rhs.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = plo[0] + static_cast<Real>(i) * dx[0];
            Real const s = std::sin(a*x);
            ex(i,j,k) = s;
            rh(i,j,k) = -a*a*s;
        });
    }

    sigma.setVal(1.0);
}

void fill_rz_eb_neumann_mms (MultiFab& exact, MultiFab& sigma,
                             MultiFab& velocity, MultiFab& rhcc,
                             MultiFab& eb_neumann,
                             Geometry const& geom, Real eb_location)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real a = pi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const ex = exact.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + static_cast<Real>(i) * dx[0];
            Real const z = plo[1] + static_cast<Real>(j) * dx[1];
            ex(i,j,k) = std::sin(a*r) * std::sin(a*z);
        });
    }

    fill_rz_sigma(sigma, geom);

    velocity.setVal(0.0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhcc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const rh = rhcc.array(mfi);
        Array4<Real> const ebneu = eb_neumann.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + (static_cast<Real>(i) + Real(0.5)) * dx[0];
            Real const z = plo[1] + (static_cast<Real>(j) + Real(0.5)) * dx[1];
            Real const sr = std::sin(a*r);
            Real const cr = std::cos(a*r);
            Real const sz = std::sin(a*z);
            rh(i,j,k) = (a*cr - Real(2.0)*r*a*a*sr) * sz;
            ebneu(i,j,k) = a * std::cos(a*eb_location) * sz;
        });
    }
}

#ifdef AMREX_USE_EB
void fill_rz_eb_neumann_curved_mms (MultiFab& exact, MultiFab& sigma,
                                    MultiFab& velocity, MultiFab& rhcc,
                                    MultiFab& eb_neumann,
                                    Geometry const& geom,
                                    EBFArrayBoxFactory const& factory)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real a = pi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const ex = exact.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + static_cast<Real>(i) * dx[0];
            Real const z = plo[1] + static_cast<Real>(j) * dx[1];
            ex(i,j,k) = std::sin(a*r) * std::sin(a*z);
        });
    }

    fill_rz_sigma(sigma, geom);

    velocity.setVal(0.0);
    rhcc.setVal(0.0);
    eb_neumann.setVal(0.0);

    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& bcent = factory.getBndryCent();
    const auto& bnorm = factory.getBndryNormal();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhcc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const rh = rhcc.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + (static_cast<Real>(i) + Real(0.5)) * dx[0];
            Real const z = plo[1] + (static_cast<Real>(j) + Real(0.5)) * dx[1];
            Real const sr = std::sin(a*r);
            Real const cr = std::cos(a*r);
            Real const sz = std::sin(a*z);
            rh(i,j,k) = (a*cr - Real(2.0)*r*a*a*sr) * sz;
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(eb_neumann, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        if (!bcent.ok(mfi) || !bnorm.ok(mfi)) { continue; }

        Box const& bx = mfi.tilebox();
        Array4<Real> const ebneu = eb_neumann.array(mfi);
        Array4<EBCellFlag const> const flag = flags.const_array(mfi);
        Array4<Real const> const bc = bcent.const_array(mfi);
        Array4<Real const> const bn = bnorm.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isSingleValued()) {
                Real const nx = bn(i,j,k,0);
                Real const ny = bn(i,j,k,1);
                Real const r = plo[0] + (static_cast<Real>(i) + Real(0.5) + bc(i,j,k,0)) * dx[0];
                Real const z = plo[1] + (static_cast<Real>(j) + Real(0.5) + bc(i,j,k,1)) * dx[1];
                Real const dphidr = a * std::cos(a*r) * std::sin(a*z);
                Real const dphidz = a * std::sin(a*r) * std::cos(a*z);
                ebneu(i,j,k) = dphidr * nx + dphidz * ny;
            }
        });
    }
}

void fill_rz_eb_inflow_velocity_curved_mms (MultiFab& eb_velocity,
                                            Geometry const& geom,
                                            EBFArrayBoxFactory const& factory)
{
    eb_velocity.setVal(0.0);

    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real a = pi;

    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& bcent = factory.getBndryCent();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(eb_velocity, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        if (!bcent.ok(mfi)) { continue; }

        Box const& bx = mfi.tilebox();
        Array4<Real> const vel = eb_velocity.array(mfi);
        Array4<EBCellFlag const> const flag = flags.const_array(mfi);
        Array4<Real const> const bc = bcent.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isSingleValued()) {
                Real const r = plo[0] + (static_cast<Real>(i) + Real(0.5) + bc(i,j,k,0)) * dx[0];
                Real const z = plo[1] + (static_cast<Real>(j) + Real(0.5) + bc(i,j,k,1)) * dx[1];
                vel(i,j,k,0) = -r * a * std::cos(a*r) * std::sin(a*z);
                vel(i,j,k,1) = -r * a * std::sin(a*r) * std::cos(a*z);
            }
        });
    }
}

void fill_rz_eb_neumann_physical_bc_mms (MultiFab& exact, MultiFab& sigma,
                                         MultiFab& velocity, MultiFab& rhcc,
                                         MultiFab& eb_neumann,
                                         Geometry const& geom,
                                         EBFArrayBoxFactory const& factory,
                                         bool all_neumann)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real ar = Real(2.0) * pi;
    Real const az = all_neumann ? Real(2.0) * pi : pi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const ex = exact.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + static_cast<Real>(i) * dx[0];
            Real const z = plo[1] + static_cast<Real>(j) * dx[1];
            Real const zfac = all_neumann ? std::cos(az*z) : std::sin(az*z);
            ex(i,j,k) = std::cos(ar*r) * zfac;
        });
    }

    fill_rz_sigma(sigma, geom);

    velocity.setVal(0.0);
    rhcc.setVal(0.0);
    eb_neumann.setVal(0.0);

    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& bcent = factory.getBndryCent();
    const auto& bnorm = factory.getBndryNormal();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhcc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const rh = rhcc.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const r = plo[0] + (static_cast<Real>(i) + Real(0.5)) * dx[0];
            Real const z = plo[1] + (static_cast<Real>(j) + Real(0.5)) * dx[1];
            Real const sr = std::sin(ar*r);
            Real const cr = std::cos(ar*r);
            Real const zfac = all_neumann ? std::cos(az*z) : std::sin(az*z);
            rh(i,j,k) = -ar*sr*zfac - r*(ar*ar + az*az)*cr*zfac;
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(eb_neumann, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        if (!bcent.ok(mfi) || !bnorm.ok(mfi)) { continue; }

        Box const& bx = mfi.tilebox();
        Array4<Real> const ebneu = eb_neumann.array(mfi);
        Array4<EBCellFlag const> const flag = flags.const_array(mfi);
        Array4<Real const> const bc = bcent.const_array(mfi);
        Array4<Real const> const bn = bnorm.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isSingleValued()) {
                Real const nx = bn(i,j,k,0);
                Real const ny = bn(i,j,k,1);
                Real const r = plo[0] + (static_cast<Real>(i) + Real(0.5) + bc(i,j,k,0)) * dx[0];
                Real const z = plo[1] + (static_cast<Real>(j) + Real(0.5) + bc(i,j,k,1)) * dx[1];
                Real const sr = std::sin(ar*r);
                Real const cr = std::cos(ar*r);
                Real const zfac = all_neumann ? std::cos(az*z) : std::sin(az*z);
                Real const dzfac = all_neumann ? -az*std::sin(az*z) : az*std::cos(az*z);
                Real const dphidr = -ar * sr * zfac;
                Real const dphidz = cr * dzfac;
                ebneu(i,j,k) = dphidr * nx + dphidz * ny;
            }
        });
    }
}

void fill_cartesian_eb_neumann_curved_mms (MultiFab& exact, MultiFab& sigma,
                                           MultiFab& velocity, MultiFab& rhcc,
                                           MultiFab& eb_neumann,
                                           Geometry const& geom,
                                           EBFArrayBoxFactory const& factory)
{
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real a = pi;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const ex = exact.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = plo[0] + static_cast<Real>(i) * dx[0];
            Real const y = plo[1] + static_cast<Real>(j) * dx[1];
            ex(i,j,k) = std::sin(a*x) * std::sin(a*y);
        });
    }

    sigma.setVal(1.0);
    velocity.setVal(0.0);
    rhcc.setVal(0.0);
    eb_neumann.setVal(0.0);

    const auto& flags = factory.getMultiEBCellFlagFab();
    const auto& bcent = factory.getBndryCent();
    const auto& bnorm = factory.getBndryNormal();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(rhcc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const rh = rhcc.array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = plo[0] + (static_cast<Real>(i) + Real(0.5)) * dx[0];
            Real const y = plo[1] + (static_cast<Real>(j) + Real(0.5)) * dx[1];
            rh(i,j,k) = -Real(2.0)*a*a*std::sin(a*x)*std::sin(a*y);
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(eb_neumann, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        if (!bcent.ok(mfi) || !bnorm.ok(mfi)) { continue; }

        Box const& bx = mfi.tilebox();
        Array4<Real> const ebneu = eb_neumann.array(mfi);
        Array4<EBCellFlag const> const flag = flags.const_array(mfi);
        Array4<Real const> const bc = bcent.const_array(mfi);
        Array4<Real const> const bn = bnorm.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (flag(i,j,k).isSingleValued()) {
                Real const nx = bn(i,j,k,0);
                Real const ny = bn(i,j,k,1);
                Real const x = plo[0] + (static_cast<Real>(i) + Real(0.5) + bc(i,j,k,0)) * dx[0];
                Real const y = plo[1] + (static_cast<Real>(j) + Real(0.5) + bc(i,j,k,1)) * dx[1];
                Real const dphidx = a * std::cos(a*x) * std::sin(a*y);
                Real const dphidy = a * std::sin(a*x) * std::cos(a*y);
                ebneu(i,j,k) = dphidx * nx + dphidy * ny;
            }
        });
    }
}
#endif

void initialize_dirichlet_guess (MultiFab& solution, MultiFab const& exact,
                                 Geometry const& geom)
{
    MultiFab::Copy(solution, exact, 0, 0, 1, 0);
    Box const interior = amrex::surroundingNodes(amrex::grow(geom.Domain(), -1));
    solution.setVal(0.0, interior, 0, 1, 0);
}

std::pair<Real,Real> compute_error_norms (MultiFab& solution, MultiFab const& exact,
                                          Geometry const& geom,
                                          iMultiFab const* active_mask = nullptr,
                                          bool remove_constant_offset = false)
{
    MultiFab error(solution.boxArray(), solution.DistributionMap(), 1, 0);
    MultiFab::Copy(error, solution, 0, 0, 1, 0);
    MultiFab::Subtract(error, exact, 0, 0, 1, 0);

    auto mask = error.OwnerMask(geom.periodicity());
    if (active_mask != nullptr) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*mask, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.tilebox();
            Array4<int> const owner = mask->array(mfi);
            Array4<int const> const active = active_mask->const_array(mfi);
            ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                owner(i,j,k) *= active(i,j,k);
            });
        }
    }

    Long const npts = mask->sum(0);
    if (npts <= 0) {
        Abort("NodeRZEBBaseline: error norm mask has no active points");
    }

    if (remove_constant_offset) {
        MultiFab ones(error.boxArray(), error.DistributionMap(), 1, 0);
        ones.setVal(1.0);
        Real const errsum = MultiFab::Dot(*mask, error, 0, ones, 0, 1, 0);
        error.plus(-errsum / static_cast<Real>(npts), 0, 1, 0);
    }

    Real const linf = error.norm0(*mask, 0, 0);
    Real const dot = MultiFab::Dot(*mask, error, 0, error, 0, 1, 0);
    Real const l2 = std::sqrt(dot / static_cast<Real>(npts));
    return {l2, linf};
}

#ifdef AMREX_USE_EB
std::unique_ptr<iMultiFab>
make_eb_active_node_mask (BoxArray const& nba, DistributionMapping const& dmap,
                          Geometry const& geom, EBFArrayBoxFactory const& factory)
{
    auto active = std::make_unique<iMultiFab>(nba, dmap, 1, 0);
    active->setVal(0);

    Box const& ccdom = geom.Domain();
    const auto domlo = amrex::lbound(ccdom);
    const auto domhi = amrex::ubound(ccdom);
    const auto& flags = factory.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*active, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<int> const act = active->array(mfi);
        Array4<EBCellFlag const> const flag = flags.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int is_active = 0;
            for (int jj = j-1; jj <= j; ++jj) {
                for (int ii = i-1; ii <= i; ++ii) {
                    bool const in_domain = (ii >= domlo.x && ii <= domhi.x &&
                                            jj >= domlo.y && jj <= domhi.y);
                    if (in_domain && !flag(ii,jj,k).isCovered()) {
                        is_active = 1;
                    }
                }
            }
            act(i,j,k) = is_active;
        });
    }

    return active;
}

std::unique_ptr<iMultiFab>
make_sphere_fluid_node_mask (BoxArray const& nba, DistributionMapping const& dmap,
                             Geometry const& geom, iMultiFab const& active_mask,
                             Real radius, Real cx, Real cy)
{
    auto fluid = std::make_unique<iMultiFab>(nba, dmap, 1, 0);
    fluid->setVal(0);
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    Real const r2 = radius * radius;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*fluid, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<int> const fl = fluid->array(mfi);
        Array4<int const> const act = active_mask.const_array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const x = plo[0] + static_cast<Real>(i) * dx[0];
            Real const y = plo[1] + static_cast<Real>(j) * dx[1];
            Real const d2 = (x-cx)*(x-cx) + (y-cy)*(y-cy);
            fl(i,j,k) = (act(i,j,k) && d2 >= r2) ? 1 : 0;
        });
    }

    return fluid;
}
#endif

LPInfo make_lp_info (Parameters const& p)
{
    LPInfo info;
    info.setMaxCoarseningLevel(p.max_coarsening_level);
    return info;
}

SolveResult run_rz_mms_case (Parameters const& p, int n_cell,
                             MLNodeLaplacian::CoarseningStrategy strategy,
                             bool full_mms)
{
    Real const t0 = second();

    Geometry geom = make_rz_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0);

    if (full_mms) {
        fill_rz_full_mms(exact, rhs, sigma, geom);
    } else {
        fill_rz_radial_mms(exact, rhs, sigma, geom);
    }
    initialize_dirichlet_guess(solution, exact, geom);

    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info);
    linop.setCoarseningStrategy(strategy);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});
    linop.setSigma(0, sigma);

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto const [l2, linf] = compute_error_norms(solution, exact, geom);
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
}

SolveResult run_cartesian_mms_case (Parameters const& p, int n_cell)
{
    Real const t0 = second();

    Geometry geom = make_cartesian_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0);

    fill_cartesian_mms(exact, rhs, sigma, geom);
    initialize_dirichlet_guess(solution, exact, geom);

    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info);
    linop.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::RAP);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});
    linop.setSigma(0, sigma);

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto const [l2, linf] = compute_error_norms(solution, exact, geom);
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
}

void report_result (char const* label, SolveResult const& r)
{
    Print() << label
            << " n_cell=" << r.n_cell
            << " mpi_ranks=" << r.nprocs
            << " mg_levels=" << r.mg_levels
            << " iterations=" << r.iterations
            << " residual=" << r.residual
            << " l2=" << r.l2
            << " linf=" << r.linf
            << " setup_s=" << r.setup_time
            << " solve_s=" << r.solve_time
            << " total_s=" << r.total_time
            << "\n";
}

void assert_mms_results (Parameters const& p, std::vector<SolveResult> const& results,
                         char const* label)
{
    if (results.empty()) {
        Abort(std::string("NodeRZEBBaseline: no ") + label + " cases were run");
    }
    for (auto const& r : results) {
        if (p.do_assert &&
            (!std::isfinite(r.residual) || !std::isfinite(r.l2) || !std::isfinite(r.linf))) {
            Abort(std::string("NodeRZEBBaseline: non-finite ") + label +
                  " residual or error norm");
        }
        if (p.do_assert && (r.iterations <= 0 || r.iterations > p.max_iter)) {
            Abort(std::string("NodeRZEBBaseline: invalid ") + label + " MLMG iteration count");
        }
        if (p.do_assert && r.mg_levels < p.min_mg_levels) {
            Abort(std::string("NodeRZEBBaseline: ") + label +
                  " MG level count below threshold");
        }
    }

    for (std::size_t i = 1; i < results.size(); ++i) {
        Real const ratio = static_cast<Real>(results[i].n_cell)
            / static_cast<Real>(results[i-1].n_cell);
        Real const l2_order = std::log(results[i-1].l2 / results[i].l2) / std::log(ratio);
        Real const linf_order = std::log(results[i-1].linf / results[i].linf) / std::log(ratio);
        Print() << label << "_ORDER"
                << " n_cell_coarse=" << results[i-1].n_cell
                << " n_cell_fine=" << results[i].n_cell
                << " l2_order=" << l2_order
                << " linf_order=" << linf_order
                << "\n";
        if (p.do_assert && (l2_order < p.min_l2_order || linf_order < p.min_linf_order)) {
            Abort(std::string("NodeRZEBBaseline: ") + label +
                  " convergence order below threshold");
        }
    }

    SolveResult const& finest = results.back();
    if (p.do_assert && (finest.l2 > p.max_finest_l2 || finest.linf > p.max_finest_linf)) {
        Abort(std::string("NodeRZEBBaseline: finest-grid ") + label +
              " error above threshold");
    }
}

void run_rz_mms (Parameters const& p)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    for (int n : n_cells) {
        if (n < 4) {
            Abort("NodeRZEBBaseline: n_cell must be at least 4");
        }
        SolveResult r = run_rz_mms_case(p, n,
                                        MLNodeLaplacian::CoarseningStrategy::Sigma,
                                        false);
        report_result("RZ_MMS", r);
        results.push_back(r);
    }

    assert_mms_results(p, results, "RZ_MMS");
    Print() << "RZ_MMS_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

void run_rz_mms_rap (Parameters const& p)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    for (int n : n_cells) {
        if (n < 4) {
            Abort("NodeRZEBBaseline: n_cell must be at least 4");
        }
        SolveResult r = run_rz_mms_case(p, n,
                                        MLNodeLaplacian::CoarseningStrategy::RAP,
                                        true);
        report_result("RZ_MMS_RAP", r);
        results.push_back(r);
    }

    assert_mms_results(p, results, "RZ_MMS_RAP");
    Print() << "RZ_MMS_RAP_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

void run_cartesian_mms (Parameters const& p)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    for (int n : n_cells) {
        if (n < 4) {
            Abort("NodeRZEBBaseline: n_cell must be at least 4");
        }
        SolveResult r = run_cartesian_mms_case(p, n);
        report_result("CART_MMS_RAP", r);
        results.push_back(r);
    }

    assert_mms_results(p, results, "CART_MMS_RAP");
    Print() << "CART_MMS_RAP_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

SolveResult run_rz_eb_rap_zero_case (Parameters const& p, int n_cell)
{
#ifdef AMREX_USE_EB
    Real const t0 = second();

    Geometry geom = make_rz_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    EB2::SphereIF sphere(0.20, {AMREX_D_DECL(0.55, 0.50, 0.0)}, false);
    auto gshop = EB2::makeShop(sphere);
    EB2::Build(gshop, geom, 0, p.max_coarsening_level);

    EB2::IndexSpace const& eb_is = EB2::IndexSpace::top();
    EB2::Level const& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2, 2, 2};
    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, EBSupport::full);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0, MFInfo(), factory);

    solution.setVal(0.0);
    exact.setVal(0.0);
    rhs.setVal(0.0);
    fill_rz_sigma(sigma, geom);

    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info,
                          Vector<EBFArrayBoxFactory const*>{&factory});
    linop.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::RAP);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});
    linop.setSigma(0, sigma);

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto const [l2, linf] = compute_error_norms(solution, exact, geom);
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
#else
    amrex::ignore_unused(p);
    amrex::ignore_unused(n_cell);
    Abort("NodeRZEBBaseline: rz_eb_rap_zero requires AMREX_USE_EB");
    return {};
#endif
}

SolveResult run_rz_eb_neumann_plane_case (Parameters const& p, int n_cell)
{
#ifdef AMREX_USE_EB
    Real const t0 = second();

    Geometry geom = make_rz_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    constexpr Real eb_location = Real(0.47);
    EB2::PlaneIF plane({AMREX_D_DECL(eb_location, 0.0, 0.0)},
                       {AMREX_D_DECL(1.0, 0.0, 0.0)}, true);
    auto gshop = EB2::makeShop(plane);
    EB2::Build(gshop, geom, 0, p.max_coarsening_level);

    EB2::IndexSpace const& eb_is = EB2::IndexSpace::top();
    EB2::Level const& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2, 2, 2};
    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, EBSupport::full);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab velocity(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
    MultiFab rhcc(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab eb_neumann(grids, dmap, 1, 1, MFInfo(), factory);

    fill_rz_eb_neumann_mms(exact, sigma, velocity, rhcc, eb_neumann, geom,
                           eb_location);
    initialize_dirichlet_guess(solution, exact, geom);
    rhs.setVal(0.0);

    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info,
                          Vector<EBFArrayBoxFactory const*>{&factory});
    linop.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::RAP);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});
    linop.setSigma(0, sigma);
    linop.setEBInhomogNeumannFlux(0, eb_neumann);
    linop.compRHS(Vector<MultiFab*>{&rhs}, Vector<MultiFab*>{&velocity},
                  Vector<const MultiFab*>(), Vector<MultiFab*>{&rhcc});

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto active_mask = make_eb_active_node_mask(nba, dmap, geom, factory);
    auto const [l2, linf] = compute_error_norms(solution, exact, geom, active_mask.get());
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
#else
    amrex::ignore_unused(p);
    amrex::ignore_unused(n_cell);
    Abort("NodeRZEBBaseline: rz_eb_neumann_plane requires AMREX_USE_EB");
    return {};
#endif
}

SolveResult run_rz_eb_neumann_sphere_case (Parameters const& p, int n_cell)
{
#ifdef AMREX_USE_EB
    Real const t0 = second();

    Geometry geom = make_rz_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    EB2::SphereIF sphere(0.33, {AMREX_D_DECL(0.92, 0.50, 0.0)}, false);
    auto gshop = EB2::makeShop(sphere);
    EB2::Build(gshop, geom, 0, p.max_coarsening_level);

    EB2::IndexSpace const& eb_is = EB2::IndexSpace::top();
    EB2::Level const& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2, 2, 2};
    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, EBSupport::full);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab velocity(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
    MultiFab rhcc(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab eb_neumann(grids, dmap, 1, 1, MFInfo(), factory);

    fill_rz_eb_neumann_curved_mms(exact, sigma, velocity, rhcc, eb_neumann, geom, factory);
    initialize_dirichlet_guess(solution, exact, geom);
    rhs.setVal(0.0);

    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info,
                          Vector<EBFArrayBoxFactory const*>{&factory});
    linop.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::RAP);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});
    linop.setSigma(0, sigma);
    linop.setEBInhomogNeumannFlux(0, eb_neumann);
    linop.compRHS(Vector<MultiFab*>{&rhs}, Vector<MultiFab*>{&velocity},
                  Vector<const MultiFab*>(), Vector<MultiFab*>{&rhcc});

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto active_mask = make_eb_active_node_mask(nba, dmap, geom, factory);
    auto fluid_node_mask = make_sphere_fluid_node_mask(nba, dmap, geom, *active_mask,
                                                       Real(0.33), Real(0.92), Real(0.50));
    auto const [l2, linf] = compute_error_norms(solution, exact, geom, fluid_node_mask.get());
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
#else
    amrex::ignore_unused(p);
    amrex::ignore_unused(n_cell);
    Abort("NodeRZEBBaseline: rz_eb_neumann_sphere requires AMREX_USE_EB");
    return {};
#endif
}

SolveResult run_rz_eb_inflow_sphere_case (Parameters const& p, int n_cell)
{
#ifdef AMREX_USE_EB
    Real const t0 = second();

    Geometry geom = make_rz_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    EB2::SphereIF sphere(0.33, {AMREX_D_DECL(0.92, 0.50, 0.0)}, false);
    auto gshop = EB2::makeShop(sphere);
    EB2::Build(gshop, geom, 0, p.max_coarsening_level);

    EB2::IndexSpace const& eb_is = EB2::IndexSpace::top();
    EB2::Level const& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2, 2, 2};
    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, EBSupport::full);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab velocity(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
    MultiFab eb_velocity(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
    MultiFab rhcc(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab eb_neumann(grids, dmap, 1, 1, MFInfo(), factory);

    fill_rz_eb_neumann_curved_mms(exact, sigma, velocity, rhcc, eb_neumann, geom, factory);
    fill_rz_eb_inflow_velocity_curved_mms(eb_velocity, geom, factory);
    initialize_dirichlet_guess(solution, exact, geom);
    rhs.setVal(0.0);

    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info,
                          Vector<EBFArrayBoxFactory const*>{&factory});
    linop.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::RAP);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});
    linop.setSigma(0, sigma);
    linop.setEBInflowVelocity(0, eb_velocity);
    linop.compRHS(Vector<MultiFab*>{&rhs}, Vector<MultiFab*>{&velocity},
                  Vector<const MultiFab*>(), Vector<MultiFab*>{&rhcc});

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto active_mask = make_eb_active_node_mask(nba, dmap, geom, factory);
    auto fluid_node_mask = make_sphere_fluid_node_mask(nba, dmap, geom, *active_mask,
                                                       Real(0.33), Real(0.92), Real(0.50));
    auto const [l2, linf] = compute_error_norms(solution, exact, geom, fluid_node_mask.get());
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
#else
    amrex::ignore_unused(p);
    amrex::ignore_unused(n_cell);
    Abort("NodeRZEBBaseline: rz_eb_inflow_sphere requires AMREX_USE_EB");
    return {};
#endif
}

SolveResult run_cartesian_eb_neumann_sphere_case (Parameters const& p, int n_cell)
{
#ifdef AMREX_USE_EB
    Real const t0 = second();

    Geometry geom = make_cartesian_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    EB2::SphereIF sphere(0.33, {AMREX_D_DECL(0.92, 0.50, 0.0)}, false);
    auto gshop = EB2::makeShop(sphere);
    EB2::Build(gshop, geom, 0, p.max_coarsening_level);

    EB2::IndexSpace const& eb_is = EB2::IndexSpace::top();
    EB2::Level const& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2, 2, 2};
    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, EBSupport::full);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab velocity(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
    MultiFab rhcc(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab eb_neumann(grids, dmap, 1, 1, MFInfo(), factory);

    fill_cartesian_eb_neumann_curved_mms(exact, sigma, velocity, rhcc, eb_neumann, geom, factory);
    initialize_dirichlet_guess(solution, exact, geom);
    rhs.setVal(0.0);

    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info,
                          Vector<EBFArrayBoxFactory const*>{&factory});
    linop.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::RAP);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet,
                                    LinOpBCType::Dirichlet)});
    linop.setSigma(0, sigma);
    linop.setEBInhomogNeumannFlux(0, eb_neumann);
    linop.compRHS(Vector<MultiFab*>{&rhs}, Vector<MultiFab*>{&velocity},
                  Vector<const MultiFab*>(), Vector<MultiFab*>{&rhcc});

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto active_mask = make_eb_active_node_mask(nba, dmap, geom, factory);
    auto fluid_node_mask = make_sphere_fluid_node_mask(nba, dmap, geom, *active_mask,
                                                       Real(0.33), Real(0.92), Real(0.50));
    auto const [l2, linf] = compute_error_norms(solution, exact, geom, fluid_node_mask.get());
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
#else
    amrex::ignore_unused(p);
    amrex::ignore_unused(n_cell);
    Abort("NodeRZEBBaseline: cart_eb_neumann_sphere requires AMREX_USE_EB");
    return {};
#endif
}

SolveResult run_rz_eb_neumann_physical_bc_case (Parameters const& p, int n_cell,
                                                bool all_neumann,
                                                bool add_incompatible_offset)
{
#ifdef AMREX_USE_EB
    Real const t0 = second();

    Geometry geom = make_rz_geometry(n_cell);
    BoxArray grids(geom.Domain());
    grids.maxSize(p.max_grid_size);
    DistributionMapping dmap(grids);

    EB2::SphereIF sphere(0.33, {AMREX_D_DECL(0.92, 0.50, 0.0)}, false);
    auto gshop = EB2::makeShop(sphere);
    EB2::Build(gshop, geom, 0, p.max_coarsening_level);

    EB2::IndexSpace const& eb_is = EB2::IndexSpace::top();
    EB2::Level const& eb_level = eb_is.getLevel(geom);
    Vector<int> ng_ebs = {2, 2, 2};
    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, EBSupport::full);

    BoxArray const nba = amrex::convert(grids, IntVect::TheNodeVector());
    MultiFab solution(nba, dmap, 1, 0);
    MultiFab exact(nba, dmap, 1, 0);
    MultiFab rhs(nba, dmap, 1, 0);
    MultiFab sigma(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab velocity(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
    MultiFab rhcc(grids, dmap, 1, 0, MFInfo(), factory);
    MultiFab eb_neumann(grids, dmap, 1, 1, MFInfo(), factory);

    fill_rz_eb_neumann_physical_bc_mms(exact, sigma, velocity, rhcc, eb_neumann,
                                       geom, factory, all_neumann);
    if (add_incompatible_offset) {
        rhcc.plus(1.0, 0, 1);
    }
    if (all_neumann) {
        solution.setVal(0.0);
    } else {
        initialize_dirichlet_guess(solution, exact, geom);
    }
    rhs.setVal(0.0);

    LinOpBCType const y_bc = all_neumann ? LinOpBCType::Neumann : LinOpBCType::Dirichlet;
    LPInfo info = make_lp_info(p);
    MLNodeLaplacian linop({geom}, {grids}, {dmap}, info,
                          Vector<EBFArrayBoxFactory const*>{&factory});
    linop.setCoarseningStrategy(MLNodeLaplacian::CoarseningStrategy::RAP);
    linop.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                    y_bc,
                                    LinOpBCType::Neumann)},
                      {AMREX_D_DECL(LinOpBCType::Neumann,
                                    y_bc,
                                    LinOpBCType::Neumann)});
    linop.setSigma(0, sigma);
    linop.setEBInhomogNeumannFlux(0, eb_neumann);
    linop.compRHS(Vector<MultiFab*>{&rhs}, Vector<MultiFab*>{&velocity},
                  Vector<const MultiFab*>(), Vector<MultiFab*>{&rhcc});

    MLMG mlmg(linop);
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.bottom_verbose);
    mlmg.setMaxIter(p.max_iter);
    mlmg.setMaxFmgIter(p.max_fmg_iter);

    Real const t1 = second();
    Vector<MultiFab*> sol{&solution};
    Vector<MultiFab const*> rh{&rhs};
    Real const residual = mlmg.solve(sol, rh, p.reltol, p.abstol);
    Real const t2 = second();

    auto active_mask = make_eb_active_node_mask(nba, dmap, geom, factory);
    auto fluid_node_mask = make_sphere_fluid_node_mask(nba, dmap, geom, *active_mask,
                                                       Real(0.33), Real(0.92), Real(0.50));
    auto const [l2, linf] = compute_error_norms(solution, exact, geom,
                                                fluid_node_mask.get(), all_neumann);
    Real const t3 = second();

    SolveResult r;
    r.n_cell = n_cell;
    r.nprocs = ParallelDescriptor::NProcs();
    r.mg_levels = linop.NMGLevels(0);
    r.iterations = mlmg.getNumIters();
    r.residual = residual;
    r.l2 = l2;
    r.linf = linf;
    r.setup_time = t1 - t0;
    r.solve_time = t2 - t1;
    r.total_time = t3 - t0;
    return r;
#else
    amrex::ignore_unused(p);
    amrex::ignore_unused(n_cell);
    amrex::ignore_unused(all_neumann);
    amrex::ignore_unused(add_incompatible_offset);
    Abort("NodeRZEBBaseline: rz_eb_neumann physical BC cases require AMREX_USE_EB");
    return {};
#endif
}

void run_rz_eb_rap_zero (Parameters const& p)
{
    int const n_cell = p.n_cell.empty() ? 32 : p.n_cell.front();
    if (n_cell < 8) {
        Abort("NodeRZEBBaseline: rz_eb_rap_zero n_cell must be at least 8");
    }

    SolveResult r = run_rz_eb_rap_zero_case(p, n_cell);
    report_result("RZ_EB_RAP_ZERO", r);

    if (p.do_assert &&
        (!std::isfinite(r.residual) || !std::isfinite(r.l2) || !std::isfinite(r.linf))) {
        Abort("NodeRZEBBaseline: non-finite RZ_EB_RAP_ZERO residual or error norm");
    }
    if (p.do_assert && r.iterations > p.max_iter) {
        Abort("NodeRZEBBaseline: invalid RZ_EB_RAP_ZERO MLMG iteration count");
    }
    if (p.do_assert && r.mg_levels < p.min_mg_levels) {
        Abort("NodeRZEBBaseline: RZ_EB_RAP_ZERO MG level count below threshold");
    }
    if (p.do_assert && (r.l2 > p.max_finest_l2 || r.linf > p.max_finest_linf)) {
        Abort("NodeRZEBBaseline: RZ_EB_RAP_ZERO zero-solution error above threshold");
    }

    Print() << "RZ_EB_RAP_ZERO_PASS"
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << " mg_levels=" << r.mg_levels
            << "\n";
}

void run_rz_eb_neumann_plane (Parameters const& p)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    for (int n : n_cells) {
        if (n < 8) {
            Abort("NodeRZEBBaseline: rz_eb_neumann_plane n_cell must be at least 8");
        }
        SolveResult r = run_rz_eb_neumann_plane_case(p, n);
        report_result("RZ_EB_NEUMANN_PLANE", r);
        results.push_back(r);
    }

    assert_mms_results(p, results, "RZ_EB_NEUMANN_PLANE");
    Print() << "RZ_EB_NEUMANN_PLANE_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

void run_rz_eb_neumann_sphere (Parameters const& p)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    for (int n : n_cells) {
        if (n < 16) {
            Abort("NodeRZEBBaseline: rz_eb_neumann_sphere n_cell must be at least 16");
        }
        SolveResult r = run_rz_eb_neumann_sphere_case(p, n);
        report_result("RZ_EB_NEUMANN_SPHERE", r);
        results.push_back(r);
    }

    assert_mms_results(p, results, "RZ_EB_NEUMANN_SPHERE");
    Print() << "RZ_EB_NEUMANN_SPHERE_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

void run_rz_eb_inflow_sphere (Parameters const& p)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    for (int n : n_cells) {
        if (n < 16) {
            Abort("NodeRZEBBaseline: rz_eb_inflow_sphere n_cell must be at least 16");
        }
        SolveResult r = run_rz_eb_inflow_sphere_case(p, n);
        report_result("RZ_EB_INFLOW_SPHERE", r);
        results.push_back(r);
    }

    assert_mms_results(p, results, "RZ_EB_INFLOW_SPHERE");
    Print() << "RZ_EB_INFLOW_SPHERE_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

void run_cartesian_eb_neumann_sphere (Parameters const& p)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    for (int n : n_cells) {
        if (n < 16) {
            Abort("NodeRZEBBaseline: cart_eb_neumann_sphere n_cell must be at least 16");
        }
        SolveResult r = run_cartesian_eb_neumann_sphere_case(p, n);
        report_result("CART_EB_NEUMANN_SPHERE", r);
        results.push_back(r);
    }

    assert_mms_results(p, results, "CART_EB_NEUMANN_SPHERE");
    Print() << "CART_EB_NEUMANN_SPHERE_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

void run_rz_eb_neumann_physical_bc (Parameters const& p, bool all_neumann,
                                    bool add_incompatible_offset = false)
{
    std::vector<int> n_cells = p.n_cell;
    std::sort(n_cells.begin(), n_cells.end());
    n_cells.erase(std::unique(n_cells.begin(), n_cells.end()), n_cells.end());

    std::vector<SolveResult> results;
    results.reserve(n_cells.size());

    char const* label = add_incompatible_offset ? "RZ_EB_NEUMANN_ALL_NEUMANN_OFFSET"
        : (all_neumann ? "RZ_EB_NEUMANN_ALL_NEUMANN"
                       : "RZ_EB_NEUMANN_MIXED_BC");
    for (int n : n_cells) {
        if (n < 16) {
            Abort(std::string("NodeRZEBBaseline: ") + label +
                  " n_cell must be at least 16");
        }
        SolveResult r = run_rz_eb_neumann_physical_bc_case(p, n, all_neumann,
                                                           add_incompatible_offset);
        report_result(label, r);
        results.push_back(r);
    }

    assert_mms_results(p, results, label);
    Print() << label << "_PASS cases=" << results.size()
            << " mpi_ranks=" << ParallelDescriptor::NProcs()
            << "\n";
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        BL_PROFILE("main");
        Parameters const p = read_parameters();
        if (p.mode == "rz_mms") {
            run_rz_mms(p);
        } else if (p.mode == "rz_mms_rap") {
            run_rz_mms_rap(p);
        } else if (p.mode == "cart_mms_rap") {
            run_cartesian_mms(p);
        } else if (p.mode == "rz_eb_rap_zero") {
            run_rz_eb_rap_zero(p);
        } else if (p.mode == "rz_eb_neumann_plane") {
            run_rz_eb_neumann_plane(p);
        } else if (p.mode == "rz_eb_neumann_sphere") {
            run_rz_eb_neumann_sphere(p);
        } else if (p.mode == "rz_eb_inflow_sphere") {
            run_rz_eb_inflow_sphere(p);
        } else if (p.mode == "cart_eb_neumann_sphere") {
            run_cartesian_eb_neumann_sphere(p);
        } else if (p.mode == "rz_eb_neumann_mixed_bc") {
            run_rz_eb_neumann_physical_bc(p, false);
        } else if (p.mode == "rz_eb_neumann_all_neumann") {
            run_rz_eb_neumann_physical_bc(p, true);
        } else if (p.mode == "rz_eb_neumann_all_neumann_offset") {
            run_rz_eb_neumann_physical_bc(p, true, true);
        } else {
            Abort("NodeRZEBBaseline: unknown mode '" + p.mode + "'");
        }
    }
    amrex::Finalize();
}
