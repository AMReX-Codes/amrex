#include "MyTest.H"

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#endif
#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <numbers>

using namespace amrex;

namespace {

struct ExactSolution
{
    AMREX_GPU_HOST_DEVICE
    Real operator() (AMREX_D_DECL(Real x, Real y, Real z)) const noexcept
    {
        constexpr Real tpi = 2.0_rt * std::numbers::pi_v<Real>;
        constexpr Real fpi = 4.0_rt * std::numbers::pi_v<Real>;

        Real phi_two_pi = std::cos(tpi*x);
        Real phi_four_pi = std::cos(fpi*x);
#if (AMREX_SPACEDIM > 1)
        phi_two_pi *= std::cos(tpi*y);
        phi_four_pi *= std::cos(fpi*y);
#endif
#if (AMREX_SPACEDIM > 2)
        phi_two_pi *= std::cos(tpi*z);
        phi_four_pi *= std::cos(fpi*z);
#endif
        return phi_two_pi + 0.25_rt*phi_four_pi;
    }
};

#ifdef AMREX_USE_EB
void unscaleEBRHS (MultiFab& rhs, EBFArrayBoxFactory const& factory)
{
    auto const& edgecent = factory.getEdgeCent();

    for (MFIter mfi(rhs, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        if (edgecent[0]->ok(mfi)) {
            Box const& bx = mfi.tilebox();
            Array4<Real> const rhsarr = rhs.array(mfi);
            GpuArray<Array4<Real const>,AMREX_SPACEDIM> const ec
                {AMREX_D_DECL(edgecent[0]->const_array(mfi),
                              edgecent[1]->const_array(mfi),
                              edgecent[2]->const_array(mfi))};
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
#if (AMREX_SPACEDIM == 2)
                amrex::ignore_unused(k);
#endif
                IntVect const iv{AMREX_D_DECL(i,j,k)};
                Real scale = 1.0_rt;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    IntVect edgem = iv;
                    --edgem[idim];
                    Real const hp = (ec[idim](iv) == 1.0_rt)
                        ? 1.0_rt : 1.0_rt+2.0_rt*ec[idim](iv);
                    Real const hm = (ec[idim](edgem) == 1.0_rt)
                        ? 1.0_rt : 1.0_rt-2.0_rt*ec[idim](edgem);
                    scale = amrex::min(scale, hm, hp);
                }
                rhsarr(iv) /= scale;
            });
        }
    }
}
#endif

}

MyTest::MyTest ()
{
    readParameters();
    initData();
}

void
MyTest::run ()
{
    solve(false);
    computeNorms("scalar sigma");

    solve(true);
    computeNorms("variable sigma");
}

void
MyTest::solve (bool variable_sigma)
{
    BL_PROFILE("FDNodalPoisson::solve()");

    LPInfo info;
    info.setAgglomeration(agglomeration);
    info.setConsolidation(consolidation);
    info.setSemicoarsening(semicoarsening);
    info.setMaxCoarseningLevel(max_coarsening_level);
    info.setMaxSemicoarseningLevel(max_semicoarsening_level);

#ifdef AMREX_USE_EB
    MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info, {factory.get()});
#else
    MLEBNodeFDLaplacian linop({geom}, {grids}, {dmap}, info);
#endif

    Array<LinOpBCType,AMREX_SPACEDIM> lobc
        {AMREX_D_DECL(LinOpBCType::Dirichlet,
                      LinOpBCType::Dirichlet,
                      LinOpBCType::Dirichlet)};
    Array<LinOpBCType,AMREX_SPACEDIM> hibc
        {AMREX_D_DECL(LinOpBCType::Dirichlet,
                      LinOpBCType::Dirichlet,
                      LinOpBCType::Dirichlet)};
#if (AMREX_SPACEDIM == 2)
    if (rz) {
        lobc[0] = LinOpBCType::Neumann;
        linop.setRZ(true);
    }
#endif
    linop.setDomainBC(lobc, hibc);
#ifdef AMREX_USE_EB
    linop.setEBDirichlet(ExactSolution{});
#endif

    if (variable_sigma) {
        linop.setSigma(0, sigma);
    } else {
        linop.setSigma(Array<Real,AMREX_SPACEDIM>
                       {{AMREX_D_DECL(1.0_rt, 1.0_rt, 1.0_rt)}});
    }

    MLMG mlmg(linop);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    // Generate a right-hand side that is exactly consistent with the
    // discretized operator used by this test.
    mlmg.apply({&rhs}, {&exact_solution});
#ifdef AMREX_USE_EB
    // Fapply includes the cut-node conditioning factor, and MLMG applies that
    // factor to the physical RHS when solve starts.  Undo it here so that this
    // manufactured RHS is conditioned exactly once.
    unscaleEBRHS(rhs, *factory);
#endif

    // The solution provides both the initial guess and the inhomogeneous
    // Dirichlet data.  Keep the exact boundary values and zero the interior.
    MultiFab::Copy(solution, exact_solution, 0, 0, 1, 0);
    Box const& interior = amrex::surroundingNodes(amrex::grow(geom.Domain(), -1));
    solution.setVal(0.0, interior, 0, 1, 0);

    mlmg.solve({&solution}, {&rhs}, reltol, 0.0);
}

void
MyTest::computeNorms (std::string const& label) const
{
    MultiFab error(solution.boxArray(), solution.DistributionMap(), 1, 0);
    MultiFab::Copy(error, solution, 0, 0, 1, 0);
    MultiFab::Subtract(error, exact_solution, 0, 0, 1, 0);

    auto mask = error.OwnerMask(geom.periodicity());
    Real const max_error = error.norm0(*mask, 0, 0);

    auto const dx = geom.CellSizeArray();
    Real const dvol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    Real const l1_error = error.norm1(0, geom.periodicity()) * dvol;

    amrex::Print() << label << "\n"
                   << "    max-norm: " << max_error << "\n"
                   << "    1-norm  : " << l1_error << "\n";

#ifdef AMREX_USE_FLOAT
    constexpr Real error_tolerance = 1.e-3_rt;
#else
    constexpr Real error_tolerance = 1.e-8_rt;
#endif

    if (!error.is_finite() ||
        !std::isfinite(max_error) || !std::isfinite(l1_error) ||
        max_error > error_tolerance) {
        amrex::Abort("FDNodalPoisson " + label + " error check failed");
    }
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("reltol", reltol);
#ifdef AMREX_USE_FLOAT
    reltol = std::max(reltol, 1.e-5_rt);
#endif

    pp.query("agglomeration", agglomeration);
    pp.query("consolidation", consolidation);
    pp.query("semicoarsening", semicoarsening);
#if (AMREX_SPACEDIM == 2)
    pp.query("rz", rz);
#endif
    pp.query("max_coarsening_level", max_coarsening_level);
    pp.query("max_semicoarsening_level", max_semicoarsening_level);
}

void
MyTest::initData ()
{
    RealBox rb({AMREX_D_DECL(0.0, 0.0, 0.0)},
               {AMREX_D_DECL(1.0, 1.0, 1.0)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
    Geometry::Setup(&rb, CoordSys::cartesian, is_periodic.data());

    Box domain(IntVect{AMREX_D_DECL(0, 0, 0)},
               IntVect{AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1)});
    geom.define(domain);

    grids.define(domain);
    grids.maxSize(max_grid_size);
    dmap.define(grids);

#ifdef AMREX_USE_EB
    EB2::SphereIF sphere(0.2_rt, {AMREX_D_DECL(0.5_rt, 0.5_rt, 0.5_rt)}, false);
    auto shop = EB2::makeShop(sphere);
    EB2::Build(shop, geom, 0, max_coarsening_level);
    factory = makeEBFabFactory(geom, grids, dmap, {2,2,2}, EBSupport::full);
#endif

    BoxArray const& nodal_grids = amrex::convert(grids, IntVect::TheNodeVector());
    solution.define(nodal_grids, dmap, 1, 0);
    rhs.define(nodal_grids, dmap, 1, 0);
    exact_solution.define(nodal_grids, dmap, 1, 0);
    sigma.define(grids, dmap, 1, 0);

    auto const dx = geom.CellSizeArray();
    auto const problo = geom.ProbLoArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact_solution, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const phi = exact_solution.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            phi(i,j,k) = ExactSolution{}
                (AMREX_D_DECL(problo[0] + i*dx[0],
                              problo[1] + j*dx[1],
                              problo[2] + k*dx[2]));
        });
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(sigma, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box const& bx = mfi.tilebox();
        Array4<Real> const sig = sigma.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            constexpr Real tpi = 2.0_rt * std::numbers::pi_v<Real>;

            Real variation = std::cos(tpi * (problo[0] + (i+0.5_rt)*dx[0]));
#if (AMREX_SPACEDIM > 1)
            variation += std::cos(tpi * (problo[1] + (j+0.5_rt)*dx[1]));
#endif
#if (AMREX_SPACEDIM > 2)
            variation += std::cos(tpi * (problo[2] + (k+0.5_rt)*dx[2]));
#endif
            sig(i,j,k) = 1.0_rt + 0.1_rt*variation;
        });
    }
}
