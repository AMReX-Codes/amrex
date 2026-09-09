#include "MyTest.H"

#include <AMReX_EB2.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_ParmParse.H>

#include <algorithm>
#include <cmath>
#include <numbers>

#if (AMREX_SPACEDIM == 1)
#error NodalPoissonEB is a 2D/3D test.
#endif

using namespace amrex;

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_phi (AMREX_D_DECL(Real x, Real y, Real z)) noexcept
{
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real tpi = Real(2.0)*pi;
    constexpr Real fpi = Real(4.0)*pi;
    Real const mode2 = AMREX_D_TERM(std::cos(tpi*x), * std::cos(tpi*y), * std::cos(tpi*z));
    Real const mode4 = AMREX_D_TERM(std::cos(fpi*x), * std::cos(fpi*y), * std::cos(fpi*z));
    return mode2 + Real(0.25) * mode4;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_rhs (AMREX_D_DECL(Real x, Real y, Real z)) noexcept
{
    constexpr Real pi = std::numbers::pi_v<Real>;
    constexpr Real tpi = Real(2.0)*pi;
    constexpr Real fpi = Real(4.0)*pi;
    constexpr Real fac = tpi*tpi*Real(AMREX_SPACEDIM);
    Real const mode2 = AMREX_D_TERM(std::cos(tpi*x), * std::cos(tpi*y), * std::cos(tpi*z));
    Real const mode4 = AMREX_D_TERM(std::cos(fpi*x), * std::cos(fpi*y), * std::cos(fpi*z));
    return -fac * mode2 - fac * mode4;
}

}

MyTest::MyTest ()
{
    readParameters();
    initData();
}

void
MyTest::solve ()
{
    BL_PROFILE("NodalPoissonEB::solve()");

    LPInfo info;
    info.setMaxCoarseningLevel(max_coarsening_level);

    MLEBNodeFDLaplacian linop(geom, grids, dmap, info,
                              {factory.get()});

    linop.setDomainBC(
        {AMREX_D_DECL(LinOpBCType::Dirichlet, LinOpBCType::Dirichlet, LinOpBCType::Dirichlet)},
        {AMREX_D_DECL(LinOpBCType::Dirichlet, LinOpBCType::Dirichlet, LinOpBCType::Dirichlet)});

    linop.setEBDirichlet([] AMREX_GPU_HOST_DEVICE (AMREX_D_DECL(Real x, Real y, Real z)) noexcept
    {
        return exact_phi(AMREX_D_DECL(x, y, z));
    });

    MLMG mlmg(linop);
    mlmg.setMaxIter(max_iter);
    mlmg.setMaxFmgIter(max_fmg_iter);
    mlmg.setVerbose(verbose);
    mlmg.setBottomVerbose(bottom_verbose);

    MultiFab::Copy(solution[0], exact_solution[0], 0, 0, 1, 0);
    const Box interior = amrex::grow(amrex::surroundingNodes(geom[0].Domain()), -1);
    solution[0].setVal(0.0, interior, 0, 1, 0);

    mlmg.solve(GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), reltol, 0.0);

    mlmg.getGradSolution({GetArrOfPtrs(grad_solution)});
}

void
MyTest::compute_norms () const
{
    MultiFab error(solution[0].boxArray(), solution[0].DistributionMap(), 1, 0);
    MultiFab::Copy(error, solution[0], 0, 0, 1, 0);
    MultiFab::Subtract(error, exact_solution[0], 0, 0, 1, 0);

    auto mask = error.OwnerMask(geom[0].periodicity());

    amrex::Print() << "max-norm: " << error.norm0(*mask, 0, 0) << "\n";
    const Real* dx = geom[0].CellSize();
    const Real dvol = AMREX_D_TERM(dx[0], * dx[1], * dx[2]);
    amrex::Print() << "1-norm  : " << error.norm1(0, geom[0].periodicity())*dvol << "\n";
}

void
MyTest::readParameters ()
{
    ParmParse pp;
    pp.query("n_cell", n_cell);
    pp.query("max_grid_size", max_grid_size);
    pp.query("max_coarsening_level", max_coarsening_level);

    pp.query("verbose", verbose);
    pp.query("bottom_verbose", bottom_verbose);
    pp.query("max_iter", max_iter);
    pp.query("max_fmg_iter", max_fmg_iter);
    pp.query("reltol", reltol);
#ifdef AMREX_USE_FLOAT
    reltol = std::max(reltol, 1.e-5F);
#endif

    pp.query("gpu_regtest", gpu_regtest);
    pp.query("do_plots", do_plots);
}

void
MyTest::initData ()
{
    geom.resize(1);
    grids.resize(1);
    dmap.resize(1);
    solution.resize(1);
    rhs.resize(1);
    exact_solution.resize(1);

    RealBox rb({AMREX_D_DECL(0.0,0.0,0.0)}, {AMREX_D_DECL(1.0,1.0,1.0)});
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0,0,0)};
    Geometry::Setup(&rb, 0, is_periodic.data());

    Box domain(IntVect{AMREX_D_DECL(0,0,0)},
               IntVect{AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)});
    geom[0].define(domain);

    grids[0].define(domain);
    grids[0].maxSize(max_grid_size);
    dmap[0].define(grids[0]);

    EB2::Build(geom[0], 0, max_coarsening_level);

    factory = makeEBFabFactory(geom[0], grids[0], dmap[0],
                               {2,2,2}, EBSupport::full);

    const BoxArray& nba = amrex::convert(grids[0], IntVect::TheNodeVector());

    solution      [0].define(nba, dmap[0], 1, 1, MFInfo(), *factory);
    rhs           [0].define(nba, dmap[0], 1, 0, MFInfo(), *factory);
    exact_solution[0].define(nba, dmap[0], 1, 0, MFInfo(), *factory);

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        IntVect typ = IntVect::TheNodeVector();
        typ[idim] = 0;
        const BoxArray& eba = amrex::convert(grids[0], typ);
        grad_solution[idim].define(eba, dmap[0], 1, 0, MFInfo(), *factory);
    }

    rhs[0].setVal(0.0);

    const auto dx = geom[0].CellSizeArray();
    const auto problo = geom[0].ProbLoArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(exact_solution[0], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<Real> const phi = exact_solution[0].array(mfi);
        Array4<Real> const rh = rhs[0].array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            AMREX_D_TERM(const Real x = problo[0] + i * dx[0];,
                         const Real y = problo[1] + j * dx[1];,
                         const Real z = problo[2] + k * dx[2];)
            phi(i,j,k) = exact_phi(AMREX_D_DECL(x, y, z));
            rh(i,j,k) = exact_rhs(AMREX_D_DECL(x, y, z));
        });
    }
}
