// Test of MLNodeCrossLaplacian: the projection u -= sigma G phi with
// D sigma G phi = D u must give D u = 0 to roundoff, which also makes it
// idempotent, and must not depend on the box layout or on sigma.

#include <AMReX.H>
#include <AMReX_GMRES_MLMG.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLNodeCrossLaplacian.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>

#include <cmath>
#include <cstdio>
#include <string>

using namespace amrex;

namespace {

struct Case {
    std::string bc;   // periodic, neumann, dirichlet, inflow_outflow, mixed
    Real aspect;      // dy/dx
};

struct Params {
    int n_cell = 64; // set per run from n_cells
    int max_grid_size = 16;
    Real sigma = Real(2.5);
    Real reltol = (sizeof(Real) == 4) ? Real(1.e-5) : Real(1.e-11);
    int verbose = 0;
    std::string multigrid_type = "algebraic";
};

struct Setup {
    Geometry geom;
    Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
};

Setup make_setup (Case const& c, int n_cell)
{
    Setup s;
    bool perx = false, pery = false;
    auto& lo = s.lobc;
    auto& hi = s.hibc;
    if (c.bc == "periodic") {
        perx = pery = true;
        lo = {LinOpBCType::Periodic, LinOpBCType::Periodic};
        hi = lo;
    } else if (c.bc == "neumann") {
        lo = {LinOpBCType::Neumann, LinOpBCType::Neumann};
        hi = lo;
    } else if (c.bc == "dirichlet") {
        lo = {LinOpBCType::Dirichlet, LinOpBCType::Dirichlet};
        hi = lo;
    } else if (c.bc == "inflow_outflow") { // periodic in x, inflow at y-lo, outflow at y-hi
        perx = true;
        lo = {LinOpBCType::Periodic, LinOpBCType::inflow};
        hi = {LinOpBCType::Periodic, LinOpBCType::Dirichlet};
    } else if (c.bc == "mixed") { // walls in x, periodic in y
        pery = true;
        lo = {LinOpBCType::Neumann, LinOpBCType::Periodic};
        hi = lo;
    } else {
        amrex::Abort("Unknown bc " + c.bc);
    }
    Box const domain(IntVect(0), IntVect(n_cell-1));
    RealBox const rb({Real(0),Real(0)}, {Real(1),c.aspect});
    s.geom.define(domain, rb, CoordSys::cartesian, {int(perx), int(pery)});
    return s;
}

// Cell-centered velocity; ghost cells hold zero at walls and outflow, and
// the inflow velocity (plus a transverse value that must be ignored) at inflow.
void fill_velocity (MultiFab& vel, Setup const& s, bool valid)
{
    Geometry const& geom = s.geom;
    auto const problo = geom.ProbLoArray();
    auto const dx = geom.CellSizeArray();
    Real const ly = geom.ProbLength(1);
    Box const domain = geom.Domain();
    auto const lobc = s.lobc;
    bool const perx = geom.isPeriodic(0);
    bool const pery = geom.isPeriodic(1);
    auto const& va = vel.arrays();
    if (valid) {
        ParallelFor(vel, [=] AMREX_GPU_DEVICE (int b, int i, int j, int) noexcept
        {
            Real const x = problo[0] + (Real(i)+Real(0.5))*dx[0];
            Real const y = (problo[1] + (Real(j)+Real(0.5))*dx[1]) / ly;
            Real const pi = Math::pi<Real>();
            Real const chk = ((i+j) & 1) ? Real(0.1) : Real(-0.1);
            va[b](i,j,0,0) = std::sin(Real(2)*pi*x)*std::cos(Real(4)*pi*y) + Real(0.5)*std::cos(Real(2)*pi*(x+y)) + chk;
            va[b](i,j,0,1) = std::cos(Real(6)*pi*x)*std::sin(Real(2)*pi*y) - Real(0.3)*std::sin(Real(4)*pi*x) - chk;
        });
    }
    vel.FillBoundary(geom.periodicity());
    ParallelFor(vel, vel.nGrowVect(), [=] AMREX_GPU_DEVICE (int b, int i, int j, int) noexcept
    {
        bool const outx = (i < domain.smallEnd(0) || i > domain.bigEnd(0)) && !perx;
        bool const outy = (j < domain.smallEnd(1) || j > domain.bigEnd(1)) && !pery;
        if (outx || outy) {
            if (j < domain.smallEnd(1) && lobc[1] == LinOpBCType::inflow && !outx) {
                Real const x = problo[0] + (Real(i)+Real(0.5))*dx[0];
                va[b](i,j,0,0) = Real(0.123); // transverse, ignored
                va[b](i,j,0,1) = Real(1.0) + Real(0.5)*std::sin(Real(2)*Math::pi<Real>()*x);
            } else {
                va[b](i,j,0,0) = Real(0.0);
                va[b](i,j,0,1) = Real(0.0);
            }
        }
    });
    Gpu::streamSynchronize();
}

struct Result {
    MultiFab u0, u1;    // before and after the projection
    Real div0, div1;    // max |D u| before and after
    Real chk;           // checkerboard component of phi relative to max |phi|
    int iters;
};

// Project the velocity.
Result project (Case const& c, Params const& p, int max_grid_size, Real sigma)
{
    Setup const s = make_setup(c, p.n_cell);
    Geometry const& geom = s.geom;
    BoxArray ba(geom.Domain());
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);
    BoxArray const nba = amrex::convert(ba, IntVect(1));

    MLNodeCrossLaplacian linop({geom}, {ba}, {dm}, LPInfo(), sigma);
    linop.setDomainBC(s.lobc, s.hibc);

    MLMG mlmg(linop);
    mlmg.setMultigridType(amrex::getEnumCaseInsensitive<MultigridType>(p.multigrid_type));
    mlmg.setVerbose(p.verbose);
    mlmg.setBottomVerbose(p.verbose > 1 ? p.verbose : 0);
    GMRESMLMG gmres(mlmg);
    gmres.setVerbose(p.verbose);
    gmres.setMaxIters(200);
    bool const use_gmres = (c.aspect != Real(1));

    Result r;
    r.u0.define(ba, dm, AMREX_SPACEDIM, 1);
    r.u1.define(ba, dm, AMREX_SPACEDIM, 1);
    MultiFab rhs(nba, dm, 1, 0);
    MultiFab phi(nba, dm, 1, 1);

    fill_velocity(r.u0, s, true);
    MultiFab::Copy(r.u1, r.u0, 0, 0, AMREX_SPACEDIM, 1);

    linop.compRHS({&rhs}, {&r.u1});
    r.div0 = rhs.norminf();
    phi.setVal(0.0);
    if (use_gmres) {
        gmres.solve(phi, rhs, p.reltol, Real(0));
        r.iters = gmres.getNumIters();
    } else {
        mlmg.solve({&phi}, {&rhs}, p.reltol, Real(0));
        r.iters = mlmg.getNumIters();
    }
    linop.updateVelocity({&r.u1}, {&phi});

    // Checkerboard component of phi (removed for singular problems).
    {
        auto const omask = phi.OwnerMask(geom.periodicity());
        auto const& pa = phi.const_arrays();
        auto const& ma = omask->const_arrays();
        auto sums = ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
            phi, IntVect(0), [=] AMREX_GPU_DEVICE (int b, int i, int j, int) noexcept
                -> GpuTuple<Real,Real>
            {
                Real const m = ma[b](i,j,0) ? Real(1) : Real(0);
                return { m * (((i+j) & 1) ? Real(-1) : Real(1)) * pa[b](i,j,0), m };
            });
        Real s1 = amrex::get<0>(sums), s2 = amrex::get<1>(sums);
        ParallelAllReduce::Sum<Real>({s1,s2}, ParallelContext::CommunicatorSub());
        r.chk = std::abs(s1/s2) / phi.norminf();
    }

    // D u1 with the same boundary ghost cells as u0.
    MultiFab u(ba, dm, AMREX_SPACEDIM, 1);
    MultiFab::Copy(u, r.u1, 0, 0, AMREX_SPACEDIM, 0);
    fill_velocity(u, s, false);
    linop.compRHS({&rhs}, {&u});
    r.div1 = rhs.norminf();
    return r;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    int nfail = 0;
    {
        Params p;
        {
            ParmParse pp;
            pp.query("max_grid_size", p.max_grid_size);
            pp.query("sigma", p.sigma);
            pp.query("reltol", p.reltol);
            pp.query("verbose", p.verbose);
            pp.query("multigrid_type", p.multigrid_type);
        }
        Vector<std::string> bcs{"periodic", "neumann", "dirichlet", "inflow_outflow", "mixed"};
        Vector<Real> aspects{Real(1), Real(2)};
        {
            ParmParse pp;
            if (pp.contains("bc")) { bcs.clear(); pp.getarr("bc", bcs); }
            if (pp.contains("aspect")) { aspects.clear(); pp.getarr("aspect", aspects); }
        }

        // An odd n_cell makes the checkerboard a non-null mode with periodic BCs.
        Vector<int> n_cells{64, 63};
        {
            ParmParse pp;
            if (pp.contains("n_cell")) { n_cells.clear(); pp.getarr("n_cell", n_cells); }
        }

        Real const tol = Real(100) * p.reltol;
        for (int n_cell : n_cells) {
            p.n_cell = n_cell;
            for (auto const& bc : bcs) {
                for (auto aspect : aspects) {
                    Case const c{.bc = bc, .aspect = aspect};
                    // Reference: one box, sigma = 1.
                    Result ref = project(c, p, p.n_cell, Real(1));
                    Result res = project(c, p, p.max_grid_size, p.sigma);

                    MultiFab uref(res.u1.boxArray(), res.u1.DistributionMap(), AMREX_SPACEDIM, 0);
                    uref.ParallelCopy(ref.u1, 0, 0, AMREX_SPACEDIM);
                    MultiFab::Subtract(uref, res.u1, 0, 0, AMREX_SPACEDIM, 0);
                    MultiFab du(res.u1.boxArray(), res.u1.DistributionMap(), AMREX_SPACEDIM, 0);
                    MultiFab::LinComb(du, Real(1), res.u1, 0, Real(-1), res.u0, 0, 0, AMREX_SPACEDIM, 0);
                    Real const layout = uref.norminf(0, AMREX_SPACEDIM, IntVect(0)) / du.norminf(0, AMREX_SPACEDIM, IntVect(0));

                    Real const div = res.div1 / res.div0;
                    bool const chk_null = (bc == "neumann") ||
                        ((bc == "periodic" || bc == "mixed") && n_cell % 2 == 0);
                    // With dx == dy the AlgMG matrix is the operator, so MLMG
                    // needs only a couple of iterations (in double precision).
                    bool const ok = div <= tol && layout <= tol
                        && (!chk_null || res.chk <= tol)
                        && (aspect != Real(1) || p.multigrid_type != "algebraic"
                            || sizeof(Real) == 4 || res.iters <= 3);
                    if (!ok) { ++nfail; }
                    amrex::Print() << (ok ? "PASS " : "FAIL ") << "n=" << n_cell << " " << bc << " dy/dx=" << aspect
                                   << (aspect != Real(1) ? " GMRES" : " MLMG") << " iters " << res.iters
                                   << " |Du1|/|Du0| " << div
                                   << " layout/sigma " << layout
                                   << " phi checkerboard " << res.chk << "\n";
                }
            }
        }
    }
    amrex::Finalize();
    if (nfail > 0) {
        std::printf("%d case(s) failed\n", nfail);
        return 1;
    }
    return 0;
}
