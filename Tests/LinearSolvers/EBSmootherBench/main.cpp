// Benchmark and reproducibility check for the MLEBABecLap Gauss-Seidel smoother.
//
// For one EB geometry and box decomposition, it reports
//   * the box census (cut/regular/covered) and whether Fsmooth takes the fused
//     MultiFab path (FabArrayBase::isFusingCandidate) or the per-box path,
//   * the time of one red+black Fsmooth pair and of one smooth() sweep
//     (median over nrepeat),
//   * the wall time, iteration count and final residual of an MLMG solve,
//   * whether two identical smooth() calls and two identical MLMG solves
//     give bitwise identical results.

#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_MLMG.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <utility>
#include <vector>

using namespace amrex;

namespace {

    Real median (std::vector<Real> v)
    {
        std::sort(v.begin(), v.end());
        return v[v.size()/2];
    }

    // Number of valid cells whose bit patterns differ, and max |a-b|.
    std::pair<Long,Real> compare (MultiFab const& a, MultiFab const& b)
    {
        auto const& ama = a.const_arrays();
        auto const& bma = b.const_arrays();
        Long ndiff = ParReduce(TypeList<ReduceOpSum>{}, TypeList<Long>{}, a, IntVect(0),
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) -> GpuTuple<Long>
        {
            Real x = ama[bno](i,j,k);
            Real y = bma[bno](i,j,k);
            std::uint64_t ux, uy;
            std::memcpy(&ux, &x, sizeof(Real));
            std::memcpy(&uy, &y, sizeof(Real));
            return { (ux != uy) ? Long(1) : Long(0) };
        });
        MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
        MultiFab::Copy(d, a, 0, 0, 1, 0);
        MultiFab::Subtract(d, b, 0, 0, 1, 0);
        return {ndiff, d.norminf(0)};
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 128;
        int max_grid_size = 32;
        std::string geom_type = "sphere";
        Real radius = 0.3;
        Vector<Real> center{AMREX_D_DECL(0.5_rt, 0.5_rt, 0.5_rt)};
        bool deterministic = false;
        bool bottom_smoother = false;
        int nsweeps = 20;
        int nwarmup = 5;
        int nrepeat = 7;
        int nsolve = 5;
        Real tol_rel = 1.e-10;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("geom_type", geom_type);
            pp.query("radius", radius);
            pp.queryarr("center", center);
            pp.query("deterministic", deterministic);
            pp.query("bottom_smoother", bottom_smoother);
            pp.query("nsweeps", nsweeps);
            pp.query("nwarmup", nwarmup);
            pp.query("nrepeat", nrepeat);
            pp.query("nsolve", nsolve);
            pp.query("tol_rel", tol_rel);
        }

        Geometry geom(Box(IntVect(0), IntVect(n_cell-1)),
                      RealBox({AMREX_D_DECL(0._rt,0._rt,0._rt)}, {AMREX_D_DECL(1._rt,1._rt,1._rt)}),
                      CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});

        const int max_coarsening_level = 30;
        RealArray c{AMREX_D_DECL(center[0], center[1], center[2])};
        if (geom_type == "cylinder") {
            EB2::CylinderIF cyl(radius, AMREX_SPACEDIM-1, c, false);
            EB2::Build(EB2::makeShop(cyl), geom, 0, max_coarsening_level);
        } else {
            EB2::SphereIF sphere(radius, c, false);
            EB2::Build(EB2::makeShop(sphere), geom, 0, max_coarsening_level);
        }

        BoxArray ba(geom.Domain());
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);
        auto factory = makeEBFabFactory(geom, ba, dm, {2,2,2}, EBSupport::full);

        int ncut = 0, nregular = 0, ncovered = 0;
        auto const& flags = factory->getMultiEBCellFlagFab();
        for (MFIter mfi(flags); mfi.isValid(); ++mfi) {
            auto t = flags[mfi].getType(mfi.validbox());
            if (t == FabType::regular) { ++nregular; }
            else if (t == FabType::covered) { ++ncovered; }
            else { ++ncut; }
        }

        MultiFab rhs(ba, dm, 1, 0, MFInfo(), *factory);
        MultiFab acoef(ba, dm, 1, 0, MFInfo(), *factory);
        MultiFab phieb(ba, dm, 1, 0, MFInfo(), *factory);
        Array<MultiFab,AMREX_SPACEDIM> bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef[idim].define(amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 0,
                               MFInfo(), *factory);
            bcoef[idim].setVal(1.0);
        }
        acoef.setVal(1.0);
        phieb.setVal(1.0);
        {
            auto const dx = geom.CellSizeArray();
            auto const& rma = rhs.arrays();
            ParallelFor(rhs, [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k)
            {
                Real x = (i+0.5_rt)*dx[0];
                Real y = (j+0.5_rt)*dx[1];
                Real z = (AMREX_SPACEDIM == 3) ? (k+0.5_rt)*dx[AMREX_SPACEDIM-1] : 0.5_rt;
                rma[bno](i,j,k) = std::sin(Real(6.283185307179586)*x)
                    *             std::cos(Real(6.283185307179586)*y)
                    *             std::cos(Real(3.141592653589793)*z);
            });
            EB_set_covered(rhs, 0.0);
        }

        LPInfo info;
        info.setDeterministic(deterministic);
        MLEBABecLap mlabec({geom}, {ba}, {dm}, info, {factory.get()});
        mlabec.setMaxOrder(2);
        mlabec.setDomainBC({AMREX_D_DECL(LinOpBCType::Dirichlet,LinOpBCType::Dirichlet,LinOpBCType::Dirichlet)},
                           {AMREX_D_DECL(LinOpBCType::Dirichlet,LinOpBCType::Dirichlet,LinOpBCType::Dirichlet)});
        MultiFab bcval(ba, dm, 1, 1, MFInfo(), *factory);
        bcval.setVal(0.0);
        mlabec.setLevelBC(0, &bcval);
        mlabec.setScalars(1.0, 1.0);
        mlabec.setACoeffs(0, acoef);
        mlabec.setBCoeffs(0, amrex::GetArrOfConstPtrs(bcoef));
        mlabec.setEBDirichlet(0, phieb, 1.0);
        mlabec.prepareForSolve();

        MultiFab x(ba, dm, 1, 1, MFInfo(), *factory);
        MultiFab x2(ba, dm, 1, 1, MFInfo(), *factory);
        const bool fused = Gpu::inLaunchRegion() && x.isFusingCandidate();

        // Fsmooth alone (one red + one black pass) and smooth() (with FillBoundary/BC).
        x.setVal(0.0);
        for (int n = 0; n < nwarmup; ++n) {
            mlabec.smooth(0, 0, x, rhs, false, 1);
        }
        std::vector<Real> t_fsmooth, t_smooth;
        for (int r = 0; r < nrepeat; ++r) {
            Gpu::streamSynchronize();
            Real t0 = amrex::second();
            for (int n = 0; n < nsweeps; ++n) {
                mlabec.Fsmooth(0, 0, x, rhs, 0);
                mlabec.Fsmooth(0, 0, x, rhs, 1);
            }
            Gpu::streamSynchronize();
            t_fsmooth.push_back((amrex::second()-t0)/nsweeps*1.e3);

            t0 = amrex::second();
            mlabec.smooth(0, 0, x, rhs, false, nsweeps);
            Gpu::streamSynchronize();
            t_smooth.push_back((amrex::second()-t0)/nsweeps*1.e3);
        }

        // Two smooth() calls on identical inputs.
        x.setVal(0.0);
        x2.setVal(0.0);
        mlabec.smooth(0, 0, x , rhs, false, nsweeps);
        mlabec.smooth(0, 0, x2, rhs, false, nsweeps);
        auto [smooth_ndiff, smooth_maxdiff] = compare(x, x2);

        // MLMG solves on identical inputs: one warm-up, then nsolve timed.
        MultiFab phi0(ba, dm, 1, 1, MFInfo(), *factory);
        MultiFab phi(ba, dm, 1, 1, MFInfo(), *factory);
        MLMG mlmg(mlabec);
        mlmg.setVerbose(0);
        mlmg.setMaxIter(200);
        if (bottom_smoother) { mlmg.setBottomSolver(MLMG::BottomSolver::smoother); }
        std::vector<Real> t_solve;
        Long solve_ndiff = 0;
        Real solve_maxdiff = 0.0;
        int iters = 0;
        bool iters_varied = false;
        Real resid = 0.0;
        for (int s = 0; s <= nsolve; ++s) {
            phi.setVal(0.0);
            Gpu::streamSynchronize();
            Real t0 = amrex::second();
            mlmg.solve({&phi}, {&rhs}, tol_rel, 0.0);
            Gpu::streamSynchronize();
            if (s == 0) {
                MultiFab::Copy(phi0, phi, 0, 0, 1, 0);
                iters = mlmg.getNumIters();
                resid = mlmg.getFinalResidual()/mlmg.getInitResidual();
            } else {
                t_solve.push_back((amrex::second()-t0)*1.e3);
                auto [nd, md] = compare(phi0, phi);
                solve_ndiff = std::max(solve_ndiff, nd);
                solve_maxdiff = std::max(solve_maxdiff, md);
                if (mlmg.getNumIters() != iters) { iters_varied = true; }
            }
        }

        amrex::Print() << "RESULT"
                       << " n_cell=" << n_cell << " max_grid_size=" << max_grid_size
                       << " geom=" << geom_type << " radius=" << radius
                       << " nboxes=" << ba.size() << " cut=" << ncut
                       << " regular=" << nregular << " covered=" << ncovered
                       << " fused=" << fused << " deterministic=" << deterministic
                       << " bottom_smoother=" << bottom_smoother
                       << " fsmooth_ms=" << median(t_fsmooth)
                       << " smooth_ms=" << median(t_smooth)
                       << " solve_ms=" << median(t_solve)
                       << " iters=" << iters << " iters_varied=" << iters_varied
                       << " resid=" << resid
                       << " smooth_ndiff=" << smooth_ndiff << " smooth_maxdiff=" << smooth_maxdiff
                       << " solve_ndiff=" << solve_ndiff << " solve_maxdiff=" << solve_maxdiff
                       << "\n";
    }
    amrex::Finalize();
}
