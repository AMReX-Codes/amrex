// MLABecLapDrift, MLABecLaplacian with the upwind drift term of setDriftCoeffs:
//
//   a_scalar a phi - b_scalar div(b grad phi) - c_scalar div(c phi) = rhs
//
// with Dirichlet and Neumann domain boundaries and face velocities that
// change sign and leave through the Neumann face. c is upwinded, which is
// stable for c_scalar of the sign opposite to b_scalar's (transport by c).
//
// - One level (max_level = 0): MLMG::apply of a known phi is compared with the
//   operator written out here, and the solve with that right-hand side must
//   give phi back.
// - Two levels (max_level = 1): the composite solve with the right-hand side
//   of MLMG::apply must give phi back (consistency of the smoother, the
//   coarse/fine fluxes and the reflux with the operator).

#include <AMReX.H>
#include <AMReX_MLABecLapDrift.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>

#include <cmath>

using namespace amrex;

static_assert(AMREX_SPACEDIM > 1, "ABecLapDrift is a 2D and 3D test");

namespace {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real exact_phi (RealVect const& x)
{
    return std::sin(Real(2.0)*x[0] + Real(0.5)) * std::cos(Real(1.5)*x[1])
#if (AMREX_SPACEDIM == 3)
        * (Real(1.0) + Real(0.3)*x[2]*x[2])
#endif
        + Real(0.2);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real velocity (int idim, RealVect const& x)
{
    constexpr Real pi = Real(3.141592653589793);
    if (idim == 0) { return std::cos(Real(2.0)*pi*x[1]) + Real(0.3); }
    // Out through the low y (Neumann) face, where it must not enter
    if (idim == 1) { return std::sin(Real(2.0)*pi*x[0]) - Real(1.2)*(Real(1.0) - x[1]); }
    return Real(0.5) - x[0];
}

// Ghost cell across a Dirichlet face (quadratic through the face value and
// the two cells inside; MLMG's max_order 3)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real dirichlet_ghost (Real bc, Real v0, Real v1)
{
    return Real(8.0/3.0)*bc - Real(2.0)*v0 + Real(1.0/3.0)*v1;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 32;
        int max_grid_size = 16;
        int max_level = 0;
        Real const ascalar = 1.0;
        Real const bscalar = 1.0;
        Real cscalar = -20.0;
        int use_gauss_seidel = 1;
        Real tol = 1.e-11;
        Real max_error = 1.e-8;
        int verbose = 1;
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("max_level", max_level);
            pp.query("cscalar", cscalar);
            pp.query("use_gauss_seidel", use_gauss_seidel);
            pp.query("tol", tol);
            pp.query("max_error", max_error);
            pp.query("verbose", verbose);
        }
        AMREX_ALWAYS_ASSERT(max_level == 0 || max_level == 1);
        int const nlevels = max_level + 1;

        // Neumann at the low y face, Dirichlet elsewhere
        Array<LinOpBCType,AMREX_SPACEDIM> lobc{AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                            LinOpBCType::Neumann,
                                                            LinOpBCType::Dirichlet)};
        Array<LinOpBCType,AMREX_SPACEDIM> hibc{AMREX_D_DECL(LinOpBCType::Dirichlet,
                                                            LinOpBCType::Dirichlet,
                                                            LinOpBCType::Dirichlet)};

        RealBox const rb(AMREX_D_DECL(0.,0.,0.), AMREX_D_DECL(1.,1.,1.));
        Vector<Geometry> geom(nlevels);
        Vector<BoxArray> grids(nlevels);
        Vector<DistributionMapping> dmap(nlevels);
        Box domain(IntVect(0), IntVect(n_cell-1));
        for (int lev = 0; lev < nlevels; ++lev) {
            geom[lev].define(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(0,0,0)});
            // On level 1, a fine patch touching the low x (Dirichlet) and low y
            // (Neumann) boundaries, with coarse/fine faces inside
            grids[lev].define(Box(IntVect(0), IntVect(n_cell-1)));
            grids[lev].maxSize(max_grid_size);
            dmap[lev].define(grids[lev]);
            domain.refine(2);
        }

        // phi_exact also gives the boundary values, in its ghost cells
        Vector<MultiFab> acoef(nlevels), phi_exact(nlevels), rhs(nlevels), sol(nlevels);
        Vector<Array<MultiFab,AMREX_SPACEDIM>> bcoef(nlevels), ccoef(nlevels);
        for (int lev = 0; lev < nlevels; ++lev) {
            auto const& ba = grids[lev];
            auto const& dm = dmap[lev];
            acoef[lev].define(ba, dm, 1, 0);
            phi_exact[lev].define(ba, dm, 1, 1);
            rhs[lev].define(ba, dm, 1, 0);
            sol[lev].define(ba, dm, 1, 1);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                BoxArray const fba = amrex::convert(ba, IntVect::TheDimensionVector(idim));
                bcoef[lev][idim].define(fba, dm, 1, 0);
                ccoef[lev][idim].define(fba, dm, 1, 0);
            }
            auto const problo = geom[lev].ProbLoArray();
            auto const probhi = geom[lev].ProbHiArray();
            auto const dx = geom[lev].CellSizeArray();
            for (MFIter mfi(phi_exact[lev]); mfi.isValid(); ++mfi) {
                auto const& pe = phi_exact[lev].array(mfi);
                auto const& a = acoef[lev].array(mfi);
                ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    IntVect const iv(AMREX_D_DECL(i,j,k));
                    RealVect x;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                        // the ghost cells outside the domain get the face values
                        x[d] = amrex::Clamp(problo[d] + (Real(iv[d])+Real(0.5))*dx[d],
                                            problo[d], probhi[d]);
                    }
                    pe(i,j,k) = exact_phi(x);
                    if (a.contains(i,j,k)) { a(i,j,k) = Real(1.0) + x[0]; }
                });
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    auto const& b = bcoef[lev][idim].array(mfi);
                    auto const& c = ccoef[lev][idim].array(mfi);
                    ParallelFor(amrex::surroundingNodes(mfi.validbox(), idim),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k)
                    {
                        IntVect const iv(AMREX_D_DECL(i,j,k));
                        RealVect x;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            x[d] = problo[d] + (Real(iv[d]) + (d == idim ? Real(0.) : Real(0.5)))*dx[d];
                        }
                        b(i,j,k) = Real(1.0) + Real(0.5)*x[idim]*x[idim];
                        c(i,j,k) = velocity(idim, x);
                    });
                }
            }
        }

        LPInfo info;
        MLABecLapDrift mlabec(geom, grids, dmap, info);
        mlabec.setMaxOrder(3);
        mlabec.setDomainBC(lobc, hibc);
        mlabec.setScalars(ascalar, bscalar, cscalar);
        for (int lev = 0; lev < nlevels; ++lev) {
            mlabec.setLevelBC(lev, &phi_exact[lev]);
            mlabec.setACoeffs(lev, acoef[lev]);
            mlabec.setBCoeffs(lev, GetArrOfConstPtrs(bcoef[lev]));
            mlabec.setDriftCoeffs(lev, GetArrOfConstPtrs(ccoef[lev]));
        }

        MLMG mlmg(mlabec);
        mlmg.setVerbose(verbose);
        if (!use_gauss_seidel) { mlabec.setGaussSeidel(false); }

        // rhs = A phi_exact, with the covered coarse cells the average of the
        // fine ones as in the composite solution
        if (nlevels > 1) {
            amrex::average_down(phi_exact[1], phi_exact[0], 0, 1, 2);
        }
        for (int lev = 0; lev < nlevels; ++lev) {
            MultiFab::Copy(sol[lev], phi_exact[lev], 0, 0, 1, 1);
        }
        mlmg.apply(GetVecOfPtrs(rhs), GetVecOfPtrs(sol));

        if (nlevels == 1) {
            // The operator written out, with the ghost cells of the boundary
            // conditions
            auto const& g = geom[0];
            Box const& dom = g.Domain();
            auto const dxinv = g.InvCellSizeArray();
            MultiFab phi(grids[0], dmap[0], 1, 1);
            MultiFab::Copy(phi, phi_exact[0], 0, 0, 1, 0);
            phi.FillBoundary(g.periodicity());
            MultiFab ref(grids[0], dmap[0], 1, 0);
            for (MFIter mfi(phi); mfi.isValid(); ++mfi) {
                Box const& vbx = mfi.validbox();
                auto const& p = phi.array(mfi);
                auto const& bd = phi_exact[0].const_array(mfi);
                for (OrientationIter oit; oit; ++oit) {
                    Orientation const ori = oit();
                    int const idim = ori.coordDir();
                    if (ori.isLow() && vbx.smallEnd(idim) != dom.smallEnd(idim)) { continue; }
                    if (ori.isHigh() && vbx.bigEnd(idim) != dom.bigEnd(idim)) { continue; }
                    bool const dirichlet = (ori.isLow() ? lobc[idim] : hibc[idim])
                        == LinOpBCType::Dirichlet;
                    IntVect const in = ori.isLow() ? IntVect::TheDimensionVector(idim)
                                                   : -IntVect::TheDimensionVector(idim);
                    ParallelFor(amrex::adjCell(vbx, ori), [=] AMREX_GPU_DEVICE (int i, int j, [[maybe_unused]] int k)
                    {
                        IntVect const gv(AMREX_D_DECL(i,j,k));
                        p(gv) = dirichlet ? dirichlet_ghost(bd(gv), p(gv+in), p(gv+2*in)) : p(gv+in);
                    });
                }
                auto const& r = ref.array(mfi);
                auto const& a = acoef[0].const_array(mfi);
                GpuArray<Array4<Real const>,AMREX_SPACEDIM> b, c;
                for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                    b[idim] = bcoef[0][idim].const_array(mfi);
                    c[idim] = ccoef[0][idim].const_array(mfi);
                }
                ParallelFor(vbx, [=] AMREX_GPU_DEVICE (int i, int j, [[maybe_unused]] int k)
                {
                    IntVect const iv(AMREX_D_DECL(i,j,k));
                    Real v = ascalar * a(iv) * p(iv);
                    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
                        IntVect const e = IntVect::TheDimensionVector(idim);
                        Real const fm = b[idim](iv)   * (p(iv)   - p(iv-e)) * dxinv[idim];
                        Real const fp = b[idim](iv+e) * (p(iv+e) - p(iv))   * dxinv[idim];
                        Real const um = c[idim](iv);
                        Real const up = c[idim](iv+e);
                        Real const dm = um * (um > Real(0.) ? p(iv-e) : p(iv));
                        Real const dp = up * (up > Real(0.) ? p(iv) : p(iv+e));
                        v -= bscalar * (fp - fm) * dxinv[idim] + cscalar * (dp - dm) * dxinv[idim];
                    }
                    r(iv) = v;
                });
            }
            MultiFab::Subtract(ref, rhs[0], 0, 0, 1, 0);
            Real const aerr = ref.norminf(0, 0) / rhs[0].norminf(0, 0);
            amrex::Print() << "MLABecLapDrift: relative max-norm difference of apply "
                           << aerr << "\n";
            if (!(aerr <= max_error)) {
                amrex::Abort("ABecLapDrift: apply differs from the operator written out");
            }
        }

        for (int lev = 0; lev < nlevels; ++lev) { sol[lev].setVal(Real(0.0)); }
        mlmg.solve(GetVecOfPtrs(sol), GetVecOfConstPtrs(rhs), tol, Real(0.0));

        Real err = 0.0;
        for (int lev = 0; lev < nlevels; ++lev) {
            MultiFab::Subtract(sol[lev], phi_exact[lev], 0, 0, 1, 0);
            err = std::max(err, sol[lev].norminf(0, 0) / phi_exact[lev].norminf(0, 0));
        }
        amrex::Print() << "MLABecLapDrift: relative max-norm difference of solve "
                       << err << " after " << mlmg.getNumIters() << " iterations\n";
        if (!(err <= max_error)) {
            amrex::Abort("ABecLapDrift: the solve does not give phi back");
        }
    }
    amrex::Finalize();
}
