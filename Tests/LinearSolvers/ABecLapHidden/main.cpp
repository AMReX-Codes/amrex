// MLABecLaplacian with a hidden direction (LPInfo::setHiddenDirection) on a
// domain one cell wide in that direction.
//
// The operator must be the two-dimensional one: the b coefficients of the
// hidden direction, set here to nonzero values, must not enter it.  The
// reference is the same operator without a hidden direction and with zero b
// coefficients in that direction, whose matrix is the two-dimensional one.
// The test checks that the two apply the same operator and give the same
// fluxes, with a zero flux in the hidden direction, and that MLMG with the
// hidden direction solves the reference problem.

#include <AMReX.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Random.H>

using namespace amrex;

namespace {

LinOpBCType parse_bc (std::string const& s)
{
    if (s == "periodic")  { return LinOpBCType::Periodic; }
    if (s == "neumann")   { return LinOpBCType::Neumann; }
    if (s == "dirichlet") { return LinOpBCType::Dirichlet; }
    amrex::Abort("Unknown bc " + s);
    return LinOpBCType::bogus;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        ParmParse pp;

        int hidden_direction = 1;
        pp.query("hidden_direction", hidden_direction);
        AMREX_ALWAYS_ASSERT(hidden_direction >= 0 && hidden_direction < AMREX_SPACEDIM);

        IntVect n_cell(32,32,32);
        pp.queryarr("n_cell", n_cell);
        n_cell[hidden_direction] = 1;
        IntVect max_grid_size(16);
        pp.queryarr("max_grid_size", max_grid_size);

        std::vector<Real> prob_hi{Real(1.0), Real(1.0), Real(1.0)};
        pp.queryarr("prob_hi", prob_hi);

        std::vector<std::string> bc_lo{"dirichlet","dirichlet","neumann"};
        std::vector<std::string> bc_hi{"dirichlet","dirichlet","dirichlet"};
        pp.queryarr("bc_lo", bc_lo);
        pp.queryarr("bc_hi", bc_hi);

        // a = 0 with no Dirichlet face makes the operator singular
        Real ascalar = Real(1.0);
        pp.query("ascalar", ascalar);
        bool semicoarsening = false;
        int semicoarsening_direction = 2;
        int max_semicoarsening_level = 3;
        pp.query("semicoarsening", semicoarsening);
        pp.query("semicoarsening_direction", semicoarsening_direction);
        pp.query("max_semicoarsening_level", max_semicoarsening_level);
        bool use_gauss_seidel = true;
        pp.query("use_gauss_seidel", use_gauss_seidel);

        int verbose = 1;
        pp.query("verbose", verbose);
        int max_iter = 200;
        pp.query("max_iter", max_iter);
        Real reltol = std::is_same_v<Real,double> ? Real(1.e-11) : Real(1.e-5);
        pp.query("reltol", reltol);
        // Apply and fluxes differ by round-off only; the solve by its tolerance
        Real const apply_tol = Real(100.)*std::numeric_limits<Real>::epsilon();
        Real check_tol = std::is_same_v<Real,double> ? Real(1.e-8) : Real(5.e-3);
        pp.query("check_tol", check_tol);

        amrex::InitRandom(42);

        Array<LinOpBCType,AMREX_SPACEDIM> lobc, hibc;
        Array<int,AMREX_SPACEDIM> is_periodic{};
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            lobc[idim] = parse_bc(bc_lo[idim]);
            hibc[idim] = parse_bc(bc_hi[idim]);
            AMREX_ALWAYS_ASSERT((lobc[idim] == LinOpBCType::Periodic) ==
                                (hibc[idim] == LinOpBCType::Periodic));
            is_periodic[idim] = (lobc[idim] == LinOpBCType::Periodic);
        }

        Box domain(IntVect(0), n_cell-1);
        RealBox rb({Real(0.),Real(0.),Real(0.)}, {prob_hi[0],prob_hi[1],prob_hi[2]});
        Geometry geom(domain, rb, CoordSys::cartesian, is_periodic);
        BoxArray ba(domain);
        ba.maxSize(max_grid_size);
        DistributionMapping dm(ba);

        auto const prob_lo = geom.ProbLoArray();
        auto const dx = geom.CellSizeArray();
        int const hd = hidden_direction;

        // Smooth, positive coefficients, periodic so that a periodic direction
        // gives a conservative operator; b of the hidden direction is nonzero
        // and must be ignored with the hidden direction.
        Real const twopi = Real(2.)*Math::pi<Real>();
        Real const Lx = prob_hi[0], Ly = prob_hi[1], Lz = prob_hi[2];
        MultiFab acoef(ba, dm, 1, 0);
        Array<MultiFab,AMREX_SPACEDIM> bcoef, bcoef_ref;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            bcoef    [idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 0);
            bcoef_ref[idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)), dm, 1, 0);
        }
        for (MFIter mfi(acoef); mfi.isValid(); ++mfi) {
            auto const& a = acoef.array(mfi);
            ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real x = prob_lo[0] + (Real(i)+Real(0.5))*dx[0];
                Real z = prob_lo[2] + (Real(k)+Real(0.5))*dx[2];
                a(i,j,k) = Real(1.0) + Real(0.5)*std::sin(twopi*(x/Lx+Real(0.5)*z/Lz));
            });
        }
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            for (MFIter mfi(bcoef[idim]); mfi.isValid(); ++mfi) {
                auto const& b = bcoef[idim].array(mfi);
                auto const& bref = bcoef_ref[idim].array(mfi);
                ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    Real x = prob_lo[0] + (Real(i)+Real(0.5)*Real(idim!=0))*dx[0];
                    Real y = prob_lo[1] + (Real(j)+Real(0.5)*Real(idim!=1))*dx[1];
                    Real z = prob_lo[2] + (Real(k)+Real(0.5)*Real(idim!=2))*dx[2];
                    b(i,j,k) = Real(1.0) + Real(0.5)*std::cos(twopi*(x/Lx-y/Ly+Real(2.)*z/Lz)
                                                              +Real(idim));
                    bref(i,j,k) = (idim == hd) ? Real(0.0) : b(i,j,k);
                });
            }
        }

        // The solution, with Dirichlet values (and ghost cells) from it
        MultiFab phi_exact(ba, dm, 1, 1);
        FillRandom(phi_exact, 0, 1);
        for (MFIter mfi(phi_exact); mfi.isValid(); ++mfi) {
            auto const& p = phi_exact.array(mfi);
            ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                Real x = prob_lo[0] + (Real(i)+Real(0.5))*dx[0];
                Real y = prob_lo[1] + (Real(j)+Real(0.5))*dx[1];
                Real z = prob_lo[2] + (Real(k)+Real(0.5))*dx[2];
                p(i,j,k) = std::sin(Real(2.)*x+Real(0.3)) * std::cos(Real(1.5)*z)
                    + Real(0.3)*std::cos(Real(3.)*y) + Real(0.01)*p(i,j,k);
            });
        }

        auto make_linop = [&] (bool hidden, Array<MultiFab,AMREX_SPACEDIM> const& b)
        {
            LPInfo info;
            if (hidden) { info.setHiddenDirection(hd); }
            if (hidden && semicoarsening) {
                info.setSemicoarsening(true)
                    .setSemicoarseningDirection(semicoarsening_direction)
                    .setMaxSemicoarseningLevel(max_semicoarsening_level);
            }
            auto op = std::make_unique<MLABecLaplacian>(Vector<Geometry>{geom}, Vector<BoxArray>{ba},
                                                        Vector<DistributionMapping>{dm}, info);
            op->setMaxOrder(2);
            op->setDomainBC(lobc, hibc);
            op->setLevelBC(0, &phi_exact);
            op->setScalars(ascalar, Real(1.0));
            op->setACoeffs(0, acoef);
            op->setBCoeffs(0, GetArrOfConstPtrs(b));
            return op;
        };

        auto linop_ref = make_linop(false, bcoef_ref);
        auto linop     = make_linop(true , bcoef);
        MLMG mlmg_ref(*linop_ref);
        MLMG mlmg(*linop);
        mlmg.setMaxIter(max_iter);
        mlmg.setVerbose(verbose);
        if (!use_gauss_seidel) { linop->setGaussSeidel(false); }

        auto rel_diff = [] (MultiFab const& a, MultiFab const& b)
        {
            MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
            MultiFab::LinComb(d, Real(1.0), a, 0, Real(-1.0), b, 0, 0, 1, 0);
            Real const scale = std::max(b.norminf(0), std::numeric_limits<Real>::min());
            return d.norminf(0) / scale;
        };
        auto failed = [] (std::string const& what, Real err, Real tol)
        {
            amrex::Print() << what << ": relative max difference " << err << "\n";
            if (!(err <= tol)) { amrex::Abort(what + " failed"); }
        };

        // The operator
        MultiFab rhs(ba, dm, 1, 0), lphi(ba, dm, 1, 0);
        mlmg_ref.apply({&rhs}, {&phi_exact});
        mlmg.apply({&lphi}, {&phi_exact});
        failed("apply", rel_diff(lphi, rhs), apply_tol);

        // The fluxes, zero in the hidden direction
        Array<MultiFab,AMREX_SPACEDIM> flux, flux_ref;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            flux    [idim].define(bcoef[idim].boxArray(), dm, 1, 0);
            flux_ref[idim].define(bcoef[idim].boxArray(), dm, 1, 0);
            flux[idim].setVal(Real(1.0));
        }
        MultiFab phi_tmp(ba, dm, 1, 1);
        MultiFab::Copy(phi_tmp, phi_exact, 0, 0, 1, 1);
        mlmg.getFluxes({GetArrOfPtrs(flux)}, {&phi_tmp}, MLMG::Location::FaceCenter);
        MultiFab::Copy(phi_tmp, phi_exact, 0, 0, 1, 1);
        mlmg_ref.getFluxes({GetArrOfPtrs(flux_ref)}, {&phi_tmp}, MLMG::Location::FaceCenter);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (idim == hd) {
                failed("flux in the hidden direction", flux[idim].norminf(0), Real(0.0));
            } else {
                failed("flux in direction " + std::to_string(idim),
                       rel_diff(flux[idim], flux_ref[idim]), apply_tol);
            }
        }

        // The fluxes on the box faces only (as in reflux), zero on both faces
        // of the hidden direction
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            flux[idim].setVal(Real(1.0));
        }
        for (MFIter mfi(phi_exact); mfi.isValid(); ++mfi) {
            MLABecLaplacian::FFlux(mfi.validbox(), geom.InvCellSize(), Real(1.0),
                                   {AMREX_D_DECL(&bcoef[0][mfi], &bcoef[1][mfi], &bcoef[2][mfi])},
                                   {AMREX_D_DECL(&flux[0][mfi], &flux[1][mfi], &flux[2][mfi])},
                                   phi_exact[mfi], 1, 1, hd);
        }
        failed("face-only flux in the hidden direction", flux[hd].norminf(0), Real(0.0));

        // The solve
        MultiFab phi(ba, dm, 1, 1);
        phi.setVal(Real(0.0));
        mlmg.solve({&phi}, {&rhs}, reltol, Real(0.0));
        if (linop->isSingular(0)) {
            Real const shift = (phi.sum(0) - phi_exact.sum(0)) / Real(domain.numPts());
            phi.plus(-shift, 0, 1);
        }
        failed("solve", rel_diff(phi, phi_exact), check_tol);
    }
    amrex::Finalize();
}
