// Regression test for MLTensorOp::setMappingFactors.
//
// The same physical problem is set up twice on index-identical grids:
//
//   ref    : a uniform grid on the physical domain [0,1]^d;
//   mapped : the uniform computational grid of a diagonal mapping
//            x_d = fac_d * xi_d with anisotropic constant factors, i.e. the
//            xi domain [0,1/fac_d], the face viscosities scaled by
//            J/fac_d^2 (J = prod fac_d) and the factors handed to the
//            operator with setMappingFactors.
//
// Both see the same velocity, viscosity and boundary data at the same
// physical positions, so the mapped operator must return J times the
// reference operator and the mapped flux on face d must be J/fac_d times
// the reference flux, to round-off (1e-12 in double, 1e-4 in single
// precision).  Without the factors the transpose and bulk terms are off by
// fac_d/fac_j, which the test also checks is visible,
// so that a regression in the ratio kernels cannot pass unnoticed.  Dirichlet
// boundaries on all sides exercise the boundary variants of the kernels;
// variable eta and kappa exercise the coefficient handling; with four boxes
// per direction both the interior and the boundary kernels run.

#include <AMReX.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLTensorOp.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <numbers>

using namespace amrex;

namespace {

// Smooth test fields on the unit cube (z is ignored in 2D).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real velfun (int n, Real x, Real y, Real z) noexcept
{
    constexpr Real pi = Real(std::numbers::pi_v<double>);
    Real const sx = std::sin(pi*x), cx = std::cos(pi*x);
    Real const sy = std::sin(Real(2.)*pi*y), cy = std::cos(Real(2.)*pi*y);
    Real const sz = std::sin(pi*z), cz = std::cos(pi*z);
    if (n == 0) { return Real(1.) + sx*cy + Real(0.3)*x*x*cz; }
    if (n == 1) { return Real(0.5) + cx*sy*sz + Real(0.2)*y*y*y; }
    return Real(-0.3) + sx*sy*cz + Real(0.1)*(x+y)*z;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real etafun (Real x, Real y, Real z) noexcept
{
    constexpr Real pi = Real(std::numbers::pi_v<double>);
    return Real(1.) + Real(0.5)*std::sin(pi*x)*std::cos(Real(2.)*pi*y)*std::cos(pi*z);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real kapfun (Real x, Real y, Real z) noexcept
{
    constexpr Real pi = Real(std::numbers::pi_v<double>);
    return Real(0.3) + Real(0.2)*std::cos(pi*x)*std::sin(pi*y)*std::cos(Real(2.)*pi*z);
}

struct Result
{
    MultiFab lap;                        // L(vel), AMREX_SPACEDIM components
    Array<MultiFab,AMREX_SPACEDIM> flux; // face fluxes, AMREX_SPACEDIM components
};

// Apply the tensor operator on the grid mapped by fac (fac = 1 gives the
// reference).  with_factors selects whether setMappingFactors is called.
Result run (Array<Real,AMREX_SPACEDIM> const& fac, bool with_factors,
            int n_cell, int max_grid_size)
{
    Real J = Real(1.);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { J *= fac[d]; }

    // Physical spacing h; computational spacing h/fac_d
    Real const h = Real(1.)/Real(n_cell);
    RealBox rb;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { rb.setLo(d, Real(0.)); rb.setHi(d, Real(1.)/fac[d]); }
    Array<int,AMREX_SPACEDIM> isper{AMREX_D_DECL(0,0,0)};
    Box domain(IntVect(0), IntVect(n_cell-1));
    Geometry geom(domain, rb, CoordSys::cartesian, isper);
    BoxArray grids(domain);
    grids.maxSize(max_grid_size);
    DistributionMapping dmap(grids);

    // Velocity with one ghost layer.  In ghost cells the Dirichlet value at
    // the domain face is stored, as the solvers expect: clamp the physical
    // coordinate to [0,1] in the direction that is outside.
    MultiFab vel(grids, dmap, AMREX_SPACEDIM, 1);
    {
        auto const& ma = vel.arrays();
        ParallelFor(vel, IntVect(1), AMREX_SPACEDIM,
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k, int n) noexcept
        {
            Real x = amrex::min(amrex::max((Real(i)+Real(0.5))*h, Real(0.)), Real(1.));
            Real y = amrex::min(amrex::max((Real(j)+Real(0.5))*h, Real(0.)), Real(1.));
#if (AMREX_SPACEDIM == 3)
            Real z = amrex::min(amrex::max((Real(k)+Real(0.5))*h, Real(0.)), Real(1.));
#else
            Real z = Real(0.25);
#endif
            ma[bno](i,j,k,n) = velfun(n, x, y, z);
        });
    }

    // Face coefficients: physical eta, kappa at the face centre, scaled by
    // J/fac_d^2 on faces of direction d; mapping factors on the same faces.
    Array<MultiFab,AMREX_SPACEDIM> eta, kappa, mapfac;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        BoxArray const& fba = amrex::convert(grids, IntVect::TheDimensionVector(d));
        eta[d].define(fba, dmap, 1, 0);
        kappa[d].define(fba, dmap, 1, 0);
        mapfac[d].define(fba, dmap, AMREX_SPACEDIM, 0);
        Real const scale = J/(fac[d]*fac[d]);
        auto const& ema = eta[d].arrays();
        auto const& kma = kappa[d].arrays();
        auto const& mma = mapfac[d].arrays();
        Real const ox = (d == 0) ? Real(0.) : Real(0.5);
        Real const oy = (d == 1) ? Real(0.) : Real(0.5);
        Real const oz = (d == 2) ? Real(0.) : Real(0.5);
        ParallelFor(eta[d],
        [=] AMREX_GPU_DEVICE (int bno, int i, int j, int k) noexcept
        {
            Real x = (Real(i)+ox)*h;
            Real y = (Real(j)+oy)*h;
#if (AMREX_SPACEDIM == 3)
            Real z = (Real(k)+oz)*h;
#else
            Real z = Real(0.25);
            amrex::ignore_unused(oz);
#endif
            ema[bno](i,j,k) = scale*etafun(x,y,z);
            kma[bno](i,j,k) = scale*kapfun(x,y,z);
            for (int c = 0; c < AMREX_SPACEDIM; ++c) { mma[bno](i,j,k,c) = fac[c]; }
        });
    }
    Gpu::streamSynchronize();

    LPInfo info;
    info.setMaxCoarseningLevel(0);
    auto owned = std::make_unique<MLTensorOp>(Vector<Geometry>{geom}, Vector<BoxArray>{grids},
                                              Vector<DistributionMapping>{dmap}, info);
    MLTensorOp* op = owned.get();
    std::array<LinOpBCType,AMREX_SPACEDIM> bc;
    bc.fill(LinOpBCType::Dirichlet);
    op->setDomainBC(bc, bc);
    op->setLevelBC(0, &vel);
    MultiFab acoef(grids, dmap, 1, 0);
    acoef.setVal(Real(0.));
    op->setACoeffs(0, acoef);
    op->setShearViscosity(0, GetArrOfConstPtrs(eta));
    op->setBulkViscosity(0, GetArrOfConstPtrs(kappa));
    if (with_factors) {
        op->setMappingFactors(0, GetArrOfConstPtrs(mapfac));
    }

    Result r;
    r.lap.define(grids, dmap, AMREX_SPACEDIM, 0);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        r.flux[d].define(amrex::convert(grids, IntVect::TheDimensionVector(d)), dmap, AMREX_SPACEDIM, 0);
    }
    MLMG mlmg(*op);
    mlmg.setVerbose(0);
    mlmg.apply({&r.lap}, {&vel});
    mlmg.getFluxes({GetArrOfPtrs(r.flux)}, {&vel}, MLMG::Location::FaceCenter);

    // Undo the mapping scaling so that the results are directly comparable
    // with the reference: L -> L/J, F_d -> F_d fac_d/J.
    r.lap.mult(Real(1.)/J);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { r.flux[d].mult(fac[d]/J); }
    return r;
}

// max |a-b| / max |b| over all components
Real reldiff (MultiFab const& a, MultiFab const& b)
{
    MultiFab d(a.boxArray(), a.DistributionMap(), a.nComp(), 0);
    MultiFab::Copy(d, a, 0, 0, a.nComp(), 0);
    MultiFab::Subtract(d, b, 0, 0, a.nComp(), 0);
    Real num = Real(0.), den = Real(0.);
    for (int n = 0; n < a.nComp(); ++n) {
        num = amrex::max(num, d.norm0(n));
        den = amrex::max(den, b.norm0(n));
    }
    return num/den;
}

} // namespace

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int n_cell = 32;
        int max_grid_size = 8;   // 4 boxes per direction: interior boxes
                                 // (interior kernels) and boundary boxes
                                 // (boundary kernels) are both present
        // Round-off on the second differences of O(1) fields at h = 1/32:
        // ~1e-14 in double, ~1e-6 in single precision.
        Real tol = (sizeof(Real) == 4) ? Real(1.e-4) : Real(1.e-12);
        Array<Real,AMREX_SPACEDIM> fac{AMREX_D_DECL(Real(2.), Real(0.5), Real(1.25))};
        {
            ParmParse pp;
            pp.query("n_cell", n_cell);
            pp.query("max_grid_size", max_grid_size);
            pp.query("tol", tol);
            Vector<Real> f;
            if (pp.queryarr("fac", f)) {
                AMREX_ALWAYS_ASSERT(f.size() == AMREX_SPACEDIM);
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { fac[d] = f[d]; }
            }
        }
        Array<Real,AMREX_SPACEDIM> const one{AMREX_D_DECL(Real(1.), Real(1.), Real(1.))};

        auto ref    = run(one, false, n_cell, max_grid_size);
        auto mapped = run(fac, true,  n_cell, max_grid_size);
        auto naive  = run(fac, false, n_cell, max_grid_size);

        Real err_lap = reldiff(mapped.lap, ref.lap);
        Real err_naive = reldiff(naive.lap, ref.lap);
        Real err_flux = Real(0.);
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            err_flux = amrex::max(err_flux, reldiff(mapped.flux[d], ref.flux[d]));
        }

        amrex::Print() << "fac = " << AMREX_D_TERM(fac[0], << " " << fac[1], << " " << fac[2]) << "\n"
                       << "  mapped vs ref, operator : " << err_lap << "\n"
                       << "  mapped vs ref, fluxes   : " << err_flux << "\n"
                       << "  without factors         : " << err_naive << "\n";

        bool pass = (err_lap < tol) && (err_flux < tol);
        // The factors must matter: without them the cross terms are off by
        // fac_d/fac_j, an O(1) relative error on the transpose part (0.67 for
        // these fields and factors).
        bool discriminates = err_naive > Real(1.e-2);

        if (!pass || !discriminates) {
            amrex::Abort("TensorMapped test FAILED");
        }
        amrex::Print() << "TensorMapped test PASSED\n";
    }
    amrex::Finalize();
}
