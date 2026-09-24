#include <AMReX_AMG.H>
#include <AMReX_MLABecLaplacian.H>
#include <AMReX_MLMG.H>
#include <AMReX_MultiFab.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <cmath>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace amrex;

namespace {

// Optional entries are only passed to the solver when given, so that the
// solver's own defaults apply otherwise.
struct Params {
    int n_cell = 32;
    int max_iter = 200;
    int repeat = 1; // run each case this many times (steady-state timing)
    int fixed_iter = 0;
    Real reltol = (sizeof(Real) == 4) ? Real(1.e-5) : Real(1.e-10);
    Real alpha = Real(1); // 1.e-6 makes the constant mode nearly null
    int variations = 1; // 1: also run variations of the options; 0: only them
    std::string bottom = "jacobi"; // jacobi, bicgstab, gmres
    std::string interp = "ext+i"; // direct, ext, ext+i
    std::string smoother = "chebyshev"; // jacobi, l1jacobi, chebyshev, l1gs (CPU)
    std::string krylov = "none"; // none, bicgstab, gmres, pcg
    // periodic, or dirichlet: homogeneous Dirichlet on the domain faces
    std::string bc = "periodic";
    // Input `problem` is a list (all four by default); each run has one.
    // constant: a*phi - lap(phi);
    // jump: coefficient `jump` in the central cube;
    // checker: blocks of `block` cells alternating 1 and `jump`;
    // aniso: coefficient `eps` in x, 1 in the other directions.
    std::string problem = "constant";
    Real jump = Real(1.e3);
    int block = 4;
    Real eps = Real(1.e-3);
    int mlmg = 1;           // also solve with geometric MLMG for comparison
    int max_grid_size = 64; // for MLMG
    std::optional<int> verbose, nu1, nu2, nu_bottom, p_max_elmts, max_levels,
                       aggressive_levels, cheby_degree, aggressive_direct;
    std::optional<Long> max_coarse_size;
    std::optional<Real> bottom_tol, trunc_factor, relax_weight, cheby_ratio, theta;
};

struct Result {
    Real error;     // max-norm error of the solution
    Real err_bound; // upper bound of the error implied by the residual
    Real rel_res;   // final residual 2-norm relative to the rhs
    int niters;
    int nlevels;
    double time;    // solve time in seconds
    bool blew_up;
};

std::string sci (Real v)
{
    std::ostringstream os;
    os << std::scientific << std::setprecision(2) << v;
    return os.str();
}

// Cell-centered diffusion coefficient of the test problems.
struct Coef
{
    int type; // 0 constant, 1 jump, 2 checker, 3 aniso
    Real jump;
    int block;
    Real eps;
    Box domain;

    AMREX_GPU_DEVICE Real operator() (IntVect const& c, int dir) const noexcept
    {
        if (type == 0) { return Real(1); }
        if (type == 3) { return (dir == 0) ? eps : Real(1); }
        if (type == 1) {
            bool inside = true;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                int const n = domain.length(d);
                inside = inside && (c[d] >= n/4) && (c[d] < 3*n/4);
            }
            return inside ? jump : Real(1);
        }
        int par = 0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { par += c[d] / block; }
        return (par % 2) ? jump : Real(1);
    }
};

// Geometric multigrid (MLMG) on the same problem, for comparison. MLMG may
// fail on some of these problems; that is reported, not fatal.
void run_mlmg (Params const& p)
{
    int const n_cell = p.n_cell;
    int const ptype = (p.problem == "constant") ? 0 : (p.problem == "jump") ? 1
                    : (p.problem == "checker") ? 2 : 3;
    Box domain(IntVect(0), IntVect(n_cell-1));
    Coef const coef{.type = ptype, .jump = p.jump, .block = p.block, .eps = p.eps, .domain = domain};
    Real const L = Real(2)*Math::pi<Real>();
    RealBox rb(AMREX_D_DECL(Real(0),Real(0),Real(0)), AMREX_D_DECL(L,L,L));
    bool const dirichlet = (p.bc == "dirichlet");
    int const per = dirichlet ? 0 : 1;
    Array<int,AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(per,per,per)};
    Geometry geom(domain, rb, CoordSys::cartesian, is_periodic);
    BoxArray ba(domain);
    ba.maxSize(p.max_grid_size);
    DistributionMapping dm(ba);
    Real const a = p.alpha;
    Real const dx = geom.CellSize(0);

    MultiFab phi(ba, dm, 1, 1), rhs(ba, dm, 1, 0), exact(ba, dm, 1, 1), res(ba, dm, 1, 0);
    auto const& exa = exact.arrays();
    ParallelFor(exact, IntVect(1), [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
    {
        amrex::ignore_unused(j, k);
        Real v = std::sin((Real(i)+Real(0.5))*dx);
#if (AMREX_SPACEDIM >= 2)
        v *= std::sin((Real(j)+Real(0.5))*dx);
#endif
#if (AMREX_SPACEDIM == 3)
        v *= std::sin((Real(k)+Real(0.5))*dx);
#endif
        exa[b](i,j,k) = Math::powi<5>(v);
    });

    // Face coefficients: harmonic mean of the two cells, as in the matrix.
    Array<MultiFab,AMREX_SPACEDIM> bcoef;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        bcoef[idim].define(amrex::convert(ba, IntVect::TheDimensionVector(idim)), dm, 1, 0);
        auto const& bfa = bcoef[idim].arrays();
        ParallelFor(bcoef[idim], [=] AMREX_GPU_DEVICE (int b, int i, int j, int k)
        {
            IntVect hi(AMREX_D_DECL(i,j,k));
            IntVect lo = hi;
            lo[idim] -= 1;
            int const n = domain.length(idim);
            if (dirichlet) { // boundary faces take the inside cell's coefficient
                if (lo[idim] < 0) { lo = hi; }
                if (hi[idim] >= n) { hi = lo; }
            } else { // periodic wrap
                lo[idim] = ((lo[idim] % n) + n) % n;
                hi[idim] = hi[idim] % n;
            }
            Real const b0 = coef(lo, idim);
            Real const b1 = coef(hi, idim);
            bfa[b](i,j,k) = Real(2)*b0*b1/(b0+b1);
        });
    }
    Gpu::streamSynchronize();

    MLABecLaplacian mlabec({geom}, {ba}, {dm});
    auto const lbc = dirichlet ? LinOpBCType::Dirichlet : LinOpBCType::Periodic;
    mlabec.setDomainBC({AMREX_D_DECL(lbc,lbc,lbc)}, {AMREX_D_DECL(lbc,lbc,lbc)});
    mlabec.setLevelBC(0, nullptr);
    mlabec.setScalars(a, Real(1));
    mlabec.setACoeffs(0, Real(1));
    mlabec.setBCoeffs(0, GetArrOfConstPtrs(bcoef));
    MLMG mlmg(mlabec);
    mlmg.setVerbose(p.verbose.value_or(0));
    mlmg.setMaxIter(p.max_iter);
    mlmg.setThrowException(true);
    mlmg.apply({&rhs}, {&exact}); // rhs = A*phi with the same operator
    phi.setVal(0);

    std::string failure;
    Gpu::streamSynchronize();
    auto const t0 = amrex::second();
    try {
        mlmg.solve({&phi}, {&rhs}, p.reltol, Real(0));
    } catch (std::exception const& e) {
        failure = e.what();
    }
    Gpu::streamSynchronize();
    auto const t1 = amrex::second();

    mlmg.compResidual({&res}, {&phi}, {&rhs});
    Real const rel_res = res.norminf(0, 0) / rhs.norminf(0, 0);
    MultiFab::Subtract(phi, exact, 0, 0, 1, 0);
    if (a == Real(0) && !dirichlet) { // solution defined up to a constant
        phi.plus(-phi.sum(0) / Real(domain.numPts()), 0, 1, 0);
    }
    amrex::Print() << "  MLMG for comparison: ";
    if (failure.empty()) {
        amrex::Print() << mlmg.getNumIters() << " iterations, "
                       << std::fixed << std::setprecision(4) << (t1-t0) << std::defaultfloat
                       << " s, rel_res " << sci(rel_res) << " (max norm), error "
                       << sci(phi.norminf(0, 0)) << "\n";
    } else {
        if (failure.back() == '.') { failure.pop_back(); }
        amrex::Print() << "did not converge after " << mlmg.getNumIters() << " iterations ("
                       << failure << "); informational only, not an AMG result\n";
    }
}

// a*phi - div(beta grad phi) with phi = prod sin^5, periodic or with
// Dirichlet on the domain faces. The rhs is the analytic discrete Laplacian
// for the periodic constant problem and A*phi otherwise.
Result run (Params const& p)
{
    int const n_cell = p.n_cell;
    std::string const& bottom = p.bottom;
    std::string const& interp = p.interp;
    Box domain(IntVect(0),IntVect(n_cell-1));
    Long n = domain.numPts();
    AlgVector<Real> xvec(n);
    AlgVector<Real> bvec(xvec.partition());
    AlgVector<Real> exact(xvec.partition());

    Real a = p.alpha;
    Real dx = Real(2)*amrex::Math::pi<Real>()/Real(domain.length(0));

    BoxIndexer box_indexer(domain);
    int const ptype = (p.problem == "constant") ? 0 : (p.problem == "jump") ? 1
                    : (p.problem == "checker") ? 2 : (p.problem == "aniso") ? 3 : -1;
    if (ptype < 0) { amrex::Abort("Unknown problem: " + p.problem); }
    Coef const coef{.type = ptype, .jump = p.jump, .block = p.block, .eps = p.eps, .domain = domain};

    {
        auto* rhs = bvec.data();
        auto* phi = exact.data();
        auto nrows = bvec.numLocalRows();
        auto ib = bvec.globalBegin();
        ParallelFor(nrows, [=] AMREX_GPU_DEVICE (Long lrow)
        {
            auto row = lrow + ib; // global row index
            IntVect cell = box_indexer.intVect(row);
#if (AMREX_SPACEDIM == 1)
            auto x = (Real(cell[0])+Real(0.5))*dx;
            auto phi0 = Math::powi<5>(std::sin(x));
            auto phixm = Math::powi<5>(std::sin(x-dx));
            auto phixp = Math::powi<5>(std::sin(x+dx));
            rhs[lrow] = a*phi0 + (Real(2)*phi0-phixm-phixp) / (dx*dx);
#elif (AMREX_SPACEDIM == 2)
            auto x = (Real(cell[0])+Real(0.5))*dx;
            auto y = (Real(cell[1])+Real(0.5))*dx;
            auto phi0 = Math::powi<5>(std::sin(x)*std::sin(y));
            auto phixm = Math::powi<5>(std::sin(x-dx)*std::sin(y));
            auto phixp = Math::powi<5>(std::sin(x+dx)*std::sin(y));
            auto phiym = Math::powi<5>(std::sin(x)*std::sin(y-dx));
            auto phiyp = Math::powi<5>(std::sin(x)*std::sin(y+dx));
            rhs[lrow] = a*phi0 + (Real(4)*phi0-phixm-phixp-phiym-phiyp) / (dx*dx);
#else
            auto x = (Real(cell[0])+Real(0.5))*dx;
            auto y = (Real(cell[1])+Real(0.5))*dx;
            auto z = (Real(cell[2])+Real(0.5))*dx;
            auto phi0 = Math::powi<5>(std::sin(x)*std::sin(y)*std::sin(z));
            auto phixm = Math::powi<5>(std::sin(x-dx)*std::sin(y)*std::sin(z));
            auto phixp = Math::powi<5>(std::sin(x+dx)*std::sin(y)*std::sin(z));
            auto phiym = Math::powi<5>(std::sin(x)*std::sin(y-dx)*std::sin(z));
            auto phiyp = Math::powi<5>(std::sin(x)*std::sin(y+dx)*std::sin(z));
            auto phizm = Math::powi<5>(std::sin(x)*std::sin(y)*std::sin(z-dx));
            auto phizp = Math::powi<5>(std::sin(x)*std::sin(y)*std::sin(z+dx));
            rhs[lrow] = a*phi0 + (Real(6)*phi0-phixm-phixp-phiym-phiyp-phizm-phizp) / (dx*dx);
#endif
            phi[lrow] = phi0;
        });
    }

    xvec.setVal(0);

    // Cross stencil with harmonic face coefficients. A Dirichlet boundary
    // face adds 2 b / dx^2 to the diagonal (as in MLMG) and has no neighbor.
    bool const dirichlet = (p.bc == "dirichlet");
    auto set_stencil = [=] AMREX_GPU_DEVICE (Long row, Long* col, Real* val)
    {
        IntVect cell = box_indexer.intVect(row);
        int i = 0;
        Real diag = a;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Real const b0 = coef(cell, idim);
            for (int side = -1; side <= 1; side += 2) {
                IntVect cell2 = cell;
                cell2[idim] += side;
                bool const outside = (cell2[idim] < domain.smallEnd(idim) ||
                                      cell2[idim] > domain.bigEnd(idim));
                if (outside && dirichlet) {
                    col[i] = -1; // dropped
                    val[i] = Real(0);
                    diag += Real(2)*b0/(dx*dx);
                } else {
                    if (outside) { // periodic wrap
                        cell2[idim] = (side < 0) ? domain.bigEnd(idim) : domain.smallEnd(idim);
                    }
                    Real const b1 = coef(cell2, idim);
                    Real const bf = Real(2)*b0*b1/(b0+b1);
                    col[i] = domain.index(cell2);
                    val[i] = -bf/(dx*dx);
                    diag += bf/(dx*dx);
                }
                ++i;
            }
        }
        col[i] = row;
        val[i] = diag;
    };

    int const nnz_row = 2*AMREX_SPACEDIM+1;
    Long const nlocal = xvec.numLocalRows();
    SpMatrix<Real> mat;
    {
        Gpu::DeviceVector<Real> vals(nlocal*nnz_row);
        Gpu::DeviceVector<Long> cols(nlocal*nnz_row);
        Gpu::DeviceVector<Long> offsets(nlocal+1);
        auto* pv = vals.data();
        auto* pc = cols.data();
        auto* po = offsets.data();
        Long const ib = xvec.globalBegin();
        ParallelFor(nlocal+1, [=] AMREX_GPU_DEVICE (Long lrow)
        {
            po[lrow] = lrow*nnz_row;
            if (lrow < nlocal) { set_stencil(lrow+ib, pc+lrow*nnz_row, pv+lrow*nnz_row); }
        });
        Gpu::streamSynchronize();
        mat.define(xvec.partition(), pv, pc, nlocal*nnz_row, po, CsrSorted{false},
                   CsrValid{!dirichlet});
    }
    if (ptype != 0 || dirichlet) { SpMV(bvec, mat, exact); } // rhs = A*phi

    AMG<Real> amg(mat);
    amg.setMaxIter(p.max_iter);
    amg.setFixedIter(p.fixed_iter);
    amg.setRelTol(p.reltol);
    if (p.verbose) { amg.setVerbose(*p.verbose); }
    if (p.nu1) { amg.setPreSmooth(*p.nu1); }
    if (p.nu2) { amg.setPostSmooth(*p.nu2); }
    if (p.nu_bottom) { amg.setBottomSmooth(*p.nu_bottom); }
    if (p.bottom_tol) { amg.setBottomTol(*p.bottom_tol); }
    if (p.p_max_elmts) { amg.setPMaxElmts(*p.p_max_elmts); }
    if (p.trunc_factor) { amg.setTruncFactor(*p.trunc_factor); }
    if (p.max_levels) { amg.setMaxLevels(*p.max_levels); }
    if (p.aggressive_levels) { amg.setAggressiveNumLevels(*p.aggressive_levels); }
    if (p.aggressive_direct) { amg.setAggressiveDirectInterp(*p.aggressive_direct != 0); }
    bool const singular = (p.alpha == Real(0) && !dirichlet);
    if (singular) { amg.setSingular(true); }
    if (p.max_coarse_size) { amg.setMaxCoarseSize(*p.max_coarse_size); }
    if (interp == "direct") {
        amg.setInterpType(AMG<Real>::InterpType::Direct);
    } else if (interp == "ext") {
        amg.setInterpType(AMG<Real>::InterpType::MMExt);
    } else if (interp == "ext+i") {
        amg.setInterpType(AMG<Real>::InterpType::MMExtI);
    } else {
        amrex::Abort("Unknown interpolation: " + interp);
    }
    if (p.smoother == "jacobi") {
        amg.setSmoother(AMG<Real>::Smoother::Jacobi);
    } else if (p.smoother == "l1jacobi") {
        amg.setSmoother(AMG<Real>::Smoother::L1Jacobi);
    } else if (p.smoother == "chebyshev") {
        amg.setSmoother(AMG<Real>::Smoother::Chebyshev);
    } else if (p.smoother == "l1gs") {
        amg.setSmoother(AMG<Real>::Smoother::L1GaussSeidel);
    } else {
        amrex::Abort("Unknown smoother: " + p.smoother);
    }
    if (p.relax_weight) { amg.setRelaxWeight(*p.relax_weight); }
    if (p.cheby_degree) { amg.setChebyshevDegree(*p.cheby_degree); }
    if (p.cheby_ratio) { amg.setChebyshevRatio(*p.cheby_ratio); }
    if (p.theta) { amg.setStrongThreshold(*p.theta); }
    if (bottom == "jacobi") {
        amg.setBottomSolver(AMG<Real>::BottomSolver::Jacobi);
    } else if (bottom == "bicgstab") {
        amg.setBottomSolver(AMG<Real>::BottomSolver::BiCGStab);
    } else if (bottom == "gmres") {
        amg.setBottomSolver(AMG<Real>::BottomSolver::GMRES);
    } else {
        amrex::Abort("Unknown bottom solver: " + bottom);
    }
    if (p.krylov == "none") {
        amg.setKrylovSolver(AMG<Real>::KrylovSolver::None);
    } else if (p.krylov == "bicgstab") {
        amg.setKrylovSolver(AMG<Real>::KrylovSolver::BiCGStab);
    } else if (p.krylov == "gmres") {
        amg.setKrylovSolver(AMG<Real>::KrylovSolver::GMRES);
    } else if (p.krylov == "pcg") {
        amg.setKrylovSolver(AMG<Real>::KrylovSolver::PCG);
    } else {
        amrex::Abort("Unknown Krylov solver: " + p.krylov);
    }

    auto bnorm = bvec.norm2();
    Gpu::streamSynchronize();
    auto const t0 = amrex::second();
    // Divergence is a failed case, not the end of the run.
    amg.setThrowException(true);
    bool blew_up = false;
    try {
        amg.solve(xvec, bvec);
    } catch (std::runtime_error const&) {
        blew_up = true;
    }
    Gpu::streamSynchronize();
    auto const t1 = amrex::second();

    // Residual computed here rather than taken from the solver.
    AlgVector<Real> rvec(xvec.partition());
    SpMV(rvec, mat, xvec);
    amrex::Axpy(rvec, Real(-1), bvec);
    auto const rnorm = rvec.norm2();
    auto rel_res = rnorm / bnorm;

    // A is symmetric and A e = -r, so |e|_inf <= |e|_2 <= |r|_2 / lambda_min.
    // Periodic: lambda_min is alpha, or for the singular problem (mean-zero
    // error) at least the smallest coefficient times the lowest nonzero
    // Laplacian eigenvalue. Dirichlet: at least alpha plus the smallest
    // coefficient times the lowest eigenvalue with the boundary one cell out.
    Real const cmin = (ptype == 0) ? Real(1) : (ptype == 3) ? std::min(Real(1), p.eps)
                                             : std::min(Real(1), p.jump);
    Real const lam_min = dirichlet
        ? a + cmin * Real(AMREX_SPACEDIM) *
              (Real(2) - Real(2)*std::cos(amrex::Math::pi<Real>()/Real(n_cell+1))) / (dx*dx)
        : singular ? cmin * (Real(2) - Real(2)*std::cos(dx)) / (dx*dx) : a;
    Real const err_bound = Real(1.01) * rnorm / lam_min
        + Real(100) * std::numeric_limits<Real>::epsilon();

    amrex::Axpy(xvec, Real(-1), exact);
    if (singular) { // the solution is defined up to a constant
        AlgVector<Real> ones(xvec.partition());
        ones.setVal(Real(1));
        auto mean = amrex::Dot(xvec, ones) / Real(xvec.partition().numGlobalRows());
        amrex::Axpy(xvec, -mean, ones);
    }
    auto error = xvec.norminf();
    return {.error = error, .err_bound = err_bound, .rel_res = rel_res,
            .niters = amg.getNumIters(), .nlevels = amg.numLevels(),
            .time = t1-t0, .blew_up = blew_up};
}

// Options of the case that differ from the base options, other than those
// with their own column, as inputs that reproduce it.
std::string variation (Params const& b, Params const& c)
{
    std::ostringstream os;
    auto sep = [&] () { if (os.tellp() > 0) { os << " "; } };
    auto opt = [&] (char const* name, auto const& x, auto const& y) {
        if (y && x != y) { sep(); os << name << "=" << *y; }
    };
    if (c.alpha != b.alpha) { sep(); os << "alpha=" << c.alpha; }
    if (c.fixed_iter != b.fixed_iter) { sep(); os << "fixed_iter=" << c.fixed_iter; }
    if (c.max_iter != b.max_iter) { sep(); os << "max_iter=" << c.max_iter; }
    opt("aggressive_levels", b.aggressive_levels, c.aggressive_levels);
    opt("aggressive_direct", b.aggressive_direct, c.aggressive_direct);
    opt("p_max_elmts", b.p_max_elmts, c.p_max_elmts);
    opt("trunc_factor", b.trunc_factor, c.trunc_factor);
    opt("max_coarse_size", b.max_coarse_size, c.max_coarse_size);
    opt("max_levels", b.max_levels, c.max_levels);
    return os.str();
}

// Options given in the inputs besides the ones with their own column.
std::string given_options (Params const& p)
{
    std::ostringstream os;
    auto opt = [&] (char const* name, auto const& x) {
        if (x) { os << " " << name << "=" << *x; }
    };
    opt("theta", p.theta);
    opt("aggressive_levels", p.aggressive_levels);
    opt("aggressive_direct", p.aggressive_direct);
    opt("p_max_elmts", p.p_max_elmts);
    opt("trunc_factor", p.trunc_factor);
    opt("max_coarse_size", p.max_coarse_size);
    opt("max_levels", p.max_levels);
    opt("nu1", p.nu1);
    opt("nu2", p.nu2);
    opt("nu_bottom", p.nu_bottom);
    opt("bottom_tol", p.bottom_tol);
    opt("relax_weight", p.relax_weight);
    opt("cheby_degree", p.cheby_degree);
    opt("cheby_ratio", p.cheby_ratio);
    return os.str();
}

std::string problem_line (Params const& p)
{
    std::ostringstream os;
    os << p.problem << ": n_cell=" << p.n_cell << " bc=" << p.bc << " alpha=" << p.alpha;
    if (p.problem == "jump" || p.problem == "checker") { os << " jump=" << p.jump; }
    if (p.problem == "checker") { os << " block=" << p.block; }
    if (p.problem == "aniso") { os << " eps=" << p.eps; }
    return os.str();
}

// Table columns: case, interp, smoother, bottom, krylov, levels, iterations,
// time, rel_res, error, bound, result, variation.
std::string table_row (std::string const& icase, std::string const& interp,
                       std::string const& smoother, std::string const& bottom,
                       std::string const& krylov, std::string const& lev,
                       std::string const& iter, std::string const& time,
                       std::string const& rel_res, std::string const& error,
                       std::string const& bound, std::string const& result,
                       std::string const& var)
{
    std::ostringstream os;
    os << std::right << std::setw(4) << icase << "   " << std::left
       << std::setw(8) << interp << std::setw(11) << smoother
       << std::setw(10) << bottom << std::setw(10) << krylov << std::right
       << std::setw(4) << lev << std::setw(6) << iter << std::setw(10) << time
       << std::setw(11) << rel_res << std::setw(11) << error << std::setw(11) << bound
       << "   " << std::left;
    if (var.empty()) {
        os << result;
    } else {
        os << std::setw(8) << result << var;
    }
    return os.str();
}

void print_row (int icase, Params const& c, Result const& r, std::string const& result,
                std::string const& var)
{
    std::ostringstream time;
    time << std::fixed << std::setprecision(4) << r.time;
    amrex::Print() << table_row(std::to_string(icase), c.interp, c.smoother, c.bottom,
                                c.krylov, std::to_string(r.nlevels), std::to_string(r.niters),
                                time.str(), sci(r.rel_res), sci(r.error), sci(r.err_bound),
                                result, var) << "\n";
}

// The cases run for one problem: the given options, or with `variations`,
// the given options with each bottom solver and variations of the other
// options with the BiCGStab bottom solver.
Vector<Params> make_cases (Params const& p)
{
    Vector<Params> cases;
    auto add = [&] (std::string const& bottom) -> Params& {
        cases.push_back(p);
        cases.back().bottom = bottom;
        return cases.back();
    };
    if (!p.variations) {
        cases.push_back(p);
        return cases;
    }
    for (auto const& b : {"jacobi", "bicgstab", "gmres"}) { add(b); }
    for (auto const& it : {"direct", "ext", "ext+i"}) {
        if (it != p.interp) { add("bicgstab").interp = it; }
    }
    if (!p.aggressive_levels) {
        add("bicgstab").aggressive_levels = 1;
        if (!p.aggressive_direct) {
            auto& c = add("bicgstab");
            c.aggressive_levels = 1;
            c.aggressive_direct = 1;
        }
    }
    if (!p.p_max_elmts && !p.trunc_factor) {
        add("bicgstab").p_max_elmts = 0; // no truncation
        add("bicgstab").trunc_factor = Real(0.2);
    }
    if (!p.max_coarse_size) { add("bicgstab").max_coarse_size = 1; }
    // PCG needs a symmetric cycle: as many pre- as post-smoothing sweeps and
    // smoother sweeps at the bottom.
    bool const symmetric = (p.nu1 == p.nu2);
    if (p.alpha != Real(0)) {
        // Singular periodic problem: plain cycles and PCG.
        add("bicgstab").alpha = Real(0);
        if (symmetric) {
            auto& c = add("jacobi");
            c.alpha = Real(0);
            c.krylov = "pcg";
        }
    }
    if (p.krylov == "none") {
        // GMRES with smoother sweeps at the bottom, so that the cycle is a
        // fixed linear operator; BiCGStab tolerates the Krylov bottom solver.
        add("bicgstab").krylov = "bicgstab";
        add("jacobi").krylov = "gmres";
        if (symmetric) { add("jacobi").krylov = "pcg"; }
    }
    Vector<std::string> smoothers = {"jacobi", "l1jacobi", "chebyshev"};
#ifndef AMREX_USE_GPU
    smoothers.push_back("l1gs");
#endif
    for (auto const& sm : smoothers) {
        if (sm != p.smoother) { add("bicgstab").smoother = sm; }
    }
#ifndef AMREX_USE_GPU
    if (p.krylov == "none") {
        // PCG with symmetric Gauss-Seidel sweeps at the bottom: one level
        // makes the bottom the whole problem.
        auto& c = add("jacobi");
        c.smoother = "l1gs";
        c.krylov = "pcg";
        c.max_levels = 1;
        c.max_iter = std::max(p.max_iter, 50*p.n_cell); // not scalable
    }
#endif
    return cases;
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    {
        ParmParse pp;
        Params p;
        pp.query("n_cell", p.n_cell);
        pp.query("max_iter", p.max_iter);
        pp.query("repeat", p.repeat);
        pp.query("fixed_iter", p.fixed_iter);
        pp.query("reltol", p.reltol);
        pp.query("alpha", p.alpha);
        pp.query("variations", p.variations);
        pp.query("bottom", p.bottom);
        pp.query("interp", p.interp);
        pp.query("smoother", p.smoother);
        pp.query("krylov", p.krylov);
        pp.query("bc", p.bc);
        if (p.bc != "periodic" && p.bc != "dirichlet") { amrex::Abort("Unknown bc: " + p.bc); }
        Vector<std::string> problems; // queryarr does not shrink a vector
        pp.queryarr("problem", problems);
        if (problems.empty()) { problems = {"constant", "jump", "checker", "aniso"}; }
        pp.query("jump", p.jump);
        pp.query("block", p.block);
        pp.query("eps", p.eps);
        pp.query("mlmg", p.mlmg);
        pp.query("max_grid_size", p.max_grid_size);
        // Seed of the random PMIS weights; each rank adds its rank.
        if (Long seed = 0; pp.query("seed", seed)) {
            auto const s = ULong(seed + ParallelDescriptor::MyProc());
            amrex::ResetRandomSeed(s, s);
        }
        auto query_opt = [&] (char const* name, auto& opt) {
            std::remove_reference_t<decltype(*opt)> v;
            if (pp.query(name, v)) { opt = v; }
        };
        query_opt("verbose", p.verbose);
        query_opt("nu1", p.nu1);
        query_opt("nu2", p.nu2);
        query_opt("nu_bottom", p.nu_bottom);
        query_opt("bottom_tol", p.bottom_tol);
        query_opt("p_max_elmts", p.p_max_elmts);
        query_opt("trunc_factor", p.trunc_factor);
        query_opt("relax_weight", p.relax_weight);
        query_opt("cheby_degree", p.cheby_degree);
        query_opt("cheby_ratio", p.cheby_ratio);
        query_opt("theta", p.theta);
        query_opt("max_levels", p.max_levels);
        query_opt("aggressive_levels", p.aggressive_levels);
        query_opt("aggressive_direct", p.aggressive_direct);
        query_opt("max_coarse_size", p.max_coarse_size);

        if (p.n_cell < 3) {
            // Periodic: the stencil would wrap onto one cell. Dirichlet: the
            // exact solution would be zero.
            amrex::Abort("n_cell must be at least 3");
        }

        int const nprocs = ParallelDescriptor::NProcs();
        amrex::Print() << "\nAMG test: " << AMREX_SPACEDIM << "D, " << nprocs
                       << (nprocs == 1 ? " MPI process" : " MPI processes")
                       << ", reltol=" << p.reltol << ", max_iter=" << p.max_iter
                       << (p.variations ? ", with variations" : "") << "\n"
                       << "Base options: interp=" << p.interp << " smoother=" << p.smoother
                       << " bottom=" << p.bottom << " krylov=" << p.krylov
                       << given_options(p) << "\n"
                       << "A case passes if rel_res < reltol (rel_res < 0.9 with fixed_iter)"
                       << " and error <= bound.\n"
                       << "Variation lists the inputs that differ from the base options.\n";
        std::string const head = table_row("#", "interp", "smoother", "bottom", "krylov",
                                           "lev", "iter", "time[s]", "rel_res", "error",
                                           "bound", "result", "variation");
        std::string const rule(head.size(), '-');

        // Report every failed case, then abort once at the end.
        int ncases = 0;
        Vector<std::string> failures;
        for (auto const& problem : problems) {
            p.problem = problem;
            amrex::Print() << "\n" << problem_line(p) << "\n" << rule << "\n" << head
                           << "\n" << rule << "\n";
            int icase = 0;
            for (auto pb : make_cases(p)) {
                // Jacobi bottom: an inexact coarse solve, so only check that
                // ten cycles reduce the residual.
                bool const inexact = (pb.bottom == "jacobi" && pb.krylov == "none"
                                      && pb.fixed_iter <= 0);
                if (inexact) { pb.fixed_iter = 10; }
                ++icase;
                for (int rep = 0; rep < (inexact ? 1 : p.repeat); ++rep) {
                    auto r = run(pb);
                    ++ncases;
                    std::string why;
                    if (r.blew_up) {
                        why = "diverged";
                    } else if (inexact && !(r.rel_res < Real(0.9))) {
                        why = "rel_res " + sci(r.rel_res) + " >= 0.9";
                    } else if (pb.fixed_iter <= 0 && !(r.rel_res < pb.reltol)) {
                        why = "rel_res " + sci(r.rel_res) + " >= reltol " + sci(pb.reltol);
                    }
                    if (!(r.error <= r.err_bound)) {
                        why += (why.empty() ? "" : ", ");
                        why += "error " + sci(r.error) + " > bound " + sci(r.err_bound);
                    }
                    print_row(icase, pb, r, why.empty() ? "pass" : "FAIL",
                              variation(p, pb));
                    if (!why.empty()) {
                        failures.push_back(problem + " #" + std::to_string(icase) + ": " + why);
                    }
                }
            }
            amrex::Print() << rule << "\n";
            if (p.mlmg) { run_mlmg(p); }
        }
        amrex::Print() << "\nSummary: " << ncases << (ncases == 1 ? " case, " : " cases, ")
                       << ncases - failures.size()
                       << " passed, " << failures.size() << " failed\n";
        for (auto const& f : failures) { amrex::Print() << "  FAILED " << f << "\n"; }
        amrex::Print() << "\n";
        if (!failures.empty()) {
            amrex::Abort(std::to_string(failures.size()) + " AMG case(s) failed");
        }
    }
    amrex::Finalize();
}
