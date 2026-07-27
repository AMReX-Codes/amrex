#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

#if (AMREX_SPACEDIM != 1)
#error "Tests/StochasticHeatEquation is a one-dimensional stochastic heat equation test."
#endif

using namespace amrex;

namespace
{
struct Inputs
{
    int npts = 32;
    int max_grid_size = 1024;
    int iper = 1;
    int seed = 43064;
    int ntherm = 0;
    int nstep = 10;
    int nout = -1;
    int nstat = 1;
    int icor = 1;
    int jmidl = -5;
    int jmidr = -5;
    int jpartl = -5;
    int jpartr = -5;
    int is_hybrid = 0;
    int ensemble = 1;
    int is_gaussian = 0;
    int irestart = 0;
    int ipdf = 0;
    int nbins = 150;
    int ensout = 1;
    int iresl = 0;
    int iresr = 0;
    int iopt = 1;
    int plot_int = -1;
    int write_ascii = 1;

    Real xlen = Real(1.0);
    Real dorand = Real(1.0);
    Real dt = Real(1.e-12);
    Real cfl = Real(-1.0);
    Real uleft = Real(300.0);
    Real uright = Real(300.0);
    Real uinit = Real(300.0);
    Real midfact = Real(1.0);
    Real binlo = Real(0.0);
    Real dbin = Real(1.0);
    Real lambda = Real(70.0);
    Real rho = Real(7870.0);
    Real cv = Real(450.0);
    Real crossA = Real(4.e-18);
};

enum StatComp : int {
    MeanAcc = 0,
    SqrAcc,
    CubAcc,
    FourthAcc,
    Cor11Acc,
    Cor12Acc,
    Cor13Acc,
    Cor21Acc,
    Cor31Acc,
    Cor22Acc,
    NumStatComps
};

enum DiagComp : int {
    Mean = 0,
    Var,
    Skew,
    Kurtosis,
    Mom3,
    Mom4,
    Cor11,
    Cor12,
    Cor13,
    Cor21,
    Cor31,
    Cor22,
    NumDiagComps
};

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real reference_state (int i, Real uinit, Real midfact, int jmidl, int jmidr) noexcept
{
    Real value = uinit;
    int const fortran_cell = i + 1;
    if (fortran_cell >= jmidl && fortran_cell <= jmidr) {
        value *= midfact;
    }
    return value;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double precision_independent_uniform (RandomEngine const& engine) noexcept
{
    constexpr unsigned int base = 65536U;
    unsigned int const hi = amrex::Random_int(base, engine);
    unsigned int const lo = amrex::Random_int(base, engine);
    return (static_cast<double>(hi) * static_cast<double>(base) +
            static_cast<double>(lo) + 0.5) * 2.3283064365386962890625e-10;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double inverse_standard_normal (double p) noexcept
{
    constexpr double a0 = -3.969683028665376e+01;
    constexpr double a1 =  2.209460984245205e+02;
    constexpr double a2 = -2.759285104469687e+02;
    constexpr double a3 =  1.383577518672690e+02;
    constexpr double a4 = -3.066479806614716e+01;
    constexpr double a5 =  2.506628277459239e+00;

    constexpr double b0 = -5.447609879822406e+01;
    constexpr double b1 =  1.615858368580409e+02;
    constexpr double b2 = -1.556989798598866e+02;
    constexpr double b3 =  6.680131188771972e+01;
    constexpr double b4 = -1.328068155288572e+01;

    constexpr double c0 = -7.784894002430293e-03;
    constexpr double c1 = -3.223964580411365e-01;
    constexpr double c2 = -2.400758277161838e+00;
    constexpr double c3 = -2.549732539343734e+00;
    constexpr double c4 =  4.374664141464968e+00;
    constexpr double c5 =  2.938163982698783e+00;

    constexpr double d0 = 7.784695709041462e-03;
    constexpr double d1 = 3.224671290700398e-01;
    constexpr double d2 = 2.445134137142996e+00;
    constexpr double d3 = 3.754408661907416e+00;

    constexpr double lo = 0.02425;
    constexpr double hi = 0.97575;

    if (p < lo) {
        double const q = std::sqrt(-2.0 * std::log(p));
        return (((((c0*q + c1)*q + c2)*q + c3)*q + c4)*q + c5) /
            ((((d0*q + d1)*q + d2)*q + d3)*q + 1.0);
    }
    if (p > hi) {
        double const q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c0*q + c1)*q + c2)*q + c3)*q + c4)*q + c5) /
            ((((d0*q + d1)*q + d2)*q + d3)*q + 1.0);
    }

    double const q = p - 0.5;
    double const r = q * q;
    return (((((a0*r + a1)*r + a2)*r + a3)*r + a4)*r + a5) * q /
        (((((b0*r + b1)*r + b2)*r + b3)*r + b4)*r + 1.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double local_random_normal (double mean, double stddev, RandomEngine const& engine) noexcept
{
    return mean + stddev * inverse_standard_normal(precision_independent_uniform(engine));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void kahan_add (Array4<Real> const& sum, Array4<Real> const& comp,
                int i, int j, int k, int n, Real value) noexcept
{
    Real const y = value - comp(i,j,k,n);
    Real const t = sum(i,j,k,n) + y;
    comp(i,j,k,n) = (t - sum(i,j,k,n)) - y;
    sum(i,j,k,n) = t;
}

struct KahanSum
{
    Real sum = Real(0.0);
    Real comp = Real(0.0);

    void add (Real value) noexcept
    {
        Real const y = value - comp;
        Real const t = sum + y;
        comp = (t - sum) - y;
        sum = t;
    }

    Real value () const noexcept { return sum; }
};

Inputs read_inputs ()
{
    Inputs p;
    ParmParse pp;

    pp.query("npts", p.npts);
    pp.query("max_grid_size", p.max_grid_size);
    pp.query("xlen", p.xlen);
    pp.query("iper", p.iper);
    pp.query("dorand", p.dorand);
    pp.query("seed", p.seed);
    pp.query("ntherm", p.ntherm);
    pp.query("nstep", p.nstep);
    pp.query("dt", p.dt);
    pp.query("nout", p.nout);
    pp.query("nstat", p.nstat);
    pp.query("icor", p.icor);
    pp.query("uleft", p.uleft);
    pp.query("uright", p.uright);
    pp.query("uinit", p.uinit);
    pp.query("cfl", p.cfl);
    pp.query("jmidl", p.jmidl);
    pp.query("jmidr", p.jmidr);
    pp.query("midfact", p.midfact);
    pp.query("jpartl", p.jpartl);
    pp.query("jpartr", p.jpartr);
    pp.query("is_hybrid", p.is_hybrid);
    pp.query("ensemble", p.ensemble);
    pp.query("is_gaussian", p.is_gaussian);
    pp.query("irestart", p.irestart);
    pp.query("ipdf", p.ipdf);
    pp.query("binlo", p.binlo);
    pp.query("nbins", p.nbins);
    pp.query("dbin", p.dbin);
    pp.query("ensout", p.ensout);
    pp.query("lambda", p.lambda);
    pp.query("rho", p.rho);
    pp.query("cv", p.cv);
    pp.query("crossA", p.crossA);
    pp.query("iresl", p.iresl);
    pp.query("iresr", p.iresr);
    pp.query("iopt", p.iopt);
    pp.query("plot_int", p.plot_int);
    pp.query("write_ascii", p.write_ascii);

    if (p.npts <= 0) {
        amrex::Abort("npts must be positive");
    }
    if (p.max_grid_size <= 0) {
        amrex::Abort("max_grid_size must be positive");
    }
    if (p.ensemble <= 0) {
        amrex::Abort("ensemble must be positive");
    }
    if (p.ensout <= 0) {
        amrex::Abort("ensout must be positive");
    }
    if (p.icor < 1 || p.icor > p.npts) {
        amrex::Abort("icor must be in [1,npts]");
    }
    if (p.is_hybrid != 0) {
        amrex::Abort("The AMReX C++ port implements the finite-volume solver path; is_hybrid must be 0");
    }
    if (p.irestart != 0) {
        amrex::Abort("Restart input is not implemented in the AMReX C++ port");
    }
    if (p.lambda <= Real(0.0) || p.rho <= Real(0.0) || p.cv <= Real(0.0) ||
        p.crossA <= Real(0.0)) {
        amrex::Abort("lambda, rho, cv and crossA must be positive");
    }
    if (p.ipdf != 0 && (p.nbins <= 0 || p.dbin <= Real(0.0))) {
        amrex::Abort("PDF output requires nbins > 0 and dbin > 0");
    }

    return p;
}

Geometry make_geometry (Inputs const& p)
{
    Box const domain(IntVect(0), IntVect(p.npts-1));
    RealBox const real_box({AMREX_D_DECL(Real(0.0), Real(0.0), Real(0.0))},
                           {AMREX_D_DECL(p.xlen, Real(1.0), Real(1.0))});
    Array<int,AMREX_SPACEDIM> const periodic{AMREX_D_DECL(p.iper, 0, 0)};

    Geometry geom;
    geom.define(domain, real_box, CoordSys::cartesian, periodic);
    return geom;
}

Vector<std::string> state_names (int ncoef)
{
    Vector<std::string> names;
    names.reserve(ncoef);
    names.push_back("u");
    for (int n = 1; n < ncoef; ++n) {
        names.push_back("coef" + std::to_string(n));
    }
    return names;
}

void fill_physical_boundary (MultiFab& u, Geometry const& geom, Inputs const& p)
{
    u.FillBoundary(geom.periodicity());

    if (p.iper == 1) {
        return;
    }

    Box const domain = geom.Domain();
    int const domlo = domain.smallEnd(0);
    int const domhi = domain.bigEnd(0);
    int const ncomp = u.nComp();
    Real const uleft = p.uleft;
    Real const uright = p.uright;
    int const iresl = p.iresl;
    int const iresr = p.iresr;

    for (MFIter mfi(u); mfi.isValid(); ++mfi) {
        Array4<Real> const ua = u.array(mfi);
        Box const fabbox = u[mfi].box();

        Box const left(IntVect(domlo-1), IntVect(domlo-1));
        if (fabbox.contains(left)) {
            ParallelFor(left, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                ua(i,j,k,n) = (iresl == 1) ? uleft : ua(domlo,j,k,n);
            });
        }

        Box const right(IntVect(domhi+1), IntVect(domhi+1));
        if (fabbox.contains(right)) {
            ParallelFor(right, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                ua(i,j,k,n) = (iresr == 1) ? uright : ua(domhi,j,k,n);
            });
        }
    }
}

void initialize_state (MultiFab& u, Inputs const& p)
{
    Real const uinit = p.uinit;
    Real const midfact = p.midfact;
    int const jmidl = p.jmidl;
    int const jmidr = p.jmidr;

    for (MFIter mfi(u); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real> const ua = u.array(mfi);

        ParallelFor(bx, u.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            ua(i,j,k,n) = reference_state(i, uinit, midfact, jmidl, jmidr);
        });
    }
}

void compute_flux (MultiFab& flux, MultiFab& ranflux, MultiFab& u, Geometry const& geom,
                   Inputs const& p, Real dx, Real dt, Real kappa, Real alpha, int ncoef)
{
    fill_physical_boundary(u, geom, p);

    for (MFIter mfi(ranflux); mfi.isValid(); ++mfi) {
        Box const& fbx = mfi.validbox();
        Array4<Real> const ra = ranflux.array(mfi);

        ParallelForRNG(fbx, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                                  RandomEngine const& engine) noexcept
        {
            double const normal = local_random_normal(0.0, 1.0, engine);
            ra(i,j,k) = (normal >= 0.5) ? Real(1.0) : Real(-1.0);
        });
    }
    ranflux.OverrideSync(geom.periodicity());

    Box const domain = geom.Domain();
    int const domlo = domain.smallEnd(0);
    int const domhi = domain.bigEnd(0);
    int const iper = p.iper;
    int const iresl = p.iresl;
    int const iresr = p.iresr;
    int const iopt = p.iopt;
    Real const uleft = p.uleft;
    Real const uright = p.uright;
    Real const dorand = p.dorand;
    Real const sqrt_two = std::sqrt(Real(2.0));
    Real const noise_scale = dorand / std::sqrt(dx*dt);
    int const coef_comp = ncoef - 1;

    for (MFIter mfi(u); mfi.isValid(); ++mfi) {
        Box const fbx = amrex::surroundingNodes(mfi.validbox(), 0);
        Array4<Real const> const ua = u.const_array(mfi);
        Array4<Real const> const ra = ranflux.const_array(mfi);
        Array4<Real> const fa = flux.array(mfi);

        ParallelFor(fbx, ncoef, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            Real f = kappa * (ua(i,j,k,n) - ua(i-1,j,k,n)) / dx;

            if (iper == 0) {
                if (i == domlo) {
                    f = (iresl == 0) ? Real(0.0) : Real(2.0)*f;
                }
                if (i == domhi+1) {
                    f = (iresr == 0) ? Real(0.0) : Real(2.0)*f;
                }
            }

            fa(i,j,k,n) = f;
        });

        ParallelFor(fbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const umin = ua(i-1,j,k,coef_comp);
            Real const uplus = ua(i,j,k,coef_comp);
            Real uave = Real(0.5) * (uplus + umin);

            if (iper == 0 && i == domlo && iresl == 1) {
                if (iopt == 1) {
                    uave = uleft;
                } else if (iopt == 2) {
                    uave = uplus;
                }
            }

            if (iper == 0 && i == domhi+1 && iresr == 1) {
                if (iopt == 1) {
                    uave = uright;
                } else if (iopt == 2) {
                    uave = umin;
                }
            }

            Real factor = Real(1.0);
            if (iper == 0 && i == domlo) {
                factor = (iresl == 1) ? sqrt_two : Real(0.0);
            }
            if (iper == 0 && i == domhi+1) {
                factor = (iresr == 1) ? sqrt_two : Real(0.0);
            }

            fa(i,j,k,0) += alpha * factor * uave * ra(i,j,k) * noise_scale;
        });
    }

    flux.OverrideSync(geom.periodicity());
}

void advance (MultiFab& u, MultiFab& unew, MultiFab const& flux, Real dx, Real dt)
{
    for (MFIter mfi(u); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real const> const ua = u.const_array(mfi);
        Array4<Real const> const fa = flux.const_array(mfi);
        Array4<Real> const unewa = unew.array(mfi);

        ParallelFor(bx, u.nComp(), [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
        {
            unewa(i,j,k,n) = ua(i,j,k,n) + dt * (fa(i+1,j,k,n) - fa(i,j,k,n)) / dx;
        });
    }

    MultiFab::Copy(u, unew, 0, 0, u.nComp(), 0);
}

void print_state (MultiFab const& u, Geometry const& geom, int step, Real time, int ens)
{
    Gpu::HostVector<Real> const line = amrex::sumToLine(u, 0, 1, geom.Domain(), 0, false);

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    Real const dx = geom.CellSize(0);
    amrex::Print() << "\n\n step,time = " << std::setw(10) << step << " "
                   << std::scientific << std::setprecision(7) << time
                   << "  ensemble = " << ens << "\n\n";
    for (int i = 0; i < geom.Domain().length(0); ++i) {
        Real const x = (Real(i) + Real(0.5)) * dx;
        amrex::Print() << std::scientific << std::setprecision(7)
                       << std::setw(15) << x << std::setw(15) << line[i] << "\n";
    }
}

void accumulate_stats (MultiFab& stats, MultiFab& stats_comp, MultiFab const& u,
                       Inputs const& p, Real uicor)
{
    Real const uinit = p.uinit;
    Real const midfact = p.midfact;
    int const jmidl = p.jmidl;
    int const jmidr = p.jmidr;
    Real const xref = reference_state(p.icor-1, uinit, midfact, jmidl, jmidr);
    Real const x = uicor - xref;
    Real const x2 = x * x;
    Real const x3 = x2 * x;

    for (MFIter mfi(stats); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real const> const ua = u.const_array(mfi);
        Array4<Real> const sa = stats.array(mfi);
        Array4<Real> const ca = stats_comp.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const y = ua(i,j,k,0) - reference_state(i, uinit, midfact, jmidl, jmidr);
            Real const y2 = y * y;
            Real const y3 = y2 * y;
            Real const y4 = y2 * y2;

            kahan_add(sa, ca, i, j, k, MeanAcc, y);
            kahan_add(sa, ca, i, j, k, SqrAcc, y2);
            kahan_add(sa, ca, i, j, k, CubAcc, y3);
            kahan_add(sa, ca, i, j, k, FourthAcc, y4);
            kahan_add(sa, ca, i, j, k, Cor11Acc, x * y);
            kahan_add(sa, ca, i, j, k, Cor12Acc, x * y2);
            kahan_add(sa, ca, i, j, k, Cor13Acc, x * y3);
            kahan_add(sa, ca, i, j, k, Cor21Acc, x2 * y);
            kahan_add(sa, ca, i, j, k, Cor31Acc, x3 * y);
            kahan_add(sa, ca, i, j, k, Cor22Acc, x2 * y2);
        });
    }
}

void build_diagnostics (MultiFab& diag, MultiFab const& stats, Inputs const& p, Long istat)
{
    if (istat <= 0) {
        amrex::Abort("Cannot build diagnostics before statistics have been accumulated");
    }

    Real const inv = Real(1.0) / Real(istat);
    Box const icor_box(IntVect(p.icor-1), IntVect(p.icor-1));
    Real const ex = stats.sum(icor_box, MeanAcc) * inv;
    Real const ex2 = stats.sum(icor_box, SqrAcc) * inv;
    Real const ex3 = stats.sum(icor_box, CubAcc) * inv;

    Real const uinit = p.uinit;
    Real const midfact = p.midfact;
    int const jmidl = p.jmidl;
    int const jmidr = p.jmidr;

    for (MFIter mfi(diag); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.validbox();
        Array4<Real const> const sa = stats.const_array(mfi);
        Array4<Real> const da = diag.array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real const ey = sa(i,j,k,MeanAcc) * inv;
            Real const ey2 = sa(i,j,k,SqrAcc) * inv;
            Real const ey3 = sa(i,j,k,CubAcc) * inv;
            Real const ey4 = sa(i,j,k,FourthAcc) * inv;

            Real const exy = sa(i,j,k,Cor11Acc) * inv;
            Real const exy2 = sa(i,j,k,Cor12Acc) * inv;
            Real const exy3 = sa(i,j,k,Cor13Acc) * inv;
            Real const ex2y = sa(i,j,k,Cor21Acc) * inv;
            Real const ex3y = sa(i,j,k,Cor31Acc) * inv;
            Real const ex2y2 = sa(i,j,k,Cor22Acc) * inv;

            Real const ey_sq = ey * ey;
            Real const ey_cu = ey_sq * ey;
            Real const ey_qu = ey_sq * ey_sq;
            Real const ex_sq = ex * ex;
            Real const ex_cu = ex_sq * ex;
            Real const mean = reference_state(i, uinit, midfact, jmidl, jmidr) + ey;
            Real const var = ey2 - ey_sq;
            Real const mom3 = ey3 - Real(3.0)*ey2*ey + Real(2.0)*ey_cu;
            Real const mom4 = ey4 - Real(4.0)*ey3*ey + Real(6.0)*ey2*ey_sq
                - Real(3.0)*ey_qu;
            Real skew = mom3;
            Real kurtosis = mom4;
            if (var > Real(0.0)) {
                skew = mom3 / (var * std::sqrt(var));
                kurtosis = mom4 / (var * var);
            }

            da(i,j,k,Mean) = mean;
            da(i,j,k,Var) = var;
            da(i,j,k,Skew) = skew;
            da(i,j,k,Kurtosis) = kurtosis;
            da(i,j,k,Mom3) = mom3;
            da(i,j,k,Mom4) = mom4;
            da(i,j,k,Cor11) = exy - ex * ey;
            da(i,j,k,Cor12) = exy2 - Real(2.0)*exy*ey - ey2*ex
                + Real(2.0)*ex*ey_sq;
            da(i,j,k,Cor21) = ex2y - Real(2.0)*exy*ex - ex2*ey
                + Real(2.0)*ex_sq*ey;
            da(i,j,k,Cor13) = exy3 - Real(3.0)*exy2*ey
                + Real(3.0)*exy*ey_sq - ey3*ex
                + Real(3.0)*ey2*ex*ey - Real(3.0)*ex*ey_cu;
            da(i,j,k,Cor31) = ex3y - Real(3.0)*ex2y*ex
                + Real(3.0)*exy*ex_sq - ex3*ey
                + Real(3.0)*ex2*ex*ey - Real(3.0)*ex_cu*ey;
            da(i,j,k,Cor22) = ex2y2 - Real(2.0)*ex2y*ey
                + ex2*ey_sq - Real(2.0)*exy2*ex
                + Real(4.0)*exy*ex*ey + ey2*ex_sq
                - Real(3.0)*ex_sq*ey_sq;
        });
    }
}

std::string output_label (int ens, int step)
{
    std::ostringstream os;
    os << std::setfill('0') << std::setw(6) << ens << "_"
       << std::setw(9) << step;
    return os.str();
}

void write_xy_file (std::string const& filename, Gpu::HostVector<Real> const& line,
                    int comp, int ncomp, Real dx)
{
    std::ofstream os(filename);
    os << std::scientific << std::setprecision(12);
    for (int i = 0; i < static_cast<int>(line.size())/ncomp; ++i) {
        Real const x = (Real(i) + Real(0.5)) * dx;
        os << std::setw(20) << x << std::setw(20) << line[comp + ncomp*i] << "\n";
    }
}

void write_correlation_file (std::string const& filename, Gpu::HostVector<Real> const& line,
                             int comp, int ncomp, Real dx, int icor)
{
    std::ofstream os(filename);
    os << std::scientific << std::setprecision(12);
    int const npts = static_cast<int>(line.size()) / ncomp;
    for (int i = 0; i < npts; ++i) {
        Real const x = (Real(i) + Real(0.5)) * dx;
        os << std::setw(20) << x << std::setw(20) << line[comp + ncomp*i] << "\n";
    }
    os << "\n";
    for (int i = 0; i < npts; ++i) {
        if (i != icor-1) {
            Real const x = (Real(i) + Real(0.5)) * dx;
            os << std::setw(20) << x << std::setw(20) << line[comp + ncomp*i] << "\n";
        }
    }
}

void write_pdf_file (std::string const& filename, Vector<Long> const& pdf,
                     Long total_pdf_points, Inputs const& p)
{
    if (total_pdf_points <= 0) {
        return;
    }

    std::ofstream os(filename);
    os << std::scientific << std::setprecision(12);
    for (int ibin = -1; ibin <= p.nbins; ++ibin) {
        Real const x = p.binlo + (Real(ibin) + Real(0.5)) * p.dbin;
        Real const prob = Real(pdf[ibin+1]) / Real(total_pdf_points);
        os << std::setw(20) << x << std::setw(20) << prob << "\n";
    }
}

void print_and_write_diagnostics (MultiFab const& diag, Geometry const& geom, Inputs const& p,
                                  Long istat, int ens, int step, std::string const& tag,
                                  Vector<Long> const& pdf, Long total_pdf_points)
{
    Gpu::HostVector<Real> const line =
        amrex::sumToLine(diag, 0, NumDiagComps, geom.Domain(), 0, false);

    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    int const npts = geom.Domain().length(0);
    Real const dx = geom.CellSize(0);

    KahanSum av_mean_sum;
    KahanSum av_var_sum;
    KahanSum av_skew_sum;
    KahanSum av_kur_sum;
    KahanSum av_mom3_sum;
    KahanSum av_mom4_sum;
    KahanSum var_mean_sum;
    KahanSum var_var_sum;
    KahanSum var_skew_sum;
    KahanSum var_kur_sum;
    KahanSum var_mom3_sum;
    KahanSum var_mom4_sum;

    KahanSum c11_m_sum;
    KahanSum c12_m_sum;
    KahanSum c13_m_sum;
    KahanSum c21_m_sum;
    KahanSum c31_m_sum;
    KahanSum c22_m_sum;
    KahanSum c11_v_sum;
    KahanSum c12_v_sum;
    KahanSum c13_v_sum;
    KahanSum c21_v_sum;
    KahanSum c31_v_sum;
    KahanSum c22_v_sum;

    for (int i = 0; i < npts; ++i) {
        Real const mean = line[Mean + NumDiagComps*i];
        Real const var = line[Var + NumDiagComps*i];
        Real const skew = line[Skew + NumDiagComps*i];
        Real const kur = line[Kurtosis + NumDiagComps*i];
        Real const mom3 = line[Mom3 + NumDiagComps*i];
        Real const mom4 = line[Mom4 + NumDiagComps*i];

        av_mean_sum.add(mean);
        av_var_sum.add(var);
        av_skew_sum.add(skew);
        av_kur_sum.add(kur);
        av_mom3_sum.add(mom3);
        av_mom4_sum.add(mom4);
        var_mean_sum.add(mean*mean);
        var_var_sum.add(var*var);
        var_skew_sum.add(skew*skew);
        var_kur_sum.add(kur*kur);
        var_mom3_sum.add(mom3*mom3);
        var_mom4_sum.add(mom4*mom4);

        if (i != p.icor-1) {
            Real const c11 = line[Cor11 + NumDiagComps*i];
            Real const c12 = line[Cor12 + NumDiagComps*i];
            Real const c13 = line[Cor13 + NumDiagComps*i];
            Real const c21 = line[Cor21 + NumDiagComps*i];
            Real const c31 = line[Cor31 + NumDiagComps*i];
            Real const c22 = line[Cor22 + NumDiagComps*i];
            c11_m_sum.add(c11);
            c12_m_sum.add(c12);
            c13_m_sum.add(c13);
            c21_m_sum.add(c21);
            c31_m_sum.add(c31);
            c22_m_sum.add(c22);
            c11_v_sum.add(c11*c11);
            c12_v_sum.add(c12*c12);
            c13_v_sum.add(c13*c13);
            c21_v_sum.add(c21*c21);
            c31_v_sum.add(c31*c31);
            c22_v_sum.add(c22*c22);
        }
    }

    Real const inv_npts = Real(1.0) / Real(npts);
    Real const av_mean = av_mean_sum.value() * inv_npts;
    Real const av_var = av_var_sum.value() * inv_npts;
    Real const av_skew = av_skew_sum.value() * inv_npts;
    Real const av_kur = av_kur_sum.value() * inv_npts;
    Real const av_mom3 = av_mom3_sum.value() * inv_npts;
    Real const av_mom4 = av_mom4_sum.value() * inv_npts;
    Real const var_mean = var_mean_sum.value()*inv_npts - av_mean*av_mean;
    Real const var_var = var_var_sum.value()*inv_npts - av_var*av_var;
    Real const var_skew = var_skew_sum.value()*inv_npts - av_skew*av_skew;
    Real const var_kur = var_kur_sum.value()*inv_npts - av_kur*av_kur;
    Real const var_mom3 = var_mom3_sum.value()*inv_npts - av_mom3*av_mom3;
    Real const var_mom4 = var_mom4_sum.value()*inv_npts - av_mom4*av_mom4;
    amrex::ignore_unused(var_mean, var_var, var_skew, var_kur, var_mom3, var_mom4);

    amrex::Print() << "istats " << istat << " " << npts << "\n";
    amrex::Print() << "aver stats " << av_mean << " " << av_var << " " << av_skew << " "
                   << av_kur << " " << av_mom3 << " " << av_mom4 << "\n";

    if (npts > 1) {
        Real const inv_cov = Real(1.0) / Real(npts-1);
        Real const c11_m = c11_m_sum.value() * inv_cov;
        Real const c12_m = c12_m_sum.value() * inv_cov;
        Real const c13_m = c13_m_sum.value() * inv_cov;
        Real const c21_m = c21_m_sum.value() * inv_cov;
        Real const c31_m = c31_m_sum.value() * inv_cov;
        Real const c22_m = c22_m_sum.value() * inv_cov;
        Real const c11_v = c11_v_sum.value()*inv_cov - c11_m*c11_m;
        Real const c12_v = c12_v_sum.value()*inv_cov - c12_m*c12_m;
        Real const c13_v = c13_v_sum.value()*inv_cov - c13_m*c13_m;
        Real const c21_v = c21_v_sum.value()*inv_cov - c21_m*c21_m;
        Real const c31_v = c31_v_sum.value()*inv_cov - c31_m*c31_m;
        Real const c22_v = c22_v_sum.value()*inv_cov - c22_m*c22_m;
        amrex::ignore_unused(c11_v, c12_v, c13_v, c21_v, c31_v, c22_v);
        amrex::Print() << "cov stats " << c11_m << " " << c12_m << " " << c13_m << " "
                       << c21_m << " " << c31_m << " " << c22_m << "\n";
    } else {
        amrex::Print() << "mean = " << line[Mean] << "\n";
        amrex::Print() << "var = " << line[Var] << "\n";
        amrex::Print() << "mom3 = " << line[Mom3] << "\n";
        amrex::Print() << "mom4 = " << line[Mom4] << "\n";
        amrex::Print() << "skew = " << line[Skew] << "\n";
        amrex::Print() << "kur = " << line[Kurtosis] << "\n";
    }

    if (p.write_ascii == 0) {
        return;
    }

    std::string const label = output_label(ens, step);
    write_xy_file("mean_" + tag + label + ".dat", line, Mean, NumDiagComps, dx);
    write_xy_file("var_" + tag + label + ".dat", line, Var, NumDiagComps, dx);
    write_xy_file("skew_" + tag + label + ".dat", line, Skew, NumDiagComps, dx);
    write_xy_file("kur_" + tag + label + ".dat", line, Kurtosis, NumDiagComps, dx);
    write_xy_file("mom3_" + tag + label + ".dat", line, Mom3, NumDiagComps, dx);
    write_xy_file("mom4_" + tag + label + ".dat", line, Mom4, NumDiagComps, dx);
    write_correlation_file("cor11_" + tag + label + ".dat", line, Cor11, NumDiagComps, dx, p.icor);
    write_correlation_file("cor12_" + tag + label + ".dat", line, Cor12, NumDiagComps, dx, p.icor);
    write_correlation_file("cor13_" + tag + label + ".dat", line, Cor13, NumDiagComps, dx, p.icor);
    write_correlation_file("cor21_" + tag + label + ".dat", line, Cor21, NumDiagComps, dx, p.icor);
    write_correlation_file("cor31_" + tag + label + ".dat", line, Cor31, NumDiagComps, dx, p.icor);
    write_correlation_file("cor22_" + tag + label + ".dat", line, Cor22, NumDiagComps, dx, p.icor);

    if (p.ipdf != 0) {
        write_pdf_file("pdf_" + tag + label + ".dat", pdf, total_pdf_points, p);
    }
}

void accumulate_pdf (Vector<Long>& pdf, Long& total_pdf_points, MultiFab const& u,
                     Geometry const& geom, Inputs const& p)
{
    if (p.ipdf == 0) {
        return;
    }

    Gpu::HostVector<Real> const line = amrex::sumToLine(u, 0, 1, geom.Domain(), 0, false);
    if (!ParallelDescriptor::IOProcessor()) {
        return;
    }

    Real const binhi = p.binlo + p.dbin * Real(p.nbins);
    for (Real const value : line) {
        int ibin = 0;
        if (value < p.binlo) {
            ibin = -1;
        } else if (value > binhi) {
            ibin = p.nbins;
        } else {
            ibin = static_cast<int>((value - p.binlo) / p.dbin);
        }
        pdf[ibin+1] += 1;
        total_pdf_points += 1;
    }
}

void maybe_write_plotfile (MultiFab const& u, Geometry const& geom, int step, Real time,
                           int ncoef, int plot_int)
{
    if (plot_int > 0 && step % plot_int == 0) {
        std::string const pltfile = amrex::Concatenate("plt", step, 5);
        WriteSingleLevelPlotfile(pltfile, u, state_names(ncoef), geom, time, step);
    }
}
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        Inputs inputs = read_inputs();

        int const ncoef = (inputs.is_gaussian == 1) ? 2 : 1;
        std::string const tag = (inputs.is_gaussian == 1) ? "ga_" : "fv_";
        Geometry geom = make_geometry(inputs);
        Real const dx = geom.CellSize(0);

        Real const kappa = inputs.lambda / (inputs.rho * inputs.cv);
        Real const kb = Real(1.380649e-23);
        Real const alpha = std::sqrt(Real(2.0)*kb*inputs.lambda/inputs.crossA)
            / (inputs.rho * inputs.cv);

        if (inputs.cfl > Real(0.0)) {
            inputs.dt = inputs.cfl * dx * dx / kappa;
        }
        if (inputs.dt <= Real(0.0)) {
            amrex::Abort("dt must be positive");
        }

        amrex::Print() << "kappa, alpha, rho, cv, lambda = "
                       << kappa << " " << alpha << " " << inputs.rho << " "
                       << inputs.cv << " " << inputs.lambda << "\n";
        amrex::Print() << "dt, dx = " << inputs.dt << " " << dx << "\n";
        amrex::Print() << "delreg, delnreg = " << Real(4.0)*dx*dx << " "
                       << Real(4.0)*dx*dx*inputs.uinit << "\n";

        BoxArray ba(geom.Domain());
        ba.maxSize(inputs.max_grid_size);
        DistributionMapping dm(ba);

        BoxArray fba = amrex::convert(ba, IntVect::TheDimensionVector(0));

        MultiFab u(ba, dm, ncoef, 1);
        MultiFab unew(ba, dm, ncoef, 0);
        MultiFab flux(fba, dm, ncoef, 0);
        MultiFab ranflux(fba, dm, 1, 0);
        MultiFab stats(ba, dm, NumStatComps, 0);
        MultiFab stats_comp(ba, dm, NumStatComps, 0);
        MultiFab diag(ba, dm, NumDiagComps, 0);
        stats.setVal(Real(0.0));
        stats_comp.setVal(Real(0.0));

        Vector<Long> pdf(inputs.nbins+2, Long(0));
        Long total_pdf_points = 0;
        Long istat = 0;
        int const total_steps = inputs.nstep + inputs.ntherm;

        amrex::InitRandom(static_cast<ULong>(inputs.seed), ParallelDescriptor::NProcs(),
                          static_cast<ULong>(inputs.seed));

        for (int ens = 1; ens <= inputs.ensemble; ++ens) {
            initialize_state(u, inputs);
            Real time = Real(0.0);
            int tout = 0;

            print_state(u, geom, 0, time, ens);
            maybe_write_plotfile(u, geom, 0, time, ncoef, inputs.plot_int);

            for (int step = 1; step <= total_steps; ++step) {
                time = Real(step) * inputs.dt;

                compute_flux(flux, ranflux, u, geom, inputs, dx, inputs.dt, kappa, alpha, ncoef);
                advance(u, unew, flux, dx, inputs.dt);

                if (step > inputs.ntherm && inputs.nout > 0 && step % inputs.nout == 0) {
                    Real const total_mass = u.sum(0) * dx;
                    Real const min_mass = u.min(0);
                    amrex::Print() << step << " time = " << time
                                   << " total mass = " << total_mass
                                   << " min state = " << min_mass << "\n";
                }

                if (step > inputs.ntherm && inputs.nstat > 0 && step % inputs.nstat == 0) {
                    Box const icor_box(IntVect(inputs.icor-1), IntVect(inputs.icor-1));
                    Real const uicor = u.sum(icor_box, 0);
                    ++istat;
                    accumulate_stats(stats, stats_comp, u, inputs, uicor);
                    accumulate_pdf(pdf, total_pdf_points, u, geom, inputs);
                    tout = 0;
                }

                if (step > inputs.ntherm && inputs.nout > 0 && step % inputs.nout == 0 &&
                    ens % inputs.ensout == 0) {
                    print_state(u, geom, step, time, ens);
                    if (istat > 0) {
                        build_diagnostics(diag, stats, inputs, istat);
                        print_and_write_diagnostics(diag, geom, inputs, istat, ens, step,
                                                    tag, pdf, total_pdf_points);
                    }
                    tout = 1;
                }

                maybe_write_plotfile(u, geom, step, time, ncoef, inputs.plot_int);
            }

            amrex::Print() << "completed " << ens << "\n";

            if (tout == 0 && ens % inputs.ensout == 0) {
                Real const time = Real(total_steps) * inputs.dt;
                print_state(u, geom, total_steps, time, ens);
            }
        }
    }
    amrex::Finalize();
    return 0;
}
