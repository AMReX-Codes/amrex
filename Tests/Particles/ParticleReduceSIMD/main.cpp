// Test and benchmark for SIMD-vectorized particle reductions:
//   amrex::ParticleReduceSIMD, amrex::ReduceSumSIMD,
//   amrex::ReduceMinSIMD, amrex::ReduceMaxSIMD
//
// The reduction kernel reproduces the "reduced beam characteristics"
// diagnostic of ImpactX: one fused reduction with 34 sums (weight +
// 9 weighted first moments + 24 weighted second moments), 6 mins and
// 6 maxes over a pure-SoA particle container with 11 real components.
//
// Three kernel variants are compared for correctness and speed:
//   A: amrex::ParticleReduce with a SuperParticle lambda (per-particle copy)
//   B: amrex::ParticleReduce with a (ConstParticleTileData, int) lambda
//   C: amrex::ParticleReduceSIMD with a single-source load_1d functor
//
// Compiles and runs with AMReX_SIMD=OFF (scalar fallback, GPU) and ON (CPU SIMD).

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParticleReduceSIMD.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_SIMD.H>
#include <AMReX_Tuple.H>
#include <AMReX_TypeList.H>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

using namespace amrex;

// SoA layout mimicking ImpactX (x, y, t double as the positions in 3D)
struct RealSoA
{
    enum { x = 0, y, t, px, py, pt, sx, sy, sz, qm, w, nattribs };
};

static constexpr int NReal = RealSoA::nattribs;  // 11
static constexpr int NInt = 0;

using PC = ParticleContainerPureSoA<NReal, NInt>;
using ConstPTDType = typename PC::ConstPTDType;
using SPType = typename PC::SuperParticleType;

// 34 sums + 6 mins + 6 maxes, exactly as in ImpactX ReducedBeamCharacteristics
static constexpr std::size_t num_sum = 34;
static constexpr std::size_t num_min = 6;
static constexpr std::size_t num_max = 6;
static constexpr std::size_t num_red = num_sum + num_min + num_max;

using ReduceOpsT = TypeMultiplier<ReduceOps,
    ReduceOpSum[num_sum], ReduceOpMin[num_min], ReduceOpMax[num_max]>;
using ReduceDataT = TypeMultiplier<ReduceData, ParticleReal[num_red]>;
using ReduceTupleT = typename ReduceDataT::Type;

//! Constant shifts subtracted before forming moments (well-conditioned squares)
struct Shifts
{
    ParticleReal x, y, t, px, py, pt, sx, sy, sz;
};

/** The per-particle moments, shared by all kernel variants.
 *
 * T is amrex::ParticleReal in the scalar variants and remainder loop,
 * and a SIMD register of WIDTH values in the vectorized main loop.
 */
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto beam_moments (T const& p_w,
                   T const& p_x, T const& p_y, T const& p_t,
                   T const& p_px, T const& p_py, T const& p_pt,
                   T const& p_sx, T const& p_sy, T const& p_sz,
                   Shifts const& s) noexcept
{
    T const dx  = p_x  - s.x;
    T const dy  = p_y  - s.y;
    T const dt  = p_t  - s.t;
    T const dpx = p_px - s.px;
    T const dpy = p_py - s.py;
    T const dpt = p_pt - s.pt;
    T const dsx = p_sx - s.sx;
    T const dsy = p_sy - s.sy;
    T const dsz = p_sz - s.sz;

    return amrex::makeTuple(
        // Sum(w)
        p_w,
        // weighted first moments (shifted): x, y, t, px, py, pt, sx, sy, sz
        dx*p_w, dy*p_w, dt*p_w, dpx*p_w, dpy*p_w, dpt*p_w,
        dsx*p_w, dsy*p_w, dsz*p_w,
        // weighted second moments (shifted): diagonal x, y, t, px, py, pt
        dx*dx*p_w, dy*dy*p_w, dt*dt*p_w, dpx*dpx*p_w, dpy*dpy*p_w, dpt*dpt*p_w,
        // same-plane correlations: xpx, ypy, tpt
        dx*dpx*p_w, dy*dpy*p_w, dt*dpt*p_w,
        // dispersive correlations: xpt, pxpt, ypt, pypt
        dx*dpt*p_w, dpx*dpt*p_w, dy*dpt*p_w, dpy*dpt*p_w,
        // cross-plane correlations: xy, xpy, xt, pxy, pxpy, pxt, yt, pyt
        dx*dy*p_w, dx*dpy*p_w, dx*dt*p_w, dpx*dy*p_w,
        dpx*dpy*p_w, dpx*dt*p_w, dy*dt*p_w, dpy*dt*p_w,
        // spin second moments (diagonal): sx, sy, sz
        dsx*dsx*p_w, dsy*dsy*p_w, dsz*dsz*p_w,
        // min of x, y, t, px, py, pt
        p_x, p_y, p_t, p_px, p_py, p_pt,
        // max of x, y, t, px, py, pt
        p_x, p_y, p_t, p_px, p_py, p_pt
    );
}

//! Variant C: single-source kernel for ParticleReduceSIMD (scalar, SIMD and GPU)
struct BeamMomentsKernel
{
    Shifts s;

    template <typename PTD, typename SI>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto operator() (PTD const& ptd, SI const si) const noexcept
    {
        auto const p_w  = simd::load_1d(ptd.rdata(RealSoA::w), si);
        auto const p_x  = simd::load_1d(ptd.rdata(RealSoA::x), si);
        auto const p_y  = simd::load_1d(ptd.rdata(RealSoA::y), si);
        auto const p_t  = simd::load_1d(ptd.rdata(RealSoA::t), si);
        auto const p_px = simd::load_1d(ptd.rdata(RealSoA::px), si);
        auto const p_py = simd::load_1d(ptd.rdata(RealSoA::py), si);
        auto const p_pt = simd::load_1d(ptd.rdata(RealSoA::pt), si);
        auto const p_sx = simd::load_1d(ptd.rdata(RealSoA::sx), si);
        auto const p_sy = simd::load_1d(ptd.rdata(RealSoA::sy), si);
        auto const p_sz = simd::load_1d(ptd.rdata(RealSoA::sz), si);

        return beam_moments(p_w, p_x, p_y, p_t, p_px, p_py, p_pt,
                            p_sx, p_sy, p_sz, s);
    }
};

//! Kernel loading a single component, for the Reduce{Sum,Min,Max}SIMD tests
template <int Comp>
struct LoadCompKernel
{
    template <typename PTD, typename SI>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto operator() (PTD const& ptd, SI const si) const noexcept
    {
        return simd::load_1d(ptd.rdata(Comp), si);
    }
};

//! Deterministic pseudo-random number in [0, 1) from an integer (splitmix64 finalizer)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
ParticleReal hash01 (std::uint64_t z) noexcept
{
    z ^= z >> 33; z *= 0xff51afd7ed558ccdULL;
    z ^= z >> 33; z *= 0xc4ceb9fe1a85ec53ULL;
    z ^= z >> 33;
    return static_cast<ParticleReal>(static_cast<double>(z >> 11) * 0x1.0p-53);
}

//! Fill one tile with np particles with deterministic pseudo-random data
void init_particles (PC& pc, int np)
{
    const int myproc = ParallelDescriptor::MyProc();

    // only the rank that owns grid 0 holds particles
    if (pc.ParticleDistributionMap(0)[0] != myproc) { return; }

    auto& ptile = pc.DefineAndReturnParticleTile(0, 0, 0);
    ptile.resize(np);
    auto ptd = ptile.getParticleTileData();

    ParallelFor(np, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        for (int k = 0; k < NReal; ++k) {
            // values in [-1, 1); positions (comps 0..AMREX_SPACEDIM-1) lie
            // inside the [-2, 2] domain
            ptd.rdata(k)[i] = ParticleReal(2.)*hash01(std::uint64_t(i)*NReal + std::uint64_t(k))
                            - ParticleReal(1.);
        }
        // positive weight in (0, 1]
        ptd.rdata(RealSoA::w)[i] = ParticleReal(1.)
            - hash01(std::uint64_t(i)*NReal + std::uint64_t(RealSoA::w));
        ptd.idcpu(i) = SetParticleIDandCPU(i+1, myproc);
    });
    Gpu::streamSynchronize();
}

ReduceTupleT run_variant_a (PC const& pc, Shifts const& s)
{
    ReduceOpsT reduce_ops;
    return ParticleReduce<ReduceDataT>(pc,
        [=] AMREX_GPU_DEVICE (const SPType& p) noexcept
        {
            return beam_moments(p.rdata(RealSoA::w),
                                p.rdata(RealSoA::x), p.rdata(RealSoA::y), p.rdata(RealSoA::t),
                                p.rdata(RealSoA::px), p.rdata(RealSoA::py), p.rdata(RealSoA::pt),
                                p.rdata(RealSoA::sx), p.rdata(RealSoA::sy), p.rdata(RealSoA::sz),
                                s);
        }, reduce_ops);
}

ReduceTupleT run_variant_b (PC const& pc, Shifts const& s)
{
    ReduceOpsT reduce_ops;
    return ParticleReduce<ReduceDataT>(pc,
        [=] AMREX_GPU_DEVICE (const ConstPTDType& ptd, const int i) noexcept
        {
            return beam_moments(ptd.rdata(RealSoA::w)[i],
                                ptd.rdata(RealSoA::x)[i], ptd.rdata(RealSoA::y)[i], ptd.rdata(RealSoA::t)[i],
                                ptd.rdata(RealSoA::px)[i], ptd.rdata(RealSoA::py)[i], ptd.rdata(RealSoA::pt)[i],
                                ptd.rdata(RealSoA::sx)[i], ptd.rdata(RealSoA::sy)[i], ptd.rdata(RealSoA::sz)[i],
                                s);
        }, reduce_ops);
}

ReduceTupleT run_variant_c (PC const& pc, Shifts const& s)
{
    ReduceOpsT reduce_ops;
    return ParticleReduceSIMD<ReduceDataT>(pc, BeamMomentsKernel{s}, reduce_ops);
}

//! like run_variant_c, but with twice the native SIMD width (more ILP)
ReduceTupleT run_variant_c2 (PC const& pc, Shifts const& s)
{
    constexpr int w2 = 2 * static_cast<int>(simd::native_simd_size_particlereal);
    ReduceOpsT reduce_ops;
    return ParticleReduceSIMD<ReduceDataT, w2>(pc, BeamMomentsKernel{s}, reduce_ops);
}

/** Compare two reduction results component-wise.
 *
 * Sums are compared with an absolute tolerance of
 * tol_scale * epsilon * max(np, 1) (the summands are O(1));
 * pass tol_scale = 0 to require bitwise identical sums.
 * Mins and maxes must always match exactly.
 */
void compare_results (ReduceTupleT const& a, ReduceTupleT const& b,
                      int np, Real tol_scale, std::string const& what)
{
    auto const ta = tupleToArray(a);
    auto const tb = tupleToArray(b);

    const auto tol = static_cast<ParticleReal>(tol_scale)
        * std::numeric_limits<ParticleReal>::epsilon()
        * static_cast<ParticleReal>(std::max(np, 1));

    for (int k = 0; k < static_cast<int>(num_red); ++k) {
        const bool is_sum = (k < static_cast<int>(num_sum));
        const auto diff = std::abs(ta[k] - tb[k]);
        if ((is_sum && diff > tol) || (!is_sum && diff != ParticleReal(0.))) {
            amrex::AllPrint() << "  mismatch (" << what << "): np=" << np
                              << " comp=" << k << " a=" << ta[k] << " b=" << tb[k]
                              << " diff=" << diff << " tol=" << tol << "\n";
            amrex::Abort("ParticleReduceSIMD result mismatch: " + what);
        }
    }
}

void correctness_tests (Geometry const& geom, DistributionMapping const& dm,
                        BoxArray const& ba, Shifts const& shifts)
{
    constexpr int W = static_cast<int>(simd::native_simd_size_particlereal);
    std::vector<int> sizes{0, 1, 3, W-1, W, W+1, 2*W+3, 1003};

    for (int np : sizes) {
        if (np < 0) { continue; }

        PC pc(geom, dm, ba);
        init_particles(pc, np);

        auto const ra = run_variant_a(pc, shifts);
        auto const rb = run_variant_b(pc, shifts);
        auto const rc = run_variant_c(pc, shifts);

        // A and B evaluate identical arithmetic in identical order
        compare_results(ra, rb, np, Real(0.), "A (SuperParticle) vs B (ptd,i)");
        // C reassociates the sums across SIMD lanes
        compare_results(rc, rb, np, Real(100.), "C (SIMD) vs B (ptd,i)");
        // non-native (2x) SIMD width
        auto const rc2 = run_variant_c2(pc, shifts);
        compare_results(rc2, rb, np, Real(100.), "C2 (SIMD 2x width) vs B (ptd,i)");

        // value helpers vs the established scalar entry points
        auto sum_simd = ReduceSumSIMD(pc, LoadCompKernel<RealSoA::w>{});
        auto sum_ref = ReduceSum(pc,
            [=] AMREX_GPU_DEVICE (const ConstPTDType& ptd, const int i) noexcept
            { return ptd.rdata(RealSoA::w)[i]; });
        const auto sum_tol = ParticleReal(100.)
            * std::numeric_limits<ParticleReal>::epsilon()
            * static_cast<ParticleReal>(std::max(np, 1));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(std::abs(sum_simd - sum_ref) <= sum_tol,
            "ReduceSumSIMD does not match ReduceSum");

        auto min_simd = ReduceMinSIMD(pc, LoadCompKernel<RealSoA::x>{});
        auto min_ref = ReduceMin(pc,
            [=] AMREX_GPU_DEVICE (const ConstPTDType& ptd, const int i) noexcept
            { return ptd.rdata(RealSoA::x)[i]; });
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(min_simd == min_ref,
            "ReduceMinSIMD does not match ReduceMin");

        auto max_simd = ReduceMaxSIMD(pc, LoadCompKernel<RealSoA::x>{});
        auto max_ref = ReduceMax(pc,
            [=] AMREX_GPU_DEVICE (const ConstPTDType& ptd, const int i) noexcept
            { return ptd.rdata(RealSoA::x)[i]; });
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(max_simd == max_ref,
            "ReduceMaxSIMD does not match ReduceMax");

        Print() << "correctness: np=" << np << " PASSED\n";
    }
}

template <typename FRun>
void bench_variant (std::string const& name, int np, long reps_per_batch,
                    int nbatch, int warmup_reps, FRun const& run, double& checksum)
{
    for (int i = 0; i < warmup_reps; ++i) {
        checksum += static_cast<double>(amrex::get<0>(run()));
    }

    double best = std::numeric_limits<double>::max();
    double total = 0.;
    for (int b = 0; b < nbatch; ++b) {
        const auto t0 = amrex::second();
        for (long r = 0; r < reps_per_batch; ++r) {
            checksum += static_cast<double>(amrex::get<0>(run()));
        }
        const auto dt = amrex::second() - t0;
        best = std::min(best, dt);
        total += dt;
    }

    const auto npd = static_cast<double>(std::max(np, 1));
    const double best_eval = best / static_cast<double>(reps_per_batch);
    const double mean_eval = total / static_cast<double>(nbatch)
                                   / static_cast<double>(reps_per_batch);
    Print() << "benchmark: np=" << np
            << " variant=" << name
            << " reps=" << reps_per_batch << "x" << nbatch
            << " best=" << best_eval / npd * 1.e9 << " ns/particle"
            << " (" << npd / best_eval / 1.e6 << " Mp/s)"
            << " mean=" << mean_eval / npd * 1.e9 << " ns/particle\n";
}

void run_benchmark (Geometry const& geom, DistributionMapping const& dm,
                    BoxArray const& ba, Shifts const& shifts)
{
    ParmParse pp("benchmark");

    int enabled = 1;
    pp.query("enabled", enabled);
    if (enabled == 0) { return; }

    std::vector<int> sizes{1000};
    if (pp.contains("sizes")) {
        sizes.clear();
        pp.queryarr("sizes", sizes);
    }
    double target_evals = 2.e6;  // total particle evaluations per measurement
    pp.query("target_evals", target_evals);
    int warmup_reps = 3;
    pp.query("warmup_reps", warmup_reps);
    const int nbatch = 5;

    Print() << "\nbenchmark setup: SIMD width=" << simd::native_simd_size_particlereal
            << " sizeof(ParticleReal)=" << sizeof(ParticleReal)
            << " OMP threads=" << OpenMP::get_max_threads()
            << " MPI ranks=" << ParallelDescriptor::NProcs() << "\n";

    double checksum = 0.;
    for (int np : sizes) {
        PC pc(geom, dm, ba);
        init_particles(pc, np);

        const long nrep = std::max(5L,
            static_cast<long>(target_evals / static_cast<double>(std::max(np, 1))));
        const long reps_per_batch = std::max(1L, nrep / nbatch);

        bench_variant("A(SuperParticle)", np, reps_per_batch, nbatch, warmup_reps,
                      [&] () { return run_variant_a(pc, shifts); }, checksum);
        bench_variant("B(ptd,i)        ", np, reps_per_batch, nbatch, warmup_reps,
                      [&] () { return run_variant_b(pc, shifts); }, checksum);
        bench_variant("C(SIMD)         ", np, reps_per_batch, nbatch, warmup_reps,
                      [&] () { return run_variant_c(pc, shifts); }, checksum);
        bench_variant("C2(SIMD 2x)     ", np, reps_per_batch, nbatch, warmup_reps,
                      [&] () { return run_variant_c2(pc, shifts); }, checksum);
    }
    Print() << "benchmark checksum: " << checksum << "\n";
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int is_per[AMREX_SPACEDIM];
        for (int& d : is_per) { d = 1; }

        RealBox real_box;
        for (int n = 0; n < AMREX_SPACEDIM; n++) {
            real_box.setLo(n, -2.0);
            real_box.setHi(n,  2.0);
        }

        IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
        IntVect domain_hi(AMREX_D_DECL(63, 63, 63));
        const Box domain(domain_lo, domain_hi);

        Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);
        BoxArray ba(domain);
        ba.maxSize(64);  // one box, one tile
        DistributionMapping dm(ba);

        const Shifts shifts{
            ParticleReal(0.1), ParticleReal(-0.2), ParticleReal(0.05),
            ParticleReal(0.3), ParticleReal(-0.15), ParticleReal(0.25),
            ParticleReal(0.01), ParticleReal(-0.02), ParticleReal(0.03)};

        correctness_tests(geom, dm, ba, shifts);
        Print() << "\nAll ParticleReduceSIMD tests PASSED.\n";

        run_benchmark(geom, dm, ba, shifts);
    }
    amrex::Finalize();
}
