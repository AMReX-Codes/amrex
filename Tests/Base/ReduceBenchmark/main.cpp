// Benchmark for ReduceOps reductions that return several values per item.
//
// Times one reduction (ReduceOps::eval plus ReduceData::value) for tuples of
// 1 to 46 values, over an index range (as used for particles) and over the
// cells of a MultiFab, and verifies each result against a serial reference
// computed on the host. The widest tuple is modeled after the beam moments
// of a particle accelerator code: first and second moments plus the extrema
// of a six-dimensional phase space, reduced in a single pass.
//
// A second part reduces tuples of mixed integer types whose results are
// known exactly, for sizes around the GPU warp and block sizes.
//
// Run with the default inputs, or override on the command line, e.g.
//   ./main3d.gnu.CUDA.ex inputs tests="sum1 moments46" nitems="100000 10000000"
//   ./main3d.gnu.CUDA.ex inputs nitems="" ncell="64 128 256"

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Reduce.H>
#include <AMReX_TypeList.H>
#include <AMReX_Utility.H>
#include <AMReX_ValLocPair.H>
#include <AMReX_Vector.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>

using namespace amrex;

namespace
{

// Deterministic pseudo-random value in [-0.5, 0.5) for item i of component
// c. The device data and the host reference are both generated from this.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real hashValue (int c, Long i) noexcept
{
    auto h = static_cast<std::uint32_t>(i) * 2654435761U
           + static_cast<std::uint32_t>(c) * 40503U;
    h ^= h >> 15;
    h *= 2246822519U;
    h ^= h >> 13;
    h *= 3266489917U;
    h ^= h >> 16;
    return static_cast<Real>(h >> 8) * Real(1.0/16777216.0) - Real(0.5);
}

// The tests. Each one reads ncomp components of an item through the getter g
// and returns nvalues values, one per reduce operation.

struct Sum1
{
    static constexpr int ncomp = 1;
    static constexpr int nvalues = 1;
    using Ops = TypeList<ReduceOpSum>;
    using Types = TypeList<Real>;

    template <typename G>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static auto eval (G const& g) noexcept
    {
        return amrex::makeTuple(g(0));
    }
};

struct MinMax2
{
    static constexpr int ncomp = 1;
    static constexpr int nvalues = 2;
    using Ops = TypeList<ReduceOpMin, ReduceOpMax>;
    using Types = TypeList<Real, Real>;

    template <typename G>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static auto eval (G const& g) noexcept
    {
        const Real x = g(0);
        return amrex::makeTuple(x, x);
    }
};

// Weighted first and second moments of three coordinates.
struct Sum8
{
    static constexpr int ncomp = 4;
    static constexpr int nvalues = 8;
    using Ops = TypeMultiplier<TypeList, ReduceOpSum[8]>;
    using Types = TypeMultiplier<TypeList, Real[8]>;

    template <typename G>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static auto eval (G const& g) noexcept
    {
        const Real w = Real(1.0) + g(0);
        const Real x = g(1);
        const Real y = g(2);
        const Real z = g(3);
        return amrex::makeTuple(w, w*x, w*y, w*z, w*x*x, w*y*y, w*z*z, w*x*y);
    }
};

// Weight, first and second moments, and extrema of a six-dimensional phase
// space with three more coordinates (e.g., a spin): 34 sums, 6 minima and 6
// maxima.
struct Moments46
{
    static constexpr int ncomp = 10;
    static constexpr int nvalues = 46;
    using Ops = TypeMultiplier<TypeList, ReduceOpSum[34], ReduceOpMin[6], ReduceOpMax[6]>;
    using Types = TypeMultiplier<TypeList, Real[46]>;

    template <typename G>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static auto eval (G const& g) noexcept
    {
        const Real w  = Real(1.0) + g(0);
        const Real x  = g(1);
        const Real y  = g(2);
        const Real t  = g(3);
        const Real px = g(4);
        const Real py = g(5);
        const Real pt = g(6);
        const Real sx = g(7);
        const Real sy = g(8);
        const Real sz = g(9);
        return amrex::makeTuple(
            w,
            x*w, y*w, t*w, px*w, py*w, pt*w, sx*w, sy*w, sz*w,
            x*x*w, y*y*w, t*t*w, px*px*w, py*py*w, pt*pt*w,
            x*px*w, y*py*w, t*pt*w,
            x*pt*w, px*pt*w, y*pt*w, py*pt*w,
            x*y*w, x*py*w, x*t*w, px*y*w, px*py*w, px*t*w, y*t*w, py*t*w,
            sx*sx*w, sy*sy*w, sz*sz*w,
            x, y, t, px, py, pt,
            x, y, t, px, py, pt);
    }
};

// Getters of the components of item i: from the device arrays of a range
// reduction, and generated on the fly for the serial host reference.
template <int ncomp>
struct ArrayGetter
{
    GpuArray<Real const*, ncomp> ptrs;
    int i;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (int c) const noexcept { return ptrs[c][i]; }
};

struct HashGetter
{
    Long i;

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real operator() (int c) const noexcept { return hashValue(c, i); }
};

enum struct Kind { sum, min, max };

template <typename Op>
constexpr Kind kindOf () noexcept
{
    return std::is_same_v<Op, ReduceOpMin> ? Kind::min
        : (std::is_same_v<Op, ReduceOpMax> ? Kind::max : Kind::sum);
}

template <typename... Ops>
constexpr std::array<Kind, sizeof...(Ops)> kindsOf (TypeList<Ops...>) noexcept
{
    return {kindOf<Ops>()...};
}

template <typename... Ts>
std::array<double, sizeof...(Ts)> toArray (GpuTuple<Ts...> const& t)
{
    return [&] <std::size_t... Is> (std::index_sequence<Is...>)
    {
        return std::array<double, sizeof...(Ts)>{static_cast<double>(amrex::get<Is>(t))...};
    }(std::index_sequence_for<Ts...>{});
}

// Serial host reference for a test over nitems items, compared with result.
// Minima and maxima must match exactly. Sums are accumulated in double
// precision and may differ by the rounding of a sum of nitems terms.
template <typename Test>
bool verify (Long nitems, std::array<double, Test::nvalues> const& result)
{
    constexpr auto kinds = kindsOf(typename Test::Ops{});
    constexpr int nvalues = Test::nvalues;
    std::array<double, nvalues> ref{};
    std::array<double, nvalues> mag{};
    for (int n = 0; n < nvalues; ++n) {
        if (kinds[n] == Kind::min) {
            ref[n] = std::numeric_limits<double>::max();
        } else if (kinds[n] == Kind::max) {
            ref[n] = std::numeric_limits<double>::lowest();
        }
    }

    for (Long i = 0; i < nitems; ++i) {
        const auto v = toArray(Test::eval(HashGetter{i}));
        for (int n = 0; n < nvalues; ++n) {
            if (kinds[n] == Kind::min) {
                ref[n] = std::min(ref[n], v[n]);
            } else if (kinds[n] == Kind::max) {
                ref[n] = std::max(ref[n], v[n]);
            } else {
                ref[n] += v[n];
                mag[n] += std::abs(v[n]);
            }
        }
    }

    const double tol = 32.0 * static_cast<double>(std::numeric_limits<Real>::epsilon())
        * std::sqrt(static_cast<double>(nitems));
    bool correct = true;
    for (int n = 0; n < nvalues; ++n) {
        if (kinds[n] == Kind::sum) {
            correct = correct && (std::abs(result[n] - ref[n]) <= tol * mag[n]);
        } else {
            correct = correct && (result[n] == ref[n]);
        }
    }
    return correct;
}

struct Timing
{
    double tmin = std::numeric_limits<double>::max();
    double tavg = 0.0;
    bool correct = true;
};

template <typename F>
Timing timeIt (int nwarmup, int nrepeat, F const& f)
{
    Timing timing;
    for (int irep = -nwarmup; irep < nrepeat; ++irep) {
        Gpu::streamSynchronize();
        const double t0 = amrex::second();
        f();
        const double dt = amrex::second() - t0;
        if (irep >= 0) {
            timing.tmin = std::min(timing.tmin, dt);
            timing.tavg += dt;
        }
    }
    timing.tavg /= nrepeat;
    return timing;
}

void report (std::string const& name, int nvalues, std::string const& layout, Long nitems,
             Timing const& timing)
{
    amrex::Print() << std::left << std::setw(12) << name << std::right
                   << std::setw(8) << nvalues
                   << std::setw(8) << layout
                   << std::setw(12) << nitems
                   << std::fixed << std::setprecision(4)
                   << std::setw(11) << timing.tmin * 1.e3
                   << std::setw(11) << timing.tavg * 1.e3
                   << std::setprecision(3)
                   << std::setw(10) << timing.tmin * 1.e9 / static_cast<double>(nitems)
                   << std::setw(8) << (timing.correct ? "ok" : "WRONG") << "\n";
}

// Reduce over an index range, reading the components from device vectors.
template <typename Test>
bool benchmarkRange (std::string const& name, int nitems, int nwarmup, int nrepeat)
{
    constexpr int ncomp = Test::ncomp;
    Array<Gpu::DeviceVector<Real>, ncomp> data;
    GpuArray<Real const*, ncomp> ptrs{};
    for (int c = 0; c < ncomp; ++c) {
        data[c].resize(nitems);
        Real* p = data[c].data();
        amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            p[i] = hashValue(c, i);
        });
        ptrs[c] = p;
    }

    auto f = [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        return Test::eval(ArrayGetter<ncomp>{ptrs, i});
    };
    auto reduce = [&] ()
    {
        Reducer<typename Test::Ops, typename Test::Types> reducer;
        reducer.eval(nitems, f);
        return reducer.getResult();
    };

    auto result = reduce();
    Timing timing = timeIt(nwarmup, nrepeat, [&] () { result = reduce(); });
    timing.correct = verify<Test>(nitems, toArray(result));
    report(name, Test::nvalues, "range", nitems, timing);
    return timing.correct;
}

// Reduce over the cells of a MultiFab on a cubic domain of ncell cells per
// direction.
template <typename Test>
bool benchmarkMultiFab (std::string const& name, int ncell, int max_grid_size,
                        int nwarmup, int nrepeat)
{
    constexpr int ncomp = Test::ncomp;
    const Box domain(IntVect(0), IntVect(ncell-1));
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    MultiFab mf(ba, DistributionMapping{ba}, ncomp, 0);

    const auto len = domain.length3d();
    const Long nx = len[0];
    const Long ny = len[1];
    auto const& ma = mf.arrays();
    amrex::ParallelFor(mf, IntVect(0), ncomp,
    [=] AMREX_GPU_DEVICE (int b, int i, int j, int k, int n) noexcept
    {
        ma[b](i,j,k,n) = hashValue(n, i + nx*(j + ny*k));
    });
    Gpu::streamSynchronize();

    auto const& cma = mf.const_arrays();
    auto f = [=] AMREX_GPU_DEVICE (int b, int i, int j, int k) noexcept
    {
        return Test::eval([=] (int c) noexcept { return cma[b](i,j,k,c); });
    };
    auto reduce = [&] ()
    {
        Reducer<typename Test::Ops, typename Test::Types> reducer;
        reducer.eval(mf, IntVect(0), f);
        return reducer.getResult();
    };

    auto result = reduce();
    Timing timing = timeIt(nwarmup, nrepeat, [&] () { result = reduce(); });
    timing.correct = verify<Test>(domain.numPts(), toArray(result));
    report(name, Test::nvalues, "cells", domain.numPts(), timing);
    return timing.correct;
}

// Reductions of mixed integer types with exactly known results.
using IntLoc = ValLocPair<Real, int>;
using LongLoc = ValLocPair<Real, Long>;
using MixedOps = TypeList<ReduceOpSum, ReduceOpSum, ReduceOpSum,
                          ReduceOpMin, ReduceOpMax, ReduceOpMin, ReduceOpMax,
                          ReduceOpLogicalAnd, ReduceOpLogicalAnd,
                          ReduceOpLogicalOr, ReduceOpLogicalOr,
                          ReduceOpMin, ReduceOpMax>;
using MixedTypes = TypeList<Long, Long, int,
                            int, int, Long, Real,
                            int, int,
                            int, int,
                            IntLoc, LongLoc>;
using MixedTuple = GpuTuple<Long, Long, int, int, int, Long, Real, int, int, int, int,
                            IntLoc, LongLoc>;

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
MixedTuple mixedValues (int i, int nitems) noexcept
{
    const auto scrambled = static_cast<int>((Long(i) * 7919) % 10007);
    return {Long(1), Long(i % 7), i % 3,
            scrambled, scrambled, Long(scrambled) - 5000, hashValue(1, i),
            int(i >= 0), int(i != nitems/2),
            int(i == nitems-1), int(i < 0),
            IntLoc{hashValue(2, i), i}, LongLoc{hashValue(3, i), Long(i)}};
}

bool mixedIsCorrect (int nitems, MixedTuple const& result)
{
    Long sum7 = 0;
    int sum3 = 0;
    int smin = std::numeric_limits<int>::max();
    int smax = std::numeric_limits<int>::lowest();
    Real hmax = std::numeric_limits<Real>::lowest();
    Real lmin = std::numeric_limits<Real>::max();
    Real lmax = std::numeric_limits<Real>::lowest();
    for (int i = 0; i < nitems; ++i) {
        const auto scrambled = static_cast<int>((Long(i) * 7919) % 10007);
        sum7 += i % 7;
        sum3 += i % 3;
        smin = std::min(smin, scrambled);
        smax = std::max(smax, scrambled);
        hmax = std::max(hmax, hashValue(1, i));
        lmin = std::min(lmin, hashValue(2, i));
        lmax = std::max(lmax, hashValue(3, i));
    }
    const auto& minloc = amrex::get<11>(result);
    const auto& maxloc = amrex::get<12>(result);
    // several items can share the extreme value: any of their indices is valid
    const bool minloc_ok = (minloc.value == lmin) && (minloc.index >= 0)
        && (minloc.index < nitems) && (hashValue(2, minloc.index) == lmin);
    const bool maxloc_ok = (maxloc.value == lmax) && (maxloc.index >= 0)
        && (maxloc.index < nitems) && (hashValue(3, maxloc.index) == lmax);
    return (amrex::get<0>(result) == Long(nitems))
        && (amrex::get<1>(result) == sum7)
        && (amrex::get<2>(result) == sum3)
        && (amrex::get<3>(result) == smin)
        && (amrex::get<4>(result) == smax)
        && (amrex::get<5>(result) == Long(smin) - 5000)
        && (amrex::get<6>(result) == hmax)
        && (amrex::get<7>(result) == 1)
        && (amrex::get<8>(result) == 0)
        && (amrex::get<9>(result) == 1)
        && (amrex::get<10>(result) == 0)
        && minloc_ok && maxloc_ok;
}

// One eval over the whole range, and three evals of unequal parts that
// accumulate into the same Reducer.
bool checkMixed (int nitems)
{
    Reducer<MixedOps, MixedTypes> whole;
    whole.eval(nitems, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        return mixedValues(i, nitems);
    });

    const Array<int,4> bounds{0, nitems/7, nitems/7 + static_cast<int>((Long(5)*nitems)/9),
                              nitems};
    Reducer<MixedOps, MixedTypes> parts;
    for (int ipart = 0; ipart < 3; ++ipart) {
        const int offset = bounds[ipart];
        parts.eval(bounds[ipart+1] - offset, [=] AMREX_GPU_DEVICE (int i) noexcept
        {
            return mixedValues(i + offset, nitems);
        });
    }

    return mixedIsCorrect(nitems, whole.getResult())
        && mixedIsCorrect(nitems, parts.getResult());
}

// The same kind of exact check for the cells of a MultiFab.
bool checkMixedMultiFab (int ncell, int max_grid_size)
{
    const Box domain(IntVect(0), IntVect(ncell-1));
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);
    const MultiFab mf(ba, DistributionMapping{ba}, 1, 0);

    const auto len = domain.length3d();
    const int nx = len[0];
    const int ny = len[1];
    Reducer<TypeList<ReduceOpSum, ReduceOpSum, ReduceOpMin, ReduceOpMax>,
            TypeList<Long, Long, int, int>> reducer;
    using Result_t = typename decltype(reducer)::Result_t;
    reducer.eval(mf, IntVect(0), [=] AMREX_GPU_DEVICE (int, int i, int j, int k) noexcept
                                 -> Result_t
    {
        const int cell = i + nx*(j + ny*k);
        return {Long(1), Long(cell % 7), cell, cell};
    });
    const auto result = reducer.getResult();

    const Long ncells = domain.numPts();
    Long sum7 = 0;
    for (Long cell = 0; cell < ncells; ++cell) { sum7 += cell % 7; }
    return (amrex::get<0>(result) == ncells)
        && (amrex::get<1>(result) == sum7)
        && (amrex::get<2>(result) == 0)
        && (amrex::get<3>(result) == static_cast<int>(ncells-1));
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        Vector<int> nitems = {1000, 100000, 1000000};
        Vector<int> ncell = {32, 64};
        int max_grid_size = 32;
        int nwarmup = 2;
        int nrepeat = 20;
        ParmParse pp;
        pp.queryarr("nitems", nitems);
        pp.queryarr("ncell", ncell);
        pp.query("max_grid_size", max_grid_size);
        pp.query("nwarmup", nwarmup);
        pp.query("nrepeat", nrepeat);
        Vector<std::string> test_list;
        pp.queryarr("tests", test_list);
        if (test_list.empty()) {
            test_list = {"sum1", "minmax2", "sum8", "moments46"};
        }

        amrex::Print() << "ReduceOps benchmark, " << nrepeat << " repetitions\n"
                       << std::left << std::setw(12) << "test" << std::right
                       << std::setw(8) << "nvalues" << std::setw(8) << "over"
                       << std::setw(12) << "nitems" << std::setw(11) << "min [ms]"
                       << std::setw(11) << "avg [ms]" << std::setw(10) << "ns/item"
                       << std::setw(8) << "check" << "\n";

        bool all_correct = true;
        for (const auto& test : test_list) {
            for (const int n : nitems) {
                bool correct = false;
                if (test == "sum1") {
                    correct = benchmarkRange<Sum1>(test, n, nwarmup, nrepeat);
                } else if (test == "minmax2") {
                    correct = benchmarkRange<MinMax2>(test, n, nwarmup, nrepeat);
                } else if (test == "sum8") {
                    correct = benchmarkRange<Sum8>(test, n, nwarmup, nrepeat);
                } else if (test == "moments46") {
                    correct = benchmarkRange<Moments46>(test, n, nwarmup, nrepeat);
                } else {
                    amrex::Abort("Unknown test " + test);
                }
                all_correct = all_correct && correct;
            }
            for (const int n : ncell) {
                bool correct = false;
                if (test == "sum1") {
                    correct = benchmarkMultiFab<Sum1>(test, n, max_grid_size, nwarmup, nrepeat);
                } else if (test == "minmax2") {
                    correct = benchmarkMultiFab<MinMax2>(test, n, max_grid_size, nwarmup, nrepeat);
                } else if (test == "sum8") {
                    correct = benchmarkMultiFab<Sum8>(test, n, max_grid_size, nwarmup, nrepeat);
                } else if (test == "moments46") {
                    correct = benchmarkMultiFab<Moments46>(test, n, max_grid_size, nwarmup, nrepeat);
                } else {
                    amrex::Abort("Unknown test " + test);
                }
                all_correct = all_correct && correct;
            }
        }
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(all_correct, "ReduceOps produced a wrong result");

        // Exact results, for sizes around the GPU warp and block sizes too.
        Vector<int> exact_sizes = {1, 2, 31, 32, 33, 63, 64, 65, 255, 256, 257, 1023, 1024,
                                   1025, 4099, 65537};
        exact_sizes.insert(exact_sizes.end(), nitems.begin(), nitems.end());
        int nwrong = 0;
        for (const int n : exact_sizes) {
            if (! checkMixed(n)) {
                ++nwrong;
                amrex::Print() << "mixed types over a range of " << n << " items: WRONG\n";
            }
        }
        for (const int n : ncell) {
            for (const int mgs : {max_grid_size, n}) {
                if (! checkMixedMultiFab(n, mgs)) {
                    ++nwrong;
                    amrex::Print() << "mixed types over cells, ncell = " << n << ", max_grid_size = "
                                   << mgs << ": WRONG\n";
                }
            }
        }
        amrex::Print() << "mixed types, exact results: "
                       << (nwrong == 0 ? "ok" : "WRONG") << "\n";
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(nwrong == 0, "ReduceOps produced a wrong result");
    }
    amrex::Finalize();
}
