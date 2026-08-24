// Test coverage for AMReX SIMD single-source design:
//   SIMDindex, ParallelForSIMD, load_1d, store_1d,
//   Vectorized, is_vectorized, is_nth_arg_non_const
//
// Compiles and runs with AMReX_SIMD=OFF (scalar fallback, GPU) and ON (CPU SIMD).

#include <AMReX.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>
#include <AMReX_ReduceSIMD.H>
#include <AMReX_SIMD.H>
#include <AMReX_Vector.H>

#include <AMReX_Algorithm.H>
#include <AMReX_Math.H>

#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <type_traits>

#include "AMReX_GpuContainers.H"

using namespace amrex;

// ---------------------------------------------------------------------------
// Helper functors / functions used by the tests
// ---------------------------------------------------------------------------

// Functor that does NOT support SIMD (no Vectorized mixin)
struct ScalarCompute
{
    template <typename T_Real>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void
    operator() (T_Real& AMREX_RESTRICT x,
                T_Real const& AMREX_RESTRICT y) const
    {
        x = x + y;
    }
};

// Functor that supports SIMD via the Vectorized mixin
struct VectorizedCompute : public simd::Vectorized<>
{
    template <typename T_Real>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    void
    operator() (T_Real& AMREX_RESTRICT x,
                T_Real const& AMREX_RESTRICT y) const
    {
        x = x + y;
    }
};

// Free functions for is_nth_arg_non_const testing
void func_mc (ParticleReal& x, ParticleReal const& y) { x += y; }
void func_cc (ParticleReal const& /*x*/, ParticleReal const& /*y*/) {}
void func_mm (ParticleReal& x, ParticleReal& y) { x += y; y *= ParticleReal(2); }

#ifdef AMREX_USE_SIMD
// Compare a SIMD math overload against its scalar counterpart over a range.
//
// A vector math library is allowed to be less accurate than scalar libm; glibc's
// libmvec documents a maximum error of 4 ULP, so this uses a tolerance of 8 ULP.
template <typename T_Simd, typename F_Simd, typename F_Scalar>
int check_simd_math (char const* name, F_Simd const& f_simd, F_Scalar const& f_scalar,
                     typename T_Simd::value_type lo, typename T_Simd::value_type hi)
{
    using T = typename T_Simd::value_type;
    constexpr std::size_t width = T_Simd::size();
    constexpr int nchunk = 16;
    constexpr int npoint = nchunk * int(width);
    constexpr T max_ulp = T(8);

    int err = 0;
    for (int c = 0; c < nchunk; ++c) {
        T in[width];
        for (std::size_t i = 0; i < width; ++i) {
            in[i] = lo + (hi - lo) * T(c * int(width) + int(i)) / T(npoint - 1);
        }
        T_Simd x;
        x.copy_from(in, simd::stdx::element_aligned);

        T_Simd const y = f_simd(x);

        for (std::size_t i = 0; i < width; ++i) {
            T const ref = f_scalar(in[i]);
            T const got = y[i];
            T const tol = max_ulp * std::numeric_limits<T>::epsilon()
                          * amrex::max(Math::abs(ref), T(1));
            if (!(Math::abs(got - ref) <= tol)) {
                ++err;
                if (err <= 2) {
                    Print() << "  " << name << " mismatch at x=" << double(in[i])
                            << ": got " << double(got) << ", expected " << double(ref) << "\n";
                }
            }
        }
    }
    if (err != 0) { Print() << "  " << name << ": FAILED (" << err << " lanes)\n"; }
    return err;
}
#endif

// ---------------------------------------------------------------------------
int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int nerrors = 0;

        // ================================================================
        // Test 1: SIMDindex struct
        // ================================================================
        {
            SIMDindex<4, int> si4{42};
            AMREX_ALWAYS_ASSERT(si4.width == 4);
            AMREX_ALWAYS_ASSERT(si4.index == 42);

            SIMDindex<1> si1{7};
            AMREX_ALWAYS_ASSERT(si1.width == 1);
            AMREX_ALWAYS_ASSERT(si1.index == 7);

            SIMDindex<8, long> si8{100L};
            AMREX_ALWAYS_ASSERT(si8.width == 8);
            AMREX_ALWAYS_ASSERT(si8.index == 100L);

            Print() << "SIMDindex: PASSED\n";
        }

        // ================================================================
        // Test 2: is_vectorized trait
        // ================================================================
        {
            static_assert(!simd::is_vectorized<ScalarCompute>,
                          "ScalarCompute must not be vectorized");
            static_assert( simd::is_vectorized<VectorizedCompute>,
                          "VectorizedCompute must be vectorized");
            static_assert(!simd::is_vectorized<int>,
                          "int must not be vectorized");
            static_assert(!simd::is_vectorized<double>,
                          "double must not be vectorized");
            Print() << "is_vectorized: PASSED\n";
        }

        // ================================================================
        // Test 3: is_nth_arg_non_const
        // ================================================================
        {
            // Free functions
            static_assert( simd::is_nth_arg_non_const(&func_mc, 0));
            static_assert(!simd::is_nth_arg_non_const(&func_mc, 1));
            static_assert(!simd::is_nth_arg_non_const(&func_cc, 0));
            static_assert(!simd::is_nth_arg_non_const(&func_cc, 1));
            static_assert( simd::is_nth_arg_non_const(&func_mm, 0));
            static_assert( simd::is_nth_arg_non_const(&func_mm, 1));

            // Functor const member (VectorizedCompute::operator() is const)
            using PR = ParticleReal;
            static_assert( simd::is_nth_arg_non_const(
                &VectorizedCompute::operator()<PR>, 0));
            static_assert(!simd::is_nth_arg_non_const(
                &VectorizedCompute::operator()<PR>, 1));

            Print() << "is_nth_arg_non_const: PASSED\n";
        }

#ifdef AMREX_USE_SIMD
        // ================================================================
        // Test 4: ParallelForSIMD<WIDTH=1> (explicit width 1, scalar)
        // ================================================================
        {
            constexpr int n = 100;
            Vector<ParticleReal> data(n, ParticleReal(0));
            auto* ptr = data.data();

            ParallelForSIMD<1>(n, [=] AMREX_GPU_DEVICE (auto si) {
                ptr[si.index] = static_cast<ParticleReal>(si.index + 1);
            });

            int err = 0;
            for (int i = 0; i < n; ++i) {
                if (data[i] != static_cast<ParticleReal>(i + 1)) { ++err; }
            }
            nerrors += err;
            Print() << "ParallelForSIMD<1>: "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 5: ParallelForSIMD<WIDTH=4> with remainder elements
        // ================================================================
        {
            constexpr int WIDTH = 4;
            constexpr int n = 11;   // 2 full lanes of 4 + 3 remainder
            Vector<ParticleReal> data(n, ParticleReal(0));
            Vector<int> widths(n, 0);
            auto* ptr  = data.data();
            auto* wptr = widths.data();

            ParallelForSIMD<WIDTH>(n, [ptr, wptr] AMREX_GPU_DEVICE (auto si) {
                for (int lane = 0; lane < si.width; ++lane) {
                    ptr[si.index + lane] =
                        static_cast<ParticleReal>(si.index + lane + 1);
                    wptr[si.index + lane] = si.width;
                }
            });

            int err = 0;
            for (int i = 0; i < n; ++i) {
                if (data[i] != static_cast<ParticleReal>(i + 1)) { ++err; }
            }
            // First 8 elements must have seen WIDTH=4, last 3 WIDTH=1
            for (int i = 0; i < 8;  ++i) { if (widths[i] != WIDTH) { ++err; } }
            for (int i = 8; i < n;  ++i) { if (widths[i] != 1)     { ++err; } }

            nerrors += err;
            Print() << "ParallelForSIMD<4> (remainder): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }
#endif // AMREX_USE_SIMD

        // ================================================================
        // Test 7: ParallelForSIMD<T> with non-Vectorized type (fallback)
        // ================================================================
        {
            constexpr int n = 50;
            Gpu::ManagedVector<ParticleReal> data(n, ParticleReal(0));
            auto* ptr = data.data();

            // ScalarCompute is NOT Vectorized → falls back to ParallelFor
            ParallelForSIMD<ScalarCompute>(n, [ptr] AMREX_GPU_DEVICE (int i) {
                ptr[i] = static_cast<ParticleReal>(i * 2);
            });
            amrex::Gpu::streamSynchronize();

            int err = 0;
            for (int i = 0; i < n; ++i) {
                if (data[i] != static_cast<ParticleReal>(i * 2)) { ++err; }
            }
            nerrors += err;
            Print() << "ParallelForSIMD<ScalarFunctor> (fallback): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 8: load_1d with plain int index
        // ================================================================
        {
            constexpr int n = 10;
            Vector<ParticleReal> data(n);
            std::iota(data.begin(), data.end(), ParticleReal(1));
            auto* ptr = data.data();

            int err = 0;
            for (int i = 0; i < n; ++i) {
                decltype(auto) val = simd::load_1d(ptr, i);
                if (val != static_cast<ParticleReal>(i + 1)) { ++err; }
            }
            nerrors += err;
            Print() << "load_1d (int index): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 9: load_1d with SIMDindex<1>
        // ================================================================
        {
            constexpr int n = 10;
            Vector<ParticleReal> data(n);
            std::iota(data.begin(), data.end(), ParticleReal(1));
            auto* ptr = data.data();

            int err = 0;
            for (int i = 0; i < n; ++i) {
                SIMDindex<1> si{i};
                decltype(auto) val = simd::load_1d(ptr, si);
                if (val != static_cast<ParticleReal>(i + 1)) { ++err; }
            }
            nerrors += err;
            Print() << "load_1d (SIMDindex<1>): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 10: store_1d no-op when ValType == T (scalar path)
        // ================================================================
        {
            constexpr int n = 10;
            Vector<ParticleReal> data(n, ParticleReal(42));
            auto* ptr = data.data();
            auto val = ParticleReal(99);

            // ValType (ParticleReal) == T (ParticleReal) → must be a no-op
            simd::store_1d<&func_mc, 0>(val, ptr, 0);

            int err = 0;
            for (int i = 0; i < n; ++i) {
                if (data[i] != ParticleReal(42)) { ++err; }
            }
            nerrors += err;
            Print() << "store_1d (same type, no-op): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 11: Combined single-source pattern
        //   ParallelForSIMD<VectorizedCompute> + load_1d + store_1d
        //   Without SIMD: falls back to ParallelFor (scalar).
        //   With SIMD: uses native SIMD lanes.
        // ================================================================
        {
            constexpr int n = 64;
            Gpu::ManagedVector<ParticleReal> x_data(n);
            Gpu::ManagedVector<ParticleReal> y_data(n);
            std::iota(x_data.begin(), x_data.end(), ParticleReal(0));
            std::ranges::fill(y_data, ParticleReal(100));

            auto* x_ptr = x_data.data();
            auto* y_ptr = y_data.data();

            VectorizedCompute vc;

            ParallelForSIMD<VectorizedCompute>(n,
                [=] AMREX_GPU_DEVICE (auto i)
            {
                decltype(auto) x = simd::load_1d(x_ptr, i);
                decltype(auto) y = simd::load_1d(y_ptr, i);

                using Val = std::remove_reference_t<decltype(x)>;
                constexpr auto method = &VectorizedCompute::operator()<Val>;

                vc(x, y);   // x += y

                simd::store_1d<method, 0>(x, x_ptr, i);
                simd::store_1d<method, 1>(y, y_ptr, i);
            });
            Gpu::streamSynchronize();

            int err = 0;
            for (int i = 0; i < n; ++i) {
                auto expected = static_cast<ParticleReal>(i) + ParticleReal(100);
                if (std::abs(x_data[i] - expected) > ParticleReal(1.e-10)) {
                    ++err;
                    Print() << "  x[" << i << "] = " << x_data[i]
                            << " expected " << expected << "\n";
                }
                // y must be unchanged (operator() arg 1 is const)
                if (y_data[i] != ParticleReal(100)) {
                    ++err;
                    Print() << "  y[" << i << "] = " << y_data[i]
                            << " expected 100\n";
                }
            }
            nerrors += err;
            Print() << "Combined ParallelForSIMD+load_1d+store_1d: "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 12 (SIMD-only): full SIMD path with native width
        //   ParallelForSIMD<VectorizedCompute> dispatches to the native
        //   SIMD width, load_1d fills SIMD registers, store_1d writes back.
        // ================================================================
        {
            constexpr int WIDTH = simd::native_simd_size_particlereal;
            constexpr int n = WIDTH * 4 + 3;  // ensure remainder
            Print() << "  (native SIMD width: " << WIDTH << ", n=" << n << ")\n";

            Gpu::ManagedVector<ParticleReal> x_data(n);
            Gpu::ManagedVector<ParticleReal> y_data(n);
            std::iota(x_data.begin(), x_data.end(), ParticleReal(0));
            std::ranges::fill(y_data, ParticleReal(10));

            auto* x_ptr = x_data.data();
            auto* y_ptr = y_data.data();

            VectorizedCompute vc;

            ParallelForSIMD<VectorizedCompute>(n,
                [=] AMREX_GPU_DEVICE (auto si)
            {
                decltype(auto) x = simd::load_1d(x_ptr, si);
                decltype(auto) y = simd::load_1d(y_ptr, si);

                using Val = std::remove_reference_t<decltype(x)>;
                constexpr auto method = &VectorizedCompute::operator()<Val>;

                vc(x, y);

                simd::store_1d<method, 0>(x, x_ptr, si);
                simd::store_1d<method, 1>(y, y_ptr, si);
            });
            Gpu::streamSynchronize();

            int err = 0;
            for (int i = 0; i < n; ++i) {
                auto expected =
                    static_cast<ParticleReal>(i) + ParticleReal(10);
                if (std::abs(x_data[i] - expected) > ParticleReal(1.e-10)) {
                    ++err;
                    Print() << "  SIMD x[" << i << "] = " << x_data[i]
                            << " expected " << expected << "\n";
                }
                if (y_data[i] != ParticleReal(10)) {
                    ++err;
                    Print() << "  SIMD y[" << i << "] = " << y_data[i]
                            << " expected 10\n";
                }
            }
            nerrors += err;
            Print() << "SIMD-path combined test: "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 13: Generic load_1d correctness with position-dependent data
        //   Verifies load_1d reads the correct array positions by using
        //   position-dependent y values: y[i] = n - i, so x[i] + y[i]
        //   always equals n + 1 regardless of SIMD lane assignment.
        //   A prime n stresses SIMD remainder handling.
        //   Works on GPU, CPU without SIMD, and CPU with SIMD.
        // ================================================================
        {
            constexpr int n = 67;  // prime, stresses remainder handling

            Gpu::ManagedVector<ParticleReal> x_data(n);
            Gpu::ManagedVector<ParticleReal> y_data(n);
            for (int i = 0; i < n; ++i) {
                x_data[i] = static_cast<ParticleReal>(i + 1);  // x[i] = i+1
                y_data[i] = static_cast<ParticleReal>(n - i);  // y[i] = n-i
            }

            auto* x_ptr = x_data.data();
            auto* y_ptr = y_data.data();

            VectorizedCompute vc;

            ParallelForSIMD<VectorizedCompute>(n,
                [=] AMREX_GPU_DEVICE (auto i)
            {
                decltype(auto) x = simd::load_1d(x_ptr, i);
                decltype(auto) y = simd::load_1d(y_ptr, i);

                using Val = std::remove_reference_t<decltype(x)>;
                constexpr auto method = &VectorizedCompute::operator()<Val>;

                vc(x, y);   // x += y

                simd::store_1d<method, 0>(x, x_ptr, i);
                simd::store_1d<method, 1>(y, y_ptr, i);
            });
            Gpu::streamSynchronize();

            int err = 0;
            for (int i = 0; i < n; ++i) {
                // x[i] = (i+1) + (n-i) = n+1 = 68 for all i
                auto expected_x = static_cast<ParticleReal>(n + 1);
                if (std::abs(x_data[i] - expected_x) > ParticleReal(1.e-10)) {
                    ++err;
                    Print() << "  x[" << i << "] = " << x_data[i]
                            << " expected " << expected_x << "\n";
                }
                // y must be unchanged (operator() arg 1 is const)
                auto expected_y = static_cast<ParticleReal>(n - i);
                if (y_data[i] != expected_y) {
                    ++err;
                    Print() << "  y[" << i << "] = " << y_data[i]
                            << " expected " << expected_y << "\n";
                }
            }
            nerrors += err;
            Print() << "Position-dependent load_1d+store_1d: "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

        // ================================================================
        // Test 14: any_of, where, select — portable single-source
        //   Uses SIMDParticleReal<>, which is a SIMD vector when AMREX_USE_SIMD=ON
        //   and a plain scalar when OFF.  The same code path exercises
        //   both the real SIMD and the scalar fallback implementations.
        // ================================================================
        {
            using PReal_t = simd::SIMDParticleReal<>;

            // safe reciprocal: 1/b where b != 0, else 0
            auto b = PReal_t(2);
            auto mask = b != PReal_t(0);
            auto safe_b = simd::stdx::select(mask, b, PReal_t(1));
            auto recip  = simd::stdx::select(mask,
                              PReal_t(1) / safe_b,
                              PReal_t(0));

            // any_of: at least one lane should be nonzero
            AMREX_ALWAYS_ASSERT(simd::stdx::any_of(mask));

            // where: masked assignment
            auto acc = PReal_t(0);
            simd::stdx::where(mask, acc) = recip;

            // verify: b=2 everywhere → recip=0.5, acc=0.5
            int err = 0;
            auto check = [&] (ParticleReal got, ParticleReal expected) {
                if (std::abs(got - expected) > ParticleReal(1.e-10)) { ++err; }
            };
#ifdef AMREX_USE_SIMD
            for (std::size_t lane = 0; lane < PReal_t::size(); ++lane) {
                check(recip[lane], ParticleReal(0.5));
                check(acc[lane],   ParticleReal(0.5));
            }
#else
            check(recip, ParticleReal(0.5));
            check(acc,   ParticleReal(0.5));
#endif

            nerrors += err;
            Print() << "any_of + where + select (portable): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }

#ifndef AMREX_USE_GPU
        // ================================================================
        // Test 15: evalReduceSIMD with a mixed-type tuple
        //   Sum(Real), Min(int), Max(Real) over the same index range.
        //   With AMReX_SIMD=ON the main loop runs vectorized (native width);
        //   without, the same code runs the scalar loop.
        // ================================================================
        {
            constexpr int n = 67;  // prime, stresses SIMD remainder handling
            Vector<Real> rdata(n);
            Vector<int> idata(n);
            for (int i = 0; i < n; ++i) {
                rdata[i] = static_cast<Real>(i + 1);  // 1 .. n
                idata[i] = (i * 7) % 23 + 1;          // in [1, 23]
            }
            auto const* rp = rdata.data();
            auto const* ip = idata.data();

            ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_ops;
            ReduceData<Real, int, Real> reduce_data(reduce_ops);
            constexpr int W = static_cast<int>(simd::native_simd_size_real);

            evalReduceSIMD<W>(n, reduce_data,
                [=] (auto si)
                {
                    auto const a = simd::load_1d(rp, si);
                    auto const b = simd::load_1d(ip, si);
                    return amrex::makeTuple(a, b, a);
                }, reduce_ops);
            auto const r = reduce_data.value(reduce_ops);

            int err = 0;
            constexpr int expected_sum = n*(n+1)/2;  // sum 1..n
            if (amrex::get<0>(r) != static_cast<Real>(expected_sum)) { ++err; }
            if (amrex::get<1>(r) != 1) { ++err; }
            if (amrex::get<2>(r) != static_cast<Real>(n)) { ++err; }

            nerrors += err;
            Print() << "evalReduceSIMD (mixed-type tuple): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
        }
#endif // !AMREX_USE_GPU

        // ================================================================
        // Test: amrex::Math transcendentals, scalar and SIMD
        // ================================================================
        {
            int err = 0;

            // The scalar overloads must exist for every build, so that a kernel
            // written once compiles with AMReX_SIMD both ON and OFF.
            {
                constexpr Real x = Real(0.5);
                if (!amrex::almostEqual(Math::sinh(x), std::sinh(x))) { ++err; }
                if (!amrex::almostEqual(Math::cosh(x), std::cosh(x))) { ++err; }
                if (!amrex::almostEqual(Math::exp(x), std::exp(x))) { ++err; }
                if (!amrex::almostEqual(Math::sqrt(x), std::sqrt(x))) { ++err; }
                if (!amrex::almostEqual(Math::pow(x, Real(3)), std::pow(x, Real(3)))) { ++err; }
                auto const [s, c] = Math::sincos(x);
                if (!amrex::almostEqual(s, std::sin(x))) { ++err; }
                if (!amrex::almostEqual(c, std::cos(x))) { ++err; }
            }

#ifdef AMREX_USE_SIMD
            using V = simd::SIMDReal<>;
            using T = Real;

#   define AMREX_CHECK_SIMD_MATH(FUNC, LO, HI)                                \
        err += check_simd_math<V>(#FUNC,                                      \
                    [] (V const& v) { return Math::FUNC(v); },                \
                    [] (T const v)  { return std::FUNC(v); },                 \
                    T(LO), T(HI))

            AMREX_CHECK_SIMD_MATH(sin,   -6.0,  6.0);
            AMREX_CHECK_SIMD_MATH(cos,   -6.0,  6.0);
            AMREX_CHECK_SIMD_MATH(tan,   -1.5,  1.5);
            AMREX_CHECK_SIMD_MATH(asin,  -1.0,  1.0);
            AMREX_CHECK_SIMD_MATH(acos,  -1.0,  1.0);
            AMREX_CHECK_SIMD_MATH(atan, -10.0, 10.0);
            AMREX_CHECK_SIMD_MATH(sinh,  -5.0,  5.0);
            AMREX_CHECK_SIMD_MATH(cosh,  -5.0,  5.0);
            AMREX_CHECK_SIMD_MATH(tanh,  -5.0,  5.0);
            AMREX_CHECK_SIMD_MATH(asinh, -5.0,  5.0);
            AMREX_CHECK_SIMD_MATH(acosh,  1.0, 10.0);
            AMREX_CHECK_SIMD_MATH(atanh, -0.9,  0.9);
            AMREX_CHECK_SIMD_MATH(exp,   -5.0,  5.0);
            AMREX_CHECK_SIMD_MATH(exp2,  -5.0,  5.0);
            AMREX_CHECK_SIMD_MATH(expm1, -1.0,  1.0);
            AMREX_CHECK_SIMD_MATH(log,    0.1, 20.0);
            AMREX_CHECK_SIMD_MATH(log2,   0.1, 20.0);
            AMREX_CHECK_SIMD_MATH(log10,  0.1, 20.0);
            AMREX_CHECK_SIMD_MATH(log1p, -0.9,  9.0);
            AMREX_CHECK_SIMD_MATH(cbrt, -20.0, 20.0);
            AMREX_CHECK_SIMD_MATH(erf,   -3.0,  3.0);
            AMREX_CHECK_SIMD_MATH(erfc,  -3.0,  3.0);
            AMREX_CHECK_SIMD_MATH(sqrt,   0.0, 20.0);
            AMREX_CHECK_SIMD_MATH(abs,  -20.0, 20.0);

#   undef AMREX_CHECK_SIMD_MATH

            // two-argument functions
            err += check_simd_math<V>("pow",
                        [] (V const& v) { return Math::pow(v, V(T(2.5))); },
                        [] (T const v)  { return std::pow(v, T(2.5)); },
                        T(0.1), T(10.0));
            err += check_simd_math<V>("atan2",
                        [] (V const& v) { return Math::atan2(v, V(T(2.0))); },
                        [] (T const v)  { return std::atan2(v, T(2.0)); },
                        T(-10.0), T(10.0));
            err += check_simd_math<V>("hypot",
                        [] (V const& v) { return Math::hypot(v, V(T(3.0))); },
                        [] (T const v)  { return std::hypot(v, T(3.0)); },
                        T(-10.0), T(10.0));

            // sincos and sincospi return both results at once
            err += check_simd_math<V>("sincos (sin)",
                        [] (V const& v) { return Math::sincos(v).first; },
                        [] (T const v)  { return std::sin(v); },
                        T(-6.0), T(6.0));
            err += check_simd_math<V>("sincos (cos)",
                        [] (V const& v) { return Math::sincos(v).second; },
                        [] (T const v)  { return std::cos(v); },
                        T(-6.0), T(6.0));
            err += check_simd_math<V>("sincospi (sin)",
                        [] (V const& v) { return Math::sincospi(v).first; },
                        [] (T const v)  { return std::sin(Math::pi<T>() * v); },
                        T(-2.0), T(2.0));
            err += check_simd_math<V>("sincospi (cos)",
                        [] (V const& v) { return Math::sincospi(v).second; },
                        [] (T const v)  { return std::cos(Math::pi<T>() * v); },
                        T(-2.0), T(2.0));

            Print() << "amrex::Math SIMD transcendentals ("
#   ifdef AMREX_SIMD_VECMATH
                    << "vector math library enabled"
#   else
                    << "vector math library not available, lane-wise fallback"
#   endif
                    << ", width " << int(V::size()) << "): "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
#else
            Print() << "amrex::Math scalar transcendentals: "
                    << (err == 0 ? "PASSED" : "FAILED") << "\n";
#endif
            nerrors += err;
        }

        // ================================================================
        // Final report
        // ================================================================
        if (nerrors > 0) {
            amrex::Finalize();
            Abort("SIMD test FAILED with "
                  + std::to_string(nerrors) + " error(s)");
        }
        Print() << "\nAll SIMD tests PASSED.\n";
    }
    amrex::Finalize();
}
