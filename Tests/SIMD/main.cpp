// Test coverage for AMReX SIMD single-source design:
//   SIMDindex, ParallelForSIMD, load_1d, store_1d,
//   Vectorized, is_vectorized, is_nth_arg_non_const
//
// Compiles and runs with AMReX_SIMD=OFF (scalar fallback, GPU) and ON (CPU SIMD).

#include <AMReX.H>
#include <AMReX_GpuLaunch.H>
#include <AMReX_Print.H>
#include <AMReX_SIMD.H>
#include <AMReX_Vector.H>

#include <cmath>
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
            std::fill(y_data.begin(), y_data.end(), ParticleReal(100));

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
            std::fill(y_data.begin(), y_data.end(), ParticleReal(10));

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
