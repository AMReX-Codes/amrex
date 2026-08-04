// Header under test first, to test transient includes
#include <AMReX_CompensatedSum.H>

#include <AMReX.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_Random.H>

#include <cmath>
#include <limits>

//
// Regression test for amrex::compensatedAdd (AMReX_CompensatedSum.H).
//
// Kahan compensation relies on (t - sum) - y evaluating to the rounding
// error dropped by sum + y, not to the exact-arithmetic value of 0.
// Compilers given license to treat floating-point addition as associative
// (e.g. -ffast-math, -fassociative-math) are free to simplify that
// expression to a literal 0, which silently turns compensatedAdd into a
// plain, uncompensated sum. AMReX enables -ffast-math by default for
// several compilers, so this test exists to catch that regression directly,
// on whatever compiler/flags CI happens to run it with.
//
// The test sums many O(1) increments onto a "big" value chosen so that the
// increments are just below the representable ULP at that magnitude: naive
// summation provably loses them, so if compensatedAdd is silently reduced
// to a plain sum, this test can tell the difference from the exact answer.
// Values derived from runtime state (ParallelDescriptor::MyProc(), a
// ParallelFor loop index) rather than compile-time constants, so the
// compiler's constant folder cannot mask the same bug by evaluating the
// whole computation ahead of time.
//

using namespace amrex;

namespace {

template <typename T>
void test_compensated_add ()
{
    int const rank = ParallelDescriptor::MyProc() + 1; // avoid a runtime-constant 0
    int const n = 12;

    // Two binary orders above the mantissa limit: the representable ULP
    // there is 2, so repeatedly adding 1 is a worst-case cancellation for
    // naive summation, while compensatedAdd must still recover it exactly.
    T const big = std::ldexp(T(1), std::numeric_limits<T>::digits + 1) * T(rank);
    T const increment = T(1) * T(rank);
    T const expected = T(n) * increment;

    // Host check
    {
        T sum = big;
        T compensation = T(0);
        for (int i = 0; i < n; ++i) {
            amrex::compensatedAdd(sum, compensation, increment);
        }
        T const result = sum - big;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(result == expected,
            "compensatedAdd lost precision on the host -- the compensation "
            "term may have been optimized away (e.g. by -ffast-math)");

        // Confirm this is actually a stringent test: naive summation at
        // this magnitude must fail, or the check above is not meaningful.
        volatile T naive_sum = big;
        for (int i = 0; i < n; ++i) { naive_sum = naive_sum + increment; }
        T const naive_result = naive_sum - big;
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(naive_result != expected,
            "test setup error: naive summation did not lose precision at "
            "this magnitude, so this test would not catch a defeated "
            "compensation term");
    }

    // Device check: same computation, run through a GPU kernel when built
    // with a GPU backend (falls back to the host path otherwise).
    {
        Gpu::DeviceScalar<T> d_result(T(0));
        T* p_result = d_result.dataPtr();

        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept
        {
            T sum = big;
            T compensation = T(0);
            for (int i = 0; i < n; ++i) {
                amrex::compensatedAdd(sum, compensation, increment);
            }
            *p_result = sum - big;
        });
        Gpu::synchronize();

        T const result = d_result.dataValue();
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(result == expected,
            "compensatedAdd lost precision on the device -- the "
            "compensation term may have been optimized away");
    }
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        test_compensated_add<float>();
        test_compensated_add<double>();
        amrex::Print() << "compensatedAdd retains precision under this compiler's default flags\n";
    }
    amrex::Finalize();
}
