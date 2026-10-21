// Header under test first, to test transient includes
#include <AMReX_Algorithm.H>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#ifdef AMREX_USE_SIMD
#   include <AMReX_SIMD.H>
#endif

#include <type_traits>

//
// Regression test for the variadic amrex::min/max (AMReX_Algorithm.H).
//
// With three or more arguments, amrex::min/max reduce pairwise through the
// two-argument min/max of the argument type. For IntVect, RealVect and SIMD
// vectors, that is an element-wise overload that returns a new object by
// value. If the variadic overload returned a reference, it would refer to
// that temporary after it was destroyed (GCC: segfault, Clang/ASan:
// stack-use-after-scope).
//

using namespace amrex::literals;

namespace
{

using IV3 = amrex::IntVectND<3>;
using RV3 = amrex::RealVectND<3>;

// Constant evaluation rejects reads of destroyed temporaries, so these fail
// to compile if the variadic overload returns a dangling reference.
static_assert(amrex::min(IV3(1,5,9), IV3(4,2,8), IV3(7,6,3)) == IV3(1,2,3));
static_assert(amrex::max(IV3(1,5,9), IV3(4,2,8), IV3(7,6,3)) == IV3(7,6,9));
static_assert(amrex::min(IV3(1,5,9), IV3(4,2,8), IV3(7,6,3), IV3(0,9,9)) == IV3(0,2,3));
static_assert(amrex::min(3, 1, 2) == 1);
static_assert(amrex::max(3, 1, 4, 2) == 4);

static_assert(!std::is_reference_v<decltype(amrex::min(std::declval<IV3>(),
                                                       std::declval<IV3>(),
                                                       std::declval<IV3>()))>);
static_assert(!std::is_reference_v<decltype(amrex::max(std::declval<RV3>(),
                                                       std::declval<RV3>(),
                                                       std::declval<RV3>()))>);

void test_scalar ()
{
    amrex::Real const a = 3.0_rt, b = -1.0_rt, c = 2.0_rt, d = 7.0_rt;

    AMREX_ALWAYS_ASSERT(amrex::min(a, b, c) == b);
    AMREX_ALWAYS_ASSERT(amrex::max(a, b, c) == a);
    AMREX_ALWAYS_ASSERT(amrex::min(a, b, c, d) == b);
    AMREX_ALWAYS_ASSERT(amrex::max(a, b, c, d) == d);

    // arguments that are temporaries in the caller
    AMREX_ALWAYS_ASSERT(amrex::min(a + 1.0_rt, b + 1.0_rt, c + 1.0_rt) == 0.0_rt);
}

void test_intvect ()
{
    IV3 const a(1,5,9), b(4,2,8), c(7,6,3), d(0,9,9);

    AMREX_ALWAYS_ASSERT(amrex::min(a, b, c) == IV3(1,2,3));
    AMREX_ALWAYS_ASSERT(amrex::max(a, b, c) == IV3(7,6,9));
    AMREX_ALWAYS_ASSERT(amrex::min(a, b, c, d) == IV3(0,2,3));
    AMREX_ALWAYS_ASSERT(amrex::max(a, b, c, d) == IV3(7,9,9));

    // binding the result to a reference extends the lifetime of the returned value
    IV3 const& r = amrex::min(a, b, c);
    AMREX_ALWAYS_ASSERT(r == IV3(1,2,3));
}

void test_realvect ()
{
    RV3 const a(1.0_rt, 5.0_rt, 9.0_rt), b(4.0_rt, 2.0_rt, 8.0_rt), c(7.0_rt, 6.0_rt, 3.0_rt);

    AMREX_ALWAYS_ASSERT(amrex::min(a, b, c) == RV3(1.0_rt, 2.0_rt, 3.0_rt));
    AMREX_ALWAYS_ASSERT(amrex::max(a, b, c) == RV3(7.0_rt, 6.0_rt, 9.0_rt));
}

#ifdef AMREX_USE_SIMD
void test_simd ()
{
    using V = amrex::simd::SIMDReal<>;
    V const a(1.0_rt), b(4.0_rt), c(-2.0_rt);

    V const mn = amrex::min(a, b, c);
    V const mx = amrex::max(a, b, c);
    for (std::size_t i = 0; i < V::size(); ++i) {
        AMREX_ALWAYS_ASSERT(mn[i] == -2.0_rt);
        AMREX_ALWAYS_ASSERT(mx[i] == 4.0_rt);
    }
}
#endif

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    test_scalar();
    test_intvect();
    test_realvect();
#ifdef AMREX_USE_SIMD
    test_simd();
#endif

    amrex::Print() << "MinMax test passed.\n";

    amrex::Finalize();
}
