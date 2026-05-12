#include <AMReX.H>
#include <AMReX_Dim3.H>
#include <AMReX_RealBox.H>
#include <AMReX_RealVect.H>

using namespace amrex;

namespace
{

void test_contains_closed_interval ()
{
    constexpr Real delta = Real(1.e-12);

    Real lo_arr[AMREX_SPACEDIM] = {AMREX_D_DECL(Real(-2.0), Real(-1.0), Real(-0.5))};
    Real hi_arr[AMREX_SPACEDIM] = {AMREX_D_DECL(Real( 3.0), Real( 4.0), Real( 5.5))};
    RealBox box(lo_arr, hi_arr);

    Real lo_outside[AMREX_SPACEDIM] =
        {AMREX_D_DECL(lo_arr[0] - delta, lo_arr[1], lo_arr[2])};
    Real hi_outside[AMREX_SPACEDIM] =
        {AMREX_D_DECL(hi_arr[0] + delta, hi_arr[1], hi_arr[2])};

    AMREX_ALWAYS_ASSERT(box.contains(lo_arr));
    AMREX_ALWAYS_ASSERT(box.contains(hi_arr));
    AMREX_ALWAYS_ASSERT(!box.contains(lo_outside));
    AMREX_ALWAYS_ASSERT(!box.contains(hi_outside));
    AMREX_ALWAYS_ASSERT(box.contains(lo_outside, delta));
    AMREX_ALWAYS_ASSERT(box.contains(hi_outside, delta));

    XDim3 lo_dim3{.x = Real(-2.0), .y = Real(-1.0), .z = Real(-0.5)};
    XDim3 hi_dim3{.x = Real( 3.0), .y = Real( 4.0), .z = Real( 5.5)};
    XDim3 lo_dim3_outside{.x = Real(-2.0) - delta, .y = Real(-1.0), .z = Real(-0.5)};
    XDim3 hi_dim3_outside{.x = Real( 3.0) + delta, .y = Real( 4.0), .z = Real( 5.5)};

    AMREX_ALWAYS_ASSERT(box.contains(lo_dim3));
    AMREX_ALWAYS_ASSERT(box.contains(hi_dim3));
    AMREX_ALWAYS_ASSERT(!box.contains(lo_dim3_outside));
    AMREX_ALWAYS_ASSERT(!box.contains(hi_dim3_outside));
    AMREX_ALWAYS_ASSERT(box.contains(lo_dim3_outside, delta));
    AMREX_ALWAYS_ASSERT(box.contains(hi_dim3_outside, delta));

    RealVect lo_rv(lo_arr);
    RealVect hi_rv(hi_arr);
    RealVect lo_rv_outside(lo_outside);
    RealVect hi_rv_outside(hi_outside);

    AMREX_ALWAYS_ASSERT(box.contains(lo_rv));
    AMREX_ALWAYS_ASSERT(box.contains(hi_rv));
    AMREX_ALWAYS_ASSERT(!box.contains(lo_rv_outside));
    AMREX_ALWAYS_ASSERT(!box.contains(hi_rv_outside));
    AMREX_ALWAYS_ASSERT(box.contains(lo_rv_outside, delta));
    AMREX_ALWAYS_ASSERT(box.contains(hi_rv_outside, delta));

    AMREX_ALWAYS_ASSERT(box.contains(box));
    AMREX_ALWAYS_ASSERT(box.contains(RealBox(lo_arr, hi_arr)));
    AMREX_ALWAYS_ASSERT(!box.contains(RealBox(lo_outside, hi_arr)));
    AMREX_ALWAYS_ASSERT(!box.contains(RealBox(lo_arr, hi_outside)));
    AMREX_ALWAYS_ASSERT(box.contains(RealBox(lo_outside, hi_arr), delta));
    AMREX_ALWAYS_ASSERT(box.contains(RealBox(lo_arr, hi_outside), delta));
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    test_contains_closed_interval();

    amrex::Finalize();
}
